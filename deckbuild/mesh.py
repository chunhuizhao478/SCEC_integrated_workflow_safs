"""mesh.py -- ingest and GATE a fault mesh.  This module does not BUILD one.

Mesh generation and improvement is an interactive, judgement-heavy craft: eighteen forked
campaign folders in the legacy tree prove it does not reduce to a function signature.  It
is handed to the `code-mesh-build-improve` skill in `skills/`, driven by Claude.  See
MESHING.md for the one-page map.

What lives HERE is the mechanical half:

  * `.msh` (Gmsh v2.2) <-> `.puml.h5` conversion, one exact format contract;
  * gates A-E, names preserved from the legacy `verify_*.py`;
  * gates F and G, two ingestion checks the legacy verifiers assumed rather than checked.

The PUML contract, reverse-engineered from a PUMGen reference:

    geometry (Nnode,3) f8 | connect (Ntet,4) u8 0-based
    boundary (Ntet,) i4   | group (Ntet,) i4
    attrs: boundary-format='i32', topology-format='geometric'
    BC packing: byte i of `boundary` holds the BC code of tet face i,
                boundary = sum_i code_i << (8*i)
    face order: f0={0,2,1} f1={0,1,3} f2={1,2,3} f3={0,3,2}

Each tet's local vertex order is normalised to a POSITIVE signed volume.  Gmsh does not
guarantee orientation, and inverted tets make SeisSol diverge with Inf/NaN bulk energy at
t ~ 1 s.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from deckbuild.config import MeshSpec, Project
from deckbuild.contract import Artifact, GateReport, HARD, Stage, WARN
from deckbuild.geometry import FACE_VERTS, GeometryError, load_fault

__all__ = ["MeshStage", "MeshError", "read_puml", "write_puml", "msh_to_puml",
           "orient_tets_positive", "pack_boundary", "unpack_boundary", "face_table"]


class MeshError(ValueError):
    """Raised for a malformed or unconvertible mesh."""


# --------------------------------------------------------------------------- format
def unpack_boundary(bnd, fmt: str = "i32") -> np.ndarray:
    """(Ntet,) packed -> (Ntet,4) per-face BC codes."""
    bnd = np.asarray(bnd)
    if bnd.ndim == 2:
        return bnd.astype(np.int64)
    bits = {"i32": 8, "i64": 16}.get(fmt)
    if bits is None:
        raise MeshError(f"unsupported boundary-format {fmt!r}; expected 'i32' or 'i64'")
    b = bnd.astype(np.int64)
    mask = (1 << bits) - 1
    return np.stack([(b >> (bits * k)) & mask for k in range(4)], axis=1)


def pack_boundary(bc, fmt: str = "i32") -> np.ndarray:
    """(Ntet,4) per-face BC codes -> the packed integer SeisSol reads."""
    bc = np.asarray(bc, np.int64)
    if bc.ndim != 2 or bc.shape[1] != 4:
        raise MeshError(f"bc must be (Ntet, 4), got {bc.shape}")
    bits = {"i32": 8, "i64": 16}.get(fmt)
    if bits is None:
        raise MeshError(f"unsupported boundary-format {fmt!r}")
    if bc.max(initial=0) >= (1 << bits):
        raise MeshError(
            f"a BC code {int(bc.max())} does not fit in {bits} bits; use boundary-format "
            f"'i64' or renumber the codes")
    out = np.zeros(len(bc), np.int64)
    for k in range(4):
        out |= bc[:, k] << (bits * k)
    return out.astype(np.int32 if fmt == "i32" else np.int64)


def orient_tets_positive(pts, tets):
    """Swap two vertices of any negative-volume tet.  Returns (tets, n_flipped).

    An inverted tet makes SeisSol blow up at t ~ 1 s with Inf/NaN bulk energy, and gmsh
    does not guarantee orientation, so this is not optional.
    """
    tets = np.array(tets, copy=True)
    p = np.asarray(pts, float)[tets.astype(np.int64)]
    vol = np.einsum("ij,ij->i", p[:, 1] - p[:, 0],
                    np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]))
    flip = vol < 0
    if flip.any():
        tets[flip] = tets[flip][:, [0, 2, 1, 3]]
    return tets, int(flip.sum())


def signed_volumes(pts, tets) -> np.ndarray:
    p = np.asarray(pts, float)[np.asarray(tets).astype(np.int64)]
    return np.einsum("ij,ij->i", p[:, 1] - p[:, 0],
                     np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0])) / 6.0


def face_table(tets):
    """(4*Ntet, 3) sorted vertex triples, one row per tet face, in face order."""
    tets = np.asarray(tets, np.int64)
    tri = np.concatenate([tets[:, list(fv)] for fv in FACE_VERTS], axis=0)
    return np.sort(tri, axis=1)


def count_face_multiplicity(tets):
    """How many tets each face belongs to, via a lexsort + searchsorted.

    NOT a packed int64 key and NOT a Counter: the legacy packed key overflows above ~2.1M
    nodes, and the Counter-based validator needed >12 GB on a production mesh.
    """
    tri = face_table(tets)
    view = np.ascontiguousarray(tri).view(
        np.dtype((np.void, tri.dtype.itemsize * tri.shape[1]))).ravel()
    uniq, inv, counts = np.unique(view, return_inverse=True, return_counts=True)
    return counts[inv], uniq.size


# --------------------------------------------------------------------------- read/write
def read_puml(path):
    """geometry, connect, per-face BC codes, group."""
    import h5py
    path = Path(path)
    if not path.is_file():
        raise MeshError(f"mesh not found: {path.resolve()}")
    with h5py.File(str(path), "r") as f:
        for req in ("geometry", "connect", "boundary"):
            if req not in f:
                raise MeshError(f"{path}: missing dataset {req!r}; have {sorted(f)}")
        geom = f["geometry"][:]
        conn = f["connect"][:].astype(np.int64)
        fmt = f.attrs.get("boundary-format", "i32")
        if isinstance(fmt, bytes):
            fmt = fmt.decode()
        bnd = f["boundary"][:]
        group = f["group"][:] if "group" in f else np.ones(len(conn), np.int32)
    return geom, conn, unpack_boundary(bnd, fmt), group


def write_puml(path, geom, conn, bc, group=None, fmt: str = "i32") -> Path:
    """Write a SeisSol PUML/HDF5 mesh, with the exact attribute and dtype contract."""
    import h5py
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    conn = np.asarray(conn)
    if group is None:
        group = np.ones(len(conn), np.int32)
    with h5py.File(str(path), "w") as f:
        f.create_dataset("geometry", data=np.asarray(geom, np.float64))
        f.create_dataset("connect", data=conn.astype(np.uint64))
        f.create_dataset("boundary", data=pack_boundary(bc, fmt))
        f.create_dataset("group", data=np.asarray(group, np.int32))
        f.attrs["boundary-format"] = fmt
        f.attrs["topology-format"] = "geometric"
    return path


def msh_to_puml(msh_path, out_path, tag_to_bc, fmt: str = "i32", validate: bool = True):
    """Gmsh v2.2 tet .msh -> PUML/HDF5.

    tag_to_bc maps the msh physical surface tags to SeisSol BC codes.  An unmapped tag is
    an ERROR: silently dropping one produces a mesh with a missing boundary condition.
    """
    import meshio
    msh_path = Path(msh_path)
    if not msh_path.is_file():
        raise MeshError(f"msh not found: {msh_path.resolve()}")
    try:
        m = meshio.read(str(msh_path))
    except (Exception, SystemExit) as exc:
        # SystemExit, not just Exception: meshio calls sys.exit() on some malformed
        # headers, and SystemExit derives from BaseException.  Letting it through would
        # kill the caller's interpreter instead of raising an actionable message -- the
        # same "a library must not exit" rule as geometry._load_puml.
        raise MeshError(
            f"{msh_path}: could not read as Gmsh.  If this is format 4.x, re-export as "
            f"v2.2 ASCII (`gmsh ... -format msh22`).  MFEM's reader mis-parses v4 and "
            f"reports the misleading 'vertices indices are not unique'.  ({exc})"
        ) from None

    tets = None
    tris, tri_tags = [], []
    for blk, tagblk in zip(m.cells, _cell_tags(m)):
        if blk.type == "tetra":
            tets = blk.data.astype(np.int64)
        elif blk.type == "triangle":
            tris.append(blk.data.astype(np.int64))
            tri_tags.append(tagblk)
    if tets is None:
        raise MeshError(f"{msh_path}: no tetrahedra")
    if not tris:
        raise MeshError(f"{msh_path}: no tagged boundary triangles")
    tris = np.concatenate(tris, axis=0)
    tri_tags = np.concatenate(tri_tags, axis=0)

    present = sorted(int(t) for t in np.unique(tri_tags))
    unmapped = [t for t in present if t not in tag_to_bc]
    if unmapped:
        raise MeshError(
            f"{msh_path}: physical tag(s) {unmapped} are not in tag_to_bc "
            f"{dict(tag_to_bc)}; an unmapped tag would silently lose a boundary condition")

    pts = np.asarray(m.points, float)
    tets, n_flipped = orient_tets_positive(pts, tets)

    # Match each tagged triangle to (tet, local face) via a sorted-triple lookup.
    ftab = face_table(tets)
    order = np.lexsort(ftab.T[::-1])
    sorted_faces = ftab[order]
    key = np.sort(tris, axis=1)
    lo, hi = _match_rows(sorted_faces, key)
    bc = np.zeros((len(tets), 4), np.int64)
    n_matched = 0
    for j in range(len(key)):
        if hi[j] <= lo[j]:
            continue
        # A tagged triangle can match MORE THAN ONE tet face: an interior fault triangle
        # is shared by exactly two tets, and BOTH sides must carry the BC.  Taking only
        # the first match leaves one side at 0, which gate C then reports as "fault face
        # not interior" -- correct, but the conversion should not produce it.
        for i in range(lo[j], hi[j]):
            flat = order[i]
            bc[flat % len(tets), flat // len(tets)] = tag_to_bc[int(tri_tags[j])]
        n_matched += 1
    if n_matched != len(tris):
        raise MeshError(
            f"{msh_path}: only {n_matched} of {len(tris)} tagged triangles matched a tet "
            f"face; the surface and volume meshes are inconsistent")

    write_puml(out_path, pts, tets, bc, fmt=fmt)
    info = {"n_nodes": len(pts), "n_tets": len(tets), "n_tagged_tris": len(tris),
            "n_flipped": n_flipped, "tags": present}
    if validate:
        geom, conn, bc2, _ = read_puml(out_path)
        if not np.array_equal(bc2, bc):
            raise MeshError("BC round-trip failed: the packed boundary did not decode "
                            "back to what was written")
        if (signed_volumes(geom, conn) <= 0).any():
            raise MeshError("inverted tets survived orientation normalisation")
    return Path(out_path), info


def _cell_tags(m):
    """meshio physical tags per cell block, or zeros when absent."""
    for key in ("gmsh:physical", "medit:ref", "cell_tags"):
        if key in m.cell_data:
            return [np.asarray(d, np.int64) for d in m.cell_data[key]]
    return [np.zeros(len(b.data), np.int64) for b in m.cells]


def _match_rows(sorted_rows, keys):
    """Half-open [lo, hi) range of each key in a lexsorted (N,3) int array.

    A RANGE, not a single index: an interior fault triangle matches two tet faces and
    both must be tagged.  np.rec (not the removed np.core.records) gives a totally
    ordered dtype that searchsorted can bisect consistently with np.lexsort.
    """
    st = np.rec.fromarrays(sorted_rows.T, names="a,b,c")
    kt = np.rec.fromarrays(keys.T, names="a,b,c")
    lo = np.searchsorted(st, kt, side="left")
    hi = np.searchsorted(st, kt, side="right")
    return lo, hi


# --------------------------------------------------------------------------- the stage
class MeshStage(Stage):
    name = "mesh"

    def build(self, cfg: Project, out_dir, *, msh: str | Path | None = None,
              mesh: str | None = None, out_name: str | None = None) -> Artifact:
        """INGEST a mesh: convert .msh -> .puml.h5, or register an existing .puml.h5.

        This stage never generates geometry -- see MESHING.md and skills/.
        """
        mspec = cfg.mesh(mesh or cfg.default_mesh)
        if msh is None:
            p = cfg.resolve_path(mspec.path)
            if not Path(p).is_file():
                raise MeshError(
                    f"no mesh at {p}.  This stage INGESTS a mesh; it does not build one. "
                    f"Build yours with the code-mesh-build-improve skill in skills/ "
                    f"(MESHING.md is the one-page map), then point "
                    f"meshes.{mesh or cfg.default_mesh}.path at it.")
            return Artifact.of(p, kind="mesh", params={"mesh": mesh or cfg.default_mesh,
                                                       "ingested": True})
        out = Path(out_dir) / (out_name or (Path(msh).stem + ".puml.h5"))
        path, info = msh_to_puml(msh, out, dict(mspec.tag_to_bc))
        return Artifact.of(path, kind="mesh",
                           params={"mesh": mesh or cfg.default_mesh,
                                   "tag_to_bc": dict(mspec.tag_to_bc), **info},
                           provenance={"source_msh": str(msh),
                                       "descriptor_sha256": cfg.sha256()})

    def verify(self, cfg: Project, artifact: Artifact, *, mesh: str | None = None,
               material_nc=None, stress_nc=None, friction_nc=None,
               fault_edge_max_m: float | None = None, f_gate_hz: float | None = None,
               **_) -> GateReport:
        """Gates A-E (legacy names) plus the two ingestion gates F and G."""
        mspec = cfg.mesh(mesh or cfg.default_mesh)
        p = Path(artifact.path)
        geom, conn, bc, _grp = read_puml(p)
        # Thresholds are PROJECT properties (see config.MeshGates), never literals here.
        if fault_edge_max_m is None:
            fault_edge_max_m = mspec.gates.fault_edge_max_m
        if f_gate_hz is None:
            f_gate_hz = mspec.gates.f_gate_hz
        rep = GateReport(f"mesh A-G ({p.name})")

        lo, hi = geom.min(axis=0), geom.max(axis=0)
        n_fault_faces = int((bc == mspec.fault_bc).sum())
        rep.add("A", True,
                f"{len(geom):,} nodes, {len(conn):,} tets, {n_fault_faces:,} BC-"
                f"{mspec.fault_bc} faces; bbox x[{lo[0]:.0f}, {hi[0]:.0f}] "
                f"y[{lo[1]:.0f}, {hi[1]:.0f}] z[{lo[2]:.0f}, {hi[2]:.0f}]",
                severity=WARN)

        try:
            fault = load_fault(mspec, cfg.data_dir, strike=cfg.strike)
        except GeometryError as exc:
            rep.add("B", False, f"cannot extract the fault: {exc}")
            return rep

        # B -- the fault edge ceiling.
        pts = fault.cent
        e = _max_edge_of_facets(geom, conn, bc, mspec.fault_bc)
        rep.add("B", e < fault_edge_max_m,
                f"fault edge max = {e:.1f} m (limit {fault_edge_max_m:.0f} m)")

        # C -- BC round trip: every tagged face is either on the geometric hull (x1) or
        # is the fault (x2, interior).
        mult, _ = count_face_multiplicity(conn)
        mult = mult.reshape(4, len(conn)).T
        tagged = bc > 0
        hull_ok = int(((mult == 1) & tagged & (bc != mspec.fault_bc)).sum())
        fault_int = int(((mult == 2) & (bc == mspec.fault_bc)).sum())
        fault_bad = int(((mult != 2) & (bc == mspec.fault_bc)).sum())
        outer_bad = int(((mult != 1) & tagged & (bc != mspec.fault_bc)).sum())
        rep.add("C", fault_bad == 0 and outer_bad == 0,
                f"BC round-trip: {hull_ok:,} boundary faces x1, {fault_int:,} fault faces "
                f"x2; {fault_bad} fault face(s) not interior, {outer_bad} boundary face(s) "
                f"not on the hull")

        # D -- orientation.
        vol = signed_volumes(geom, conn)
        n_inv = int((vol <= 0).sum())
        rep.add("D", n_inv == 0,
                f"inverted tets: {n_inv} (an inverted tet makes SeisSol diverge with "
                f"Inf/NaN bulk energy at t ~ 1 s)")

        # E -- the volume resolution gate.  Needs a material nc; SKIP is not PASS.
        if material_nc is None:
            rep.skip("E", "no material nc supplied; the f = Vs/dx resolution gate did "
                          "not run")
        else:
            from deckbuild.asagi import trilinear_sample
            bary = geom[conn.astype(np.int64)].mean(axis=1)
            smp = trilinear_sample(material_nc, bary[:, 0], bary[:, 1], bary[:, 2],
                                   fields=["rho", "mu"])
            with np.errstate(divide="ignore", invalid="ignore"):
                vs = np.sqrt(np.where(smp["rho"] > 0, smp["mu"] / smp["rho"], np.nan))
            dx = _max_edge_per_tet(geom, conn)
            f = vs / dx
            n_below = int(np.nansum(f < f_gate_hz))
            rep.add("E", n_below == 0,
                    f"resolution f = Vs/dx: worst {np.nanmin(f):.4f} Hz, median "
                    f"{np.nanmedian(f):.3f}; {n_below:,} of {len(f):,} cells below "
                    f"{f_gate_hz}")

        # F -- every fault facet inside the field hulls.
        from deckbuild.asagi import hull_containment
        any_field = False
        for label, nc in (("stress", stress_nc), ("friction", friction_nc),
                          ("material", material_nc)):
            if nc is None:
                continue
            any_field = True
            hull_containment(nc, pts, label=f"fault vs {label} nc", gate=f"F-{label}",
                             severity=HARD, report=rep)
        if not any_field:
            rep.skip("F", "no field nc supplied; hull containment not checked")

        # G -- the pickpoint partition heuristic.  A WARN by construction: the partition
        # is not known until SeisSol runs.
        rep.skip("G", "pickpoint/receiver partition-boundary screen needs a receiver "
                      "list; run it from DeckStage once receivers exist")
        return rep


def _fault_tris(geom, conn, bc, fault_bc):
    tris = np.concatenate(
        [conn[bc[:, k] == fault_bc][:, list(FACE_VERTS[k])] for k in range(4)], axis=0)
    _, keep = np.unique(np.sort(tris, axis=1), axis=0, return_index=True)
    return tris[np.sort(keep)]


def _max_edge_of_facets(geom, conn, bc, fault_bc) -> float:
    tris = _fault_tris(geom, conn, bc, fault_bc)
    if len(tris) == 0:
        return float("inf")
    p = geom[tris.astype(np.int64)]
    e = np.stack([np.linalg.norm(p[:, (i + 1) % 3] - p[:, i], axis=1) for i in range(3)],
                 axis=1)
    return float(e.max())


def _max_edge_per_tet(geom, conn) -> np.ndarray:
    p = geom[conn.astype(np.int64)]
    pairs = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
    return np.max(np.stack([np.linalg.norm(p[:, b] - p[:, a], axis=1)
                            for a, b in pairs], axis=1), axis=1)
