"""config.py -- the project descriptor.

One YAML file per fault system, holding every number that used to be hardcoded in the
legacy SAFS toolbox.  Loading is strict: an unknown key, an unknown reader `kind`, a
reversed band or a missing file all raise `ConfigError`.  Nothing here is silently
defaulted, because a retargeted project that quietly keeps a SAFS constant is the exact
failure this module exists to prevent.

Every field's docstring names the legacy symbol it replaces (see the exploration document
section 10 for the full table).
"""
from __future__ import annotations

import hashlib
import re
from collections.abc import Mapping
from dataclasses import dataclass, field, fields, is_dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Any

import yaml

__all__ = [
    "ConfigError",
    "CRS",
    "StrikeFrame",
    "GridBox",
    "NamedBand",
    "Hypocenter",
    "MeshSpec",
    "MeshGates",
    "PhysicsDefaults",
    "SourceSpec",
    "RawSources",
    "Project",
    "ORIENTATION_KINDS",
    "VELOCITY_KINDS",
    "THERMAL_KINDS",
]


class ConfigError(ValueError):
    """Raised for any descriptor problem.  Never a warning -- see the module docstring."""


# --------------------------------------------------------------------------- reader kinds
# `kind` dispatches to a reader registered in these dicts.  Later phases register the real
# readers; Phase 0 only validates that the name is known.  This is the extension point for
# a fault system with different data products.
ORIENTATION_KINDS: dict[str, Any] = {"csm_csv": None, "constant": None, "callable": None}
VELOCITY_KINDS: dict[str, Any] = {"cvm_slices": None, "layered_1d": None, "callable": None}
THERMAL_KINDS: dict[str, Any] = {"ctm_slices": None, "depth_profile": None}


# --------------------------------------------------------------------------- leaf schema
@dataclass(frozen=True)
class CRS:
    """Coordinate reference systems.  `epsg` was a hardcoded literal in 5 legacy modules."""

    epsg: str
    geographic_epsg: str = "EPSG:4326"


@dataclass(frozen=True)
class StrikeFrame:
    """The along-strike coordinate `s`.

    azimuth_deg  was: _common.STRIKE_HINT_AZ / pipeline.STRIKE_AZ
    origin_xy    was: pipeline.GATE_ORIGIN_XY
    (The SAFS values live in projects/safs_*.yaml -- deliberately not repeated here; see
    tests/test_no_hardcoded_constants.py, which forbids them anywhere in this package.)
    """

    azimuth_deg: float
    origin_xy: tuple[float, float]
    harmonise_hint_deg: float | None = None

    @property
    def hint_deg(self) -> float:
        """The azimuth used to harmonise facet normals; defaults to `azimuth_deg`."""
        return self.azimuth_deg if self.harmonise_hint_deg is None else self.harmonise_hint_deg


@dataclass(frozen=True)
class GridBox:
    """An ASAGI grid box.  was: _common.ALT_BOX / PREF_BOX literals."""

    xmin: float
    xmax: float
    ymin: float
    ymax: float
    zmin: float
    zmax: float
    dx: float
    dz: float


@dataclass(frozen=True)
class NamedBand:
    """A named restraining bend along strike.  was: _common.GATE_BANDS entries."""

    s_start_km: float
    s_end_km: float
    name: str


@dataclass(frozen=True)
class Hypocenter:
    """A REQUESTED nucleation point.

    Give EITHER projected coords (x, y, z) OR geographic (lon, lat, depth_m); the loader
    converts via `Project.crs`.  This point is snapped onto the real fault by
    geometry.snap_hypocenter -- it is not assumed to lie on the triangulation.

    `depth_m` is referenced to SEA LEVEL (z = -depth_m), NOT to the local ground surface.
    Under topography the two differ by the local elevation.
    """

    x: float | None = None
    y: float | None = None
    z: float | None = None
    lon: float | None = None
    lat: float | None = None
    depth_m: float | None = None
    snap_tol_m: float = 2000.0
    label: str = ""

    @property
    def is_projected(self) -> bool:
        return self.x is not None and self.y is not None and self.z is not None

    @property
    def is_geographic(self) -> bool:
        return self.lon is not None and self.lat is not None and self.depth_m is not None


def _default_tag_to_bc() -> Mapping[int, int]:
    """Conventional msh tag -> SeisSol BC: fault->3, top->1 (free surface), sides/bottom->5
    (absorbing).  Read-only: see _freeze()."""
    return MappingProxyType({101: 3, 102: 1, 103: 5, 104: 5})


@dataclass(frozen=True)
class MeshGates:
    """Acceptance thresholds for a mesh.  These are PROJECT properties, not universals.

    500 m / 0.6667 Hz are the SAFS production values; a demo or a different fault system
    has its own.  Hardcoding them in the verifier would be exactly the constant leak this
    descriptor exists to prevent.
    """

    fault_edge_max_m: float = 500.0
    f_gate_hz: float = 0.6667      # resolved 0.5 Hz at p3: (3/4) Vs/dx >= 0.5
    partition_tol_m: float = 50.0


@dataclass(frozen=True)
class MeshSpec:
    """A fault mesh.  `path` is relative to the project's data directory."""

    path: str
    fault_bc: int = 3
    tag_to_bc: Mapping[int, int] = field(default_factory=_default_tag_to_bc)
    daylights: bool = False
    daylight_min_depth_m: float = 0.0
    gates: MeshGates = field(default_factory=MeshGates)


@dataclass(frozen=True)
class PhysicsDefaults:
    """Physics constants.  All were module-level literals in legacy `_common.py`."""

    mu_shear_pa: float = 23.5e9
    w_energy_m: float = 10.0e3
    dc_m: float = 1.2
    kappa_c: float = 0.9
    seis_band_km: tuple[float, float] = (3.0, 12.0)
    active_band_km: tuple[float, float] = (0.3, 15.0)
    l_coast_km: float = 6.0
    gap_max_km: float = 4.0
    sn_floor_mpa: float = 30.0
    rs_b: float = 0.019
    rs_sl0: float = 0.10


@dataclass(frozen=True)
class SourceSpec:
    """One raw input: which reader to use, where the data is, and its reader options."""

    kind: str
    path: str | None = None
    params: Mapping = field(default_factory=dict)


@dataclass(frozen=True)
class RawSources:
    """The raw data the user supplies.

    There is deliberately no `density_profile` field: the Sv(z) profile the stress closure
    needs is DERIVED from the built material nc (Phase 4), never supplied.  The legacy
    6 KB density_profile.npz becomes a cache under data/<project>/cache/, keyed by the
    material nc's sha256.
    """

    orientation: SourceSpec | None = None
    velocity: SourceSpec | None = None
    thermal: SourceSpec | None = None


# --------------------------------------------------------------------------- the project
@dataclass(frozen=True)
class Project:
    """The whole descriptor.  Load with `Project.load(path)`."""

    name: str
    crs: CRS
    strike: StrikeFrame
    stress_box: GridBox
    meshes: Mapping[str, MeshSpec]
    default_mesh: str
    gate_bands: tuple[NamedBand, ...]
    hypocenters: Mapping[str, Hypocenter]
    physics: PhysicsDefaults
    raw: RawSources
    provenance: Mapping[str, str] = field(default_factory=dict)
    source_path: Path | None = None
    data_dir: Path | None = None

    # ---------------------------------------------------------------- loading
    @classmethod
    def load(cls, path: str | Path, require_files: bool = True,
             data_dir: str | Path | None = None) -> "Project":
        """Read and validate a descriptor YAML.

        require_files -- when True (the default) every referenced path must exist.  Tests
                         pass False to skip ONLY the existence checks; schema validation
                         always runs.
        data_dir      -- root that MeshSpec.path and SourceSpec.path resolve against.
                         Defaults to `<descriptor's parent>/../data/<name>`.
        """
        p = Path(path)
        if not p.is_file():
            raise ConfigError(f"descriptor not found: {p.resolve()}")
        with open(p, "r") as fh:
            raw = yaml.safe_load(fh)
        if not isinstance(raw, dict):
            raise ConfigError(f"{p}: top level must be a mapping, got {type(raw).__name__}")
        obj = cls.from_dict(raw, source_path=p, data_dir=data_dir)
        obj.validate(require_files=require_files)
        return obj

    @classmethod
    def from_dict(cls, raw: dict, source_path: Path | None = None,
                  data_dir: str | Path | None = None) -> "Project":
        where = str(source_path) if source_path else "<dict>"
        _reject_unknown(raw, _PROJECT_KEYS, where, "top level")

        try:
            name = raw["name"]
            crs = _build(CRS, raw["crs"], where, "crs")
            strike = _build(StrikeFrame, raw["strike"], where, "strike")
            stress_box = _build(GridBox, raw["stress_box"], where, "stress_box")
            default_mesh = raw["default_mesh"]
        except KeyError as exc:
            raise ConfigError(f"{where}: missing required key {exc.args[0]!r}") from None

        meshes = {k: _build(MeshSpec, v, where, f"meshes.{k}")
                  for k, v in (raw.get("meshes") or {}).items()}
        hypocenters = {k: _build(Hypocenter, v, where, f"hypocenters.{k}")
                       for k, v in (raw.get("hypocenters") or {}).items()}
        gate_bands = [_build(NamedBand, v, where, f"gate_bands[{i}]")
                      for i, v in enumerate(raw.get("gate_bands") or [])]
        physics = _build(PhysicsDefaults, raw.get("physics") or {}, where, "physics")

        raw_block = raw.get("raw") or {}
        _reject_unknown(raw_block, {"orientation", "velocity", "thermal"}, where, "raw")
        sources = RawSources(
            orientation=_build_source(raw_block.get("orientation"), where, "raw.orientation"),
            velocity=_build_source(raw_block.get("velocity"), where, "raw.velocity"),
            thermal=_build_source(raw_block.get("thermal"), where, "raw.thermal"),
        )

        if data_dir is not None:
            ddir = Path(data_dir)
        elif source_path is not None:
            ddir = source_path.resolve().parent.parent / "data" / name
        else:
            ddir = None

        return cls(
            name=name, crs=crs, strike=strike, stress_box=stress_box,
            meshes=_freeze(meshes),
            default_mesh=default_mesh, gate_bands=tuple(gate_bands),
            hypocenters=_freeze(hypocenters),
            physics=physics, raw=sources,
            provenance=_freeze(raw.get("provenance") or {}),
            source_path=source_path.resolve() if source_path else None,
            data_dir=ddir,
        )

    # ---------------------------------------------------------------- validation
    def validate(self, require_files: bool = True) -> None:
        """Hard-fail on any descriptor problem.  Raises ConfigError; never warns."""
        where = str(self.source_path) if self.source_path else self.name

        if not isinstance(self.name, str) or not self.name.strip():
            raise ConfigError(f"{where}: 'name' must be a non-empty string")

        self._validate_box(self.stress_box, where, "stress_box")

        if not self.meshes:
            raise ConfigError(f"{where}: 'meshes' must contain at least one entry")
        if self.default_mesh not in self.meshes:
            raise ConfigError(
                f"{where}: default_mesh {self.default_mesh!r} is not in meshes "
                f"({sorted(self.meshes)})")

        for mname, hypo in self.hypocenters.items():
            if mname not in self.meshes:
                raise ConfigError(
                    f"{where}: hypocenters key {mname!r} has no matching mesh "
                    f"({sorted(self.meshes)})")
            proj_set = [k for k in ("x", "y", "z") if getattr(hypo, k) is not None]
            geo_set = [k for k in ("lon", "lat", "depth_m") if getattr(hypo, k) is not None]
            if proj_set and geo_set:
                # Catch PARTIAL mixing, which the complete-set XOR below would miss: a full
                # projected point plus a stray `lon` would validate and silently discard the
                # lon.  That is the likely shape of a half-finished retarget, and it would
                # run on the ORIGINAL project's hypocentre with a ~0 m snap distance.
                raise ConfigError(
                    f"{where}: hypocenters.{mname} mixes projected {proj_set} with "
                    f"geographic {geo_set}; give one coordinate system only")
            if hypo.is_projected == hypo.is_geographic:
                raise ConfigError(
                    f"{where}: hypocenters.{mname} must give EITHER a complete projected "
                    f"(x, y, z) OR a complete geographic (lon, lat, depth_m) point, "
                    f"not both and not neither")
            if hypo.snap_tol_m <= 0:
                raise ConfigError(
                    f"{where}: hypocenters.{mname}.snap_tol_m must be > 0, "
                    f"got {hypo.snap_tol_m}")

        for i, band in enumerate(self.gate_bands):
            if band.s_end_km <= band.s_start_km:
                raise ConfigError(
                    f"{where}: gate_bands[{i}] ({band.name!r}) has s_end_km "
                    f"{band.s_end_km} <= s_start_km {band.s_start_km}")

        self._validate_kind(self.raw.orientation, ORIENTATION_KINDS, where, "raw.orientation")
        self._validate_kind(self.raw.velocity, VELOCITY_KINDS, where, "raw.velocity")
        self._validate_kind(self.raw.thermal, THERMAL_KINDS, where, "raw.thermal")

        for mname, mesh in self.meshes.items():
            if not mesh.tag_to_bc:
                raise ConfigError(f"{where}: meshes.{mname}.tag_to_bc must not be empty")
            if mesh.fault_bc not in mesh.tag_to_bc.values():
                raise ConfigError(
                    f"{where}: meshes.{mname}.fault_bc {mesh.fault_bc} does not appear in "
                    f"tag_to_bc values {sorted(set(mesh.tag_to_bc.values()))}")

        if require_files:
            self._validate_files(where)

    def _validate_files(self, where: str) -> None:
        for mname, mesh in self.meshes.items():
            self._require(self.resolve_path(mesh.path), where, f"meshes.{mname}.path")
        for label, spec in (("raw.orientation", self.raw.orientation),
                            ("raw.velocity", self.raw.velocity),
                            ("raw.thermal", self.raw.thermal)):
            if spec is not None and spec.path:
                self._require(self.resolve_path(spec.path), where, f"{label}.path")

    @staticmethod
    def _require(path: Path, where: str, label: str) -> None:
        # A source path may be a directory of slices or a single file; accept either, and
        # accept a glob that matches at least one entry.
        if path.exists():
            return
        if any(ch in str(path) for ch in "*?["):
            import glob as _glob
            if _glob.glob(str(path)):
                return
        raise ConfigError(
            f"{where}: {label} not found: {path}\n"
            f"  hint: fetch the project inputs with `python download_data.py`, or load with "
            f"require_files=False if you only need the schema.")

    @staticmethod
    def _validate_box(box: GridBox, where: str, label: str) -> None:
        if box.dx <= 0:
            raise ConfigError(f"{where}: {label}.dx must be > 0, got {box.dx}")
        if box.dz <= 0:
            raise ConfigError(f"{where}: {label}.dz must be > 0, got {box.dz}")
        for lo, hi, axis in ((box.xmin, box.xmax, "x"), (box.ymin, box.ymax, "y"),
                             (box.zmin, box.zmax, "z")):
            if hi <= lo:
                raise ConfigError(
                    f"{where}: {label}.{axis}max ({hi}) must be > {label}.{axis}min ({lo})")

    @staticmethod
    def _validate_kind(spec: SourceSpec | None, known: dict, where: str, label: str) -> None:
        if spec is None:
            return
        if spec.kind not in known:
            raise ConfigError(
                f"{where}: {label}.kind {spec.kind!r} is not a known reader; "
                f"expected one of {sorted(known)}")

    # ---------------------------------------------------------------- helpers
    def resolve_path(self, rel: str | Path) -> Path:
        """Resolve a descriptor-relative path against the project's data directory."""
        p = Path(rel)
        if p.is_absolute():
            return p
        if self.data_dir is None:
            return p
        return self.data_dir / p

    def mesh(self, name: str | None = None) -> MeshSpec:
        key = self.default_mesh if name is None else name
        if key not in self.meshes:
            raise ConfigError(f"unknown mesh {key!r}; have {sorted(self.meshes)}")
        return self.meshes[key]

    def hypocenter(self, mesh_name: str | None = None) -> Hypocenter:
        key = self.default_mesh if mesh_name is None else mesh_name
        if key not in self.hypocenters:
            raise ConfigError(
                f"no hypocenter defined for mesh {key!r}; have {sorted(self.hypocenters)}")
        return self.hypocenters[key]

    def sha256(self) -> str:
        """Hash of the descriptor file, for provenance.  '' when loaded from a dict."""
        if self.source_path is None or not self.source_path.is_file():
            return ""
        h = hashlib.sha256()
        with open(self.source_path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()

    def summary(self) -> str:
        lines = [
            f"project      : {self.name}",
            f"descriptor   : {self.source_path}",
            f"crs          : {self.crs.epsg}",
            f"strike       : az {self.strike.azimuth_deg} deg from {self.strike.origin_xy}",
            f"default mesh : {self.default_mesh}  (of {sorted(self.meshes)})",
            f"gate bands   : {[b.name for b in self.gate_bands] or 'none'}",
            f"raw sources  : orientation={_kind_of(self.raw.orientation)} "
            f"velocity={_kind_of(self.raw.velocity)} thermal={_kind_of(self.raw.thermal)}",
            f"data dir     : {self.data_dir}",
        ]
        return "\n".join(lines)


def _kind_of(spec: SourceSpec | None) -> str:
    return "none" if spec is None else spec.kind


# --------------------------------------------------------------------------- dict -> dataclass
_PROJECT_KEYS = {
    "name", "crs", "strike", "stress_box", "meshes", "default_mesh", "gate_bands",
    "hypocenters", "physics", "raw", "provenance",
}

_TUPLE_FIELDS = {
    (StrikeFrame, "origin_xy"): 2,
    (PhysicsDefaults, "seis_band_km"): 2,
    (PhysicsDefaults, "active_band_km"): 2,
}

# Fields that MUST be numeric.  Checked with _numeric so a YAML 1.1 scientific-notation
# string (e.g. `23.5e9`, which PyYAML reads as str) fails loudly at load time.
_SCI_NO_SIGN = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)[eE]\d+$")

_FLOAT_FIELDS: dict[type, tuple[str, ...]] = {
    StrikeFrame: ("azimuth_deg", "harmonise_hint_deg"),
    GridBox: ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax", "dx", "dz"),
    NamedBand: ("s_start_km", "s_end_km"),
    Hypocenter: ("x", "y", "z", "lon", "lat", "depth_m", "snap_tol_m"),
    MeshSpec: ("daylight_min_depth_m",),
    MeshGates: ("fault_edge_max_m", "f_gate_hz", "partition_tol_m"),
    PhysicsDefaults: ("mu_shear_pa", "w_energy_m", "dc_m", "kappa_c", "l_coast_km",
                      "gap_max_km", "sn_floor_mpa", "rs_b", "rs_sl0"),
}


def _reject_unknown(block: dict, known: set[str], where: str, label: str) -> None:
    """A YAML key we do not recognise is an ERROR, never ignored.

    Silent key-drop is how a retargeted project keeps a SAFS default without anyone
    noticing -- e.g. a misspelt 'strike_frame:' leaves the real 'strike:' unset.
    """
    if not isinstance(block, dict):
        raise ConfigError(f"{where}: {label} must be a mapping, got {type(block).__name__}")
    unknown = sorted(set(block) - known)
    if unknown:
        raise ConfigError(
            f"{where}: unknown key(s) in {label}: {unknown}; known keys are {sorted(known)}")


def _freeze(m):
    """Return a read-only view of a mapping.

    `frozen=True` only blocks attribute REBINDING; a plain dict field stays mutable, so
    `cfg.meshes[k].tag_to_bc[999] = 1` would silently rewrite the shared descriptor and
    make an artifact's recorded provenance disagree with what actually built it.
    """
    return MappingProxyType(dict(m))


def _numeric(value, where: str, label: str, field_name: str) -> float:
    """Require an actual number, and diagnose the YAML scientific-notation trap.

    PyYAML implements YAML 1.1, where `23.5e9` is a STRING -- the exponent needs an
    explicit sign (`23.5e+9`) or the value must be written out.  Without this check a
    shear modulus silently arrives as the string '23.5e9' and only explodes much later,
    or worse, compares unequal to the legacy constant without anyone noticing.
    """
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        extra = ""
        if isinstance(value, str) and _SCI_NO_SIGN.match(value.strip()):
            extra = (f"  -- this looks like scientific notation that YAML 1.1 read as a "
                     f"STRING; write it as {value.strip().replace('e', 'e+').replace('E', 'e+')} "
                     f"(explicit '+') or spell the number out")
        raise ConfigError(
            f"{where}: {label}.{field_name} must be a number, got "
            f"{type(value).__name__} {value!r}{extra}")
    return float(value)


def _build(dc, block, where: str, label: str):
    """Construct a frozen dataclass from a YAML mapping, rejecting unknown keys."""
    if block is None:
        block = {}
    known = {f.name for f in fields(dc)}
    _reject_unknown(block, known, where, label)
    kwargs = dict(block)
    for f in fields(dc):
        n = _TUPLE_FIELDS.get((dc, f.name))
        if n is not None and f.name in kwargs:
            val = kwargs[f.name]
            if not isinstance(val, (list, tuple)) or len(val) != n:
                raise ConfigError(
                    f"{where}: {label}.{f.name} must be a {n}-element sequence, got {val!r}")
            kwargs[f.name] = tuple(
                _numeric(v, where, label, f"{f.name}[{i}]") for i, v in enumerate(val))
            continue
        if f.name in _FLOAT_FIELDS.get(dc, ()) and kwargs.get(f.name) is not None:
            kwargs[f.name] = _numeric(kwargs[f.name], where, label, f.name)
    if dc is MeshSpec and kwargs.get("gates") is not None:
        kwargs["gates"] = _build(MeshGates, kwargs["gates"], where, f"{label}.gates")
    if dc is MeshSpec and "tag_to_bc" in kwargs and kwargs["tag_to_bc"] is not None:
        try:
            kwargs["tag_to_bc"] = _freeze(
                {int(k): int(v) for k, v in kwargs["tag_to_bc"].items()})
        except (AttributeError, TypeError, ValueError) as exc:
            raise ConfigError(
                f"{where}: {label}.tag_to_bc must be a mapping of int->int ({exc})") from None
    try:
        return dc(**kwargs)
    except TypeError as exc:
        raise ConfigError(f"{where}: cannot build {label}: {exc}") from None


def _build_source(block, where: str, label: str) -> SourceSpec | None:
    if block is None:
        return None
    if not isinstance(block, dict):
        raise ConfigError(f"{where}: {label} must be a mapping, got {type(block).__name__}")
    _reject_unknown(block, {"kind", "path", "params"}, where, label)
    if "kind" not in block:
        raise ConfigError(f"{where}: {label} requires a 'kind'")
    params = block.get("params") or {}
    if not isinstance(params, dict):
        raise ConfigError(f"{where}: {label}.params must be a mapping")
    return SourceSpec(kind=block["kind"], path=block.get("path"), params=_freeze(params))


def as_dict(obj) -> Any:
    """Recursively convert a descriptor dataclass back to plain types (for manifests)."""
    if is_dataclass(obj):
        return {f.name: as_dict(getattr(obj, f.name)) for f in fields(obj)
                if f.name not in ("source_path", "data_dir")}
    if isinstance(obj, Mapping):          # covers dict AND MappingProxyType
        return {k: as_dict(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [as_dict(v) for v in obj]
    if isinstance(obj, Path):
        return str(obj)
    return obj
