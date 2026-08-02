"""deckbuild -- build SeisSol dynamic-rupture deck inputs from raw data.

Four ingredients, one contract.  Three of them (material, stress, friction) are built here
from raw data; the fourth (the mesh) is ingested and gated here, and BUILT with the
`code-mesh-build-improve` skill in skills/.  See docs/PLAN_*.md.

    from deckbuild.bootstrap import init
    B = init(project="safs_alt")
"""
from __future__ import annotations

__version__ = "0.1.0.dev0"
__all__ = [
    "__version__",
    # config
    "Project", "ConfigError", "CRS", "StrikeFrame", "GridBox", "NamedBand",
    "Hypocenter", "MeshSpec", "PhysicsDefaults", "SourceSpec", "RawSources",
    # contract
    "Stage", "Gate", "GateReport", "GateFailure", "Artifact", "Manifest",
    "HARD", "WARN", "SKIP", "sha256_file",
    # geometry
    "Fault", "SnappedHypocenter", "GeometryError", "strike_s_km", "load_fault",
    "fault_trace", "resolve_tractions", "build_grid", "snap_hypocenter",
    # asagi
    "AsagiError", "write_asagi", "read_asagi", "trilinear_sample",
    "hull_containment", "roundtrip_selfcheck",
    # stress
    "StressStage", "KDesign", "StressError",
    "OrientationField", "read_orientation", "OrientationError",
    # friction / material
    "FrictionStage", "FwDesign", "FrictionError",
    "MaterialStage", "MaterialArtifacts", "PlasticitySpec", "AttenuationSpec",
    "MaterialError",
    # bootstrap
    "init", "Bootstrap",
]

from deckbuild.config import (  # noqa: E402
    ConfigError, CRS, GridBox, Hypocenter, MeshSpec, NamedBand, PhysicsDefaults,
    Project, RawSources, SourceSpec, StrikeFrame,
)
from deckbuild.contract import (  # noqa: E402
    Artifact, Gate, GateFailure, GateReport, HARD, Manifest, SKIP, Stage, WARN,
    sha256_file,
)
from deckbuild.geometry import (  # noqa: E402
    Fault, GeometryError, SnappedHypocenter, build_grid, fault_trace, load_fault,
    resolve_tractions, snap_hypocenter, strike_s_km,
)
from deckbuild.asagi import (  # noqa: E402
    AsagiError, hull_containment, read_asagi, roundtrip_selfcheck, trilinear_sample,
    write_asagi,
)
from deckbuild.orientation import (  # noqa: E402
    OrientationError, OrientationField, read_orientation,
)
from deckbuild.stress import KDesign, StressError, StressStage  # noqa: E402
from deckbuild.friction import (  # noqa: E402
    FrictionError, FrictionStage, FwDesign,
)
from deckbuild.material import (  # noqa: E402
    AttenuationSpec, MaterialArtifacts, MaterialError, MaterialStage, PlasticitySpec,
)
from deckbuild.bootstrap import Bootstrap, init  # noqa: E402
