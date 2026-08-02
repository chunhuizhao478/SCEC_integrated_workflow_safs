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
from deckbuild.bootstrap import Bootstrap, init  # noqa: E402
