"""contract.py -- the three-verb stage contract and its supporting types.

Every stage of the workflow exposes the same three verbs:

    build(cfg, out_dir, **params) -> Artifact     writes a file
    verify(cfg, artifact, **params) -> GateReport re-opens THAT FILE and gates it
    plot(cfg, artifact, **params)                 draws it

`verify` re-opening the written file rather than checking an in-memory array is the whole
point: several legacy defects were verify-on-memory / ship-a-different-array.  The gate
names (V1-V4, F0-F6, G1-G4, M1-M5, A-G, P1-P8) are preserved verbatim from the legacy
verifiers so existing run notes stay readable.
"""
from __future__ import annotations

import hashlib
import json
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

__all__ = [
    "Gate",
    "GateReport",
    "GateFailure",
    "Artifact",
    "Stage",
    "Manifest",
    "sha256_file",
    "HARD",
    "WARN",
    "SKIP",
]

HARD = "hard"
WARN = "warn"
SKIP = "skip"
_SEVERITIES = (HARD, WARN, SKIP)

# Printed markers.  Asserted against a golden string in tests/test_contract.py -- do not
# change these without updating that test, because run notes quote them.
_MARK = {
    (HARD, True): "[PASS]",
    (HARD, False): "[FAIL]",
    (WARN, True): "[PASS]",
    (WARN, False): "[WARN]",
    (SKIP, True): "[ -- ]",
    (SKIP, False): "[ -- ]",
}


class GateFailure(RuntimeError):
    """Raised by `GateReport.raise_if_failed` when a hard gate did not pass."""


@dataclass
class Gate:
    """One named pass/fail check.

    severity: "hard" -- a failure blocks;
              "warn" -- reported, never blocks (e.g. the kappa corridor, which is
                        necessary-not-sufficient: a design can pass it and still arrest);
              "skip" -- not evaluated (a missing optional input); NEVER report a skipped
                        gate as a pass.
    """

    name: str
    passed: bool
    detail: str
    severity: str = HARD

    def __post_init__(self) -> None:
        if self.severity not in _SEVERITIES:
            raise ValueError(
                f"gate {self.name!r}: severity must be one of {_SEVERITIES}, "
                f"got {self.severity!r}")
        if not self.name:
            raise ValueError("a gate must have a name")

    @property
    def marker(self) -> str:
        return _MARK[(self.severity, bool(self.passed))]

    def line(self) -> str:
        return f"{self.marker} {self.name}: {self.detail}"


class GateReport:
    """An ordered collection of gates, with a single pass/fail verdict."""

    def __init__(self, title: str = "") -> None:
        self.title = title
        self.gates: list[Gate] = []

    # ------------------------------------------------------------------ building
    def add(self, name: str, passed: bool, detail: str, severity: str = HARD) -> Gate:
        gate = Gate(name=name, passed=bool(passed), detail=detail, severity=severity)
        self.gates.append(gate)
        return gate

    def skip(self, name: str, detail: str) -> Gate:
        """Record a gate that was not evaluated.  A skip is not a pass."""
        return self.add(name, False, detail, severity=SKIP)

    def extend(self, other: "GateReport", prefix: str = "") -> "GateReport":
        """Absorb another report's gates, optionally prefixing their names."""
        for g in other.gates:
            self.gates.append(Gate(name=f"{prefix}{g.name}", passed=g.passed,
                                   detail=g.detail, severity=g.severity))
        return self

    # ------------------------------------------------------------------ verdict
    @property
    def ok(self) -> bool:
        """True when every HARD gate passed.  Warns and skips never block."""
        return all(g.passed for g in self.gates if g.severity == HARD)

    @property
    def failed(self) -> list[Gate]:
        return [g for g in self.gates if g.severity == HARD and not g.passed]

    @property
    def warnings(self) -> list[Gate]:
        return [g for g in self.gates if g.severity == WARN and not g.passed]

    @property
    def skipped(self) -> list[Gate]:
        return [g for g in self.gates if g.severity == SKIP]

    def raise_if_failed(self) -> None:
        if self.ok:
            return
        names = ", ".join(g.name for g in self.failed)
        raise GateFailure(
            f"{self.title or 'gate battery'}: {len(self.failed)} hard gate(s) failed "
            f"({names})\n" + "\n".join(g.line() for g in self.failed))

    # ------------------------------------------------------------------ output
    def lines(self) -> list[str]:
        out = []
        if self.title:
            out.append(f"=== {self.title} ===")
        out.extend(g.line() for g in self.gates)
        n_hard = sum(1 for g in self.gates if g.severity == HARD)
        n_pass = sum(1 for g in self.gates if g.severity == HARD and g.passed)
        out.append(
            f"{'PASS' if self.ok else 'FAIL'}: {n_pass}/{n_hard} hard, "
            f"{len(self.warnings)} warn, {len(self.skipped)} skipped")
        return out

    def print(self) -> None:
        print("\n".join(self.lines()))

    def to_dict(self) -> dict:
        return {
            "title": self.title,
            "ok": self.ok,
            "gates": [
                {"name": g.name, "passed": g.passed, "detail": g.detail,
                 "severity": g.severity}
                for g in self.gates
            ],
        }

    def __len__(self) -> int:
        return len(self.gates)

    def __repr__(self) -> str:
        return (f"GateReport({self.title!r}, {len(self.gates)} gates, "
                f"ok={self.ok})")


# --------------------------------------------------------------------------- artifacts
def sha256_file(path: str | Path) -> str:
    """Streaming sha256 of a file.  Artifacts can be hundreds of MB."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


@dataclass
class Artifact:
    """One file a stage produced, with everything needed to reproduce it."""

    path: Path
    kind: str
    sha256: str = ""
    params: dict = field(default_factory=dict)
    provenance: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.path = Path(self.path)

    @classmethod
    def of(cls, path: str | Path, kind: str, params: dict | None = None,
           provenance: dict | None = None) -> "Artifact":
        """Build an Artifact for an existing file, hashing it now."""
        p = Path(path)
        if not p.is_file():
            raise FileNotFoundError(
                f"cannot register artifact {kind!r}: no such file: {p.resolve()}")
        return cls(path=p, kind=kind, sha256=sha256_file(p),
                   params=dict(params or {}), provenance=dict(provenance or {}))

    def to_dict(self) -> dict:
        return {"path": str(self.path), "kind": self.kind, "sha256": self.sha256,
                "params": self.params, "provenance": self.provenance}


# --------------------------------------------------------------------------- the stage
class Stage(ABC):
    """Base class for every workflow stage.

    Subclasses set `name` and implement build/verify.  `plot` and `check_outputs` have
    usable defaults so a stage is never forced to grow a figure it does not have.
    """

    name: str = "stage"

    @abstractmethod
    def build(self, cfg, out_dir: Path, **params) -> Artifact:
        """Write the stage's artifact into out_dir and return it."""

    @abstractmethod
    def verify(self, cfg, artifact: Artifact, **params) -> GateReport:
        """Re-open the artifact FROM DISK and run this stage's gate battery."""

    def plot(self, cfg, artifact: Artifact, **params):
        raise NotImplementedError(f"{self.name}: no plot() implemented")

    def check_outputs(self, out_dir: Path, keep: Iterable[str | Path]) -> dict:
        """Report which files in out_dir are deliverables and which are scratch.

        Some stages build intermediate files (the graded-k V1 references, for example).
        Exactly the files in `keep` are deliverables; everything else is listed as extra so
        a stale artifact cannot be shipped unnoticed.
        """
        out_dir = Path(out_dir)
        keep_resolved = {Path(k).resolve() for k in keep}
        present = {p.resolve() for p in out_dir.glob("*") if p.is_file()}
        missing = sorted(str(p) for p in keep_resolved - present)
        extra = sorted(str(p) for p in present - keep_resolved)
        return {"stage": self.name, "dir": str(out_dir),
                "deliverables": sorted(str(p) for p in keep_resolved),
                "missing": missing, "extra": extra,
                "ok": not missing}


# --------------------------------------------------------------------------- manifest
class Manifest:
    """Accumulates artifacts and gate reports, and writes deck_manifest.json.

    Records the descriptor's name and sha256 alongside every artifact, so a deck can be
    traced back to the exact configuration that produced it -- and so two artifacts built
    from different descriptors can be detected before they are assembled together.
    """

    FILENAME = "deck_manifest.json"

    def __init__(self, project_name: str = "", descriptor_sha256: str = "") -> None:
        self.project_name = project_name
        self.descriptor_sha256 = descriptor_sha256
        self.stages: dict[str, dict] = {}

    @classmethod
    def for_project(cls, cfg) -> "Manifest":
        return cls(project_name=cfg.name, descriptor_sha256=cfg.sha256())

    def record(self, stage: str, artifact: Artifact | None = None,
               report: GateReport | None = None, **extra) -> None:
        entry = self.stages.setdefault(stage, {})
        if artifact is not None:
            entry.setdefault("artifacts", []).append(artifact.to_dict())
        if report is not None:
            entry["gates"] = report.to_dict()
        if extra:
            entry.update(extra)

    def to_dict(self) -> dict:
        return {
            "project": self.project_name,
            "descriptor_sha256": self.descriptor_sha256,
            "stages": self.stages,
        }

    def write(self, out_dir: str | Path) -> Path:
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        path = out / self.FILENAME
        with open(path, "w") as fh:
            json.dump(self.to_dict(), fh, indent=2, sort_keys=True)
            fh.write("\n")
        return path

    @classmethod
    def read(cls, path: str | Path) -> dict:
        with open(path, "r") as fh:
            return json.load(fh)

    def artifacts(self) -> list[dict]:
        out: list[dict] = []
        for entry in self.stages.values():
            out.extend(entry.get("artifacts", []))
        return out

    def __repr__(self) -> str:
        return (f"Manifest({self.project_name!r}, {len(self.stages)} stages, "
                f"{len(self.artifacts())} artifacts)")
