"""Acceptance criterion 5: no SAFS constant may be a literal inside deckbuild/.

The whole premise of the workflow is that retargeting to a new fault system means editing
one YAML.  A constant left in the package silently keeps the San Andreas value for every
project.  This test greps the package for the four load-bearing literals.

It scans SOURCE ONLY: the descriptors in projects/ are exactly where these numbers belong.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
PKG = ROOT / "deckbuild"

# (regex, human name).  Word boundaries so 3146 or 06971 do not false-positive.
FORBIDDEN = [
    (re.compile(r"\b32611\b"), "EPSG:32611 (the SAFS UTM zone)"),
    (re.compile(r"\b314\.0\b"), "314.0 (the SAFS strike azimuth)"),
    (re.compile(r"\b606971\b"), "606971 (the SAFS s-origin easting)"),
    (re.compile(r"\b3707270\b"), "3707270 (the SAFS s-origin northing)"),
]


def _sources() -> list[Path]:
    return sorted(p for p in PKG.rglob("*.py") if "__pycache__" not in p.parts)


def test_package_has_sources():
    assert _sources(), f"no python sources found under {PKG}"


@pytest.mark.parametrize("path", _sources(), ids=lambda p: p.name)
def test_no_safs_literal_in_source(path: Path):
    text = path.read_text()
    hits = []
    for lineno, line in enumerate(text.splitlines(), 1):
        for pattern, label in FORBIDDEN:
            if pattern.search(line):
                hits.append(f"  {path.relative_to(ROOT)}:{lineno}: {label}\n      {line.strip()}")
    assert not hits, (
        "SAFS constants must live in projects/*.yaml, not in the package:\n" + "\n".join(hits))


def test_the_guard_itself_would_catch_a_leak(tmp_path):
    """Guard against the guard silently matching nothing (e.g. a broken regex)."""
    sample = "epsg = 'EPSG:32611'  # leaked\naz = 314.0\nox = 606971\noy = 3707270\n"
    found = {label for pattern, label in FORBIDDEN if pattern.search(sample)}
    assert len(found) == len(FORBIDDEN), f"guard missed: {[l for _, l in FORBIDDEN]}"


def test_safs_constants_ARE_present_in_the_descriptors():
    """The mirror image: the numbers must exist somewhere, or they were lost in transit."""
    text = (ROOT / "projects" / "safs_alt.yaml").read_text()
    for pattern, label in FORBIDDEN:
        assert pattern.search(text), f"{label} missing from projects/safs_alt.yaml"
