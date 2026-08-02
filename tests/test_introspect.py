"""Phase 8 acceptance tests: deck introspection and the reproduction ladder."""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from deckbuild.introspect import UNKNOWN, introspect_deck  # noqa: E402

DECKS = Path(os.environ.get("DECKBUILD_DECKS",
                            str(Path.home() / "Downloads" / "seisol_quakeworx")))
TARGET = next(DECKS.glob("*ALT_THERMAL_CASE1_intermediate*deep40km"), None) \
    if DECKS.is_dir() else None
needs_deck = pytest.mark.skipif(TARGET is None, reason="the shipped SAFS deck is absent")


def test_missing_deck_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="deck not found"):
        introspect_deck(tmp_path / "nope")


def test_introspect_a_synthetic_deck(tmp_path):
    """Anything not inferable is UNKNOWN, never guessed."""
    d = tmp_path / "deck"
    d.mkdir()
    f = introspect_deck(d)
    for k in ("stress_nc", "material_nc", "mesh"):
        assert f[k] == UNKNOWN
    assert "stress_nc" in f.unknown


@needs_deck
def test_recovers_the_stress_design():
    f = introspect_deck(TARGET)
    assert f["k_from_filename"] == 1.7
    # The Phase 2 naming rule: no _freeze tag means the freeze is OFF.
    assert f["freeze_from_filename"] == 0.0
    b = f["stress_box"]
    assert (b["dx"], b["dz"]) == (1000.0, 250.0)


@needs_deck
def test_recovers_the_strike_frame_from_the_lua():
    """The frame is only present as two coefficients inside the emitted Lua."""
    f = introspect_deck(TARGET)
    assert f["strike_azimuth_deg"] == pytest.approx(314.0, abs=1e-6)
    assert f["strike_origin_xy"] == (606971.0, 3707270.0)


@needs_deck
def test_recovers_the_seven_region_fw_design():
    f = introspect_deck(TARGET)
    assert f["fw_n_regions"] == 7
    assert f["fw_values"] == [0.0, 0.0, 0.045, 0.03, 0.06, 0.0175, 0.05]
    assert f["fw_boundaries_s_km"][0] == (20.0, 28.0)
    assert f["fw_boundaries_s_km"][-1] == (244.0, 252.0)
    # the design must agree with the rs_muw yaml's own FILENAME
    y = next(TARGET.glob("*rs_muw*.yaml"), None)
    if y:
        for v in ("0.045", "0.03", "0.06", "0.0175", "0.05"):
            assert v in y.name


@needs_deck
def test_recovers_the_hypocentre_and_nucleation():
    f = introspect_deck(TARGET)
    assert f["hypocenter_xyz"] == pytest.approx(
        (604446.944, 3704576.3853, -10067.9819), abs=1e-6)
    assert f["nucleation_radius_m"] == 2000.0
    assert f["nucleation_amplitude"] == pytest.approx(75.0e6)


@needs_deck
def test_reads_phi_from_the_field_not_the_prose():
    """A deck's yaml comment has been wrong before; bulkFriction is ground truth."""
    f = introspect_deck(TARGET)
    assert f["phi_deg_from_field"] == pytest.approx([30.0, 40.0], abs=1e-3)


@needs_deck
def test_recovers_the_scalars_and_par_settings():
    f = introspect_deck(TARGET)
    assert f["rs_b_constant"] == 0.019
    assert f["rs_sl0_constant"] == 0.10
    assert f["Plasticity"] == 1 and f["FL"] == 103
    assert f["FreqCentral"] == 0.5


@needs_deck
def test_nothing_is_left_unknown_on_a_real_deck():
    """A field the introspector cannot infer must be MARKED, and there should be none
    here -- if that changes, the exercise's E2 rung degrades and must say so."""
    f = introspect_deck(TARGET)
    assert f.unknown == [], f"could not infer: {f.unknown}"


@needs_deck
def test_the_four_grids_are_independent():
    """Measured on the real deck: material 1500/250, stress 1000/250, friction 1500/200."""
    f = introspect_deck(TARGET)
    assert "dx=1500 dz=250" in f["material_grid"]
    assert "dx=1000 dz=250" in f["stress_grid"]
    assert "dx=1500 dz=200" in f["friction_grid"]
