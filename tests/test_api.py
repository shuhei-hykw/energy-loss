"""Tests for the top-level scope.md-style API and the model registry."""

from __future__ import annotations

import math

import pytest

from energy_loss import (
  compare_energy_from_range,
  compare_range_from_energy,
  energy_from_range,
  list_models,
  load_table,
  range_from_energy,
)
from energy_loss.models.registry import resolve_model

# -------- model registry / auto policy -------------------------------------

def test_registry_lists_nist_pstar_for_proton_emulsion():
  assert "nist_pstar" in list_models(
    particle="proton", material="nuclear_emulsion"
  )


def test_registry_lists_nist_astar_for_alpha_emulsion():
  assert "nist_astar" in list_models(
    particle="alpha", material="nuclear_emulsion"
  )


def test_auto_resolves_for_proton():
  assert resolve_model("proton", "nuclear_emulsion", "auto") == "nist_pstar"


def test_auto_resolves_for_alpha():
  assert resolve_model("alpha", "nuclear_emulsion", "auto") == "nist_astar"


def test_auto_fails_for_unknown_particle():
  with pytest.raises(ValueError, match="auto-policy"):
    resolve_model("muon-", "nuclear_emulsion", "auto")


def test_unknown_model_raises_with_available_list():
  with pytest.raises(ValueError, match="Available models"):
    resolve_model("proton", "nuclear_emulsion", "fictitious_model")


def test_material_synonyms_resolve():
  # "emulsion" -> "nuclear_emulsion" -> NIST PSTAR.
  assert resolve_model("proton", "emulsion", "auto") == "nist_pstar"
  assert resolve_model("proton", "photographic_emulsion", "auto") == "nist_pstar"


def test_particle_synonyms_resolve():
  assert resolve_model("p", "nuclear_emulsion", "auto") == "nist_pstar"
  assert resolve_model("4He", "nuclear_emulsion", "auto") == "nist_astar"


def test_load_table_returns_consistent_object():
  t1 = load_table("proton", "nuclear_emulsion")
  t2 = load_table("proton", "nuclear_emulsion")
  assert t1 is t2  # cached


# -------- energy_from_range / range_from_energy ----------------------------

def test_energy_from_range_um_shorthand():
  t = energy_from_range("alpha", "nuclear_emulsion", range_um=28.0)
  assert 4.0 < t < 8.0


def test_energy_from_range_via_range_value_unit():
  # 10 um == 0.001 cm; should agree.
  via_um = energy_from_range("proton", "nuclear_emulsion", range_um=10.0)
  via_cm = energy_from_range(
    "proton", "nuclear_emulsion",
    range_value=1.0e-3, range_unit="cm",
  )
  assert math.isclose(via_um, via_cm, rel_tol=1e-12)


def test_energy_from_range_unit_conversion_on_output():
  t_mev = energy_from_range("proton", "nuclear_emulsion", range_um=10.0)
  t_kev = energy_from_range(
    "proton", "nuclear_emulsion", range_um=10.0, energy_unit="keV"
  )
  assert math.isclose(t_mev * 1000.0, t_kev, rel_tol=1e-12)


def test_range_from_energy_round_trip():
  r_um = range_from_energy(
    "proton", "nuclear_emulsion", energy_mev=2.0
  )
  t_back = energy_from_range(
    "proton", "nuclear_emulsion", range_um=r_um
  )
  assert math.isclose(t_back, 2.0, rel_tol=1e-3)


def test_range_from_energy_unit_conversion():
  r_um = range_from_energy("proton", "nuclear_emulsion", energy_mev=10.0)
  r_mm = range_from_energy(
    "proton", "nuclear_emulsion", energy_mev=10.0, range_unit="mm"
  )
  assert math.isclose(r_um, r_mm * 1.0e3, rel_tol=1e-12)


def test_energy_from_range_requires_one_of_two_inputs():
  with pytest.raises(ValueError, match="exactly one"):
    energy_from_range("proton", "nuclear_emulsion")
  with pytest.raises(ValueError, match="exactly one"):
    energy_from_range(
      "proton", "nuclear_emulsion", range_um=1.0, range_value=2.0
    )


# -------- compare_* helpers ------------------------------------------------

def test_compare_energy_from_range_reports_unit_and_value():
  out = compare_energy_from_range(
    "proton", "nuclear_emulsion", range_um=10.0
  )
  assert out["unit"] == "MeV"
  assert "nist_pstar" in out
  assert 0.7 < out["nist_pstar"] < 0.9


def test_compare_range_from_energy_reports_unit_and_value():
  out = compare_range_from_energy(
    "alpha", "nuclear_emulsion", energy_mev=5.486
  )
  assert out["unit"] == "um"
  assert "nist_astar" in out
  # 5.486 MeV alpha range in NIST emulsion is around 28 um.
  assert 20.0 < out["nist_astar"] < 35.0


def test_compare_models_with_explicit_list_includes_auto():
  out = compare_energy_from_range(
    "proton", "nuclear_emulsion", range_um=10.0,
    models=["auto", "nist_pstar"],
  )
  # Both should resolve to the same concrete model, so we get one key.
  assert "nist_pstar" in out
  assert len([k for k in out if k != "unit"]) == 1


def test_compare_with_unknown_particle_raises():
  with pytest.raises(ValueError, match="No models registered"):
    compare_energy_from_range(
      "fictitious_particle", "nuclear_emulsion", range_um=10.0
    )
