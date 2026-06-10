"""Tests for v0.6 model-validity tracking."""

from __future__ import annotations

import math
import warnings

import pytest

from energy_loss import (
  compare_energy_from_range,
  energy_from_range,
  list_models,
  range_from_energy,
)
from energy_loss.models.registry import (
  get_model_entry,
  resolve_model,
)

# -------- list_models with energy filter -----------------------------------

def test_list_models_filters_by_valid_at():
  # ATIMA recommended range is [2, 1000] MeV.
  models_lo = list_models(
    particle="proton", material="nuclear_emulsion", valid_at_mev=0.5,
  )
  assert "geant4_atima_11_4_1" not in models_lo
  models_mid = list_models(
    particle="proton", material="nuclear_emulsion", valid_at_mev=100.0,
  )
  assert "geant4_atima_11_4_1" in models_mid
  models_hi = list_models(
    particle="proton", material="nuclear_emulsion", valid_at_mev=5000.0,
  )
  assert "geant4_atima_11_4_1" not in models_hi


# -------- ModelEntry exposes the recommended band --------------------------

def test_atima_entry_advertises_2_to_1000_mev():
  entry = get_model_entry("proton", "nuclear_emulsion", "geant4_atima")
  assert entry.valid_min_mev == pytest.approx(2.0)
  assert entry.valid_max_mev == pytest.approx(1000.0)
  assert entry.covers(100.0)
  assert not entry.covers(0.5)


def test_pstar_entry_advertises_full_table():
  entry = get_model_entry("proton", "nuclear_emulsion", "nist_pstar")
  assert entry.valid_min_mev <= 0.001
  assert entry.valid_max_mev >= 1.0e4


# -------- energy-aware auto resolution -------------------------------------

def test_auto_resolves_to_pstar_at_low_energy_for_proton():
  # 0.5 MeV is below ATIMA min — PSTAR is the only NIST-tagged
  # candidate; auto should pick it.
  assert (
    resolve_model("proton", "nuclear_emulsion", "auto", energy_mev=0.5)
    == "nist_pstar"
  )


def test_auto_resolves_to_pstar_at_mid_energy_for_proton():
  # At 50 MeV all models cover; auto returns the first preference
  # which is PSTAR.
  assert (
    resolve_model("proton", "nuclear_emulsion", "auto", energy_mev=50.0)
    == "nist_pstar"
  )


def test_auto_resolves_to_pstar_above_atima_emax():
  # 5 GeV is above ATIMA emax; PSTAR still covers.
  assert (
    resolve_model("proton", "nuclear_emulsion", "auto", energy_mev=5000.0)
    == "nist_pstar"
  )


# -------- out-of-range warnings --------------------------------------------

def test_atima_below_band_warns():
  # 5 um proton track -> ~0.5 MeV via ATIMA, below the 2 MeV floor.
  with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    energy_from_range(
      "proton", "nuclear_emulsion", range_um=5.0, model="geant4_atima"
    )
  messages = [str(rec.message) for rec in w]
  assert any("outside its recommended" in m for m in messages)


def test_in_range_query_does_not_warn():
  # 50 MeV proton via PSTAR (well within its range) — no warning.
  with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    range_from_energy(
      "proton", "nuclear_emulsion", energy_mev=50.0, model="nist_pstar"
    )
  messages = [str(rec.message) for rec in w]
  assert not any("outside its recommended" in m for m in messages)


def test_warning_message_carries_recommended_range_and_note():
  with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    energy_from_range(
      "alpha", "nuclear_emulsion", range_um=1.0, model="geant4_atima"
    )
  matching = [
    str(rec.message) for rec in w if "geant4_atima" in str(rec.message)
  ]
  assert matching
  msg = matching[0]
  assert "[2" in msg and "1000" in msg
  assert "ATIMA" in msg or "Atima" in msg


# -------- compare skip_out_of_range ----------------------------------------

def test_compare_skip_out_of_range_drops_atima_at_low_E():
  out = compare_energy_from_range(
    "proton", "nuclear_emulsion", range_um=5.0,
    skip_out_of_range=True,
  )
  assert "geant4_atima_11_4_1" not in out
  assert "nist_pstar" in out


def test_compare_without_skip_includes_atima_even_when_out_of_band():
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    out = compare_energy_from_range(
      "proton", "nuclear_emulsion", range_um=5.0,
      skip_out_of_range=False,
    )
  assert "geant4_atima_11_4_1" in out


def test_compare_inverse_skip_out_of_range_uses_input_energy():
  out = compare_energy_from_range(
    "proton", "nuclear_emulsion", range_um=5.0,
    skip_out_of_range=True,
  )
  # PSTAR/option4 cover 5 MeV/30 keV range; ATIMA does not.
  assert "geant4_atima_11_4_1" not in out


# -------- legacy back-compat -----------------------------------------------

def test_explicit_model_still_resolvable_without_energy():
  # resolve_model with model="nist_pstar" should work regardless of
  # whether energy_mev is supplied.
  assert (
    resolve_model("proton", "nuclear_emulsion", "nist_pstar")
    == "nist_pstar"
  )
  assert (
    resolve_model(
      "proton", "nuclear_emulsion", "nist_pstar", energy_mev=math.inf,
    )
    == "nist_pstar"
  )
