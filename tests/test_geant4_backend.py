"""Tests for the Geant4 backend.

The full backend test (subprocess + CSV round trip) only runs when the
C++ generator is available; otherwise the test is skipped so the suite
stays green on environments where Geant4 has not been built. Pure
Python unit tests (path resolution, error paths) run unconditionally.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from energy_loss import (
  compare_energy_from_range,
  energy_from_range,
  list_models,
  load_table,
)
from energy_loss.backends.geant4_runner import (
  Geant4TableSpec,
  Geant4Unavailable,
  detect_geant4_executable,
)


def _geant4_available() -> bool:
  try:
    detect_geant4_executable()
  except Geant4Unavailable:
    return False
  return True


requires_geant4 = pytest.mark.skipif(
  not _geant4_available(),
  reason="Geant4 generator executable not built",
)


# -------- pure unit tests --------------------------------------------------

def test_spec_cache_filename_is_deterministic():
  s1 = Geant4TableSpec(
    particle="proton", material="nuclear_emulsion",
    emin_mev=0.01, emax_mev=1000.0, n_points=100,
  )
  s2 = Geant4TableSpec(
    particle="proton", material="nuclear_emulsion",
    emin_mev=0.01, emax_mev=1000.0, n_points=100,
  )
  assert s1.cache_filename() == s2.cache_filename()


def test_spec_cache_filename_changes_with_physics_list():
  s_o4 = Geant4TableSpec(
    particle="proton", material="nuclear_emulsion",
    emin_mev=0.01, emax_mev=100.0, n_points=20,
    physics_list="option4",
  )
  s_at = Geant4TableSpec(
    particle="proton", material="nuclear_emulsion",
    emin_mev=0.01, emax_mev=100.0, n_points=20,
    physics_list="atima",
  )
  assert s_o4.cache_filename() != s_at.cache_filename()
  # Filename carries the physics-list label for human readability.
  assert "atima" in s_at.cache_filename()
  assert "option4" in s_o4.cache_filename()


def test_spec_changes_with_n_points():
  s1 = Geant4TableSpec(
    particle="proton", material="nuclear_emulsion",
    emin_mev=0.01, emax_mev=1000.0, n_points=100,
  )
  s2 = Geant4TableSpec(
    particle="proton", material="nuclear_emulsion",
    emin_mev=0.01, emax_mev=1000.0, n_points=200,
  )
  assert s1.cache_filename() != s2.cache_filename()


def test_detect_falls_back_to_explicit_path(tmp_path: Path):
  exe = tmp_path / "fake_g4"
  exe.write_text("")
  exe.chmod(0o755)
  resolved = detect_geant4_executable(exe)
  assert resolved == exe


def test_detect_raises_on_missing(tmp_path: Path, monkeypatch):
  # Force every resolution slot to point at non-existent paths.
  monkeypatch.delenv("ENERGY_LOSS_G4_EXECUTABLE", raising=False)
  bogus = tmp_path / "does_not_exist"
  # Patch the module-level default and PATH search.
  from energy_loss.backends import geant4_runner

  monkeypatch.setattr(geant4_runner, "_DEFAULT_BUILT", bogus)
  monkeypatch.setattr(geant4_runner.shutil, "which", lambda _: None)
  with pytest.raises(Geant4Unavailable):
    detect_geant4_executable()


def test_env_variable_overrides_default(tmp_path: Path, monkeypatch):
  exe = tmp_path / "via_env"
  exe.write_text("")
  exe.chmod(0o755)
  monkeypatch.setenv("ENERGY_LOSS_G4_EXECUTABLE", str(exe))
  assert detect_geant4_executable() == exe


# -------- registry: model name and listing ---------------------------------

def test_registry_advertises_geant4_for_proton_emulsion():
  models = list_models(particle="proton", material="nuclear_emulsion")
  assert "geant4_11_4_1" in models


def test_geant4_keyword_resolves_to_version_specific_name():
  # The "geant4" alias should resolve regardless of whether the
  # executable is present; resolve_model only checks the registration.
  from energy_loss.models.registry import resolve_model

  assert resolve_model("proton", "nuclear_emulsion", "geant4") == "geant4_11_4_1"


def test_geant4_atima_alias_resolves():
  from energy_loss.models.registry import resolve_model

  assert (
    resolve_model("proton", "nuclear_emulsion", "geant4_atima")
    == "geant4_atima_11_4_1"
  )
  assert (
    resolve_model("alpha", "nuclear_emulsion", "geant4_atima")
    == "geant4_atima_11_4_1"
  )


def test_registry_lists_both_geant4_models_for_emulsion():
  models = list_models(particle="proton", material="nuclear_emulsion")
  assert "geant4_11_4_1" in models
  assert "geant4_atima_11_4_1" in models


# -------- integration: runs the actual generator ---------------------------

@requires_geant4
def test_generate_proton_table_round_trip(tmp_path: Path):
  from energy_loss.backends import generate_geant4_table
  from energy_loss.range import RangeEnergyTable

  spec = Geant4TableSpec(
    particle="proton", material="nuclear_emulsion",
    emin_mev=0.1, emax_mev=100.0, n_points=20,
  )
  csv = generate_geant4_table(spec, output_path=tmp_path / "g4.csv")
  table = RangeEnergyTable.from_nist_csv(csv)
  assert table.metadata.density_g_per_cm3 == pytest.approx(3.815, rel=1e-4)
  assert "11.4.1" in table.metadata.source
  # Range monotonic / positive over the requested range.
  assert table.range_csda_g_per_cm2[0] < table.range_csda_g_per_cm2[-1]


@requires_geant4
def test_load_table_resolves_geant4_model():
  t = load_table("proton", "nuclear_emulsion", "geant4")
  assert "11.4.1" in t.metadata.source


@requires_geant4
def test_geant4_agrees_with_pstar_within_a_few_percent():
  out = compare_energy_from_range(
    "proton", "nuclear_emulsion", range_um=35.2,
    models=["nist_pstar", "geant4"],
  )
  rel = abs(out["nist_pstar"] - out["geant4_11_4_1"]) / out["nist_pstar"]
  assert rel < 0.05, f"PSTAR vs Geant4 diverged by {rel:.3%}"


@requires_geant4
def test_geant4_alpha_agrees_with_astar_within_a_few_percent():
  t_g4 = energy_from_range(
    "alpha", "nuclear_emulsion", range_um=28.0, model="geant4"
  )
  t_nist = energy_from_range(
    "alpha", "nuclear_emulsion", range_um=28.0, model="nist_astar"
  )
  assert abs(t_g4 - t_nist) / t_nist < 0.05


# -------- cache behaviour --------------------------------------------------

@requires_geant4
def test_atima_table_has_atima_metadata():
  t = load_table("proton", "nuclear_emulsion", "geant4_atima")
  assert "atima" in t.metadata.source.lower()


@requires_geant4
def test_atima_disagrees_with_pstar_at_low_energy_emulsion():
  # The whole point of registering ATIMA is to see a difference vs
  # PSTAR / option4 in the emulsion regime. At 10 um proton range
  # PSTAR / option4 agree to <1% and ATIMA differs by >10%.
  out = compare_energy_from_range(
    "proton", "nuclear_emulsion", range_um=10.0,
    models=["nist_pstar", "geant4", "geant4_atima"],
  )
  ratio_o4 = abs(out["nist_pstar"] - out["geant4_11_4_1"]) / out["nist_pstar"]
  ratio_at = abs(out["nist_pstar"] - out["geant4_atima_11_4_1"]) / out["nist_pstar"]
  assert ratio_o4 < 0.02
  assert ratio_at > 0.10


@requires_geant4
def test_cached_csv_is_reused(tmp_path: Path):
  from energy_loss.backends import generate_geant4_table

  spec = Geant4TableSpec(
    particle="proton", material="nuclear_emulsion",
    emin_mev=0.5, emax_mev=50.0, n_points=10,
  )
  out = tmp_path / "cache.csv"
  generate_geant4_table(spec, output_path=out)
  mtime1 = os.stat(out).st_mtime_ns
  generate_geant4_table(spec, output_path=out)
  mtime2 = os.stat(out).st_mtime_ns
  assert mtime1 == mtime2  # second call short-circuits on the cached file.
