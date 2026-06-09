"""Tests for RangeEnergyTable and the bundled emulsion tables."""

from __future__ import annotations

import math
import textwrap
from pathlib import Path

import numpy as np
import pytest

from energy_loss import (
  RangeEnergyTable,
  RangeEnergyTableMetadata,
  energy_from_emulsion_range,
  get_emulsion_range_energy,
  list_bundled_emulsion_tables,
)

# -------- generic table behaviour ------------------------------------------

def test_from_arrays_round_trip():
  # Pick a simple monotone curve: R = T^2
  t = np.array([1.0, 2.0, 5.0, 10.0, 50.0, 100.0])
  r = t**2
  tab = RangeEnergyTable.from_arrays(t, r)
  for ti in (2.0, 5.0, 10.0, 50.0):
    assert math.isclose(tab.range_from_energy(ti), ti**2, rel_tol=1e-6)
    assert math.isclose(tab.energy_from_range(ti**2), ti, rel_tol=1e-6)


def test_log_interp_is_monotone():
  t = np.geomspace(1.0, 1.0e4, 50)
  r = t**1.7
  tab = RangeEnergyTable.from_arrays(t, r)
  prev = -math.inf
  for ti in np.geomspace(2.0, 5.0e3, 25):
    ri = tab.range_from_energy(float(ti))
    assert ri > prev
    prev = ri


def test_out_of_range_raises():
  t = np.array([1.0, 2.0, 5.0])
  r = np.array([1.0, 4.0, 25.0])
  tab = RangeEnergyTable.from_arrays(t, r)
  with pytest.raises(ValueError, match="outside the tabulated"):
    tab.range_from_energy(0.5)
  with pytest.raises(ValueError, match="outside the tabulated"):
    tab.range_from_energy(10.0)
  with pytest.raises(ValueError, match="outside the tabulated"):
    tab.energy_from_range(0.5)
  with pytest.raises(ValueError, match="outside the tabulated"):
    tab.energy_from_range(100.0)


def test_non_monotone_rejected():
  t = np.array([1.0, 5.0, 2.0])  # not increasing
  r = np.array([1.0, 4.0, 25.0])
  with pytest.raises(ValueError, match="strictly increasing"):
    RangeEnergyTable.from_arrays(t, r)


def test_non_positive_rejected():
  t = np.array([1.0, 2.0, 5.0])
  r = np.array([0.0, 4.0, 25.0])
  with pytest.raises(ValueError, match="strictly positive"):
    RangeEnergyTable.from_arrays(t, r)


def test_linear_range_helper_uses_metadata_density(tmp_path: Path):
  meta = RangeEnergyTableMetadata(density_g_per_cm3=2.0)
  tab = RangeEnergyTable.from_arrays(
    np.array([1.0, 10.0]), np.array([1.0, 10.0]), metadata=meta,
  )
  # range = 2 g/cm^2 -> energy 2 MeV. linear 1 cm * 2 g/cm^3 = 2 g/cm^2.
  assert math.isclose(tab.energy_from_linear_range(1.0), 2.0, rel_tol=1e-12)


def test_linear_range_helper_needs_density():
  tab = RangeEnergyTable.from_arrays(
    np.array([1.0, 10.0]), np.array([1.0, 10.0])
  )
  with pytest.raises(ValueError, match="needs a density"):
    tab.energy_from_linear_range(1.0)


# -------- bundled NIST tables ----------------------------------------------

def test_list_bundled_tables_includes_proton_and_alpha():
  bundled = list_bundled_emulsion_tables()
  assert "proton" in bundled
  assert "alpha" in bundled


def test_bundled_proton_table_metadata():
  tab = get_emulsion_range_energy("proton")
  assert "Photographic Emulsion" in tab.metadata.material_name
  assert tab.metadata.density_g_per_cm3 == pytest.approx(3.815, rel=1e-4)
  assert tab.metadata.mean_excitation_energy_ev == pytest.approx(331.0)
  # NIST PSTAR composition: weight fraction of Ag should be dominant.
  assert tab.metadata.composition["Ag"] > 0.4
  assert "PSTAR" in tab.metadata.source
  # Version trace: fetched date is recorded.
  assert len(tab.metadata.fetched) == len("YYYY-MM-DD")


def test_bundled_alpha_table_metadata():
  tab = get_emulsion_range_energy("alpha")
  assert "ASTAR" in tab.metadata.source
  assert tab.metadata.particle == "alpha"


def test_proton_28um_in_emulsion_gives_few_mev_range():
  # Standard sanity check: 28 um proton range in NIST emulsion.
  t = energy_from_emulsion_range("proton", 28.0, "um")
  assert 1.0 < t < 3.0


def test_alpha_28um_in_emulsion_matches_po214_alpha():
  # 5.3-7.7 MeV alpha range in nuclear emulsion is ~25-40 um. A 28 um
  # track should give ~5-7 MeV.
  t = energy_from_emulsion_range("alpha", 28.0, "um")
  assert 4.0 < t < 8.0


def test_emulsion_helper_accepts_unit_conversion():
  # 10 um == 0.001 cm; both should give the same energy.
  t_um = energy_from_emulsion_range("proton", 10.0, "um")
  t_cm = energy_from_emulsion_range("proton", 1.0e-3, "cm")
  assert math.isclose(t_um, t_cm, rel_tol=1e-12)


def test_unknown_particle_emulsion_table():
  with pytest.raises(ValueError, match="No bundled emulsion"):
    get_emulsion_range_energy("muon")


# -------- bundled alpha + new built-in particles ---------------------------

def test_alpha_built_in_particle():
  from energy_loss import get_particle

  alpha = get_particle("alpha")
  assert alpha.charge == 2
  assert 3727.0 < alpha.mass_mev < 3728.0
  # alias
  assert get_particle("4He").name == "alpha"


def test_deuteron_and_triton_built_in():
  from energy_loss import get_particle

  d = get_particle("d")
  t = get_particle("t")
  assert d.name == "deuteron" and d.charge == 1
  assert t.name == "triton" and t.charge == 1


# -------- NIST CSV custom load ---------------------------------------------

def test_nist_csv_loads_with_metadata(tmp_path: Path):
  # Round-trip a fake mini NIST-style CSV through from_nist_csv.
  csv = textwrap.dedent(
    """\
    # NIST PSTAR table for FakeMaterial (matno=999)
    # Particle: proton
    # Source:   https://example.invalid/PSTAR.html
    # Density [g/cm^3]: 1.50000
    # Mean excitation excitation energy [eV]: 100.000
    # Composition (weight fraction): H=0.500000, C=0.500000
    # Fetched: 2026-06-09 via scripts/fetch_nist_emulsion.py
    # Columns: T_MeV, S_elec, S_nuc, S_total, R_csda, R_proj, detour
    T_MeV,S_elec,S_nuc,S_total,R_csda,R_proj,detour
    1.0,10.0,1.0,11.0,0.1,0.09,0.9
    10.0,5.0,0.1,5.1,1.0,0.95,0.95
    100.0,2.0,0.05,2.05,10.0,9.8,0.98
    """
  )
  p = tmp_path / "fake.csv"
  p.write_text(csv)
  tab = RangeEnergyTable.from_nist_csv(p)
  assert tab.metadata.particle == "proton"
  assert tab.metadata.density_g_per_cm3 == 1.5
  assert tab.metadata.composition == {"H": 0.5, "C": 0.5}
  assert math.isclose(tab.range_from_energy(10.0), 1.0, rel_tol=1e-12)
  assert math.isclose(tab.energy_from_range(1.0), 10.0, rel_tol=1e-12)
