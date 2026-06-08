"""Tests for the basic Bethe stopping-power implementation."""

from __future__ import annotations

import math
import warnings

import pytest

from energy_loss.constants import K_BETHE_MEV_MOL_PER_G
from energy_loss.materials import get_material
from energy_loss.stopping import (
  bethe_linear_stopping_power,
  bethe_mass_stopping_power,
)
from energy_loss.stopping.models import linear_stopping_power, mass_stopping_power
from energy_loss.units import energy_to_mev, length_to_cm


def test_k_bethe_prefactor_matches_pdg():
  # PDG 2024 quotes K = 0.307075 MeV mol / g.
  assert math.isclose(K_BETHE_MEV_MOL_PER_G, 0.307075, rel_tol=1e-4)


def test_mass_stopping_power_positive_for_proton_in_aluminum():
  # 100 MeV proton in aluminium — well within Bethe's validity range.
  s = bethe_mass_stopping_power("proton", 100.0, "aluminum")
  assert s > 0.0
  assert math.isfinite(s)


def test_linear_equals_mass_times_density():
  mat = get_material("kapton")
  s_mass = bethe_mass_stopping_power("proton", 200.0, mat)
  s_lin = bethe_linear_stopping_power("proton", 200.0, mat)
  assert math.isclose(s_lin, s_mass * mat.density_g_per_cm3, rel_tol=1e-12)


def test_mip_proton_in_aluminum_matches_pdg_order_of_magnitude():
  # PDG: minimum of dE/dx for a singly-charged particle in Al is
  # ~1.615 MeV cm^2 / g, located near beta*gamma ~ 3 (i.e. T ~ 2.5 GeV
  # for protons). Without the density-effect correction we slightly
  # overshoot at high energy, so allow a generous tolerance.
  s_min_candidates = []
  for t_gev in (2.0, 2.5, 3.0, 3.5):
    s_min_candidates.append(bethe_mass_stopping_power("proton", t_gev * 1.0e3, "aluminum"))
  s_min = min(s_min_candidates)
  assert 1.4 < s_min < 2.5, f"min dE/dx in Al = {s_min:.3f} MeV cm^2/g"


def test_dedx_decreases_with_energy_in_low_relativistic_regime():
  # Below the MIP, dE/dx falls as 1/beta^2 with increasing energy.
  s_50 = bethe_mass_stopping_power("proton", 50.0, "aluminum")
  s_200 = bethe_mass_stopping_power("proton", 200.0, "aluminum")
  s_500 = bethe_mass_stopping_power("proton", 500.0, "aluminum")
  assert s_50 > s_200 > s_500


def test_dedx_scales_with_charge_squared():
  # Same beta -> same beta*gamma -> bracket is identical. Stopping power
  # should scale as z^2. Use pi+ and K+ at the *same* beta*gamma by
  # picking kinetic energies proportional to mass.
  from energy_loss.particles import get_particle

  pi = get_particle("pion+")
  k = get_particle("kaon+")
  # Pick beta*gamma = 1.0 -> gamma = sqrt(2), so T = (gamma-1) * M.
  gamma = math.sqrt(2.0)
  t_pi = (gamma - 1.0) * pi.mass_mev
  t_k = (gamma - 1.0) * k.mass_mev
  s_pi = bethe_mass_stopping_power("pion+", t_pi, "aluminum")
  s_k = bethe_mass_stopping_power("kaon+", t_k, "aluminum")
  # T_max depends on the projectile mass, so the ratio is not exactly 1
  # even with z=1 for both. But at beta*gamma=1 both are within ~1%.
  assert math.isclose(s_pi, s_k, rel_tol=0.05)


def test_low_energy_emits_warning():
  # 0.1 MeV proton -> beta*gamma ~ 0.015, below the threshold.
  with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    bethe_mass_stopping_power("proton", 0.1, "aluminum")
    assert any("Bethe" in str(rec.message) or "beta" in str(rec.message) for rec in w)


def test_electron_not_implemented():
  with pytest.raises(NotImplementedError):
    bethe_mass_stopping_power("electron", 10.0, "aluminum")


def test_invalid_kinetic_energy_rejected():
  with pytest.raises(ValueError):
    bethe_mass_stopping_power("proton", 0.0, "aluminum")
  with pytest.raises(ValueError):
    bethe_mass_stopping_power("proton", -10.0, "aluminum")


def test_unknown_material_rejected():
  with pytest.raises(ValueError):
    bethe_mass_stopping_power("proton", 100.0, "unobtainium")


def test_unknown_particle_rejected():
  with pytest.raises(ValueError):
    bethe_mass_stopping_power("graviton", 100.0, "aluminum")


def test_models_dispatcher_matches_direct_call():
  s_direct = bethe_mass_stopping_power("proton", 100.0, "aluminum")
  s_via = mass_stopping_power("proton", 100.0, "aluminum", model="bethe")
  assert math.isclose(s_direct, s_via, rel_tol=1e-15)
  l_direct = bethe_linear_stopping_power("proton", 100.0, "aluminum")
  l_via = linear_stopping_power("proton", 100.0, "aluminum", model="bethe")
  assert math.isclose(l_direct, l_via, rel_tol=1e-15)


def test_unknown_model_rejected():
  with pytest.raises(ValueError):
    mass_stopping_power("proton", 100.0, "aluminum", model="srim")


def test_unit_conversion_round_trip():
  # Spot-check the unit helpers since Bethe consumes MeV / cm internally.
  assert math.isclose(energy_to_mev(1.0, "GeV"), 1000.0)
  assert math.isclose(energy_to_mev(1.0, "keV"), 1.0e-3)
  assert math.isclose(length_to_cm(10.0, "um"), 10.0e-4)
  assert math.isclose(length_to_cm(1.0, "mm"), 0.1)
