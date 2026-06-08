"""Tests for the step-wise RK4 transport integrator."""

from __future__ import annotations

import math
import warnings
from pathlib import Path

import pytest

from energy_loss import Layer, load_config, propagate, propagate_config
from energy_loss.stopping import bethe_mass_stopping_power

REPO = Path(__file__).resolve().parent.parent
CFG_DIR = REPO / "examples" / "configs"


# -------- basic invariants --------------------------------------------------

def test_zero_thickness_layer_is_a_no_op():
  layers = [Layer.from_mass_thickness("Be", 0.0)]
  r = propagate("proton", 100.0, layers)
  assert math.isclose(r.total_energy_loss_mev, 0.0)
  assert math.isclose(r.exit_kinetic_energy_mev, 100.0)


def test_empty_layer_list_raises():
  with pytest.raises(ValueError, match="at least one"):
    propagate("proton", 100.0, [])


def test_propagation_loses_energy_in_thick_target():
  layers = [Layer.from_mass_thickness("Be", 3.5)]
  r = propagate("pion-", 1068.5, layers)
  assert 0.0 < r.total_energy_loss_mev < r.initial_kinetic_energy_mev
  assert r.exit_kinetic_energy_mev > 0.0
  assert not r.stopped
  # close to single-point linear estimate at this energy
  s = bethe_mass_stopping_power("pion-", 1068.5, "Be")
  assert math.isclose(r.total_energy_loss_mev, s * 3.5, rel_tol=0.01)


def test_multi_layer_dE_sums_to_total():
  layers = [
    Layer.from_thickness("kapton", 5.0e-3),
    Layer.from_thickness("air", 5.0),
    Layer.from_mass_thickness("Be", 3.5),
  ]
  r = propagate("pion-", 1068.5, layers)
  per_layer_sum = sum(pl.energy_loss_mev for pl in r.per_layer)
  assert math.isclose(r.total_energy_loss_mev, per_layer_sum, rel_tol=1e-12)


def test_layer_chains_energy_through_correctly():
  layers = [
    Layer.from_mass_thickness("Be", 1.0),
    Layer.from_mass_thickness("Be", 1.0),
    Layer.from_mass_thickness("Be", 1.5),
  ]
  r = propagate("pion-", 1068.5, layers)
  # exit of layer i must equal entry of layer i+1
  for i in range(len(r.per_layer) - 1):
    a = r.per_layer[i].exit_kinetic_energy_mev
    b = r.per_layer[i + 1].entry_kinetic_energy_mev
    assert math.isclose(a, b, rel_tol=1e-12)


# -------- grammage invariance / step convergence ---------------------------

def test_thickness_form_invariance():
  """A single layer specified as mass thickness or via density and linear
  thickness should produce identical results."""
  by_mt = [Layer.from_mass_thickness("Be", 3.5)]
  # equivalent linear thickness:
  rho = by_mt[0].material.density_g_per_cm3
  by_lin = [Layer.from_thickness("Be", 3.5 / rho)]
  r_mt = propagate("pion-", 1068.5, by_mt)
  r_lin = propagate("pion-", 1068.5, by_lin)
  assert math.isclose(
    r_mt.exit_kinetic_energy_mev,
    r_lin.exit_kinetic_energy_mev,
    rel_tol=1e-12,
  )


def test_step_size_convergence_under_refinement():
  """Halving the step size should change the integrated dE by less than
  a part-per-thousand at high energy."""
  layers = [Layer.from_mass_thickness("Be", 3.5)]
  r_coarse = propagate("pion-", 1068.5, layers, step_g_per_cm2=3.5 / 100)
  r_fine = propagate("pion-", 1068.5, layers, step_g_per_cm2=3.5 / 2000)
  rel = abs(r_coarse.total_energy_loss_mev - r_fine.total_energy_loss_mev) / r_fine.total_energy_loss_mev
  assert rel < 1.0e-3, f"step refinement diverges: rel diff = {rel}"


# -------- stopping ----------------------------------------------------------

def test_low_energy_particle_stops_in_target():
  # 5 MeV proton has range ~0.04 g/cm^2 in Be (rough). A 1 g/cm^2 layer
  # is well beyond that, so the integrator should flag `stopped`.
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    r = propagate("proton", 5.0, [Layer.from_mass_thickness("Be", 1.0)])
  assert r.stopped
  assert r.exit_kinetic_energy_mev == 0.0
  assert r.per_layer[0].stopped_in_layer


def test_stopped_layers_after_stop_are_no_op():
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    layers = [
      Layer.from_mass_thickness("Be", 1.0),   # stops here
      Layer.from_mass_thickness("aluminum", 1.0),
    ]
    r = propagate("proton", 5.0, layers)
  assert r.stopped
  assert r.per_layer[1].energy_loss_mev == 0.0


# -------- config integration -----------------------------------------------

def test_propagate_config_on_jparc_e10_pim():
  cfg = load_config(CFG_DIR / "jparc_e10_pim.yaml")
  r = propagate_config(cfg)
  assert not r.stopped
  # ~6 MeV expected
  assert 4.0 < r.total_energy_loss_mev < 9.0


def test_propagate_config_on_jparc_stack_sums_correctly():
  cfg = load_config(CFG_DIR / "jparc_e10_stack.yaml")
  r = propagate_config(cfg)
  assert len(r.per_layer) == 3
  # the Be layer should dominate
  be_loss = r.per_layer[-1].energy_loss_mev
  assert be_loss / r.total_energy_loss_mev > 0.95


def test_propagate_alpha_in_emulsion_runs_and_loses_energy():
  # 5 MeV alpha sits below the basic-Bethe validity edge so the formula
  # clamps to zero before T reaches the stopping threshold. We only
  # assert that the integrator runs and reduces T (the physics is
  # documented as not trustworthy in this regime).
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    cfg = load_config(CFG_DIR / "alpha_in_emulsion.yaml")
    r = propagate_config(cfg)
  assert r.total_energy_loss_mev > 0.0
  assert r.exit_kinetic_energy_mev < cfg.beam.kinetic_energy_mev
