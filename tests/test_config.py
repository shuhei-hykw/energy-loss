"""Tests for the YAML setup loader (energy_loss.config)."""

from __future__ import annotations

import math
import textwrap
from pathlib import Path

import pytest

from energy_loss import (
  compute_linear_stopping_power,
  compute_mass_stopping_power,
  get_material,
  get_particle,
  load_config,
)
from energy_loss.config import Beam, Config, Target

REPO = Path(__file__).resolve().parent.parent
CFG_DIR = REPO / "examples" / "configs"


def test_jparc_e10_pim_config_loads():
  cfg = load_config(CFG_DIR / "jparc_e10_pim.yaml")
  assert isinstance(cfg, Config)
  assert cfg.beam.particle.name == "pion-"
  assert math.isclose(cfg.beam.momentum_mev_c, 1200.0, rel_tol=1e-12)
  # T_kin = sqrt(p^2 + m^2) - m for m_pi = 139.57 MeV
  m_pi = cfg.beam.particle.mass_mev
  expected_t = math.sqrt(1200.0**2 + m_pi**2) - m_pi
  assert math.isclose(cfg.beam.kinetic_energy_mev, expected_t, rel_tol=1e-12)
  assert cfg.target.material.name == "Be"
  assert math.isclose(cfg.target.mass_thickness_g_per_cm2, 3.5, rel_tol=1e-12)
  # thickness must be grammage / density
  expected_x = 3.5 / cfg.target.material.density_g_per_cm3
  assert math.isclose(cfg.target.thickness_cm, expected_x, rel_tol=1e-12)


def test_jparc_e10_km_config_loads():
  cfg = load_config(CFG_DIR / "jparc_e10_km.yaml")
  assert cfg.beam.particle.name == "kaon-"
  assert math.isclose(cfg.beam.momentum_mev_c, 1500.0, rel_tol=1e-12)


def test_jparc_pi_and_k_stopping_powers_are_close():
  # Both beams at ~1 GeV kinetic energy on the same target -> dE/dx values
  # are within ~10% of each other and well above MIP (Be).
  cfg_pi = load_config(CFG_DIR / "jparc_e10_pim.yaml")
  cfg_k = load_config(CFG_DIR / "jparc_e10_km.yaml")
  s_pi = compute_mass_stopping_power(cfg_pi)
  s_k = compute_mass_stopping_power(cfg_k)
  assert 1.2 < s_pi < 2.5
  assert 1.2 < s_k < 2.5
  assert abs(s_pi - s_k) / s_pi < 0.2


def test_mean_energy_loss_jparc_pim_is_reasonable():
  cfg = load_config(CFG_DIR / "jparc_e10_pim.yaml")
  s_mass = compute_mass_stopping_power(cfg)
  de = s_mass * cfg.target.mass_thickness_g_per_cm2
  # ~6 MeV for pi- 1.2 GeV/c through 3.5 g/cm^2 Be.
  assert 4.0 < de < 9.0, f"unexpected dE = {de} MeV"


def test_alpha_inline_particle_registered_and_used():
  cfg = load_config(CFG_DIR / "alpha_in_emulsion.yaml")
  assert cfg.beam.particle.name == "alpha"
  assert cfg.beam.particle.charge == 2
  # The custom material should also be findable post-load.
  m = get_material("my_emulsion_2024")
  assert m.density_g_per_cm3 == pytest.approx(3.80)
  # alpha was registered globally
  assert get_particle("alpha").mass_mev == pytest.approx(3727.379)


def test_beam_kinetic_energy_path(tmp_path: Path):
  yaml_text = textwrap.dedent(
    """
    beam:
      particle: proton
      kinetic_energy: 0.2
      energy_unit: GeV
    target:
      material: aluminum
      thickness: 1.0
      thickness_unit: mm
    """
  )
  p = tmp_path / "ke.yaml"
  p.write_text(yaml_text)
  cfg = load_config(p)
  assert math.isclose(cfg.beam.kinetic_energy_mev, 200.0, rel_tol=1e-12)
  # momentum derived: sqrt((T+m)^2 - m^2)
  m = cfg.beam.particle.mass_mev
  expected_p = math.sqrt((200.0 + m) ** 2 - m * m)
  assert math.isclose(cfg.beam.momentum_mev_c, expected_p, rel_tol=1e-12)


def test_beam_must_specify_exactly_one_kinematic_input(tmp_path: Path):
  yaml_text = textwrap.dedent(
    """
    beam:
      particle: proton
      kinetic_energy: 100.0
      momentum: 1.0
      momentum_unit: GeV/c
    target:
      material: aluminum
      thickness: 1.0
      thickness_unit: mm
    """
  )
  p = tmp_path / "bad.yaml"
  p.write_text(yaml_text)
  with pytest.raises(ValueError, match="exactly one"):
    load_config(p)


def test_target_thickness_and_mass_thickness_are_mutually_exclusive(tmp_path: Path):
  yaml_text = textwrap.dedent(
    """
    beam:
      particle: proton
      kinetic_energy: 100.0
    target:
      material: aluminum
      thickness: 1.0
      thickness_unit: mm
      mass_thickness: 1.0
    """
  )
  p = tmp_path / "bad.yaml"
  p.write_text(yaml_text)
  with pytest.raises(ValueError, match="exactly one"):
    load_config(p)


def test_missing_beam_section(tmp_path: Path):
  yaml_text = textwrap.dedent(
    """
    target:
      material: aluminum
      thickness: 1.0
      thickness_unit: mm
    """
  )
  p = tmp_path / "bad.yaml"
  p.write_text(yaml_text)
  with pytest.raises(ValueError, match="'beam' section"):
    load_config(p)


def test_missing_target_section(tmp_path: Path):
  yaml_text = textwrap.dedent(
    """
    beam:
      particle: proton
      kinetic_energy: 100.0
    """
  )
  p = tmp_path / "bad.yaml"
  p.write_text(yaml_text)
  with pytest.raises(ValueError, match="'target' section"):
    load_config(p)


def test_compute_linear_matches_mass_times_density():
  cfg = load_config(CFG_DIR / "jparc_e10_pim.yaml")
  s_mass = compute_mass_stopping_power(cfg)
  s_lin = compute_linear_stopping_power(cfg)
  assert math.isclose(
    s_lin, s_mass * cfg.target.material.density_g_per_cm3, rel_tol=1e-12
  )


def test_beam_and_target_are_dataclasses():
  cfg = load_config(CFG_DIR / "jparc_e10_pim.yaml")
  assert isinstance(cfg.beam, Beam)
  assert isinstance(cfg.target, Target)
