"""Tests for the YAML-driven material registry."""

from __future__ import annotations

import math
import textwrap
from pathlib import Path

import periodictable as pt
import pytest

from energy_loss import get_material, list_materials, load_materials_from_yaml
from energy_loss.materials import Material
from energy_loss.stopping import bethe_mass_stopping_power


def test_known_materials_include_target_elements():
  names = set(list_materials())
  for name in ("Be", "carbon", "aluminum", "H2", "He", "air", "P10",
               "kapton", "mylar", "plastic_scintillator",
               "nuclear_emulsion"):
    assert name in names, f"{name!r} missing from registry: {names}"


def test_alias_lookup_for_carbon():
  via_alias = get_material("C")
  direct = get_material("carbon")
  assert via_alias is direct


def test_be_uses_periodictable_9be_isotope():
  # Be is configured as 9Be(Enriched): Z/A from the isotope mass.
  be = get_material("Be")
  iso = pt.elements.symbol("Be")[9]
  assert math.isclose(be.z_over_a, iso.element.number / iso.mass, rel_tol=1e-12)
  # Density is the natural Be density (Be is mononuclidic in nature).
  assert math.isclose(be.density_g_per_cm3, 1.848, rel_tol=1e-3)


def test_9be_alias_resolves_to_same_material():
  assert get_material("9Be") is get_material("Be")
  assert get_material("Be9") is get_material("Be")


def test_carbon_uses_periodictable_density():
  c = get_material("carbon")
  el = pt.elements.symbol("C")
  assert math.isclose(c.density_g_per_cm3, el.density, rel_tol=1e-12)
  assert math.isclose(c.z_over_a, el.number / el.mass, rel_tol=1e-12)


def test_h2_gas_density_overrides_periodictable_liquid():
  # periodictable returns ~0.0708 g/cm^3 for H (= LH2). YAML overrides
  # for H2 gas with the NTP density.
  h2 = get_material("H2")
  assert h2.density_g_per_cm3 < 1.0e-3  # gas, not liquid
  lh2 = get_material("LH2")
  assert lh2.density_g_per_cm3 > 0.05  # liquid


def test_bethe_proton_in_be_positive_and_finite():
  s = bethe_mass_stopping_power("proton", 100.0, "Be")
  assert s > 0.0 and math.isfinite(s)


def test_bethe_proton_in_carbon_falls_with_energy():
  s_50 = bethe_mass_stopping_power("proton", 50.0, "C")
  s_200 = bethe_mass_stopping_power("proton", 200.0, "C")
  s_500 = bethe_mass_stopping_power("proton", 500.0, "C")
  assert s_50 > s_200 > s_500


def test_bethe_relative_order_be_c_al():
  # At 100 MeV, mass stopping power tends to fall slowly with Z because
  # Z/A drops (Be ~0.444, C ~0.499, Al ~0.482). Order is not strictly
  # monotonic but Be should be below C.
  s_be = bethe_mass_stopping_power("proton", 100.0, "Be")
  s_c = bethe_mass_stopping_power("proton", 100.0, "C")
  assert s_be < s_c


def test_unknown_element_in_yaml_rejected(tmp_path: Path):
  yaml_text = textwrap.dedent(
    """
    materials:
      unobtainium:
        element: Xx
        mean_excitation_energy_ev: 100.0
    """
  )
  p = tmp_path / "bad.yaml"
  p.write_text(yaml_text)
  with pytest.raises(ValueError, match="Unknown element"):
    load_materials_from_yaml(p)


def test_missing_required_field_rejected(tmp_path: Path):
  yaml_text = textwrap.dedent(
    """
    materials:
      no_i:
        z_over_a: 0.5
        density_g_per_cm3: 1.0
    """
  )
  p = tmp_path / "bad.yaml"
  p.write_text(yaml_text)
  with pytest.raises(ValueError, match="mean_excitation_energy_ev"):
    load_materials_from_yaml(p)


def test_user_yaml_can_add_target(tmp_path: Path):
  yaml_text = textwrap.dedent(
    """
    materials:
      custom_target:
        element: Cu
        mean_excitation_energy_ev: 322.0
        reference: "test"
    """
  )
  p = tmp_path / "extra.yaml"
  p.write_text(yaml_text)
  load_materials_from_yaml(p)
  cu = get_material("custom_target")
  assert isinstance(cu, Material)
  s = bethe_mass_stopping_power("proton", 100.0, "custom_target")
  assert s > 0.0
