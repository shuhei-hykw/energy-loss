"""Material definitions.

A :class:`Material` is represented by an *effective* (Z, A, I, density),
which is the form the Bethe formula needs. For compounds and mixtures
this is the effective <Z/A> times A together with a Bragg-additivity-
based mean excitation energy I.

The values here are reasonable defaults to get v0.1 working; they are
NOT a substitute for experiment-specific calibration. In particular
``nuclear_emulsion`` density and composition vary between emulsion
types and processing — see README.

Units
-----
- density : g / cm^3
- mean excitation energy I : eV

References:
- PDG "Atomic and nuclear properties of materials" 2024 table
- NIST PSTAR/ESTAR documentation for I values
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Material:
  """A material treated as a uniform stopping medium.

  Attributes
  ----------
  name : str
    Canonical name.
  z_over_a : float
    Effective <Z/A>. Dimensionless.
  mean_excitation_energy_ev : float
    Mean excitation energy I [eV].
  density_g_per_cm3 : float
    Density [g/cm^3].
  """

  name: str
  z_over_a: float
  mean_excitation_energy_ev: float
  density_g_per_cm3: float


# PDG-style values. I values for compounds are approximate effective
# values; for high-precision work use NIST or experiment calibration.
_MATERIALS: dict[str, Material] = {
  # Hydrogen gas, NTP.
  "H2": Material("H2", z_over_a=0.99212, mean_excitation_energy_ev=19.2,
                 density_g_per_cm3=8.376e-5),
  # Liquid hydrogen.
  "LH2": Material("LH2", z_over_a=0.99212, mean_excitation_energy_ev=21.8,
                  density_g_per_cm3=0.0708),
  # Helium gas, NTP.
  "He": Material("He", z_over_a=0.49967, mean_excitation_energy_ev=41.8,
                 density_g_per_cm3=1.663e-4),
  # Dry air, NTP (PDG).
  "air": Material("air", z_over_a=0.49919, mean_excitation_energy_ev=85.7,
                  density_g_per_cm3=1.205e-3),
  # P10: 90% Ar + 10% CH4 by volume, NTP, density ~ 1.561e-3 g/cm^3.
  "P10": Material("P10", z_over_a=0.4505, mean_excitation_energy_ev=174.0,
                  density_g_per_cm3=1.561e-3),
  # Generic polystyrene-based plastic scintillator (e.g. NE-102 / BC-408).
  "plastic_scintillator": Material(
    "plastic_scintillator", z_over_a=0.5414,
    mean_excitation_energy_ev=64.7, density_g_per_cm3=1.032,
  ),
  # Kapton polyimide.
  "kapton": Material("kapton", z_over_a=0.5126,
                     mean_excitation_energy_ev=79.6,
                     density_g_per_cm3=1.42),
  # Mylar (PET).
  "mylar": Material("mylar", z_over_a=0.5197,
                    mean_excitation_energy_ev=78.7,
                    density_g_per_cm3=1.40),
  # Aluminium.
  "aluminum": Material("aluminum", z_over_a=13.0 / 26.9815385,
                       mean_excitation_energy_ev=166.0,
                       density_g_per_cm3=2.699),
  # Amorphous carbon (graphite ~2.21).
  "carbon": Material("carbon", z_over_a=6.0 / 12.0107,
                     mean_excitation_energy_ev=78.0,
                     density_g_per_cm3=2.21),
  # Nuclear emulsion (rough effective values; CALIBRATE per emulsion).
  "nuclear_emulsion": Material(
    "nuclear_emulsion", z_over_a=0.4255,
    mean_excitation_energy_ev=331.0, density_g_per_cm3=3.815,
  ),
}

_ALIASES: dict[str, str] = {
  "Al": "aluminum",
  "C": "carbon",
  "emulsion": "nuclear_emulsion",
}


def get_material(name: str | Material) -> Material:
  """Look up a :class:`Material` by name (or pass through if already one)."""
  if isinstance(name, Material):
    return name
  key = _ALIASES.get(name, name)
  try:
    return _MATERIALS[key]
  except KeyError as exc:
    known = ", ".join(sorted(_MATERIALS))
    raise ValueError(
      f"Unknown material {name!r}. Known materials: {known}"
    ) from exc
