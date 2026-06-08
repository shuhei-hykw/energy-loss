"""Unit conversion helpers.

The package's *internal* units are MeV for energy and cm for length.
These helpers convert user-supplied values into internal units, and back.
Keep the conversion table minimal and explicit; do not silently accept
unknown unit names.
"""

from __future__ import annotations

_ENERGY_TO_MEV: dict[str, float] = {
  "eV": 1.0e-6,
  "keV": 1.0e-3,
  "MeV": 1.0,
  "GeV": 1.0e3,
  "TeV": 1.0e6,
}

_LENGTH_TO_CM: dict[str, float] = {
  "nm": 1.0e-7,
  "um": 1.0e-4,
  "micron": 1.0e-4,
  "mm": 1.0e-1,
  "cm": 1.0,
  "m": 1.0e2,
}

# Momentum is internally stored in MeV/c. Numerically p[MeV/c] = p*c[MeV],
# so the factors here are the same as the energy table.
_MOMENTUM_TO_MEV_C: dict[str, float] = {
  "eV/c": 1.0e-6,
  "keV/c": 1.0e-3,
  "MeV/c": 1.0,
  "GeV/c": 1.0e3,
  "TeV/c": 1.0e6,
}

_MASS_THICKNESS_TO_G_PER_CM2: dict[str, float] = {
  "mg/cm^2": 1.0e-3,
  "g/cm^2": 1.0,
  "kg/m^2": 0.1,        # 1 kg/m^2 = 1e3 g / 1e4 cm^2 = 0.1 g/cm^2
}


def energy_to_mev(value: float, unit: str) -> float:
  """Convert ``value`` in ``unit`` to MeV."""
  try:
    factor = _ENERGY_TO_MEV[unit]
  except KeyError as exc:
    raise ValueError(f"Unknown energy unit: {unit!r}") from exc
  return value * factor


def mev_to_energy(value_mev: float, unit: str) -> float:
  """Convert a value in MeV to ``unit``."""
  try:
    factor = _ENERGY_TO_MEV[unit]
  except KeyError as exc:
    raise ValueError(f"Unknown energy unit: {unit!r}") from exc
  return value_mev / factor


def length_to_cm(value: float, unit: str) -> float:
  """Convert ``value`` in ``unit`` to cm."""
  try:
    factor = _LENGTH_TO_CM[unit]
  except KeyError as exc:
    raise ValueError(f"Unknown length unit: {unit!r}") from exc
  return value * factor


def cm_to_length(value_cm: float, unit: str) -> float:
  """Convert a value in cm to ``unit``."""
  try:
    factor = _LENGTH_TO_CM[unit]
  except KeyError as exc:
    raise ValueError(f"Unknown length unit: {unit!r}") from exc
  return value_cm / factor


def momentum_to_mev_c(value: float, unit: str) -> float:
  """Convert ``value`` in momentum ``unit`` to MeV/c."""
  try:
    factor = _MOMENTUM_TO_MEV_C[unit]
  except KeyError as exc:
    raise ValueError(f"Unknown momentum unit: {unit!r}") from exc
  return value * factor


def mass_thickness_to_g_per_cm2(value: float, unit: str) -> float:
  """Convert ``value`` in mass-thickness ``unit`` to g/cm^2."""
  try:
    factor = _MASS_THICKNESS_TO_G_PER_CM2[unit]
  except KeyError as exc:
    raise ValueError(f"Unknown mass-thickness unit: {unit!r}") from exc
  return value * factor


def kinetic_energy_from_momentum(p_mev_c: float, mass_mev: float) -> float:
  """Convert a relativistic momentum [MeV/c] to kinetic energy [MeV].

  Uses ``T = sqrt(p^2 + m^2) - m`` with c = 1.
  """
  import math as _math

  if p_mev_c < 0.0:
    raise ValueError(f"momentum must be non-negative, got {p_mev_c}")
  if mass_mev <= 0.0:
    raise ValueError(f"mass must be positive, got {mass_mev}")
  return _math.sqrt(p_mev_c * p_mev_c + mass_mev * mass_mev) - mass_mev
