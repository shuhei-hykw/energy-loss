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
