"""Thin shim around the model registry for emulsion convenience access.

The v0.3 API (``energy_from_emulsion_range`` / ``get_emulsion_range_energy``)
is preserved for back-compatibility, but everything now routes through
:mod:`energy_loss.models.registry`. New code should prefer the general
:func:`energy_loss.api.energy_from_range`.
"""

from __future__ import annotations

from energy_loss.models.registry import (
  AUTO_POLICY,
  list_models,
  load_table,
)
from energy_loss.range.table import RangeEnergyTable
from energy_loss.units import length_to_cm


def list_bundled_emulsion_tables() -> list[str]:
  """Particles for which a bundled NIST table exists for nuclear_emulsion."""
  out: set[str] = set()
  for p, preferences in AUTO_POLICY.items():
    available = list_models(particle=p, material="nuclear_emulsion")
    # The first NIST-tagged preference that is actually registered.
    nist_match = next(
      (m for m in preferences if m.startswith("nist_") and m in available),
      None,
    )
    if nist_match is not None:
      out.add(p)
  return sorted(out)


def get_emulsion_range_energy(particle: str) -> RangeEnergyTable:
  """Return the default nuclear-emulsion table for ``particle``."""
  try:
    return load_table(particle=particle, material="nuclear_emulsion", model="auto")
  except ValueError as exc:
    raise ValueError(
      f"No bundled emulsion range-energy table for particle {particle!r}. "
      f"Available: {list_bundled_emulsion_tables()}."
    ) from exc


def energy_from_emulsion_range(
  particle: str, range_value: float, range_unit: str = "um",
  table: RangeEnergyTable | None = None,
) -> float:
  """Kinetic energy [MeV] for an emulsion track of ``range_value`` ``range_unit``.

  By default uses the bundled NIST table for ``particle``; pass an
  explicit ``table`` to use a different calibration.
  """
  tab = table if table is not None else get_emulsion_range_energy(particle)
  range_cm = length_to_cm(float(range_value), range_unit)
  return tab.energy_from_linear_range(range_cm)
