"""Default range-energy tables for nuclear emulsion.

The bundled tables come from NIST PSTAR/ASTAR for "Photographic
Emulsion" (matno=215). See the CSV headers under
``energy_loss/data/nist/`` for source URLs, density, composition and
the date the data was fetched. Future v0.x can replace these defaults
with experiment-specific calibrations without changing the API.
"""

from __future__ import annotations

from importlib import resources
from pathlib import Path

from energy_loss.range.table import RangeEnergyTable
from energy_loss.units import length_to_cm

_BUNDLED: dict[str, str] = {
  "proton": "pstar_photographic_emulsion.csv",
  "alpha": "astar_photographic_emulsion.csv",
}

_CACHE: dict[str, RangeEnergyTable] = {}


def list_bundled_emulsion_tables() -> list[str]:
  """Particle keys for which v0.3 bundles a default emulsion table."""
  return sorted(_BUNDLED)


def get_emulsion_range_energy(particle: str) -> RangeEnergyTable:
  """Return the default :class:`RangeEnergyTable` for ``particle``.

  Currently supported: ``"proton"`` (PSTAR) and ``"alpha"`` (ASTAR).
  """
  key = particle.lower()
  if key in _CACHE:
    return _CACHE[key]
  try:
    fname = _BUNDLED[key]
  except KeyError as exc:
    raise ValueError(
      f"No bundled emulsion range-energy table for particle {particle!r}. "
      f"Available: {list_bundled_emulsion_tables()}."
    ) from exc
  res = resources.files("energy_loss.data").joinpath("nist", fname)
  with resources.as_file(res) as path:
    table = RangeEnergyTable.from_nist_csv(Path(path))
  _CACHE[key] = table
  return table


def energy_from_emulsion_range(
  particle: str, range_value: float, range_unit: str = "um",
  table: RangeEnergyTable | None = None,
) -> float:
  """Estimate kinetic energy [MeV] from a measured emulsion track length.

  By default uses the bundled NIST table for ``particle``; pass an
  explicit ``table`` to use a different calibration.
  """
  tab = table if table is not None else get_emulsion_range_energy(particle)
  range_cm = length_to_cm(float(range_value), range_unit)
  return tab.energy_from_linear_range(range_cm)
