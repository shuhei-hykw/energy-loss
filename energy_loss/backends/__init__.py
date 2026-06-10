"""External-tool backends (Geant4, SRIM, ATIMA, ...).

Each backend module exposes a thin Python wrapper around an external
table generator. The wrappers run the generator via :mod:`subprocess`
and return the path of the produced CSV; loading the CSV back into a
:class:`energy_loss.range.RangeEnergyTable` is the job of the model
registry layer.
"""

from energy_loss.backends.geant4_runner import (
  Geant4TableSpec,
  detect_geant4_executable,
  generate_geant4_table,
)

__all__ = [
  "Geant4TableSpec",
  "detect_geant4_executable",
  "generate_geant4_table",
]
