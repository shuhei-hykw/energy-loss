"""Convenience access to range-energy data for nuclear emulsion.

v0.3 ships with the NIST PSTAR/ASTAR tables for *Photographic Emulsion*
(``matno=215``) — density 3.815 g/cm^3, mean excitation energy
331.0 eV, with the standard Ag/Br-dominated composition used by ICRU 49.
That table is the default for proton (PSTAR) and alpha (ASTAR) tracks;
for experiment-specific emulsion batches, load your own calibration
CSV with :meth:`energy_loss.range.RangeEnergyTable.from_nist_csv` (or
``from_arrays``) and pass it explicitly.
"""

from energy_loss.emulsion.range_energy import (
  energy_from_emulsion_range,
  get_emulsion_range_energy,
  list_bundled_emulsion_tables,
)

__all__ = [
  "energy_from_emulsion_range",
  "get_emulsion_range_energy",
  "list_bundled_emulsion_tables",
]
