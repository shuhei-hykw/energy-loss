"""Range and range-energy tables.

The submodule provides a generic :class:`RangeEnergyTable` that
interpolates between kinetic energy and range (CSDA range in grammage,
optionally combined with a density to give a linear range). The
bundled NIST PSTAR/ASTAR tables for Photographic Emulsion (matno=215)
are loaded via :func:`energy_loss.emulsion.get_emulsion_range_energy`.
"""

from energy_loss.range.table import RangeEnergyTable, RangeEnergyTableMetadata

__all__ = ["RangeEnergyTable", "RangeEnergyTableMetadata"]
