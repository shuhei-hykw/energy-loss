"""Model registry for range-energy tables.

The package treats every stopping-power / range-energy source — NIST
PSTAR/ASTAR, SRIM, ATIMA, Geant4-generated CSV, the basic Bethe formula
— as a *model* identified by a short string. The registry maps
``(particle, material, model_name)`` triples to factory functions that
return a :class:`energy_loss.range.RangeEnergyTable`.

This is the layer that the high-level
:func:`energy_loss.api.energy_from_range` API talks to, so adding a new
backend (e.g. SRIM) is just a matter of registering tables here without
touching the user-facing entry points.

v0.3.1 bundles:

* ``nist_pstar`` for ``proton`` in ``nuclear_emulsion``
* ``nist_astar`` for ``alpha`` in ``nuclear_emulsion``

Future v0.4+ will add SRIM / ATIMA / Geant4 sources.
"""

from energy_loss.models.registry import (
  AUTO_POLICY,
  ModelKey,
  list_models,
  load_table,
  register_table_factory,
  resolve_model,
)

__all__ = [
  "AUTO_POLICY",
  "ModelKey",
  "list_models",
  "load_table",
  "register_table_factory",
  "resolve_model",
]
