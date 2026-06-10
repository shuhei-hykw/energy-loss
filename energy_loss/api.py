"""High-level user-facing functions for range-energy conversion.

These wrap the model registry (:mod:`energy_loss.models`) with the
signatures suggested by ``scope.md``: ``energy_from_range``,
``range_from_energy``, ``compare_energy_from_range`` and ``load_table``.

Naming is uniform: arguments are ``particle``, ``material``, the
quantity in its natural unit (``range_um`` / ``energy_mev``) and an
optional ``model`` identifier. Units can be overridden via ``range_unit``
and ``energy_unit``.

v0.6 makes the model selection energy-aware:

* ``range_from_energy(..., model="auto")`` knows the kinetic energy at
  resolve time and picks the first preferred model whose recommended
  range covers it.
* ``energy_from_range(..., model="auto")`` only learns the energy after
  the table lookup; if the result falls outside the model's
  recommended range the call emits a :class:`UserWarning` but still
  returns the value.
* ``compare_*`` gain a ``skip_out_of_range`` keyword that drops models
  outside their recommended band from the report.
"""

from __future__ import annotations

import warnings

from energy_loss.models.registry import get_model_entry, resolve_model
from energy_loss.models.registry import load_table as _load_table
from energy_loss.range.table import RangeEnergyTable
from energy_loss.units import energy_to_mev, length_to_cm, mev_to_energy

__all__ = [
  "compare_energy_from_range",
  "compare_range_from_energy",
  "energy_from_range",
  "load_table",
  "range_from_energy",
]


def _warn_out_of_range(
  particle: str, material: str, model: str, t_mev: float,
) -> None:
  entry = get_model_entry(particle, material, model)
  if not entry.covers(t_mev):
    warnings.warn(
      f"Model '{entry.key.model}' used outside its recommended energy "
      f"range [{entry.valid_min_mev:g}, {entry.valid_max_mev:g}] MeV "
      f"for ({entry.key.particle}, {entry.key.material}); "
      f"queried T = {t_mev:g} MeV. {entry.note}".rstrip(),
      UserWarning,
      stacklevel=3,
    )


def load_table(
  particle: str, material: str, model: str = "auto"
) -> RangeEnergyTable:
  """Return the :class:`RangeEnergyTable` for ``(particle, material, model)``."""
  return _load_table(particle, material, model)


def energy_from_range(
  particle: str,
  material: str,
  range_um: float | None = None,
  *,
  range_value: float | None = None,
  range_unit: str = "um",
  energy_unit: str = "MeV",
  model: str = "auto",
) -> float:
  """Convert a measured range to kinetic energy.

  The kinetic energy is not known until after the table lookup, so
  ``model="auto"`` picks the registry's first preferred model and the
  function emits a :class:`UserWarning` if the resulting energy falls
  outside that model's recommended range.
  """
  if (range_um is None) == (range_value is None):
    raise ValueError(
      "Specify exactly one of 'range_um' or 'range_value'."
    )
  if range_um is not None:
    range_cm = length_to_cm(float(range_um), "um")
  else:
    range_cm = length_to_cm(float(range_value), range_unit)
  table = _load_table(particle, material, model)
  t_mev = table.energy_from_linear_range(range_cm)
  _warn_out_of_range(particle, material, model, t_mev)
  return mev_to_energy(t_mev, energy_unit)


def range_from_energy(
  particle: str,
  material: str,
  energy_mev: float | None = None,
  *,
  energy_value: float | None = None,
  energy_unit: str = "MeV",
  range_unit: str = "um",
  model: str = "auto",
) -> float:
  """Convert a kinetic energy to range.

  ``model="auto"`` walks the per-particle preference list and picks the
  first registered model whose recommended range covers ``energy_mev``;
  if every preferred model is out of range, the first registered one is
  used and a :class:`UserWarning` is emitted.
  """
  if (energy_mev is None) == (energy_value is None):
    raise ValueError(
      "Specify exactly one of 'energy_mev' or 'energy_value'."
    )
  if energy_mev is not None:
    t_mev = float(energy_mev)
  else:
    t_mev = energy_to_mev(float(energy_value), energy_unit)
  resolved = resolve_model(particle, material, model, energy_mev=t_mev)
  _warn_out_of_range(particle, material, resolved, t_mev)
  table = _load_table(particle, material, resolved)
  r_g = table.range_from_energy(t_mev)
  rho = table.metadata.density_g_per_cm3
  if rho is None or rho <= 0.0:
    raise ValueError(
      "Table is missing a density; cannot convert grammage to linear range."
    )
  range_cm = r_g / rho
  from energy_loss.units import cm_to_length

  return cm_to_length(range_cm, range_unit)


def compare_energy_from_range(
  particle: str,
  material: str,
  range_um: float | None = None,
  *,
  range_value: float | None = None,
  range_unit: str = "um",
  energy_unit: str = "MeV",
  models: list[str] | None = None,
  skip_out_of_range: bool = False,
) -> dict[str, float | str]:
  """Run :func:`energy_from_range` across several models and report them.

  Returns a mapping with ``"unit": energy_unit`` plus one numeric entry
  per *resolved* model name. ``models=None`` means every registered
  model for the (particle, material) pair (the report is reproducible).

  Set ``skip_out_of_range=True`` to drop models whose recommended range
  does not cover the resulting kinetic energy. Out-of-range warnings
  are silenced in that case (the explicit skip is treated as the user
  acknowledging the limit).
  """
  from energy_loss.models.registry import list_models

  if models is None:
    models = list_models(particle=particle, material=material)
    if not models:
      raise ValueError(
        f"No models registered for ({particle}, {material})."
      )
  results: dict[str, float | str] = {"unit": energy_unit}
  for raw_name in models:
    concrete = resolve_model(particle, material, raw_name)
    with warnings.catch_warnings():
      warnings.simplefilter("ignore", UserWarning)
      t = energy_from_range(
        particle=particle, material=material,
        range_um=range_um, range_value=range_value,
        range_unit=range_unit, energy_unit=energy_unit, model=concrete,
      )
    entry = get_model_entry(particle, material, concrete)
    t_mev_for_check = (
      t if energy_unit == "MeV" else energy_to_mev(t, energy_unit)
    )
    if skip_out_of_range and not entry.covers(t_mev_for_check):
      continue
    results[concrete] = t
  return results


def compare_range_from_energy(
  particle: str,
  material: str,
  energy_mev: float | None = None,
  *,
  energy_value: float | None = None,
  energy_unit: str = "MeV",
  range_unit: str = "um",
  models: list[str] | None = None,
  skip_out_of_range: bool = False,
) -> dict[str, float | str]:
  """Mirror of :func:`compare_energy_from_range` for the inverse direction."""
  from energy_loss.models.registry import list_models

  if energy_mev is not None:
    t_mev = float(energy_mev)
  elif energy_value is not None:
    t_mev = energy_to_mev(float(energy_value), energy_unit)
  else:
    t_mev = None  # let range_from_energy raise the usual error.

  if models is None:
    models = list_models(particle=particle, material=material)
    if not models:
      raise ValueError(
        f"No models registered for ({particle}, {material})."
      )
  results: dict[str, float | str] = {"unit": range_unit}
  for raw_name in models:
    concrete = resolve_model(particle, material, raw_name)
    entry = get_model_entry(particle, material, concrete)
    if skip_out_of_range and t_mev is not None and not entry.covers(t_mev):
      continue
    with warnings.catch_warnings():
      warnings.simplefilter("ignore", UserWarning)
      results[concrete] = range_from_energy(
        particle=particle, material=material,
        energy_mev=energy_mev, energy_value=energy_value,
        energy_unit=energy_unit, range_unit=range_unit, model=concrete,
      )
  return results
