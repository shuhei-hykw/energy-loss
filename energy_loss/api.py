"""High-level user-facing functions for range-energy conversion.

These wrap the model registry (:mod:`energy_loss.models`) with the
signatures suggested by ``scope.md``: ``energy_from_range``,
``range_from_energy``, ``compare_energy_from_range`` and ``load_table``.

Naming is uniform: arguments are ``particle``, ``material``, the
quantity in its natural unit (``range_um`` / ``energy_mev``) and an
optional ``model`` identifier. Units can be overridden via ``range_unit``
and ``energy_unit``.
"""

from __future__ import annotations

from energy_loss.models.registry import load_table as _load_table
from energy_loss.models.registry import resolve_model
from energy_loss.range.table import RangeEnergyTable
from energy_loss.units import energy_to_mev, length_to_cm, mev_to_energy

__all__ = [
  "compare_energy_from_range",
  "compare_range_from_energy",
  "energy_from_range",
  "load_table",
  "range_from_energy",
]


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

  Parameters
  ----------
  particle, material : str
    Looked up in the model registry. ``material`` accepts
    ``nuclear_emulsion`` and its synonyms.
  range_um : float, optional
    Range in micrometres. Shorthand for ``range_value`` with
    ``range_unit="um"``.
  range_value, range_unit : float, str, optional
    Range expressed in any unit understood by
    :mod:`energy_loss.units`.
  energy_unit : str
    Unit of the returned energy (default MeV).
  model : str
    Backend model name (default ``"auto"`` per scope.md policy).

  Returns
  -------
  float
    Kinetic energy in ``energy_unit``.
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

  Mirror of :func:`energy_from_range`. Returns the range in the
  requested ``range_unit`` (default micrometres).
  """
  if (energy_mev is None) == (energy_value is None):
    raise ValueError(
      "Specify exactly one of 'energy_mev' or 'energy_value'."
    )
  if energy_mev is not None:
    t_mev = float(energy_mev)
  else:
    t_mev = energy_to_mev(float(energy_value), energy_unit)
  table = _load_table(particle, material, model)
  r_g = table.range_from_energy(t_mev)
  rho = table.metadata.density_g_per_cm3
  if rho is None or rho <= 0.0:
    raise ValueError(
      "Table is missing a density; cannot convert grammage to linear range."
    )
  range_cm = r_g / rho
  # length conversions are linear: cm -> requested unit via inverse factor.
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
) -> dict[str, float | str]:
  """Run :func:`energy_from_range` across several models and report them.

  The returned dict carries ``"unit": energy_unit`` plus one numeric
  entry per resolved model name. ``models=None`` means *every*
  registered model for the (particle, material) pair (auto-policy
  is silently expanded so the result is reproducible).
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
    results[concrete] = energy_from_range(
      particle=particle, material=material,
      range_um=range_um, range_value=range_value,
      range_unit=range_unit, energy_unit=energy_unit, model=concrete,
    )
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
) -> dict[str, float | str]:
  """Mirror of :func:`compare_energy_from_range` for the inverse direction."""
  from energy_loss.models.registry import list_models

  if models is None:
    models = list_models(particle=particle, material=material)
    if not models:
      raise ValueError(
        f"No models registered for ({particle}, {material})."
      )
  results: dict[str, float | str] = {"unit": range_unit}
  for raw_name in models:
    concrete = resolve_model(particle, material, raw_name)
    results[concrete] = range_from_energy(
      particle=particle, material=material,
      energy_mev=energy_mev, energy_value=energy_value,
      energy_unit=energy_unit, range_unit=range_unit, model=concrete,
    )
  return results
