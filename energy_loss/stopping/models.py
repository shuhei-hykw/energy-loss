"""Stopping-power model dispatcher.

v0.1 only knows ``"bethe"``. The function exists so downstream modules
(:mod:`energy_loss.transport`, :mod:`energy_loss.range`) can be written
in a model-agnostic way once more models are added.
"""

from __future__ import annotations

from collections.abc import Callable

from energy_loss.materials import Material
from energy_loss.particles import Particle
from energy_loss.stopping.bethe import (
  bethe_linear_stopping_power,
  bethe_mass_stopping_power,
)

_LINEAR_MODELS: dict[str, Callable[[str | Particle, float, str | Material], float]] = {
  "bethe": bethe_linear_stopping_power,
}

_MASS_MODELS: dict[str, Callable[[str | Particle, float, str | Material], float]] = {
  "bethe": bethe_mass_stopping_power,
}


def linear_stopping_power(
  particle: str | Particle,
  kinetic_energy_mev: float,
  material: str | Material,
  model: str = "bethe",
) -> float:
  """Linear stopping power -dE/dx [MeV/cm] for the given model."""
  try:
    fn = _LINEAR_MODELS[model]
  except KeyError as exc:
    known = ", ".join(sorted(_LINEAR_MODELS))
    raise ValueError(
      f"Unknown stopping model {model!r}. Known models: {known}"
    ) from exc
  return fn(particle, kinetic_energy_mev, material)


def mass_stopping_power(
  particle: str | Particle,
  kinetic_energy_mev: float,
  material: str | Material,
  model: str = "bethe",
) -> float:
  """Mass stopping power -dE/(rho dx) [MeV cm^2/g] for the given model."""
  try:
    fn = _MASS_MODELS[model]
  except KeyError as exc:
    known = ", ".join(sorted(_MASS_MODELS))
    raise ValueError(
      f"Unknown stopping model {model!r}. Known models: {known}"
    ) from exc
  return fn(particle, kinetic_energy_mev, material)
