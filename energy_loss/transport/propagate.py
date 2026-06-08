"""RK4 step-wise integration of dT/d(rho x) = -S_mass(T).

The independent variable is grammage xi = rho * x [g/cm^2], so the
result is independent of how layer thicknesses are expressed and
multi-layer stacks compose cleanly. Step size defaults to total
grammage / 1000 and is bounded by ``max_steps_per_layer``; the user can
override either.

Stopping
--------
If the kinetic energy drops to (or below) ``min_kinetic_energy_mev``
the particle is considered to have stopped. The default is 0.1 MeV
because the basic Bethe formula is already past its validity threshold
well above that; the warning from
:mod:`energy_loss.stopping.bethe` will fire before the integrator gets
there. Callers that care about precise range / residual at very low
energy should switch to a table-based stopping model (planned v0.3).
"""

from __future__ import annotations

import warnings
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from energy_loss.particles import Particle, get_particle
from energy_loss.stopping.models import mass_stopping_power
from energy_loss.transport.layer import Layer

_DEFAULT_STEPS_PER_LAYER: int = 1000
_DEFAULT_MIN_KE_MEV: float = 0.1


@dataclass(frozen=True)
class LayerResult:
  """Energy bookkeeping for a single layer."""

  layer: Layer
  entry_kinetic_energy_mev: float
  exit_kinetic_energy_mev: float
  energy_loss_mev: float
  stopped_in_layer: bool


@dataclass(frozen=True)
class PropagationResult:
  """Aggregate result of propagating one particle through a layer stack."""

  initial_kinetic_energy_mev: float
  exit_kinetic_energy_mev: float
  total_energy_loss_mev: float
  stopped: bool
  per_layer: tuple[LayerResult, ...]

  @property
  def total_mass_thickness_g_per_cm2(self) -> float:
    return sum(r.layer.mass_thickness_g_per_cm2 for r in self.per_layer)


def _rk4_step(
  t_mev: float, h_g_per_cm2: float, particle: Particle,
  material_name: str, model: str,
) -> float:
  """One RK4 step on dT/d(xi) = -S_mass(T). Returns the new T [MeV]."""
  def s(t: float) -> float:
    if t <= 0.0:
      return 0.0
    return mass_stopping_power(particle, t, material_name, model=model)

  k1 = -s(t_mev)
  k2 = -s(t_mev + 0.5 * h_g_per_cm2 * k1)
  k3 = -s(t_mev + 0.5 * h_g_per_cm2 * k2)
  k4 = -s(t_mev + h_g_per_cm2 * k3)
  return t_mev + h_g_per_cm2 * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0


def _propagate_layer(
  t_in_mev: float, layer: Layer, particle: Particle, model: str,
  step_g_per_cm2: float | None, max_steps: int,
  min_kinetic_energy_mev: float,
) -> LayerResult:
  total_xi = layer.mass_thickness_g_per_cm2
  if total_xi <= 0.0:
    return LayerResult(
      layer=layer,
      entry_kinetic_energy_mev=t_in_mev,
      exit_kinetic_energy_mev=t_in_mev,
      energy_loss_mev=0.0,
      stopped_in_layer=False,
    )
  if step_g_per_cm2 is None:
    n = _DEFAULT_STEPS_PER_LAYER
  else:
    n = max(1, int(round(total_xi / step_g_per_cm2)))
  if n > max_steps:
    raise ValueError(
      f"layer requires {n} integration steps, which exceeds "
      f"max_steps_per_layer={max_steps}. Increase max_steps_per_layer "
      "or coarsen step_g_per_cm2."
    )
  h = total_xi / n

  t = t_in_mev
  stopped = False
  # Silence the low-beta-gamma warning from individual sub-Bethe steps
  # during the integration: the integrator emits its own once-per-layer
  # warning when the particle stops or enters the unreliable region.
  with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=UserWarning)
    for _ in range(n):
      t_next = _rk4_step(t, h, particle, layer.material.name, model)
      if t_next <= min_kinetic_energy_mev:
        t = 0.0
        stopped = True
        break
      t = t_next
  if stopped:
    warnings.warn(
      f"Particle stopped inside layer of {layer.material.name}.",
      UserWarning,
      stacklevel=3,
    )
  return LayerResult(
    layer=layer,
    entry_kinetic_energy_mev=t_in_mev,
    exit_kinetic_energy_mev=t,
    energy_loss_mev=t_in_mev - t,
    stopped_in_layer=stopped,
  )


def propagate(
  particle: str | Particle,
  initial_kinetic_energy_mev: float,
  layers: Iterable[Layer],
  model: str = "bethe",
  step_g_per_cm2: float | None = None,
  max_steps_per_layer: int = 1_000_000,
  min_kinetic_energy_mev: float = _DEFAULT_MIN_KE_MEV,
) -> PropagationResult:
  """Integrate kinetic energy through a stack of layers.

  Parameters
  ----------
  particle : str or Particle
    Projectile.
  initial_kinetic_energy_mev : float
    Kinetic energy at the entrance to the first layer [MeV].
  layers : iterable of Layer
    Material layers, in the order the particle traverses them.
  model : str
    Stopping-power model name (only ``"bethe"`` in v0.2).
  step_g_per_cm2 : float, optional
    Fixed step size in grammage. Default: ``layer.grammage / 1000``.
  max_steps_per_layer : int
    Safety cap on integration steps per layer.
  min_kinetic_energy_mev : float
    Treat the particle as stopped when ``T`` falls below this value.
  """
  p = get_particle(particle)
  if initial_kinetic_energy_mev <= 0.0:
    raise ValueError(
      f"initial kinetic energy must be positive, "
      f"got {initial_kinetic_energy_mev}"
    )
  layers_tuple: Sequence[Layer] = tuple(layers)
  if not layers_tuple:
    raise ValueError("propagate(): at least one layer is required.")

  t = initial_kinetic_energy_mev
  per_layer: list[LayerResult] = []
  stopped = False
  for layer in layers_tuple:
    if stopped:
      per_layer.append(
        LayerResult(
          layer=layer,
          entry_kinetic_energy_mev=0.0,
          exit_kinetic_energy_mev=0.0,
          energy_loss_mev=0.0,
          stopped_in_layer=False,
        )
      )
      continue
    result = _propagate_layer(
      t, layer, p, model, step_g_per_cm2,
      max_steps_per_layer, min_kinetic_energy_mev,
    )
    per_layer.append(result)
    t = result.exit_kinetic_energy_mev
    if result.stopped_in_layer:
      stopped = True

  return PropagationResult(
    initial_kinetic_energy_mev=initial_kinetic_energy_mev,
    exit_kinetic_energy_mev=t,
    total_energy_loss_mev=initial_kinetic_energy_mev - t,
    stopped=stopped,
    per_layer=tuple(per_layer),
  )
