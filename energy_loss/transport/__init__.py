"""Transport: integrate a particle's kinetic energy through layered matter.

The integrator solves ``dT/d(rho x) = -S_mass(T)`` along grammage so the
result is independent of how a thickness is expressed (linear vs.
mass-thickness). Particles that lose all their energy in a layer are
flagged as ``stopped`` and the propagation halts.
"""

from energy_loss.transport.layer import Layer
from energy_loss.transport.propagate import (
  LayerResult,
  PropagationResult,
  propagate,
)

__all__ = [
  "Layer",
  "LayerResult",
  "PropagationResult",
  "propagate",
]
