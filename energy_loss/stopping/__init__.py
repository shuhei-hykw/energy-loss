"""Stopping power models.

Currently provides the Bethe formula for heavy charged particles.
"""

from energy_loss.stopping.bethe import (
  bethe_linear_stopping_power,
  bethe_mass_stopping_power,
)

__all__ = [
  "bethe_linear_stopping_power",
  "bethe_mass_stopping_power",
]
