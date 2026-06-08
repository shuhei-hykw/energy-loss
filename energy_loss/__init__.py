"""energy_loss: charged-particle energy loss, stopping power, range-energy relations.

Internal unit convention (used everywhere unless stated otherwise):
  - energy:     MeV
  - mass:       MeV/c^2
  - length:     cm
  - density:    g/cm^3
  - stopping power (mass):    MeV cm^2 / g
  - stopping power (linear):  MeV / cm
  - mean excitation energy I: eV  (kept in eV because that's how it's tabulated)

Public API re-exports the most commonly used objects.
"""

from energy_loss.materials import (
  Material,
  get_material,
  list_materials,
  load_materials_from_yaml,
)
from energy_loss.particles import Particle, get_particle

__all__ = [
  "Material",
  "Particle",
  "get_material",
  "get_particle",
  "list_materials",
  "load_materials_from_yaml",
]

__version__ = "0.1.0"
