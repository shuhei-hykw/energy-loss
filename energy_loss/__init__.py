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

from energy_loss.config import (
  Beam,
  Config,
  Target,
  compute_linear_stopping_power,
  compute_mass_stopping_power,
  load_config,
  propagate_config,
)
from energy_loss.emulsion import (
  energy_from_emulsion_range,
  get_emulsion_range_energy,
  list_bundled_emulsion_tables,
)
from energy_loss.materials import (
  Material,
  get_material,
  list_materials,
  load_materials_from_yaml,
)
from energy_loss.particles import (
  Particle,
  get_particle,
  list_particles,
  register_particle,
)
from energy_loss.range import RangeEnergyTable, RangeEnergyTableMetadata
from energy_loss.transport import (
  Layer,
  LayerResult,
  PropagationResult,
  propagate,
)

__all__ = [
  "Beam",
  "Config",
  "Layer",
  "LayerResult",
  "Material",
  "Particle",
  "PropagationResult",
  "RangeEnergyTable",
  "RangeEnergyTableMetadata",
  "Target",
  "compute_linear_stopping_power",
  "compute_mass_stopping_power",
  "energy_from_emulsion_range",
  "get_emulsion_range_energy",
  "get_material",
  "get_particle",
  "list_bundled_emulsion_tables",
  "list_materials",
  "list_particles",
  "load_config",
  "load_materials_from_yaml",
  "propagate",
  "propagate_config",
  "register_particle",
]

__version__ = "0.1.0"
