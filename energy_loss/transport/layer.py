"""A single uniform layer of material with a thickness."""

from __future__ import annotations

from dataclasses import dataclass

from energy_loss.materials import Material, get_material


@dataclass(frozen=True)
class Layer:
  """A single uniform layer of a material.

  Stores both linear thickness [cm] and mass thickness (grammage)
  [g/cm^2] so the integrator can use whichever is more natural without
  recomputing. Either can be ``None`` only if the other is provided.
  """

  material: Material
  thickness_cm: float
  mass_thickness_g_per_cm2: float

  @staticmethod
  def from_thickness(
    material: str | Material, thickness_cm: float
  ) -> Layer:
    """Construct a layer from a linear thickness."""
    m = get_material(material)
    return Layer(
      material=m,
      thickness_cm=float(thickness_cm),
      mass_thickness_g_per_cm2=float(thickness_cm) * m.density_g_per_cm3,
    )

  @staticmethod
  def from_mass_thickness(
    material: str | Material, mass_thickness_g_per_cm2: float
  ) -> Layer:
    """Construct a layer from a grammage."""
    m = get_material(material)
    return Layer(
      material=m,
      thickness_cm=float(mass_thickness_g_per_cm2) / m.density_g_per_cm3,
      mass_thickness_g_per_cm2=float(mass_thickness_g_per_cm2),
    )
