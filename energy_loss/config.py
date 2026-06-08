"""Single-file YAML setup describing a beam + target calculation.

The intent is that one YAML file fully captures a stopping-power /
energy-loss problem: which beam particle at what energy, hitting which
target. Optional ``particles:`` and ``materials:`` blocks let a user
define custom entries inline so a configuration is self-contained and
reproducible.

Schema
------
::

    # Optional: register new particles (added to the global registry).
    particles:
      <name>:
        mass_mev: <float>
        charge:   <int>
        aliases:  [<str>, ...]   # optional

    # Optional: register new materials (same schema as materials.yaml).
    materials:
      <name>:
        element: <symbol>                 # OR z_over_a + density_g_per_cm3
        mean_excitation_energy_ev: <float>
        density_g_per_cm3: <float>        # optional override
        aliases:  [<str>, ...]            # optional

    # Required: beam definition.
    beam:
      particle: <name>          # looked up via get_particle
      # Specify either kinetic energy or momentum; not both.
      kinetic_energy: <float>
      energy_unit: MeV          # default MeV
      # ... or ...
      momentum: <float>
      momentum_unit: GeV/c      # default MeV/c

    # Required: target definition.
    target:
      material: <name>
      # Optional thickness, specified either linearly or as grammage.
      thickness: <float>
      thickness_unit: cm        # default cm
      # ... or ...
      mass_thickness: <float>
      mass_thickness_unit: g/cm^2   # default g/cm^2
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml as _yaml

from energy_loss.materials import (
  Material,
  get_material,
)
from energy_loss.materials import (
  _ingest as _ingest_materials,
)
from energy_loss.particles import Particle, get_particle, register_particle
from energy_loss.stopping.models import (
  linear_stopping_power,
  mass_stopping_power,
)
from energy_loss.units import (
  energy_to_mev,
  kinetic_energy_from_momentum,
  length_to_cm,
  mass_thickness_to_g_per_cm2,
  momentum_to_mev_c,
)


@dataclass(frozen=True)
class Beam:
  """A beam: a :class:`Particle` plus its initial kinematics.

  Both ``kinetic_energy_mev`` and ``momentum_mev_c`` are populated so
  downstream code can use whichever is more convenient.
  """

  particle: Particle
  kinetic_energy_mev: float
  momentum_mev_c: float


@dataclass(frozen=True)
class Target:
  """A target: a :class:`Material` and, optionally, its thickness.

  Both linear thickness [cm] and mass thickness [g/cm^2] are populated
  whenever a thickness is given in the config, regardless of which form
  the user wrote.
  """

  material: Material
  thickness_cm: float | None = None
  mass_thickness_g_per_cm2: float | None = None


@dataclass(frozen=True)
class Config:
  """Fully resolved beam + target configuration loaded from a YAML file."""

  beam: Beam
  target: Target
  source_path: Path | None = None


def _register_particles_from_doc(doc: dict[str, Any]) -> None:
  particles = doc.get("particles") or {}
  if not isinstance(particles, dict):
    raise ValueError("'particles' must be a mapping (name -> entry).")
  for name, entry in particles.items():
    if not isinstance(entry, dict):
      raise ValueError(f"Particle {name!r}: entry must be a mapping.")
    try:
      mass = entry["mass_mev"]
      charge = entry["charge"]
    except KeyError as exc:
      raise ValueError(
        f"Particle {name!r}: missing required field {exc.args[0]!r}."
      ) from exc
    register_particle(
      name=name,
      mass_mev=float(mass),
      charge=int(charge),
      aliases=list(entry.get("aliases") or []),
    )


def _register_materials_from_doc(doc: dict[str, Any]) -> None:
  if doc.get("materials"):
    _ingest_materials({"materials": doc["materials"]})


def _build_beam(doc: dict[str, Any]) -> Beam:
  try:
    section = doc["beam"]
  except KeyError as exc:
    raise ValueError("Config YAML is missing required 'beam' section.") from exc
  if not isinstance(section, dict):
    raise ValueError("'beam' must be a mapping.")
  try:
    particle_name = section["particle"]
  except KeyError as exc:
    raise ValueError(
      f"beam: missing required field {exc.args[0]!r}."
    ) from exc
  particle = get_particle(str(particle_name))

  has_ke = "kinetic_energy" in section
  has_mom = "momentum" in section
  if has_ke == has_mom:
    raise ValueError(
      "beam: specify exactly one of 'kinetic_energy' or 'momentum'."
    )
  if has_ke:
    unit = section.get("energy_unit", "MeV")
    ke_mev = energy_to_mev(float(section["kinetic_energy"]), str(unit))
    import math as _math
    total_e = ke_mev + particle.mass_mev
    p_mev_c = _math.sqrt(max(total_e * total_e - particle.mass_mev**2, 0.0))
  else:
    unit = section.get("momentum_unit", "MeV/c")
    p_mev_c = momentum_to_mev_c(float(section["momentum"]), str(unit))
    ke_mev = kinetic_energy_from_momentum(p_mev_c, particle.mass_mev)

  return Beam(
    particle=particle,
    kinetic_energy_mev=ke_mev,
    momentum_mev_c=p_mev_c,
  )


def _build_target(doc: dict[str, Any]) -> Target:
  try:
    section = doc["target"]
  except KeyError as exc:
    raise ValueError("Config YAML is missing required 'target' section.") from exc
  if not isinstance(section, dict):
    raise ValueError("'target' must be a mapping.")
  try:
    material_name = section["material"]
  except KeyError as exc:
    raise ValueError(
      f"target: missing required field {exc.args[0]!r}."
    ) from exc
  material = get_material(str(material_name))

  has_thick = "thickness" in section
  has_mt = "mass_thickness" in section
  if has_thick and has_mt:
    raise ValueError(
      "target: specify at most one of 'thickness' or 'mass_thickness'."
    )

  thickness_cm: float | None = None
  mass_thickness: float | None = None
  if has_thick:
    unit = section.get("thickness_unit", "cm")
    thickness_cm = length_to_cm(float(section["thickness"]), str(unit))
    mass_thickness = thickness_cm * material.density_g_per_cm3
  elif has_mt:
    unit = section.get("mass_thickness_unit", "g/cm^2")
    mass_thickness = mass_thickness_to_g_per_cm2(
      float(section["mass_thickness"]), str(unit)
    )
    thickness_cm = mass_thickness / material.density_g_per_cm3

  return Target(
    material=material,
    thickness_cm=thickness_cm,
    mass_thickness_g_per_cm2=mass_thickness,
  )


def load_config(path: str | Path) -> Config:
  """Load a setup YAML and return the fully resolved :class:`Config`.

  Side effect: any ``particles:`` or ``materials:`` block in the file is
  merged into the global registry, so subsequent calls to
  :func:`get_particle` / :func:`get_material` can also see them.
  """
  path = Path(path)
  with open(path, encoding="utf-8") as f:
    doc = _yaml.safe_load(f) or {}
  if not isinstance(doc, dict):
    raise ValueError(f"{path}: top-level YAML must be a mapping.")

  _register_particles_from_doc(doc)
  _register_materials_from_doc(doc)
  beam = _build_beam(doc)
  target = _build_target(doc)
  return Config(beam=beam, target=target, source_path=path)


def compute_mass_stopping_power(config: Config, model: str = "bethe") -> float:
  """Convenience: mass stopping power [MeV cm^2/g] for ``config``."""
  return mass_stopping_power(
    config.beam.particle,
    config.beam.kinetic_energy_mev,
    config.target.material,
    model=model,
  )


def compute_linear_stopping_power(config: Config, model: str = "bethe") -> float:
  """Convenience: linear stopping power [MeV/cm] for ``config``."""
  return linear_stopping_power(
    config.beam.particle,
    config.beam.kinetic_energy_mev,
    config.target.material,
    model=model,
  )
