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

    # Required: target definition. Two equivalent forms.
    #
    # (a) Single layer (back-compat with v0.1):
    target:
      material: <name>
      # Optional thickness, specified either linearly or as grammage.
      thickness: <float>
      thickness_unit: cm        # default cm
      # ... or ...
      mass_thickness: <float>
      mass_thickness_unit: g/cm^2   # default g/cm^2
    #
    # (b) Layer stack (v0.2+):
    target:
      layers:
        - material: kapton
          thickness: 50
          thickness_unit: um
        - material: Be
          mass_thickness: 3.5
          mass_thickness_unit: g/cm^2
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
from energy_loss.transport import Layer, PropagationResult, propagate
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
  """A target as an ordered stack of :class:`Layer`.

  Single-material targets are represented as a single-element tuple.
  Convenience properties forward to that single layer; multi-layer
  consumers should iterate :attr:`layers` directly.
  """

  layers: tuple[Layer, ...]

  @property
  def is_single_layer(self) -> bool:
    return len(self.layers) == 1

  @property
  def material(self) -> Material:
    """Material of the (first) layer.

    For multi-layer stacks this returns the first layer's material; the
    caller is responsible for handling the rest.
    """
    return self.layers[0].material

  @property
  def thickness_cm(self) -> float:
    """Linear thickness of the first layer [cm]."""
    return self.layers[0].thickness_cm

  @property
  def mass_thickness_g_per_cm2(self) -> float:
    """Grammage of the first layer [g/cm^2]."""
    return self.layers[0].mass_thickness_g_per_cm2

  @property
  def total_mass_thickness_g_per_cm2(self) -> float:
    """Sum of grammage over all layers [g/cm^2]."""
    return sum(layer.mass_thickness_g_per_cm2 for layer in self.layers)


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


def _layer_from_section(section: dict[str, Any], where: str) -> Layer:
  try:
    material_name = section["material"]
  except KeyError as exc:
    raise ValueError(
      f"{where}: missing required field {exc.args[0]!r}."
    ) from exc
  material = get_material(str(material_name))

  has_thick = "thickness" in section
  has_mt = "mass_thickness" in section
  if has_thick == has_mt:
    raise ValueError(
      f"{where}: specify exactly one of 'thickness' or 'mass_thickness'."
    )
  if has_thick:
    unit = section.get("thickness_unit", "cm")
    thickness_cm = length_to_cm(float(section["thickness"]), str(unit))
    return Layer.from_thickness(material, thickness_cm)
  unit = section.get("mass_thickness_unit", "g/cm^2")
  mass_thickness = mass_thickness_to_g_per_cm2(
    float(section["mass_thickness"]), str(unit)
  )
  return Layer.from_mass_thickness(material, mass_thickness)


def _build_target(doc: dict[str, Any]) -> Target:
  try:
    section = doc["target"]
  except KeyError as exc:
    raise ValueError("Config YAML is missing required 'target' section.") from exc
  if not isinstance(section, dict):
    raise ValueError("'target' must be a mapping.")

  if "layers" in section:
    if "material" in section or "thickness" in section or "mass_thickness" in section:
      raise ValueError(
        "target: 'layers' cannot be combined with 'material' / "
        "'thickness' / 'mass_thickness' at the same level."
      )
    raw_layers = section["layers"]
    if not isinstance(raw_layers, list) or not raw_layers:
      raise ValueError("target: 'layers' must be a non-empty list.")
    layers = tuple(
      _layer_from_section(entry, where=f"target.layers[{i}]")
      for i, entry in enumerate(raw_layers)
    )
  else:
    layers = (_layer_from_section(section, where="target"),)

  return Target(layers=layers)


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
  """Convenience: mass stopping power [MeV cm^2/g] at the beam's initial
  kinetic energy, in the *first* layer of the target. This is the v0.1
  single-point evaluator; for accurate energy loss through a layer
  stack use :func:`propagate_config`.
  """
  return mass_stopping_power(
    config.beam.particle,
    config.beam.kinetic_energy_mev,
    config.target.layers[0].material,
    model=model,
  )


def compute_linear_stopping_power(config: Config, model: str = "bethe") -> float:
  """Linear stopping power [MeV/cm] counterpart to
  :func:`compute_mass_stopping_power`.
  """
  return linear_stopping_power(
    config.beam.particle,
    config.beam.kinetic_energy_mev,
    config.target.layers[0].material,
    model=model,
  )


def propagate_config(
  config: Config,
  model: str = "bethe",
  step_g_per_cm2: float | None = None,
) -> PropagationResult:
  """Run the step-wise integrator on ``config``'s beam and target stack."""
  return propagate(
    config.beam.particle,
    config.beam.kinetic_energy_mev,
    config.target.layers,
    model=model,
    step_g_per_cm2=step_g_per_cm2,
  )
