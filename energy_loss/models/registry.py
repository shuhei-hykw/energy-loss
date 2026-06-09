"""Particle/material/model registry.

A *model* is a short string such as ``"nist_pstar"`` or ``"srim"`` that
identifies the source of a range-energy table. The registry is keyed by
``(particle_name, material_name, model_name)`` and stores a zero-arg
factory function so tables are only loaded on first use.

The ``"auto"`` policy resolves to a concrete model name per scope.md:

* proton -> ``nist_pstar``
* alpha  -> ``nist_astar``
* deuteron / triton / heavier ions -> ``srim`` (when available)

If no concrete model is registered for the requested triple, a clear
:class:`ValueError` is raised instead of silently falling back to an
unreliable default.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from importlib import resources
from pathlib import Path

from energy_loss.range.table import RangeEnergyTable

# Canonical particle alias map for the registry (so users can type
# "p" or "4He" and still land on a registered table). This is independent
# of the particle dataclass aliases; we want to keep the registry keys
# canonical without forcing every particle alias through it.
_PARTICLE_CANON: dict[str, str] = {
  "p": "proton",
  "proton": "proton",
  "4He": "alpha",
  "alpha": "alpha",
  "d": "deuteron",
  "deuteron": "deuteron",
  "t": "triton",
  "triton": "triton",
  "3He": "helion",
  "helion": "helion",
}

_MATERIAL_CANON: dict[str, str] = {
  "nuclear_emulsion": "nuclear_emulsion",
  "emulsion": "nuclear_emulsion",
  "photographic_emulsion": "nuclear_emulsion",
}

# scope.md model-selection policy: particle -> preferred model name.
AUTO_POLICY: dict[str, str] = {
  "proton": "nist_pstar",
  "alpha": "nist_astar",
}


@dataclass(frozen=True)
class ModelKey:
  """The canonical ``(particle, material, model)`` triple of a registry entry."""

  particle: str
  material: str
  model: str


_REGISTRY: dict[ModelKey, Callable[[], RangeEnergyTable]] = {}
_CACHE: dict[ModelKey, RangeEnergyTable] = {}


def _canonical_particle(name: str) -> str:
  return _PARTICLE_CANON.get(name, name)


def _canonical_material(name: str) -> str:
  return _MATERIAL_CANON.get(name, name)


def register_table_factory(
  particle: str, material: str, model: str,
  factory: Callable[[], RangeEnergyTable],
) -> None:
  """Register a lazy ``RangeEnergyTable`` factory for the given triple."""
  key = ModelKey(
    particle=_canonical_particle(particle),
    material=_canonical_material(material),
    model=model,
  )
  _REGISTRY[key] = factory


def resolve_model(particle: str, material: str, model: str) -> str:
  """Map a possibly-``auto`` ``model`` to a concrete registered name.

  Raises :class:`ValueError` if no entry is registered for the triple.
  """
  particle_c = _canonical_particle(particle)
  material_c = _canonical_material(material)
  if model == "auto":
    try:
      model = AUTO_POLICY[particle_c]
    except KeyError as exc:
      available = list_models(particle=particle_c, material=material_c)
      raise ValueError(
        f"No auto-policy default for particle {particle!r}; "
        f"registered models for ({particle_c}, {material_c}): {available}. "
        "Pass an explicit model name."
      ) from exc
  key = ModelKey(particle=particle_c, material=material_c, model=model)
  if key not in _REGISTRY:
    available = list_models(particle=particle_c, material=material_c)
    raise ValueError(
      f"No table registered for ({particle_c}, {material_c}, {model}). "
      f"Available models for this (particle, material): {available}."
    )
  return model


def load_table(
  particle: str, material: str, model: str = "auto",
) -> RangeEnergyTable:
  """Load and cache the :class:`RangeEnergyTable` for the triple."""
  concrete = resolve_model(particle, material, model)
  key = ModelKey(
    particle=_canonical_particle(particle),
    material=_canonical_material(material),
    model=concrete,
  )
  if key not in _CACHE:
    _CACHE[key] = _REGISTRY[key]()
  return _CACHE[key]


def list_models(
  particle: str | None = None, material: str | None = None,
) -> list[str]:
  """List registered model names, optionally filtered by particle/material."""
  particle_c = _canonical_particle(particle) if particle else None
  material_c = _canonical_material(material) if material else None
  names = sorted({
    key.model for key in _REGISTRY
    if (particle_c is None or key.particle == particle_c)
    and (material_c is None or key.material == material_c)
  })
  return names


# ---------------------------------------------------------------------------
# Built-in registrations.
# ---------------------------------------------------------------------------


def _bundled_nist_factory(filename: str) -> Callable[[], RangeEnergyTable]:
  def factory() -> RangeEnergyTable:
    res = resources.files("energy_loss.data").joinpath("nist", filename)
    with resources.as_file(res) as path:
      return RangeEnergyTable.from_nist_csv(Path(path))
  return factory


register_table_factory(
  "proton", "nuclear_emulsion", "nist_pstar",
  _bundled_nist_factory("pstar_photographic_emulsion.csv"),
)
register_table_factory(
  "alpha", "nuclear_emulsion", "nist_astar",
  _bundled_nist_factory("astar_photographic_emulsion.csv"),
)
