"""Material registry.

Materials are defined declaratively in
``energy_loss/data/materials.yaml`` (and any user YAML files registered
via :func:`load_materials_from_yaml`). The loader resolves pure-element
entries via the :mod:`periodictable` package so atomic number, atomic
mass and (where applicable) density don't have to be repeated by hand.
"""

from __future__ import annotations

from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Any

import periodictable as _pt
import yaml as _yaml


@dataclass(frozen=True)
class Material:
  """A material treated as a uniform stopping medium.

  Attributes
  ----------
  name : str
    Canonical name.
  z_over_a : float
    Effective <Z/A>. Dimensionless.
  mean_excitation_energy_ev : float
    Mean excitation energy I [eV].
  density_g_per_cm3 : float
    Density [g/cm^3].
  reference : str
    Short note on the source of the values (free-form).
  """

  name: str
  z_over_a: float
  mean_excitation_energy_ev: float
  density_g_per_cm3: float
  reference: str = ""


_MATERIALS: dict[str, Material] = {}
_ALIASES: dict[str, str] = {}


def _resolve_element_defaults(symbol: str) -> tuple[float, float]:
  """Return ``(z_over_a, density_g_per_cm3)`` for an element symbol.

  Uses :mod:`periodictable`. Raises ``ValueError`` if the symbol is
  unknown or the element has no density listed.
  """
  try:
    el = _pt.elements.symbol(symbol)
  except ValueError as exc:
    raise ValueError(f"Unknown element symbol {symbol!r}") from exc
  z_over_a = el.number / el.mass
  density = el.density
  if density is None:
    raise ValueError(
      f"periodictable does not provide a density for element {symbol!r}; "
      "specify density_g_per_cm3 explicitly in the YAML entry."
    )
  return float(z_over_a), float(density)


def _build_material(name: str, entry: dict[str, Any]) -> Material:
  """Construct a :class:`Material` from one YAML entry."""
  element = entry.get("element")
  if element is not None:
    z_over_a_default, density_default = _resolve_element_defaults(element)
  else:
    z_over_a_default, density_default = None, None

  z_over_a = entry.get("z_over_a", z_over_a_default)
  density = entry.get("density_g_per_cm3", density_default)

  if z_over_a is None:
    raise ValueError(
      f"Material {name!r}: missing 'z_over_a' (and no 'element' fallback)."
    )
  if density is None:
    raise ValueError(
      f"Material {name!r}: missing 'density_g_per_cm3' "
      "(and no 'element' fallback)."
    )
  try:
    i_ev = entry["mean_excitation_energy_ev"]
  except KeyError as exc:
    raise ValueError(
      f"Material {name!r}: 'mean_excitation_energy_ev' is required."
    ) from exc

  return Material(
    name=name,
    z_over_a=float(z_over_a),
    mean_excitation_energy_ev=float(i_ev),
    density_g_per_cm3=float(density),
    reference=str(entry.get("reference", "")),
  )


def _ingest(doc: dict[str, Any]) -> None:
  """Add all materials from a parsed YAML document into the registry."""
  try:
    materials = doc["materials"]
  except KeyError as exc:
    raise ValueError(
      "YAML document must have a top-level 'materials' mapping."
    ) from exc
  if not isinstance(materials, dict):
    raise ValueError("'materials' must be a mapping (name -> entry).")

  for name, entry in materials.items():
    if not isinstance(entry, dict):
      raise ValueError(f"Material {name!r}: entry must be a mapping.")
    mat = _build_material(name, entry)
    _MATERIALS[name] = mat
    for alias in entry.get("aliases", []) or []:
      _ALIASES[alias] = name


def load_materials_from_yaml(path: str | Path) -> None:
  """Register additional materials from a user-supplied YAML file.

  Later entries with the same name overwrite earlier ones, so users can
  override built-in defaults (e.g. for an experiment-specific emulsion
  calibration).
  """
  with open(path, encoding="utf-8") as f:
    doc = _yaml.safe_load(f)
  _ingest(doc)


def _load_builtin() -> None:
  """Load the bundled ``materials.yaml`` once at import time."""
  with resources.files("energy_loss.data").joinpath("materials.yaml").open(
    "r", encoding="utf-8"
  ) as f:
    doc = _yaml.safe_load(f)
  _ingest(doc)


_load_builtin()


def get_material(name: str | Material) -> Material:
  """Look up a :class:`Material` by name (or pass through if already one)."""
  if isinstance(name, Material):
    return name
  key = _ALIASES.get(name, name)
  try:
    return _MATERIALS[key]
  except KeyError as exc:
    known = ", ".join(sorted(_MATERIALS))
    raise ValueError(
      f"Unknown material {name!r}. Known materials: {known}"
    ) from exc


def list_materials() -> list[str]:
  """Return the canonical names of all currently registered materials."""
  return sorted(_MATERIALS)
