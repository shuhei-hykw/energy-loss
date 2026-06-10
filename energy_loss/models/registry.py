"""Particle/material/model registry.

A *model* is a short string such as ``"nist_pstar"`` or ``"srim"`` that
identifies the source of a range-energy table. The registry is keyed by
``(particle_name, material_name, model_name)`` and stores a
:class:`RegistryEntry` carrying a lazy factory plus the model's
recommended energy range.

The ``"auto"`` policy resolves to a concrete model name. v0.6 makes the
policy *energy-aware*: when a kinetic energy is known at lookup time
(as it is for ``range_from_energy``), ``auto`` walks a preferred-order
list per particle and picks the first entry whose recommended range
covers that energy. When the energy is not yet known (as in
``energy_from_range``), the resolver returns the canonical preferred
model and the caller surfaces a warning if the result lands outside
the model's recommended range.

If no concrete model is registered for the requested triple, a clear
:class:`ValueError` is raised instead of silently falling back to an
unreliable default.
"""

from __future__ import annotations

import math
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

# Per-particle preferred-order list. ``auto`` picks the first entry
# that is registered for the (particle, material) pair *and* whose
# recommended energy range covers the requested kinetic energy. If the
# energy is not provided, ``auto`` returns the first registered entry
# regardless of range.
AUTO_POLICY: dict[str, list[str]] = {
  "proton": ["nist_pstar", "geant4_11_4_1", "geant4_atima_11_4_1"],
  "alpha": ["nist_astar", "geant4_11_4_1", "geant4_atima_11_4_1"],
  "deuteron": ["geant4_11_4_1", "geant4_atima_11_4_1"],
  "triton": ["geant4_11_4_1", "geant4_atima_11_4_1"],
  "helion": ["geant4_11_4_1"],
}

# Geant4-derived models live under names like ``geant4_11_4_1``. The
# tag the registry uses is set from this constant; we keep it as a name
# so user code can request ``model="geant4"`` and the resolver maps it
# to the bundled Geant4 backend automatically (similar to ``auto``).
GEANT4_MODEL_NAME: str = "geant4_11_4_1"
GEANT4_ATIMA_MODEL_NAME: str = "geant4_atima_11_4_1"


@dataclass(frozen=True)
class ModelKey:
  """The canonical ``(particle, material, model)`` triple of a registry entry."""

  particle: str
  material: str
  model: str


@dataclass(frozen=True)
class RegistryEntry:
  """One row of the model registry."""

  key: ModelKey
  factory: Callable[[], RangeEnergyTable]
  # Recommended energy range [MeV]. Outside this band the model still
  # returns a value (the underlying table will interpolate) but the
  # high-level API emits a UserWarning. Defaults to "no restriction".
  valid_min_mev: float = 0.0
  valid_max_mev: float = math.inf
  note: str = ""

  def covers(self, t_mev: float) -> bool:
    return self.valid_min_mev <= t_mev <= self.valid_max_mev


_REGISTRY: dict[ModelKey, RegistryEntry] = {}
_CACHE: dict[ModelKey, RangeEnergyTable] = {}


def _canonical_particle(name: str) -> str:
  return _PARTICLE_CANON.get(name, name)


def _canonical_material(name: str) -> str:
  return _MATERIAL_CANON.get(name, name)


def register_table_factory(
  particle: str, material: str, model: str,
  factory: Callable[[], RangeEnergyTable],
  valid_min_mev: float = 0.0,
  valid_max_mev: float = math.inf,
  note: str = "",
) -> None:
  """Register a lazy ``RangeEnergyTable`` factory for the given triple.

  ``valid_min_mev`` and ``valid_max_mev`` advertise the energy band in
  which the underlying model is considered reliable. They are advisory
  — out-of-range usage produces a :class:`UserWarning` from the
  high-level API, never an outright failure.
  """
  key = ModelKey(
    particle=_canonical_particle(particle),
    material=_canonical_material(material),
    model=model,
  )
  _REGISTRY[key] = RegistryEntry(
    key=key, factory=factory,
    valid_min_mev=float(valid_min_mev),
    valid_max_mev=float(valid_max_mev),
    note=note,
  )


def get_model_entry(
  particle: str, material: str, model: str,
) -> RegistryEntry:
  """Return the registered :class:`RegistryEntry` for the resolved triple."""
  concrete = resolve_model(particle, material, model)
  key = ModelKey(
    particle=_canonical_particle(particle),
    material=_canonical_material(material),
    model=concrete,
  )
  return _REGISTRY[key]


def resolve_model(
  particle: str, material: str, model: str,
  energy_mev: float | None = None,
) -> str:
  """Map a possibly-``auto`` ``model`` to a concrete registered name.

  When ``energy_mev`` is provided, the ``auto`` policy walks the
  per-particle preference list and picks the first entry whose
  recommended range covers ``energy_mev``. When it isn't, ``auto``
  returns the first preference that is registered, regardless of range.

  Raises :class:`ValueError` if no entry is registered for the triple.
  """
  particle_c = _canonical_particle(particle)
  material_c = _canonical_material(material)
  if model == "auto":
    preferences = AUTO_POLICY.get(particle_c, [])
    if not preferences:
      available = list_models(particle=particle_c, material=material_c)
      raise ValueError(
        f"No auto-policy default for particle {particle!r}; "
        f"registered models for ({particle_c}, {material_c}): {available}. "
        "Pass an explicit model name."
      )
    fallback: str | None = None
    for candidate in preferences:
      key = ModelKey(
        particle=particle_c, material=material_c, model=candidate,
      )
      if key not in _REGISTRY:
        continue
      if fallback is None:
        fallback = candidate
      if energy_mev is None or _REGISTRY[key].covers(float(energy_mev)):
        return candidate
    if fallback is not None:
      return fallback
    available = list_models(particle=particle_c, material=material_c)
    raise ValueError(
      f"No registered model for ({particle_c}, {material_c}); "
      f"preferences {preferences} are all missing. Available: {available}."
    )
  if model == "geant4":
    model = GEANT4_MODEL_NAME
  elif model == "geant4_atima":
    model = GEANT4_ATIMA_MODEL_NAME
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
  energy_mev: float | None = None,
) -> RangeEnergyTable:
  """Load and cache the :class:`RangeEnergyTable` for the triple."""
  concrete = resolve_model(particle, material, model, energy_mev=energy_mev)
  key = ModelKey(
    particle=_canonical_particle(particle),
    material=_canonical_material(material),
    model=concrete,
  )
  if key not in _CACHE:
    _CACHE[key] = _REGISTRY[key].factory()
  return _CACHE[key]


def list_models(
  particle: str | None = None, material: str | None = None,
  valid_at_mev: float | None = None,
) -> list[str]:
  """List registered model names, optionally filtered by particle/material.

  When ``valid_at_mev`` is provided, only models whose recommended
  range covers that energy are returned.
  """
  particle_c = _canonical_particle(particle) if particle else None
  material_c = _canonical_material(material) if material else None
  names = sorted({
    key.model for key, entry in _REGISTRY.items()
    if (particle_c is None or key.particle == particle_c)
    and (material_c is None or key.material == material_c)
    and (valid_at_mev is None or entry.covers(float(valid_at_mev)))
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


# NIST PSTAR proton table covers 1 keV to 10 GeV — full tabulated range.
register_table_factory(
  "proton", "nuclear_emulsion", "nist_pstar",
  _bundled_nist_factory("pstar_photographic_emulsion.csv"),
  valid_min_mev=1.0e-3, valid_max_mev=1.0e4,
  note="NIST PSTAR (ICRU 49); reference for proton in matter.",
)
# NIST ASTAR alpha table goes to 1 GeV; beyond that the data isn't tabulated.
register_table_factory(
  "alpha", "nuclear_emulsion", "nist_astar",
  _bundled_nist_factory("astar_photographic_emulsion.csv"),
  valid_min_mev=1.0e-3, valid_max_mev=1.0e3,
  note="NIST ASTAR (ICRU 49/73); reference for alpha in matter.",
)


# ---------------------------------------------------------------------------
# Geant4 backend registration (lazy — does not call subprocess at import).
# ---------------------------------------------------------------------------


def _geant4_factory(
  particle: str, material: str, emin_mev: float = 0.001,
  emax_mev: float = 1.0e4, n_points: int = 200,
  physics_list: str = "option4",
) -> Callable[[], RangeEnergyTable]:
  """Build a lazy factory that runs the Geant4 generator on first call."""

  def factory() -> RangeEnergyTable:
    # Imported here so importing the registry doesn't drag the backend
    # dependency tree in. Raises Geant4Unavailable if the executable
    # isn't built; resolve_model() reports that as a clear "not
    # registered" once the registration is skipped.
    from energy_loss.backends.geant4_runner import (
      Geant4TableSpec,
      generate_geant4_table,
    )

    spec = Geant4TableSpec(
      particle=particle, material=material,
      emin_mev=emin_mev, emax_mev=emax_mev,
      n_points=n_points, physics_list=physics_list,
    )
    csv_path = generate_geant4_table(spec)
    return RangeEnergyTable.from_nist_csv(csv_path)

  return factory


# Geant4 option4 — same ICRU 49/73 base as NIST. Reliable from a few
# keV to >10 GeV; effectively as good as PSTAR/ASTAR for these
# (particle, material) combinations.
_OPTION4_RANGE = dict(valid_min_mev=1.0e-3, valid_max_mev=1.0e4)
_OPTION4_NOTE = (
  "Geant4 G4EmStandardPhysics_option4 (ICRU 49/73-aligned)."
)

for _particle in ("proton", "alpha", "deuteron", "triton", "helion"):
  register_table_factory(
    _particle, "nuclear_emulsion", GEANT4_MODEL_NAME,
    _geant4_factory(_particle, "nuclear_emulsion"),
    **_OPTION4_RANGE, note=_OPTION4_NOTE,
  )
for _material in ("aluminum", "carbon", "Be"):
  for _particle in ("proton", "alpha"):
    register_table_factory(
      _particle, _material, GEANT4_MODEL_NAME,
      _geant4_factory(_particle, _material),
      **_OPTION4_RANGE, note=_OPTION4_NOTE,
    )

# Hadron beams for the J-PARC E10 9Be target comparison. Pi/K/mu in Be
# are valid down to ~1 MeV (below that EM transport breaks down for
# rapidly decaying particles) and up to >10 GeV.
_E10_RANGE = dict(valid_min_mev=1.0, valid_max_mev=1.0e4)
for _particle in ("pion-", "pion+", "kaon-", "kaon+", "muon-", "muon+"):
  register_table_factory(
    _particle, "Be", GEANT4_MODEL_NAME,
    _geant4_factory(_particle, "Be"),
    **_E10_RANGE, note=_OPTION4_NOTE,
  )

# Geant4 ATIMA (G4AtimaEnergyLossModel). The model header sets a low-
# energy cutoff at 2 MeV (sezi_dedx_e branch below 30 MeV/u, transition
# 10-30 MeV/u, Bethek_dedx_e above). Empirically the model disagrees
# with PSTAR/option4 by 15-25% in nuclear emulsion even at relativistic
# energies (different density-effect data; see Geant4 source). We
# advertise [2 MeV, 1 GeV] as the recommended range — it is the band
# the model is most defensible in. Geant4's range table also plateaus
# around 1 GeV so the factory caps emax there.
_ATIMA_EMAX_MEV = 1000.0
_ATIMA_RANGE = dict(valid_min_mev=2.0, valid_max_mev=_ATIMA_EMAX_MEV)
_ATIMA_NOTE = (
  "Geant4 G4AtimaEnergyLossModel; designed for heavy charged "
  "particles at intermediate energies; differs from PSTAR/option4."
)

for _particle in ("proton", "alpha", "deuteron", "triton"):
  register_table_factory(
    _particle, "nuclear_emulsion", GEANT4_ATIMA_MODEL_NAME,
    _geant4_factory(
      _particle, "nuclear_emulsion",
      physics_list="atima", emax_mev=_ATIMA_EMAX_MEV,
    ),
    **_ATIMA_RANGE, note=_ATIMA_NOTE,
  )
for _material in ("aluminum", "carbon", "Be"):
  for _particle in ("proton", "alpha"):
    register_table_factory(
      _particle, _material, GEANT4_ATIMA_MODEL_NAME,
      _geant4_factory(
        _particle, _material,
        physics_list="atima", emax_mev=_ATIMA_EMAX_MEV,
      ),
      **_ATIMA_RANGE, note=_ATIMA_NOTE,
    )
