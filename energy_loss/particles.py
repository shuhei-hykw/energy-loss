"""Particle definitions.

Where :mod:`scipy.constants` provides a CODATA mass for a particle
(electron, muon, proton), that value is used directly. For mesons that
scipy does not carry (pions, kaons), PDG 2024 values are used; those
are the only hand-typed masses in this module.
"""

from __future__ import annotations

from dataclasses import dataclass

from scipy import constants as _sc


def _mev(key: str) -> float:
  """Helper: pull '<particle> mass energy equivalent in MeV' from scipy."""
  return _sc.physical_constants[key][0]


_ELECTRON_MEV = _mev("electron mass energy equivalent in MeV")
_MUON_MEV = _mev("muon mass energy equivalent in MeV")
_PROTON_MEV = _mev("proton mass energy equivalent in MeV")
_DEUTERON_MEV = _mev("deuteron mass energy equivalent in MeV")
_TRITON_MEV = _mev("triton mass energy equivalent in MeV")
_HELION_MEV = _mev("helion mass energy equivalent in MeV")
_ALPHA_MEV = _mev("alpha particle mass energy equivalent in MeV")

# PDG 2024 (not in scipy.constants):
#   pi+- : 139.57039 MeV/c^2
#   K+-  : 493.677  MeV/c^2
_PION_MEV = 139.57039
_KAON_MEV = 493.677


@dataclass(frozen=True)
class Particle:
  """A charged particle.

  Attributes
  ----------
  name : str
    Canonical name (e.g. ``"proton"``).
  mass_mev : float
    Rest-mass energy [MeV].
  charge : int
    Charge in units of the elementary charge (signed).
  """

  name: str
  mass_mev: float
  charge: int


_PARTICLES: dict[str, Particle] = {
  "proton": Particle("proton", _PROTON_MEV, +1),
  "pion+": Particle("pion+", _PION_MEV, +1),
  "pion-": Particle("pion-", _PION_MEV, -1),
  "kaon+": Particle("kaon+", _KAON_MEV, +1),
  "kaon-": Particle("kaon-", _KAON_MEV, -1),
  "muon+": Particle("muon+", _MUON_MEV, +1),
  "muon-": Particle("muon-", _MUON_MEV, -1),
  "electron": Particle("electron", _ELECTRON_MEV, -1),
  "positron": Particle("positron", _ELECTRON_MEV, +1),
  "deuteron": Particle("deuteron", _DEUTERON_MEV, +1),
  "triton": Particle("triton", _TRITON_MEV, +1),
  "helion": Particle("helion", _HELION_MEV, +2),       # 3He nucleus
  "alpha": Particle("alpha", _ALPHA_MEV, +2),           # 4He nucleus
}

_ALIASES: dict[str, str] = {
  "p": "proton",
  "pi+": "pion+",
  "pi-": "pion-",
  "K+": "kaon+",
  "K-": "kaon-",
  "mu+": "muon+",
  "mu-": "muon-",
  "e-": "electron",
  "e+": "positron",
  "d": "deuteron",
  "t": "triton",
  "3He": "helion",
  "4He": "alpha",
}


def get_particle(name: str | Particle) -> Particle:
  """Look up a :class:`Particle` by name (or pass through if already one)."""
  if isinstance(name, Particle):
    return name
  key = _ALIASES.get(name, name)
  try:
    return _PARTICLES[key]
  except KeyError as exc:
    known = ", ".join(sorted(_PARTICLES))
    raise ValueError(
      f"Unknown particle {name!r}. Known particles: {known}"
    ) from exc


def register_particle(
  name: str,
  mass_mev: float,
  charge: int,
  aliases: list[str] | None = None,
) -> Particle:
  """Register a new particle so :func:`get_particle` can find it.

  Used by :func:`energy_loss.config.load_config` to handle inline
  ``particles:`` blocks in a setup YAML. Later calls with the same name
  overwrite earlier ones, so users can override built-ins.
  """
  if mass_mev <= 0.0:
    raise ValueError(f"mass_mev must be positive, got {mass_mev}")
  p = Particle(name=name, mass_mev=float(mass_mev), charge=int(charge))
  _PARTICLES[name] = p
  for alias in aliases or []:
    _ALIASES[alias] = name
  return p


def list_particles() -> list[str]:
  """Return the canonical names of all currently registered particles."""
  return sorted(_PARTICLES)
