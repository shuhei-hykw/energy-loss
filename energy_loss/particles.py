"""Particle definitions.

v0.1 covers heavy charged particles (proton, pi, K, mu) plus electron
as a placeholder. Electron stopping power is NOT implemented in v0.1 —
Bethe formula for electrons differs and is out of scope.

Masses are from PDG 2024 (rest mass energy in MeV).
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Particle:
  """A charged particle.

  Attributes
  ----------
  name : str
    Canonical name (e.g. ``"proton"``).
  mass_mev : float
    Rest mass energy [MeV].
  charge : int
    Charge in units of the elementary charge (signed).
  """

  name: str
  mass_mev: float
  charge: int


_PARTICLES: dict[str, Particle] = {
  "proton": Particle("proton", 938.27208816, +1),
  "pion+": Particle("pion+", 139.57039, +1),
  "pion-": Particle("pion-", 139.57039, -1),
  "kaon+": Particle("kaon+", 493.677, +1),
  "kaon-": Particle("kaon-", 493.677, -1),
  "muon+": Particle("muon+", 105.6583755, +1),
  "muon-": Particle("muon-", 105.6583755, -1),
  "electron": Particle("electron", 0.51099895069, -1),
  "positron": Particle("positron", 0.51099895069, +1),
}

# Common aliases.
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
