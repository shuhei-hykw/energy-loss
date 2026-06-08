"""Bethe formula for heavy charged particles.

The implementation here is the *basic* Bethe-Bloch formula without
density-effect correction, shell correction, Barkas, Bloch, or
effective-charge corrections. Those are explicitly out of scope for
v0.1 but the API leaves placeholders for them.

Formula (PDG eq. 34.5, simplified — no delta/2 term):

  -dE/dx [MeV cm^2 / g]
    = K * z^2 * (Z/A) / beta^2
        * ( 0.5 * ln( 2 m_e c^2 beta^2 gamma^2 T_max / I^2 ) - beta^2 )

with

  K           = 4 pi N_A r_e^2 m_e c^2   ~= 0.307075 MeV mol / g
  z           = projectile charge (in units of e)
  Z, A        = atomic number / mass number of medium
  I           = mean excitation energy of medium
  beta, gamma = projectile velocity / Lorentz factor
  T_max       = maximum kinetic energy transferable to an electron
                in a single collision

  T_max = 2 m_e c^2 beta^2 gamma^2
          / ( 1 + 2 gamma m_e / M + (m_e / M)^2 )

where M is the projectile rest mass.

Validity & caveats
------------------
- Intended for *heavy* charged particles (M >> m_e). Electrons and
  positrons are NOT handled by this function.
- The formula breaks down at low velocity (typically below
  beta * gamma ~ 0.05, i.e. roughly < 1 MeV/nucleon protons). In
  that regime a warning is issued; for accurate range / dE/dx at
  low energy, use a table-based approach (SRIM, ATIMA, NIST PSTAR)
  via :mod:`energy_loss.range`.
- No density-effect correction. At high energy (beta*gamma ≳ 10)
  this leads to an over-estimate of dE/dx; that correction is a
  planned v0.2 addition.

Internal units (see :mod:`energy_loss`):
  energy: MeV, length: cm, density: g/cm^3, I: eV.
"""

from __future__ import annotations

import math
import warnings

from energy_loss.constants import ELECTRON_MASS_MEV, K_BETHE_MEV_MOL_PER_G
from energy_loss.materials import Material, get_material
from energy_loss.particles import Particle, get_particle

# Below this beta*gamma the basic Bethe formula is not trustworthy and
# we emit a warning. The exact threshold is not sharp; 0.05 corresponds
# roughly to T ~ 1 MeV/nucleon for protons.
_LOW_BETA_GAMMA_WARN: float = 0.05


def _beta_gamma_from_kinetic(kinetic_energy_mev: float, mass_mev: float) -> tuple[float, float]:
  """Return (beta, gamma) for a particle of mass M and kinetic energy T.

  Both inputs are in MeV (mass means rest mass energy M c^2).
  """
  if kinetic_energy_mev <= 0.0:
    raise ValueError(
      f"kinetic energy must be positive, got {kinetic_energy_mev}"
    )
  if mass_mev <= 0.0:
    raise ValueError(f"mass must be positive, got {mass_mev}")
  gamma = 1.0 + kinetic_energy_mev / mass_mev
  beta2 = 1.0 - 1.0 / (gamma * gamma)
  # Guard against tiny negative beta2 from round-off.
  beta = math.sqrt(max(beta2, 0.0))
  return beta, gamma


def _t_max_mev(beta: float, gamma: float, mass_mev: float) -> float:
  """Maximum kinetic energy transferable to an electron in one collision."""
  me = ELECTRON_MASS_MEV
  ratio = me / mass_mev
  numer = 2.0 * me * (beta * gamma) ** 2
  denom = 1.0 + 2.0 * gamma * ratio + ratio * ratio
  return numer / denom


def bethe_mass_stopping_power(
  particle: str | Particle,
  kinetic_energy_mev: float,
  material: str | Material,
) -> float:
  """Mass stopping power -dE/(rho dx) [MeV cm^2 / g] (basic Bethe).

  Parameters
  ----------
  particle : str or Particle
    Projectile. Must be a heavy charged particle (not electron / positron).
  kinetic_energy_mev : float
    Projectile kinetic energy [MeV].
  material : str or Material
    Stopping medium.

  Returns
  -------
  float
    Mass stopping power [MeV cm^2 / g]. Always positive (the minus sign
    in -dE/dx is absorbed; we return the energy *loss* rate per unit
    grammage).

  Warns
  -----
  UserWarning
    If beta*gamma is below the threshold where the basic Bethe formula
    is reliable, or if the logarithm argument is non-positive.
  """
  p = get_particle(particle)
  m = get_material(material)

  if p.name in ("electron", "positron"):
    raise NotImplementedError(
      "Bethe formula for electrons/positrons is not implemented in v0.1; "
      "the basic heavy-particle formula does not apply."
    )

  beta, gamma = _beta_gamma_from_kinetic(kinetic_energy_mev, p.mass_mev)
  beta_gamma = beta * gamma

  if beta_gamma < _LOW_BETA_GAMMA_WARN:
    warnings.warn(
      f"beta*gamma = {beta_gamma:.4g} is below {_LOW_BETA_GAMMA_WARN}; "
      "basic Bethe formula is unreliable in this regime. Use a "
      "table-based approach (SRIM/ATIMA/NIST PSTAR) instead.",
      UserWarning,
      stacklevel=2,
    )

  beta2 = beta * beta
  if beta2 <= 0.0:
    # Particle is at rest; infinite/ill-defined. Return inf to signal
    # the integrator that it should stop.
    return float("inf")

  t_max = _t_max_mev(beta, gamma, p.mass_mev)
  # I in MeV (stored in eV).
  i_mev = m.mean_excitation_energy_ev * 1.0e-6
  arg = 2.0 * ELECTRON_MASS_MEV * beta2 * gamma * gamma * t_max / (i_mev * i_mev)
  if arg <= 0.0:
    warnings.warn(
      "Bethe logarithm argument is non-positive; returning 0.",
      UserWarning,
      stacklevel=2,
    )
    return 0.0

  z = float(p.charge)
  prefactor = K_BETHE_MEV_MOL_PER_G * z * z * m.z_over_a / beta2
  bracket = 0.5 * math.log(arg) - beta2
  value = prefactor * bracket
  # In the deep non-relativistic limit the simplified formula can go
  # negative due to the missing shell/Barkas terms; clamp to zero and
  # warn, since negative stopping power is unphysical.
  if value < 0.0:
    warnings.warn(
      "Bethe formula returned a negative value (likely outside its "
      "validity range); clamping to 0.",
      UserWarning,
      stacklevel=2,
    )
    return 0.0
  return value


def bethe_linear_stopping_power(
  particle: str | Particle,
  kinetic_energy_mev: float,
  material: str | Material,
) -> float:
  """Linear stopping power -dE/dx [MeV / cm] (basic Bethe).

  Equivalent to :func:`bethe_mass_stopping_power` multiplied by the
  material density.
  """
  m = get_material(material)
  s_mass = bethe_mass_stopping_power(particle, kinetic_energy_mev, m)
  return s_mass * m.density_g_per_cm3
