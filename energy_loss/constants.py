"""Physical constants used by the package.

All values are taken from :mod:`scipy.constants` (CODATA recommended
values) so this file does not contain hand-typed magic numbers for any
quantity that an established scientific library already provides.

Internal units throughout the package:
  - energy  : MeV
  - length  : cm
  - density : g / cm^3
"""

from __future__ import annotations

import math

from scipy import constants as _sc

# Electron rest-mass energy [MeV]. CODATA via scipy.
ELECTRON_MASS_MEV: float = _sc.physical_constants[
  "electron mass energy equivalent in MeV"
][0]

# Classical electron radius. scipy returns metres; convert to cm.
R_E_CM: float = _sc.physical_constants["classical electron radius"][0] * 100.0

# Avogadro number [1/mol]. Defined exactly in SI.
N_A: float = _sc.Avogadro

# Speed of light [cm/s]. Rarely needed directly here.
C_CM_PER_S: float = _sc.c * 100.0

# Bethe prefactor K = 4 pi N_A r_e^2 m_e c^2  [MeV mol / g].
# PDG quotes ~0.307075 MeV mol/g; this is derived, not typed.
K_BETHE_MEV_MOL_PER_G: float = (
  4.0 * math.pi * N_A * (R_E_CM**2) * ELECTRON_MASS_MEV
)
