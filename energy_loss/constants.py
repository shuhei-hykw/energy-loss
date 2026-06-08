"""Physical constants.

All values are in SI-derived "particle physics" units used throughout the
package. Sources are PDG 2024 / CODATA where applicable. Keep this file
as the single source of truth for constants.
"""

from __future__ import annotations

import math

# Electron rest mass energy [MeV]
ELECTRON_MASS_MEV: float = 0.51099895069

# Classical electron radius [cm]
R_E_CM: float = 2.8179403262e-13

# Avogadro number [1/mol]
N_A: float = 6.02214076e23

# Speed of light [cm/s] (rarely needed directly in this package)
C_CM_PER_S: float = 2.99792458e10

# K = 4 * pi * N_A * r_e^2 * m_e * c^2   [MeV mol / g]
# This is the prefactor of the Bethe formula (Z/A times this gives the
# coefficient with the right units). PDG quotes K = 0.307075 MeV mol / g.
K_BETHE_MEV_MOL_PER_G: float = (
  4.0 * math.pi * N_A * (R_E_CM**2) * ELECTRON_MASS_MEV
)
# K_BETHE_MEV_MOL_PER_G ~= 0.307075 MeV mol / g
