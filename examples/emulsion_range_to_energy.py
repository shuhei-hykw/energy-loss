"""Convert measured emulsion track lengths to kinetic energies.

Usage:
    python examples/emulsion_range_to_energy.py PARTICLE LENGTH [UNIT]

Defaults: UNIT = "um". PARTICLE is one of {proton, alpha}.

Example:
    python examples/emulsion_range_to_energy.py alpha 28
    python examples/emulsion_range_to_energy.py proton 10 um
"""

from __future__ import annotations

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
  sys.path.insert(0, str(_REPO_ROOT))

from energy_loss import (  # noqa: E402
  energy_from_emulsion_range,
  get_emulsion_range_energy,
)


def main(argv: list[str]) -> int:
  if len(argv) < 3 or len(argv) > 4:
    print(__doc__, file=sys.stderr)
    return 2
  particle = argv[1]
  length = float(argv[2])
  unit = argv[3] if len(argv) == 4 else "um"
  t_mev = energy_from_emulsion_range(particle, length, unit)
  table = get_emulsion_range_energy(particle)
  meta = table.metadata
  print(f"{particle} track {length} {unit} in nuclear emulsion:")
  print(f"  -> T_kin = {t_mev:.4f} MeV")
  print()
  print(f"  table source : {meta.source}")
  print(
    f"  material     : {meta.material_name}  "
    f"(rho={meta.density_g_per_cm3} g/cm^3, I={meta.mean_excitation_energy_ev} eV)"
  )
  print(f"  fetched      : {meta.fetched}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main(sys.argv))
