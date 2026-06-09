"""Convert measured track lengths to kinetic energies via the top-level API.

Usage:
    python examples/emulsion_range_to_energy.py PARTICLE LENGTH [UNIT] [MODEL]

Defaults: UNIT = "um", MODEL = "auto".

Examples:
    python examples/emulsion_range_to_energy.py alpha 28
    python examples/emulsion_range_to_energy.py proton 10 um
    python examples/emulsion_range_to_energy.py proton 10 um nist_pstar
"""

from __future__ import annotations

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
  sys.path.insert(0, str(_REPO_ROOT))

from energy_loss import (  # noqa: E402
  energy_from_range,
  list_models,
  load_table,
)


def main(argv: list[str]) -> int:
  if not 3 <= len(argv) <= 5:
    print(__doc__, file=sys.stderr)
    return 2
  particle = argv[1]
  length = float(argv[2])
  unit = argv[3] if len(argv) >= 4 else "um"
  model = argv[4] if len(argv) >= 5 else "auto"

  t_mev = energy_from_range(
    particle=particle,
    material="nuclear_emulsion",
    range_value=length,
    range_unit=unit,
    model=model,
  )
  table = load_table(particle, "nuclear_emulsion", model)
  meta = table.metadata

  print(f"{particle} track {length} {unit} in nuclear emulsion (model={model}):")
  print(f"  -> T_kin = {t_mev:.4f} MeV")
  print()
  print(f"  table source : {meta.source}")
  print(
    f"  material     : {meta.material_name}  "
    f"(rho={meta.density_g_per_cm3} g/cm^3, I={meta.mean_excitation_energy_ev} eV)"
  )
  print(f"  fetched      : {meta.fetched}")
  print(
    "  registered models for this (particle, nuclear_emulsion): "
    f"{list_models(particle=particle, material='nuclear_emulsion')}"
  )
  return 0


if __name__ == "__main__":
  raise SystemExit(main(sys.argv))
