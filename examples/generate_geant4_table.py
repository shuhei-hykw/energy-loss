"""Drive the Geant4 table generator from a YAML config.

Usage:
    python examples/generate_geant4_table.py <config.yaml>

The config schema is one top-level mapping ``geant4_table`` carrying:
particle, material, emin_mev, emax_mev, n_points, optional grid.
"""

from __future__ import annotations

import sys
from pathlib import Path

import yaml

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
  sys.path.insert(0, str(_REPO_ROOT))

from energy_loss import RangeEnergyTable  # noqa: E402
from energy_loss.backends import (  # noqa: E402
  Geant4TableSpec,
  generate_geant4_table,
)


def main(argv: list[str]) -> int:
  if len(argv) != 2:
    print(__doc__, file=sys.stderr)
    return 2
  with open(argv[1], encoding="utf-8") as f:
    doc = yaml.safe_load(f) or {}
  try:
    section = doc["geant4_table"]
  except KeyError:
    print("YAML must have a top-level 'geant4_table' mapping.", file=sys.stderr)
    return 2

  spec = Geant4TableSpec(
    particle=str(section["particle"]),
    material=str(section["material"]),
    emin_mev=float(section["emin_mev"]),
    emax_mev=float(section["emax_mev"]),
    n_points=int(section["n_points"]),
    grid=str(section.get("grid", "log")),
  )
  print(f"Running Geant4 generator for {spec.particle} in {spec.material}")
  print(
    f"  energy grid : {spec.n_points} points ({spec.grid}) "
    f"from {spec.emin_mev} to {spec.emax_mev} MeV"
  )
  csv = generate_geant4_table(spec)
  print(f"  -> {csv}")

  table = RangeEnergyTable.from_nist_csv(csv)
  m = table.metadata
  print(
    f"  table : {len(table.kinetic_energy_mev)} rows, "
    f"density={m.density_g_per_cm3} g/cm^3, I={m.mean_excitation_energy_ev} eV"
  )
  print(f"  source : {m.source}")
  print(f"  fetched: {m.fetched}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main(sys.argv))
