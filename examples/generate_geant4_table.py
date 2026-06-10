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
from energy_loss.models.registry import (  # noqa: E402
  _REGISTRY,
  GEANT4_ATIMA_MODEL_NAME,
  GEANT4_MODEL_NAME,
  ModelKey,
  _canonical_material,
  _canonical_particle,
)


def _validity_check(spec: Geant4TableSpec) -> None:
  """Warn if the requested energy span sticks out of the model's band."""
  model = (
    GEANT4_ATIMA_MODEL_NAME if spec.physics_list == "atima"
    else GEANT4_MODEL_NAME
  )
  key = ModelKey(
    particle=_canonical_particle(spec.particle),
    material=_canonical_material(spec.material),
    model=model,
  )
  entry = _REGISTRY.get(key)
  if entry is None:
    return
  lo, hi = entry.valid_min_mev, entry.valid_max_mev
  if spec.emin_mev < lo or spec.emax_mev > hi:
    print(
      f"  NOTE: requested [{spec.emin_mev}, {spec.emax_mev}] MeV "
      f"extends beyond {model}'s recommended band [{lo:g}, {hi:g}] MeV. "
      f"{entry.note}",
      file=sys.stderr,
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
    physics_list=str(section.get("physics_list", "option4")),
  )

  # Warn if the requested energy span pokes outside the physics-list's
  # recommended band (per the model registry advisory).
  _validity_check(spec)
  print(
    f"Running Geant4 generator for {spec.particle} in {spec.material} "
    f"(physics={spec.physics_list})"
  )
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
