"""Run a stopping-power calculation from a YAML config.

Usage:
    python examples/spectroscopy_loss.py <config.yaml> [<config.yaml> ...]
"""

from __future__ import annotations

import sys
from pathlib import Path

# Allow running the script even if the package isn't installed (e.g. when
# `pip install -e .` is shadowed by an iCloud-Drive-hidden .pth file).
_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
  sys.path.insert(0, str(_REPO_ROOT))

from energy_loss import (  # noqa: E402
  Config,
  compute_linear_stopping_power,
  compute_mass_stopping_power,
  load_config,
)


def _report(cfg: Config) -> None:
  beam = cfg.beam
  target = cfg.target
  mat = target.material
  s_mass = compute_mass_stopping_power(cfg)
  s_lin = compute_linear_stopping_power(cfg)

  print(f"=== {cfg.source_path.name if cfg.source_path else 'config'} ===")
  print(
    f"  beam   : {beam.particle.name}, p = {beam.momentum_mev_c / 1000.0:.4f} GeV/c, "
    f"T = {beam.kinetic_energy_mev:.2f} MeV"
  )
  print(
    f"  target : {mat.name}, rho = {mat.density_g_per_cm3:.4f} g/cm^3"
  )
  if target.thickness_cm is not None:
    print(
      f"           thickness = {target.thickness_cm:.4f} cm, "
      f"grammage = {target.mass_thickness_g_per_cm2:.4f} g/cm^2"
    )
  print(f"  mass dE/dx   : {s_mass:.4f} MeV cm^2 / g")
  print(f"  linear dE/dx : {s_lin:.4f} MeV / cm")
  if target.mass_thickness_g_per_cm2 is not None:
    de = s_mass * target.mass_thickness_g_per_cm2
    print(
      f"  ~mean dE     : {de:.4f} MeV  (linear approx., grammage * mass dE/dx)"
    )
  print()


def main(argv: list[str]) -> int:
  if len(argv) < 2:
    print(__doc__, file=sys.stderr)
    return 2
  for path in argv[1:]:
    cfg = load_config(Path(path))
    _report(cfg)
  return 0


if __name__ == "__main__":
  raise SystemExit(main(sys.argv))
