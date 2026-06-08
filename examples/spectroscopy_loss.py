"""Run a stopping-power / energy-loss calculation from a YAML config.

Reports the step-wise RK4 integrated energy loss for each layer and the
single-point linear approximation, so the difference between the two is
visible.

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
  load_config,
  propagate_config,
)
from energy_loss.stopping.models import mass_stopping_power  # noqa: E402


def _report(cfg: Config) -> None:
  beam = cfg.beam
  target = cfg.target
  result = propagate_config(cfg)

  print(f"=== {cfg.source_path.name if cfg.source_path else 'config'} ===")
  print(
    f"  beam   : {beam.particle.name}, p = {beam.momentum_mev_c / 1000.0:.4f} GeV/c, "
    f"T = {beam.kinetic_energy_mev:.3f} MeV"
  )
  print(
    f"  target : {len(target.layers)} layer(s), "
    f"total grammage = {target.total_mass_thickness_g_per_cm2:.4f} g/cm^2"
  )
  for i, lr in enumerate(result.per_layer):
    layer = lr.layer
    print(
      f"    [{i}] {layer.material.name:<22s} "
      f"x={layer.thickness_cm:.4e} cm  "
      f"rho*x={layer.mass_thickness_g_per_cm2:.4e} g/cm^2  "
      f"T: {lr.entry_kinetic_energy_mev:.3f} -> {lr.exit_kinetic_energy_mev:.3f}  "
      f"dE={lr.energy_loss_mev:.4f} MeV"
      + ("  [STOPPED]" if lr.stopped_in_layer else "")
    )

  # Per-layer single-point approximation: sum_i S_i(T_in) * (rho*x)_i.
  # Compares well with the integrator only when dE/T_in is small.
  linear = sum(
    mass_stopping_power(beam.particle, beam.kinetic_energy_mev, layer.material)
    * layer.mass_thickness_g_per_cm2
    for layer in target.layers
  )

  total = result.total_energy_loss_mev
  print(f"  integrated total dE  : {total:.4f} MeV "
        f"(T_exit = {result.exit_kinetic_energy_mev:.3f} MeV)")
  print(f"  linear approx (sum S_i(T_in)*rho*x_i): {linear:.4f} MeV "
        f"(diff = {linear - total:+.4f} MeV)")
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
