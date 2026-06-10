"""Bethe (this library) vs Geant4 backend comparison plots.

Two figures are produced under ``examples/figures/``:

1. ``bethe_vs_geant4_e10.png``
   J-PARC E10 setup: mass stopping power vs kinetic energy for pi-
   and K- on 9Be. The basic Bethe formula (energy_loss.stopping) is
   overlaid with the Geant4-derived dE/dx (computed numerically from
   the Geant4 CSDA range table by finite differences), with markers
   at the two beam working points 1.2 GeV/c and 1.5 GeV/c.

2. ``bethe_vs_geant4_emulsion.png``
   Nuclear emulsion: CSDA range vs kinetic energy for proton and
   alpha. The library's basic-Bethe CSDA (numerically integrated) is
   compared against the Geant4 backend. The deviation at low energy
   is the same physics behind preferring tabulated stopping powers
   in the emulsion range-energy workflow.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
  sys.path.insert(0, str(_REPO_ROOT))

from energy_loss import (  # noqa: E402
  load_config,
  load_table,
  propagate_config,
)
from energy_loss.stopping import bethe_mass_stopping_power  # noqa: E402

FIG_DIR = _REPO_ROOT / "examples" / "figures"
CFG_DIR = _REPO_ROOT / "examples" / "configs"


def _bethe_curve(particle: str, material: str, t_grid: np.ndarray) -> np.ndarray:
  import warnings as _w

  with _w.catch_warnings():
    _w.simplefilter("ignore", category=UserWarning)
    return np.array(
      [bethe_mass_stopping_power(particle, float(t), material) for t in t_grid]
    )


def _bethe_csda_grammage(
  particle: str, material: str, t_grid: np.ndarray,
  t_min_mev: float = 1.0,
) -> np.ndarray:
  """Cumulative CSDA range [g/cm^2] vs T grid using basic Bethe.

  Uses trapezoidal integration of ``1 / S(T)`` on a fine log grid
  between ``t_min_mev`` and each grid point. ``t_min_mev`` should sit
  above the basic-Bethe validity edge; below it the formula clamps to
  zero and the integral becomes meaningless.
  """
  import warnings as _w

  # Build a dense log grid covering [t_min, max(t_grid)] once and
  # compute the cumulative integral; sample at t_grid via interp.
  hi = float(t_grid.max())
  dense = np.geomspace(t_min_mev, hi, num=2000)
  with _w.catch_warnings():
    _w.simplefilter("ignore", category=UserWarning)
    inv_s = np.array(
      [1.0 / bethe_mass_stopping_power(particle, float(t), material)
       for t in dense]
    )
  # Cumulative trapezoidal.
  cum = np.concatenate(([0.0], np.cumsum(0.5 * (inv_s[:-1] + inv_s[1:]) * np.diff(dense))))
  # Anything below t_min has no Bethe support — fill with NaN.
  out = np.full_like(t_grid, np.nan, dtype=float)
  mask = t_grid >= t_min_mev
  out[mask] = np.interp(np.log(t_grid[mask]), np.log(dense), cum)
  return out


def _g4_dedx_from_table(particle: str, material: str) -> tuple[np.ndarray, np.ndarray]:
  """Return (T_MeV, dE/dx [MeV cm^2/g]) from the Geant4 backend CSV.

  Re-runs / re-uses the same generator that ``load_table`` would have
  invoked but reads dE/dx (column 3) directly from the CSV body.
  """
  from energy_loss.backends.geant4_runner import (
    Geant4TableSpec,
    generate_geant4_table,
  )

  spec = Geant4TableSpec(
    particle=particle, material=material,
    emin_mev=0.001, emax_mev=1.0e4, n_points=200,
  )
  csv_path = generate_geant4_table(spec)
  rows: list[list[float]] = []
  with open(csv_path, encoding="utf-8") as f:
    for line in f:
      if line.startswith("#") or "," not in line:
        continue
      try:
        rows.append([float(v) for v in line.split(",")])
      except ValueError:
        continue
  data = np.asarray(rows, dtype=float)
  return data[:, 0], data[:, 3]


def plot_e10_comparison() -> Path:
  fig, ax = plt.subplots(figsize=(8.5, 5.5))
  for cfg_name, particle, color, marker in [
    ("jparc_e10_pim.yaml", "pion-", "C0", "o"),
    ("jparc_e10_km.yaml", "kaon-", "C1", "s"),
  ]:
    # Curves: dE/dx vs T for the species in 9Be.
    t_grid = np.geomspace(1.0, 1.0e4, num=400)
    s_bethe = _bethe_curve(particle, "Be", t_grid)
    g4_t, g4_s = _g4_dedx_from_table(particle, "Be")
    nice = {"pion-": r"$\pi^{-}$", "kaon-": r"$K^{-}$"}[particle]
    ax.loglog(t_grid, s_bethe, "--", color=color, alpha=0.85,
              label=f"{nice} basic Bethe (this lib)")
    ax.loglog(g4_t, g4_s, "-", color=color,
              label=f"{nice} Geant4 11.4.1")
    # E10 working point marker — read T from the YAML.
    cfg = load_config(CFG_DIR / cfg_name)
    t_op = cfg.beam.kinetic_energy_mev
    s_op_g4 = float(np.interp(np.log(t_op), np.log(g4_t), g4_s))
    ax.loglog(
      t_op, s_op_g4, marker, markersize=11, markeredgewidth=1.5,
      markeredgecolor=color, markerfacecolor="white",
    )
  ax.set_xlabel("kinetic energy [MeV]")
  ax.set_ylabel(r"mass stopping power $-dE/(\rho\,dx)$ [MeV cm$^2$/g]")
  ax.set_title(
    r"J-PARC E10 on 9Be: basic Bethe vs Geant4  ($\pi^{-}$ 1.2 GeV/c, "
    r"$K^{-}$ 1.5 GeV/c)"
  )
  ax.grid(True, which="both", linestyle=":", alpha=0.6)
  ax.legend(fontsize=9, loc="upper right")

  # Annotate the ΔE through 3.5 g/cm^2 9Be using the propagate path.
  lines = []
  for cfg_name in ("jparc_e10_pim.yaml", "jparc_e10_km.yaml"):
    cfg = load_config(CFG_DIR / cfg_name)
    r = propagate_config(cfg)
    lines.append(
      f"{cfg.beam.particle.name} @ {cfg.beam.momentum_mev_c/1000:.2f} GeV/c: "
      f"Bethe $\\Delta E$={r.total_energy_loss_mev:.3f} MeV"
    )
  ax.text(
    0.02, 0.04, "\n".join(lines), transform=ax.transAxes,
    fontsize=9, ha="left", va="bottom",
    bbox={"boxstyle": "round,pad=0.4", "fc": "white",
          "ec": "0.6", "lw": 0.8, "alpha": 0.9},
  )

  fig.tight_layout()
  out = FIG_DIR / "bethe_vs_geant4_e10.png"
  fig.savefig(out, dpi=150)
  plt.close(fig)
  return out


def plot_emulsion_comparison() -> Path:
  fig, ax = plt.subplots(figsize=(8.5, 5.5))
  rho_emul = 3.815
  for particle, color in [("proton", "C0"), ("alpha", "C1")]:
    g4 = load_table(particle, "nuclear_emulsion", "geant4")
    t_g4 = g4.kinetic_energy_mev
    r_g4 = g4.range_csda_g_per_cm2
    ax.loglog(
      t_g4, r_g4 / rho_emul * 1.0e4, "-", color=color,
      label=f"{particle} Geant4 11.4.1",
    )
    # Bethe CSDA over the same grid (start at 1 MeV).
    r_bethe = _bethe_csda_grammage(particle, "nuclear_emulsion", t_g4)
    mask = np.isfinite(r_bethe)
    ax.loglog(
      t_g4[mask], r_bethe[mask] / rho_emul * 1.0e4, "--", color=color,
      alpha=0.85, label=f"{particle} basic Bethe CSDA (this lib)",
    )

  ax.set_xlabel("kinetic energy [MeV]")
  ax.set_ylabel(r"range in nuclear emulsion [$\mu$m]")
  ax.set_title("Nuclear emulsion: basic Bethe CSDA vs Geant4")
  ax.grid(True, which="both", linestyle=":", alpha=0.6)
  ax.legend(fontsize=9, loc="upper left")
  ax.set_xlim(1.0e-2, 1.0e4)
  fig.tight_layout()
  out = FIG_DIR / "bethe_vs_geant4_emulsion.png"
  fig.savefig(out, dpi=150)
  plt.close(fig)
  return out


def main() -> int:
  FIG_DIR.mkdir(exist_ok=True)
  paths = [plot_e10_comparison(), plot_emulsion_comparison()]
  print("Wrote:")
  for p in paths:
    print(f"  {p.relative_to(_REPO_ROOT)}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
