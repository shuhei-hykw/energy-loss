"""Comparison plots for the stopping-power / range models in this package.

Figures produced under ``examples/figures/``:

1. ``bethe_vs_geant4_e10.png``
   J-PARC E10 setup: mass stopping power vs kinetic energy for pi-
   and K- on 9Be (basic Bethe vs Geant4 option4). Markers at the
   1.2 GeV/c pi- and 1.5 GeV/c K- working points.

2. ``bethe_vs_geant4_emulsion.png``
   Nuclear emulsion CSDA range for proton and alpha — basic Bethe
   vs Geant4 option4.

3. ``models_emulsion_low_energy.png`` (v0.5)
   Four-model emulsion comparison zoomed into the regime where
   ATIMA differs most from PSTAR/option4: NIST PSTAR/ASTAR, basic
   Bethe CSDA, Geant4 option4 and Geant4 ATIMA, for proton and
   alpha. Ratios to NIST are annotated for the relevant track
   lengths.
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
  get_emulsion_range_energy,
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


def plot_emulsion_four_models() -> Path:
  """4-model emulsion comparison emphasising the ATIMA differences.

  Solid: NIST table (PSTAR/ASTAR).
  Dotted: Geant4 option4.
  Dash-dot: Geant4 ATIMA.
  Dashed: basic Bethe CSDA (this library).

  Annotations report the ratio of each model to NIST at a couple of
  emulsion-relevant track lengths.
  """
  fig, ax = plt.subplots(figsize=(8.5, 6.0))
  rho_emul = 3.815

  for particle, color, ref_um in [
    ("proton", "C0", 10.0),
    ("alpha", "C1", 28.0),
  ]:
    nist = get_emulsion_range_energy(particle)
    g4o4 = load_table(particle, "nuclear_emulsion", "geant4")
    g4at = load_table(particle, "nuclear_emulsion", "geant4_atima")
    bethe_t = nist.kinetic_energy_mev[nist.kinetic_energy_mev >= 1.0]

    ax.loglog(
      nist.kinetic_energy_mev, nist.range_csda_g_per_cm2 / rho_emul * 1e4,
      "-", color=color, linewidth=1.8,
      label=f"{particle} NIST "
            f"{'PSTAR' if particle == 'proton' else 'ASTAR'}",
    )
    ax.loglog(
      g4o4.kinetic_energy_mev, g4o4.range_csda_g_per_cm2 / rho_emul * 1e4,
      ":", color=color, linewidth=1.6,
      label=f"{particle} Geant4 option4",
    )
    ax.loglog(
      g4at.kinetic_energy_mev, g4at.range_csda_g_per_cm2 / rho_emul * 1e4,
      "-.", color=color, linewidth=1.6,
      label=f"{particle} Geant4 ATIMA",
    )
    r_bethe = _bethe_csda_grammage(particle, "nuclear_emulsion", bethe_t)
    mask = np.isfinite(r_bethe)
    ax.loglog(
      bethe_t[mask], r_bethe[mask] / rho_emul * 1e4,
      "--", color=color, alpha=0.7, linewidth=1.2,
      label=f"{particle} basic Bethe CSDA",
    )

    # Ratios at a representative track length.
    ref_g = ref_um * 1e-4 * rho_emul
    try:
      t_nist = nist.energy_from_range(ref_g)
      t_o4 = g4o4.energy_from_range(ref_g)
      t_at = g4at.energy_from_range(ref_g)
      ax.text(
        0.62, 0.05 + (0.0 if particle == "proton" else 0.12),
        f"{particle} @ {ref_um} um (T from R):  "
        f"NIST={t_nist:.3f} MeV,  option4={t_o4:.3f} ({100*t_o4/t_nist:.1f}%),  "
        f"ATIMA={t_at:.3f} ({100*t_at/t_nist:.1f}%)",
        transform=ax.transAxes, fontsize=8.5, ha="left", va="bottom",
        color=color,
        bbox={"boxstyle": "round,pad=0.3", "fc": "white",
              "ec": color, "lw": 0.8, "alpha": 0.9},
      )
    except ValueError:
      pass

  ax.set_xlabel("kinetic energy [MeV]")
  ax.set_ylabel(r"range in nuclear emulsion [$\mu$m]")
  ax.set_title(
    "Nuclear emulsion range-energy: NIST / Geant4 option4 / Geant4 ATIMA "
    "/ basic Bethe"
  )
  ax.grid(True, which="both", linestyle=":", alpha=0.6)
  ax.legend(fontsize=8, loc="upper left", ncol=2)
  ax.set_xlim(1.0e-1, 1.0e3)
  ax.set_ylim(1.0e-1, 1.0e5)
  fig.tight_layout()
  out = FIG_DIR / "models_emulsion_low_energy.png"
  fig.savefig(out, dpi=150)
  plt.close(fig)
  return out


def main() -> int:
  FIG_DIR.mkdir(exist_ok=True)
  paths = [
    plot_e10_comparison(),
    plot_emulsion_comparison(),
    plot_emulsion_four_models(),
  ]
  print("Wrote:")
  for p in paths:
    print(f"  {p.relative_to(_REPO_ROOT)}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
