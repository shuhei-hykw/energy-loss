"""Verification plot for the Bethe stopping-power implementation.

Renders two figures into ``examples/figures/``:

1. ``stopping_power_curves.png``
   Mass stopping power vs kinetic energy for proton / pi- / K- on a
   9Be target. Compares the qualitative shape (1/beta^2 rise at low T,
   broad MIP near beta*gamma ~ 3) against expectation.

2. ``jparc_e10_marker.png``
   Same curves but with the J-PARC E10 beam working points
   (pi- 1.2 GeV/c, K- 1.5 GeV/c) marked, plus the mean energy loss
   through 3.5 g/cm^2 of 9Be annotated.

Run from the repo root:

    python examples/plot_verification.py
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

from energy_loss import get_particle, load_config  # noqa: E402
from energy_loss.config import compute_mass_stopping_power  # noqa: E402
from energy_loss.stopping import bethe_mass_stopping_power  # noqa: E402

FIG_DIR = _REPO_ROOT / "examples" / "figures"
CFG_DIR = _REPO_ROOT / "examples" / "configs"


def _kinetic_energy_grid(mass_mev: float) -> np.ndarray:
  # Logarithmic grid in kinetic energy from ~ 1 MeV/nucleon (Bethe lower
  # edge) up to ~ 100 GeV. The lower edge intentionally enters the
  # warning region so the curve still gets drawn, but we filter NaNs.
  return np.geomspace(1.0, 1.0e5, num=400)


def _curve(particle_name: str, material_name: str) -> tuple[np.ndarray, np.ndarray]:
  import warnings as _w

  p = get_particle(particle_name)
  t = _kinetic_energy_grid(p.mass_mev)
  # Silence the low-beta-gamma warning; the curve intentionally crosses
  # the basic-Bethe validity edge so the reader can see it.
  with _w.catch_warnings():
    _w.simplefilter("ignore", category=UserWarning)
    s = np.array(
      [bethe_mass_stopping_power(p, float(ti), material_name) for ti in t]
    )
  return t, s


def plot_stopping_curves() -> Path:
  fig, ax = plt.subplots(figsize=(7.0, 5.0))
  for name, label, style in [
    ("proton", r"proton", "-"),
    ("pion-", r"$\pi^{-}$", "--"),
    ("kaon-", r"$K^{-}$", "-."),
  ]:
    t, s = _curve(name, "Be")
    ax.loglog(t, s, style, label=label)
  ax.set_xlabel("kinetic energy [MeV]")
  ax.set_ylabel(r"mass stopping power $-dE/(\rho\,dx)$ [MeV cm$^2$/g]")
  ax.set_title("Basic Bethe stopping power on 9Be (Enriched)")
  ax.grid(True, which="both", linestyle=":", alpha=0.6)
  ax.legend()
  fig.tight_layout()
  out = FIG_DIR / "stopping_power_curves.png"
  fig.savefig(out, dpi=150)
  plt.close(fig)
  return out


def _marker_for_config(yaml_path: Path) -> tuple[str, float, float, float]:
  cfg = load_config(yaml_path)
  s_mass = compute_mass_stopping_power(cfg)
  de = s_mass * (cfg.target.mass_thickness_g_per_cm2 or 0.0)
  return (cfg.beam.particle.name, cfg.beam.kinetic_energy_mev, s_mass, de)


def plot_jparc_e10_marker() -> Path:
  fig, ax = plt.subplots(figsize=(8.0, 5.0))
  curves = {
    "pion-": (_curve("pion-", "Be"), r"$\pi^{-}$", "C0"),
    "kaon-": (_curve("kaon-", "Be"), r"$K^{-}$", "C1"),
  }
  for _, ((t, s), label, color) in curves.items():
    ax.loglog(t, s, "-", color=color, label=label)

  marker_specs = [
    ("jparc_e10_pim.yaml", "o", "C0"),
    ("jparc_e10_km.yaml", "s", "C1"),
  ]
  marker_labels: list[str] = []
  for cfg_name, marker, color in marker_specs:
    name, t, s, de = _marker_for_config(CFG_DIR / cfg_name)
    ax.loglog(
      t, s, marker,
      markersize=10, markeredgewidth=1.5,
      markeredgecolor=color, markerfacecolor="white",
    )
    marker_labels.append(
      f"{name} @ J-PARC E10: T={t / 1000.0:.3f} GeV, "
      f"dE/dx={s:.3f} MeV cm$^2$/g, $\\Delta E\\approx${de:.2f} MeV"
    )

  # Single text box in the empty upper-left corner; no leader arrows.
  ax.text(
    0.02, 0.97,
    "\n".join(marker_labels),
    transform=ax.transAxes,
    fontsize=9,
    ha="left", va="top",
    bbox={"boxstyle": "round,pad=0.4", "fc": "white",
          "ec": "0.6", "lw": 0.8, "alpha": 0.9},
  )

  ax.set_xlabel("kinetic energy [MeV]")
  ax.set_ylabel(r"mass stopping power $-dE/(\rho\,dx)$ [MeV cm$^2$/g]")
  ax.set_title(r"J-PARC E10 working points on 9Be, 3.5 g/cm$^2$")
  ax.grid(True, which="both", linestyle=":", alpha=0.6)
  ax.legend(fontsize=10, loc="upper right")
  fig.tight_layout()
  out = FIG_DIR / "jparc_e10_marker.png"
  fig.savefig(out, dpi=150)
  plt.close(fig)
  return out


def main() -> int:
  FIG_DIR.mkdir(exist_ok=True)
  paths = [plot_stopping_curves(), plot_jparc_e10_marker()]
  print("Wrote:")
  for p in paths:
    print(f"  {p.relative_to(_REPO_ROOT)}")
  # Print numeric summary so the test log is self-contained.
  print()
  for cfg in ["jparc_e10_pim.yaml", "jparc_e10_km.yaml"]:
    name, t, s, de = _marker_for_config(CFG_DIR / cfg)
    print(
      f"{name:6s} : T={t / 1000.0:.4f} GeV, dE/dx={s:.4f} MeV cm^2/g, "
      f"<dE>~{de:.4f} MeV through 3.5 g/cm^2 9Be"
    )
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
