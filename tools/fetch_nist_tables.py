"""Fetch NIST PSTAR/ASTAR tables and convert them into the package CSV format.

The script POSTs to the NIST STAR form and writes one CSV per
``(program, material)`` pair into ``energy_loss/data/nist/``. Each CSV
carries a metadata header (source URL, density, mean excitation
energy, composition, fetch date) so the data version is reproducible.

By default it regenerates the bundled v0.3 tables (Photographic
Emulsion, matno=215, for both PSTAR and ASTAR). To add a new material,
extend ``_DEFAULT_JOBS`` below or invoke the module functions directly.

This script is *not* required to use the package; ``pip install`` users
get the bundled CSVs untouched. It lives in the repo so the bundled
data's provenance is fully reproducible.

Sources (public domain, NIST PML)
- PSTAR: https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html
- ASTAR: https://physics.nist.gov/PhysRefData/Star/Text/ASTAR.html

ICRU Reports 49 and 73 underlie PSTAR and ASTAR respectively.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import re
import subprocess
import sys
from pathlib import Path

_FORM_URL = "https://physics.nist.gov/cgi-bin/Star/ap_table.pl"
_COMPOS_URL = "https://physics.nist.gov/cgi-bin/Star/compos.pl"

# Mapping of NIST element Z -> symbol for the composition table.
_Z_TO_SYMBOL: dict[int, str] = {
  1: "H", 6: "C", 7: "N", 8: "O", 16: "S", 35: "Br", 47: "Ag", 53: "I",
}


def _post_table(prog: str, matno: int) -> str:
  """Submit the NIST PSTAR/ASTAR form via curl and return the HTML."""
  cmd = [
    "curl", "-s", "-L", "-X", "POST", _FORM_URL,
    "-F", f"prog={prog}",
    "-F", f"matno={matno:03d}",
    "-F", "GraphType=None",
    "-F", "ShowDefault=on",
    "-F", "Energies=",
    "-F", "userfile=",
  ]
  out = subprocess.run(cmd, check=True, capture_output=True, timeout=60)
  return out.stdout.decode("utf-8", errors="replace")


def _parse_table_rows(html: str) -> list[tuple[float, ...]]:
  """Pull the numeric rows out of the NIST results page."""
  pattern = re.compile(
    r'<tr BGColor="[A-Fa-f0-9]+">\s*'
    r'((?:<td[^>]*>[^<]+</td>\s*)+)'
  )
  cell = re.compile(r'<td[^>]*>([^<]+)</td>')
  rows: list[tuple[float, ...]] = []
  for m in pattern.finditer(html):
    values = [c.strip() for c in cell.findall(m.group(1))]
    if len(values) < 7:
      continue
    try:
      rows.append(tuple(float(v) for v in values[:7]))
    except ValueError:
      continue
  return rows


def _fetch_composition(matno: int) -> tuple[float, float, dict[str, float]]:
  """Density [g/cm^3], mean excitation energy I [eV], weight fractions."""
  url = f"{_COMPOS_URL}?ap{matno:d}"
  out = subprocess.run(
    ["curl", "-s", "-L", url], check=True, capture_output=True, timeout=60
  )
  html = out.stdout.decode("utf-8", errors="replace")
  density = float(
    re.search(r"Density.*?>([0-9.+\-Ee]+)</td>", html, re.S).group(1)
  )
  i_ev = float(
    re.search(r"Mean Excitation Energy.*?>([0-9.+\-Ee]+)</td>",
              html, re.S).group(1)
  )
  # Composition rows: <td>Z</td><td>fraction</td>
  pairs = re.findall(
    r"<td>\s*([0-9]+)\s*</td>\s*<td>\s*([0-9.+\-Ee]+)\s*</td>", html
  )
  composition = {
    _Z_TO_SYMBOL.get(int(z), f"Z{z}"): float(frac)
    for z, frac in pairs if _Z_TO_SYMBOL.get(int(z)) is not None
  }
  return density, i_ev, composition


def _write_csv(
  out_path: Path, prog: str, matno: int, particle: str,
  density: float, i_ev: float, composition: dict[str, float],
  rows: list[tuple[float, ...]],
) -> None:
  today = _dt.date.today().isoformat()
  header_lines = [
    f"# NIST {prog} table for Photographic Emulsion (matno={matno})",
    f"# Particle: {particle}",
    f"# Source:   https://physics.nist.gov/PhysRefData/Star/Text/{prog}.html",
    "# Underlying physics: ICRU Report 49 (proton) / ICRU 49 + 73 (alpha)",
    f"# Density [g/cm^3]: {density:.5f}",
    f"# Mean excitation energy I [eV]: {i_ev:.3f}",
    "# Composition (weight fraction): "
    + ", ".join(f"{k}={v:.6f}" for k, v in composition.items()),
    f"# Fetched: {today} via tools/fetch_nist_tables.py",
    "# Columns: T_MeV, S_elec_MeV_cm2_per_g, S_nuc_MeV_cm2_per_g, "
    "S_total_MeV_cm2_per_g, R_csda_g_per_cm2, R_proj_g_per_cm2, detour_factor",
    "T_MeV,S_elec_MeV_cm2_per_g,S_nuc_MeV_cm2_per_g,"
    "S_total_MeV_cm2_per_g,R_csda_g_per_cm2,R_proj_g_per_cm2,detour_factor",
  ]
  body = "\n".join(",".join(f"{v:.6e}" for v in r) for r in rows)
  out_path.write_text("\n".join(header_lines) + "\n" + body + "\n")


# One row of _DEFAULT_JOBS = (prog, matno, particle, material_label, output_filename).
# Extend this to add new materials when going past v0.3.
_DEFAULT_JOBS: list[tuple[str, int, str, str, str]] = [
  ("PSTAR", 215, "proton", "Photographic Emulsion",
   "pstar_photographic_emulsion.csv"),
  ("ASTAR", 215, "alpha", "Photographic Emulsion",
   "astar_photographic_emulsion.csv"),
]


def main() -> int:
  ap = argparse.ArgumentParser(description=__doc__)
  ap.add_argument(
    "--out-dir",
    default=str(Path(__file__).resolve().parent.parent
                / "energy_loss" / "data" / "nist"),
  )
  args = ap.parse_args()
  out_dir = Path(args.out_dir)
  out_dir.mkdir(parents=True, exist_ok=True)

  # Group jobs by matno so we only fetch the composition once per material.
  matno_to_meta: dict[int, tuple[float, float, dict[str, float]]] = {}
  for _, matno, _, _, _ in _DEFAULT_JOBS:
    if matno not in matno_to_meta:
      matno_to_meta[matno] = _fetch_composition(matno)

  for prog, matno, particle, _material_label, fname in _DEFAULT_JOBS:
    density, i_ev, composition = matno_to_meta[matno]
    print(f"Fetching {prog} for matno={matno} ({particle}) ...", file=sys.stderr)
    html = _post_table(prog, matno)
    rows = _parse_table_rows(html)
    if not rows:
      print(f"  no rows extracted; aborting {prog}.", file=sys.stderr)
      return 1
    out = out_dir / fname
    _write_csv(out, prog, matno, particle, density, i_ev, composition, rows)
    print(f"  wrote {out} ({len(rows)} rows)", file=sys.stderr)
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
