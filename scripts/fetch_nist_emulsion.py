"""Fetch PSTAR and ASTAR tables for Photographic Emulsion (matno=215)
from NIST and convert them into CSV with a documented header.

This script is run only when the bundled tables need to be regenerated
or extended. The resulting CSV files are committed to
``energy_loss/data/nist/``; ``pip`` users do *not* need to run this
script. The reason we keep it in the repo is to make the provenance of
the bundled data fully reproducible.

Sources (public domain, NIST PML)
- PSTAR: https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html
- ASTAR: https://physics.nist.gov/PhysRefData/Star/Text/ASTAR.html

ICRU Reports 49 and 73 underlie PSTAR and ASTAR respectively. The
"Photographic Emulsion" predefined material is matno=215 in both.
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
    f"# Fetched: {today} via scripts/fetch_nist_emulsion.py",
    "# Columns: T_MeV, S_elec_MeV_cm2_per_g, S_nuc_MeV_cm2_per_g, "
    "S_total_MeV_cm2_per_g, R_csda_g_per_cm2, R_proj_g_per_cm2, detour_factor",
    "T_MeV,S_elec_MeV_cm2_per_g,S_nuc_MeV_cm2_per_g,"
    "S_total_MeV_cm2_per_g,R_csda_g_per_cm2,R_proj_g_per_cm2,detour_factor",
  ]
  body = "\n".join(",".join(f"{v:.6e}" for v in r) for r in rows)
  out_path.write_text("\n".join(header_lines) + "\n" + body + "\n")


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

  matno = 215  # Photographic Emulsion
  density, i_ev, composition = _fetch_composition(matno)
  for prog, particle, fname in [
    ("PSTAR", "proton", "pstar_photographic_emulsion.csv"),
    ("ASTAR", "alpha", "astar_photographic_emulsion.csv"),
  ]:
    print(f"Fetching {prog} for matno={matno} ...", file=sys.stderr)
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
