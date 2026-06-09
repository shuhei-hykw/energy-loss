"""Range-energy tables and their interpolation.

A :class:`RangeEnergyTable` stores parallel arrays of kinetic energy
[MeV] and CSDA range [g/cm^2], plus optional metadata that documents
where the numbers came from. Both directions of the relation
(``range_from_energy`` and ``energy_from_range``) are exposed; the
interpolation is log-log linear because range varies over many decades.

The bundled NIST PSTAR/ASTAR CSV format is parsed by
:meth:`RangeEnergyTable.from_nist_csv`. The header lines (``# ...``) in
those files record the source URL, density, mean excitation energy,
composition and fetch date so the provenance of every bundled table is
explicit.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class RangeEnergyTableMetadata:
  """Provenance / version information for a :class:`RangeEnergyTable`.

  Populated automatically by :meth:`RangeEnergyTable.from_nist_csv`
  from the ``# ...`` header lines of the bundled NIST CSV files.
  """

  source: str = ""
  particle: str = ""
  material_name: str = ""
  density_g_per_cm3: float | None = None
  mean_excitation_energy_ev: float | None = None
  composition: dict[str, float] = field(default_factory=dict)
  fetched: str = ""
  notes: str = ""


@dataclass(frozen=True)
class RangeEnergyTable:
  """A tabulated range-energy relation with log-log linear interpolation.

  ``kinetic_energy_mev`` and ``range_csda_g_per_cm2`` must be sorted by
  kinetic energy and strictly positive.
  """

  kinetic_energy_mev: np.ndarray
  range_csda_g_per_cm2: np.ndarray
  metadata: RangeEnergyTableMetadata = field(
    default_factory=RangeEnergyTableMetadata
  )

  def __post_init__(self) -> None:
    t = np.asarray(self.kinetic_energy_mev, dtype=float)
    r = np.asarray(self.range_csda_g_per_cm2, dtype=float)
    if t.shape != r.shape or t.ndim != 1:
      raise ValueError(
        "kinetic_energy_mev and range_csda_g_per_cm2 must be 1-D arrays "
        f"of the same shape; got {t.shape} vs {r.shape}"
      )
    if np.any(t <= 0.0) or np.any(r <= 0.0):
      raise ValueError("table values must be strictly positive")
    if np.any(np.diff(t) <= 0.0):
      raise ValueError("kinetic_energy_mev must be strictly increasing")
    if np.any(np.diff(r) <= 0.0):
      raise ValueError(
        "range_csda_g_per_cm2 must be strictly increasing with energy"
      )
    # Store as plain float arrays (frozen dataclass — use object.__setattr__).
    object.__setattr__(self, "kinetic_energy_mev", t)
    object.__setattr__(self, "range_csda_g_per_cm2", r)

  # ----- public lookup ------------------------------------------------------

  def range_from_energy(self, kinetic_energy_mev: float) -> float:
    """Return CSDA range [g/cm^2] for kinetic energy ``T`` [MeV]."""
    self._check_in_energy_range(kinetic_energy_mev)
    log_t = np.log(self.kinetic_energy_mev)
    log_r = np.log(self.range_csda_g_per_cm2)
    return float(np.exp(np.interp(np.log(kinetic_energy_mev), log_t, log_r)))

  def energy_from_range(self, range_g_per_cm2: float) -> float:
    """Return kinetic energy [MeV] for CSDA range ``R`` [g/cm^2]."""
    self._check_in_range_range(range_g_per_cm2)
    log_t = np.log(self.kinetic_energy_mev)
    log_r = np.log(self.range_csda_g_per_cm2)
    return float(np.exp(np.interp(np.log(range_g_per_cm2), log_r, log_t)))

  def energy_from_linear_range(
    self, range_cm: float, density_g_per_cm3: float | None = None
  ) -> float:
    """Same as :meth:`energy_from_range` but takes a linear range.

    If ``density_g_per_cm3`` is omitted, the metadata density is used.
    """
    rho = density_g_per_cm3 or self.metadata.density_g_per_cm3
    if rho is None or rho <= 0.0:
      raise ValueError(
        "energy_from_linear_range needs a density (either explicit or "
        "stored in metadata)."
      )
    return self.energy_from_range(range_cm * rho)

  # ----- factories ----------------------------------------------------------

  @staticmethod
  def from_arrays(
    kinetic_energy_mev: np.ndarray,
    range_csda_g_per_cm2: np.ndarray,
    metadata: RangeEnergyTableMetadata | None = None,
  ) -> RangeEnergyTable:
    return RangeEnergyTable(
      kinetic_energy_mev=np.asarray(kinetic_energy_mev, dtype=float),
      range_csda_g_per_cm2=np.asarray(range_csda_g_per_cm2, dtype=float),
      metadata=metadata or RangeEnergyTableMetadata(),
    )

  @staticmethod
  def from_nist_csv(path: str | Path) -> RangeEnergyTable:
    """Parse a NIST PSTAR/ASTAR-format CSV bundled with the package.

    The CSV is the format produced by ``scripts/fetch_nist_emulsion.py``.
    Header lines starting with ``#`` carry the source, density,
    composition and fetch date; the table body is comma-separated.
    """
    path = Path(path)
    meta = _parse_nist_header(path)
    rows: list[list[float]] = []
    with open(path, encoding="utf-8") as f:
      for line in f:
        if line.startswith("#") or "," not in line:
          continue
        try:
          rows.append([float(v) for v in line.split(",")])
        except ValueError:
          # Column-name row — skip it.
          continue
    data = np.asarray(rows, dtype=float)
    if data.ndim != 2 or data.shape[1] < 5:
      raise ValueError(
        f"{path}: expected at least 5 numeric columns "
        "(T, S_elec, S_nuc, S_total, R_csda)."
      )
    return RangeEnergyTable(
      kinetic_energy_mev=data[:, 0],
      range_csda_g_per_cm2=data[:, 4],
      metadata=meta,
    )

  # ----- internals ----------------------------------------------------------

  def _check_in_energy_range(self, t: float) -> None:
    lo, hi = self.kinetic_energy_mev[0], self.kinetic_energy_mev[-1]
    if not (lo <= t <= hi):
      raise ValueError(
        f"kinetic energy {t} MeV is outside the tabulated range "
        f"[{lo}, {hi}] MeV"
      )

  def _check_in_range_range(self, r: float) -> None:
    lo, hi = self.range_csda_g_per_cm2[0], self.range_csda_g_per_cm2[-1]
    if not (lo <= r <= hi):
      raise ValueError(
        f"range {r} g/cm^2 is outside the tabulated range "
        f"[{lo}, {hi}] g/cm^2"
      )


# ---------- NIST CSV header parser -------------------------------------------

_HEADER_PATTERNS: dict[str, re.Pattern[str]] = {
  "source": re.compile(r"#\s*Source:\s*(\S+)"),
  "particle": re.compile(r"#\s*Particle:\s*(\S+)"),
  "material_name": re.compile(r"#\s*NIST\s+\S+\s+table\s+for\s+(.+)\s*$"),
  "density": re.compile(r"#\s*Density.*?:\s*([0-9.+\-Ee]+)"),
  "i_ev": re.compile(r"#\s*Mean excitation energy.*?:\s*([0-9.+\-Ee]+)"),
  "composition": re.compile(r"#\s*Composition.*?:\s*(.+?)\s*$"),
  "fetched": re.compile(r"#\s*Fetched:\s*(\S+)"),
}


def _parse_nist_header(path: Path) -> RangeEnergyTableMetadata:
  source = ""
  particle = ""
  material_name = ""
  density: float | None = None
  i_ev: float | None = None
  composition: dict[str, float] = {}
  fetched = ""
  notes_lines: list[str] = []
  with open(path, encoding="utf-8") as f:
    for line in f:
      if not line.startswith("#"):
        break
      stripped = line.rstrip("\n")
      m = _HEADER_PATTERNS["source"].match(stripped)
      if m:
        source = m.group(1)
        continue
      m = _HEADER_PATTERNS["particle"].match(stripped)
      if m:
        particle = m.group(1)
        continue
      m = _HEADER_PATTERNS["material_name"].match(stripped)
      if m:
        material_name = m.group(1).split("(")[0].strip()
        continue
      m = _HEADER_PATTERNS["density"].match(stripped)
      if m:
        density = float(m.group(1))
        continue
      m = _HEADER_PATTERNS["i_ev"].match(stripped)
      if m:
        i_ev = float(m.group(1))
        continue
      m = _HEADER_PATTERNS["composition"].match(stripped)
      if m:
        for token in m.group(1).split(","):
          if "=" in token:
            k, v = token.split("=", 1)
            try:
              composition[k.strip()] = float(v)
            except ValueError:
              pass
        continue
      m = _HEADER_PATTERNS["fetched"].match(stripped)
      if m:
        fetched = m.group(1)
        continue
      notes_lines.append(stripped.lstrip("# "))
  return RangeEnergyTableMetadata(
    source=source,
    particle=particle,
    material_name=material_name,
    density_g_per_cm3=density,
    mean_excitation_energy_ev=i_ev,
    composition=composition,
    fetched=fetched,
    notes="\n".join(notes_lines),
  )
