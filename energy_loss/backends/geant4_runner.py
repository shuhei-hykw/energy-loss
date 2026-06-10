"""Subprocess-based wrapper around the in-tree Geant4 table generator.

scope.md's principle is that the Python package does not link Geant4 at
runtime: a separate C++ executable is built once under ``geant4/`` and
the Python side talks to it via :mod:`subprocess`. This module finds
that executable, runs it with a parameter set described by
:class:`Geant4TableSpec`, caches the resulting CSV under
``energy_loss/data/geant4/`` and returns the path.

The executable's path is resolved, in order, from:

1. The ``ENERGY_LOSS_G4_EXECUTABLE`` environment variable.
2. An explicit argument passed to :func:`generate_geant4_table`.
3. ``geant4/build/g4_table_generator`` relative to the repo root.
4. ``g4_table_generator`` on ``PATH``.

If none of these resolves to an executable file,
:class:`Geant4Unavailable` is raised — callers (notably the model
registry) treat that as "this leg of the trio is not available right
now" and surface it as a clear message.
"""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


class Geant4Unavailable(RuntimeError):
  """Raised when the Geant4 backend executable cannot be located."""


# Map this library's canonical particle / material names to the
# Geant4 / NIST identifiers the C++ side expects.
_PARTICLE_TO_G4: dict[str, str] = {
  "proton": "proton",
  "alpha": "alpha",
  "deuteron": "deuteron",
  "triton": "triton",
  "helion": "He3",
  "electron": "e-",
  "positron": "e+",
  "muon-": "mu-",
  "muon+": "mu+",
  "pion-": "pi-",
  "pion+": "pi+",
  "kaon-": "kaon-",
  "kaon+": "kaon+",
}

_MATERIAL_TO_G4: dict[str, str] = {
  "nuclear_emulsion": "G4_PHOTO_EMULSION",
  "aluminum": "G4_Al",
  "carbon": "G4_C",
  "Be": "G4_Be",
  "air": "G4_AIR",
  "kapton": "G4_KAPTON",
  "mylar": "G4_MYLAR",
  "LH2": "G4_lH2",
}


_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
_DEFAULT_BUILT = _REPO_ROOT / "geant4" / "build" / "g4_table_generator"
_CACHE_DIR = _REPO_ROOT / "energy_loss" / "data" / "geant4"


@dataclass(frozen=True)
class Geant4TableSpec:
  """Parameters that uniquely determine a generated table file."""

  particle: str
  material: str
  emin_mev: float
  emax_mev: float
  n_points: int
  grid: str = "log"            # "log" or "linear"
  physics_list: str = "option4"  # "option4" (PSTAR-aligned) or "atima"

  def cache_filename(self) -> str:
    """Stable filename derived from the spec — runs are deterministic."""
    key = (
      f"{self.particle}|{self.material}|{self.emin_mev:.6g}|"
      f"{self.emax_mev:.6g}|{self.n_points}|{self.grid}|"
      f"{self.physics_list}"
    )
    digest = hashlib.sha1(key.encode("utf-8")).hexdigest()[:10]
    safe_mat = self.material.replace("/", "_")
    return f"g4_{self.physics_list}_{self.particle}_{safe_mat}_{digest}.csv"


def detect_geant4_executable(explicit: str | Path | None = None) -> Path:
  """Resolve the Geant4 table-generator executable path.

  See module docstring for the resolution order. Raises
  :class:`Geant4Unavailable` if nothing is found.
  """
  candidates: list[Path] = []
  env = os.environ.get("ENERGY_LOSS_G4_EXECUTABLE")
  if env:
    candidates.append(Path(env))
  if explicit:
    candidates.append(Path(explicit))
  candidates.append(_DEFAULT_BUILT)
  on_path = shutil.which("g4_table_generator")
  if on_path:
    candidates.append(Path(on_path))

  for p in candidates:
    if p.is_file() and os.access(p, os.X_OK):
      return p

  raise Geant4Unavailable(
    "Could not find the Geant4 table generator executable. Build it "
    "with `cmake -S geant4 -B geant4/build && cmake --build geant4/build`, "
    "or set ENERGY_LOSS_G4_EXECUTABLE to its full path."
  )


def generate_geant4_table(
  spec: Geant4TableSpec,
  output_path: str | Path | None = None,
  executable: str | Path | None = None,
  overwrite: bool = False,
) -> Path:
  """Run the Geant4 generator and return the CSV path.

  By default the CSV is cached under ``energy_loss/data/geant4/`` using a
  deterministic filename derived from ``spec``. Existing files are
  reused unless ``overwrite=True``.
  """
  exe = detect_geant4_executable(executable)
  particle_g4 = _PARTICLE_TO_G4.get(spec.particle, spec.particle)
  material_g4 = _MATERIAL_TO_G4.get(spec.material, spec.material)

  if output_path is None:
    _CACHE_DIR.mkdir(parents=True, exist_ok=True)
    output_path = _CACHE_DIR / spec.cache_filename()
  output_path = Path(output_path)

  if output_path.exists() and not overwrite:
    return output_path

  cmd = [
    str(exe),
    "--particle", particle_g4,
    "--material", material_g4,
    "--emin", f"{spec.emin_mev}",
    "--emax", f"{spec.emax_mev}",
    "--n", str(int(spec.n_points)),
    "--grid", spec.grid,
    "--physics", spec.physics_list,
    "--output", str(output_path),
  ]
  try:
    subprocess.run(cmd, check=True, capture_output=True, timeout=300)
  except subprocess.CalledProcessError as exc:
    raise RuntimeError(
      f"Geant4 generator failed (exit {exc.returncode}): "
      f"{exc.stderr.decode(errors='replace')}"
    ) from exc
  return output_path
