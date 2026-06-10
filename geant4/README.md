# geant4/ — table-generator backend

This is the optional Geant4 backend for the `energy_loss` Python
package. It is a small standalone C++ CLI that prints stopping-power
and CSDA-range tables for one `(particle, material)` combination to a
CSV file. The CSV format matches the NIST PSTAR/ASTAR header layout
that `energy_loss.range.RangeEnergyTable.from_nist_csv` already
parses, so the Python side reads Geant4-generated tables through the
same interpolation code path as the bundled NIST data.

The Python package does **not** link Geant4 — it talks to the
executable via `subprocess`. `pip install energy-loss` works without
Geant4; this subdirectory only gets built when you want the 3rd model
leg available.

## Build

Requirements: a Geant4 install (this repo was developed against
11.4.1), CMake ≥ 3.16, a C++17 compiler. With the build tree under
`geant4/build/`:

```bash
cmake -S geant4 -B geant4/build \
      -DGeant4_DIR=$HOME/software/geant4/11.4.1/lib/Geant4-11.4.1
cmake --build geant4/build -j
```

The resulting binary is `geant4/build/g4_table_generator`. The Python
runner finds it via:

1. `ENERGY_LOSS_G4_EXECUTABLE` environment variable
2. an explicit `executable=` argument
3. `geant4/build/g4_table_generator` (this default)
4. `g4_table_generator` on `$PATH`

## CLI

```text
g4_table_generator
  --particle <name>     Geant4 particle name (e.g. "proton", "alpha")
  --material <G4 NIST>  Geant4 NIST material (e.g. "G4_PHOTO_EMULSION")
  --emin <MeV>          lower bound of the kinetic-energy grid
  --emax <MeV>          upper bound
  --n <int>             number of grid points (>= 2)
  --grid <log|linear>   default: log
  --output <path>       CSV path to write
```

The output CSV contains comment-prefixed metadata (Geant4 version,
material density, mean excitation energy, composition, fetch date)
followed by columns `T_MeV, S_elec, S_nuc, S_total, R_csda, R_proj,
detour_factor`. Units mirror the NIST PSTAR format: stopping power
in MeV cm²/g, range in g/cm².

## Physics

The physics list is `G4EmStandardPhysics_option4`, the same EM line
that PSTAR/ASTAR are aligned with (ICRU 49 / 73). `BuildCSDARange` is
enabled before initialisation; when the CSDA table is missing for a
particle/material combination the runner falls back to the standard
energy-loss range table (`G4EmCalculator::GetRange`).
