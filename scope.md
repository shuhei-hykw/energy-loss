# Scope of `energy-loss`

## Purpose

`energy-loss` is a Python-based project for calculating stopping power, range-energy relations, and range-to-energy conversion for charged particles in materials.

The initial motivation is the analysis of nuclear emulsion data in the J-PARC E07 experiment, where the kinetic energy of charged particles must be determined from their measured ranges. However, the package should be designed as a general-purpose, reusable tool for detector materials, targets, absorbers, and nuclear emulsion studies.

The main use cases are:

* calculate or tabulate stopping power;
* convert kinetic energy to range;
* convert measured range to kinetic energy;
* compare different stopping-power/range models;
* estimate model-dependent systematic uncertainties;
* generate and use range-energy tables from external tools such as NIST PSTAR/ASTAR, SRIM, ATIMA, and Geant4.

## Basic design philosophy

The project should follow a Geant4-like strategy.

It should not rely on a single Bethe-Bloch formula over the full energy range. Instead, the appropriate model or tabulated data should be selected depending on:

* particle species;
* kinetic-energy range;
* material;
* available reference data;
* intended accuracy.

For example:

* NIST PSTAR should be used as a reference for protons when available.
* NIST ASTAR should be used as a reference for alpha particles when available.
* SRIM-based tables should be supported for deuterons, tritons, and heavier charged fragments.
* ATIMA and Geant4-derived tables should be supported as additional models or cross-checks.
* A simple Bethe-Bloch implementation may be included for educational use, quick estimates, or high-energy regions, but it should not be the only model.

The package should provide a common interface for all models.

## Repository layout

The repository should use a Python `src` layout. Python is the main product of this repository. Geant4 C++ code is included in the same top-level repository, but it is separated as an optional table-generation backend.

Expected layout:

```text
energy-loss/
  README.md
  scope.md
  pyproject.toml

  src/
    energy_loss/
      __init__.py
      tables.py
      interpolation.py
      particles.py
      materials.py
      converters.py

      models/
        __init__.py
        nist.py
        srim.py
        atima.py
        geant4.py
        bethe.py

      backends/
        __init__.py
        geant4_runner.py

  data/
    nist/
    srim/
    atima/
    geant4/

  tools/
    convert_nist_tables.py
    convert_srim_tables.py
    inspect_table.py

  geant4/
    README.md
    CMakeLists.txt
    src/
      main.cc
      range_table_generator.cc
    include/
      range_table_generator.hh
    macros/

  examples/
    range_to_energy.py
    energy_to_range.py
    compare_models.py
    generate_geant4_table.py

  docs/
    design.md
    models.md
    nuclear_emulsion.md
```

The top-level `src/` directory is for the Python package.
The `geant4/` directory is an independent CMake-based C++ subproject.

## Python package role

The Python package is the main user-facing component.

It should provide lightweight functions such as:

```python
from energy_loss import energy_from_range, range_from_energy

energy = energy_from_range(
    particle="proton",
    material="nuclear_emulsion",
    range_um=35.2,
    model="nist_pstar",
)

range_um = range_from_energy(
    particle="alpha",
    material="nuclear_emulsion",
    energy_mev=5.486,
    model="nist_astar",
)
```

The Python package should mainly do the following:

* read range-energy tables;
* read stopping-power tables;
* validate table metadata;
* interpolate tables;
* invert range-energy relations;
* compare models;
* provide uncertainties or model differences when possible.

The normal Python package must not require Geant4 at install time.

In other words:

```bash
pip install energy-loss
```

should install the Python package only. It should not build Geant4.

## Geant4 role

Geant4 is an optional backend.

The Geant4 C++ code should be used as a table generator, not as a required runtime dependency of the normal Python functions.

The Geant4 backend should generate stopping-power and range-energy tables for specified particles, materials, and energy ranges.

The basic workflow should be:

```text
Python
  -> calls Geant4 C++ executable through subprocess
  -> Geant4 generates CSV/JSON table
  -> Python reads the generated table
  -> Python uses the same interpolation/inversion routines as other models
```

The first implementation should use a CLI backend.

That means Python should call a compiled Geant4 executable, for example:

```python
from energy_loss.backends.geant4_runner import generate_table

generate_table(
    particle="proton",
    material="nuclear_emulsion",
    energy_min_mev=0.01,
    energy_max_mev=100.0,
    output="data/geant4/geant4_11_3_nuclear_emulsion_proton.csv",
)
```

Internally, this function should use `subprocess` to run the C++ executable built under `geant4/`.

Do not use Geant4 Python bindings in the first implementation.
A direct binding backend, such as one using `geant4-pybind`, may be considered later, but the initial design should use the CLI backend.

## Geant4 dependency policy

The normal Python package must not depend on Geant4.

Geant4 should be optional and external.

The user who wants to generate Geant4 tables should manually build the Geant4 subproject:

```bash
cd geant4
cmake -S . -B build
cmake --build build
```

Then Python may call the built executable.

The Geant4 executable path should be configurable by:

* function argument;
* environment variable;
* configuration file;
* or automatic search under `geant4/build/`.

## Table format

All models should be converted into a common table format.

A CSV table should include metadata comments followed by data columns.

Example:

```text
# model: geant4
# geant4_version: 11.3.0
# physics_list: G4EmStandardPhysics_option4
# material: nuclear_emulsion
# density_g_cm3: ...
# particle: proton
# range_type: csda
# energy_unit: MeV
# range_unit: um
# stopping_power_unit: MeV/um
energy_mev,range_um,stopping_power_mev_per_um
0.010, ...
0.011, ...
```

Metadata is important. It should be preserved when reading, writing, or converting tables.

The table reader should check at least:

* model name;
* particle;
* material;
* density;
* units;
* range type;
* energy grid;
* monotonicity of range vs energy.

## API goals

The high-level API should be simple.

Target examples:

```python
energy = energy_from_range(
    particle="proton",
    material="nuclear_emulsion",
    range_um=35.2,
    model="auto",
)
```

```python
result = compare_energy_from_range(
    particle="proton",
    material="nuclear_emulsion",
    range_um=35.2,
    models=["nist_pstar", "srim", "geant4_11_3"],
)
```

```python
table = load_table(
    model="geant4_11_3",
    particle="proton",
    material="nuclear_emulsion",
)
```

Possible return value for comparison:

```python
{
    "nist_pstar": 2.31,
    "srim": 2.28,
    "geant4_11_3": 2.33,
    "unit": "MeV",
}
```

## Model selection policy

The model `"auto"` should choose a reasonable default.

Initial policy:

* proton: prefer NIST PSTAR if available;
* alpha: prefer NIST ASTAR if available;
* deuteron/triton: prefer SRIM or another prepared table;
* heavier ions: prefer SRIM, ATIMA, or Geant4-derived tables depending on availability;
* unknown case: raise a clear error instead of silently choosing an unreliable model.

The package should make it easy to compare models rather than hide model dependence.

## Interpolation policy

Range-energy relations should be treated as monotonic tables.

The basic conversion should use interpolation, not repeated numerical simulation.

Recommended approach:

* interpolate energy-to-range relation;
* invert the monotonic relation for range-to-energy;
* use log-log interpolation where appropriate;
* avoid extrapolation by default;
* raise a clear error when the requested range or energy is outside the table.

Extrapolation may be allowed only with an explicit option.

## Nuclear emulsion support

The first important material is nuclear emulsion.

The package should support at least:

* NIST photographic emulsion, if using NIST PSTAR/ASTAR;
* user-defined nuclear emulsion composition and density;
* Geant4-defined nuclear emulsion material;
* SRIM-compatible material definition.

The documentation should clearly distinguish:

* NIST photographic emulsion;
* actual E07 emulsion material;
* Geant4 material definition;
* SRIM material definition.

These may not be exactly identical. The difference should be treated as a possible systematic uncertainty.

## Out of scope for the initial version

The initial version should not try to do everything.

Out of scope for the first implementation:

* full Geant4 detector simulation;
* event-by-event particle transport in Python;
* direct Geant4 Python binding;
* GUI;
* automatic SRIM execution;
* Docker-based workflows;
* complex uncertainty propagation;
* detector response simulation;
* emulsion scanning or track reconstruction.

These may be considered later, but the first goal is a reliable range-energy conversion library.

## Development priorities

Suggested implementation order:

1. Create the Python package skeleton using `src/energy_loss/`.
2. Implement table loading and metadata parsing.
3. Implement interpolation and range-energy inversion.
4. Add NIST PSTAR/ASTAR table support.
5. Add SRIM table import support.
6. Add comparison utilities.
7. Add Geant4 table format support.
8. Add the Geant4 C++ table generator under `geant4/`.
9. Add Python CLI backend to call the Geant4 executable.
10. Add examples and documentation.

## Naming conventions

Use `snake_case` for Python variables, function names, and file names.

Use clear names such as:

* `energy_from_range`
* `range_from_energy`
* `load_table`
* `compare_models`
* `generate_table`
* `nuclear_emulsion`

Avoid unclear abbreviations unless they are standard names such as `srim`, `nist`, `atima`, or `geant4`.

## Main principle

The main principle is:

The Python package should be lightweight and easy to install.
Geant4 should be available as an optional CLI backend for generating tables.
All generated or imported tables should be converted to a common format and used through the same Python interpolation and inversion routines.
