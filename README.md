# energy-loss

Energy loss, stopping power, and range-energy relations for charged
particles in matter. Small, focused, and easy to test.

## Purpose

Two main use cases:

1. **Spectroscopy experiments** — compute the average energy loss of a
   beam or reaction product as it traverses targets, detector gas,
   windows, scintillators, air, He bags, etc. Multi-layer material
   stacks are supported.
2. **Nuclear emulsion** — convert a measured short range (e.g. 10 µm)
   to a kinetic energy via a tabulated range-energy table. The Bethe
   formula alone is *not* trusted in this short-range regime; tables
   are the primary route.

## Scope

### v0.1

- Heavy-charged-particle Bethe formula (no density / shell / Barkas /
  Bloch corrections yet).
- YAML-driven material registry (see `energy_loss/data/materials.yaml`).
- Pure elements: Z, A and standard density come from `periodictable`;
  isotope-specific entries (e.g. `9Be`) are supported via `isotope: <N>`.
  Only the empirical mean excitation energy `I` lives in YAML.
- Particles: proton, π±, K±, μ±, e±. Masses for proton/muon/electron
  come from `scipy.constants` (CODATA); π and K from PDG 2024.
- Single-file setup YAML (`load_config(path)`) capturing beam
  (particle + kinetic energy *or* momentum) and target (material +
  optional linear thickness *or* grammage). Optional inline
  `particles:` / `materials:` blocks register custom entries.
- Stopping power evaluated at the beam's *initial* kinetic energy
  (single-point); see v0.2 for integrated energy loss.

### v0.2

- **`transport` submodule**: step-wise RK4 integration of
  `dT/d(rho x) = -S(T)` along grammage. Independent of how a thickness
  is expressed (linear vs grammage). Particles that drop below a
  threshold are flagged as `stopped`.
- **Layer stacks** via YAML `target.layers: [...]` — the loader keeps
  back-compat with the v0.1 single-target form.
- `propagate()` / `propagate_config()` return a `PropagationResult`
  with per-layer breakdown, exit energy, total `dE` and `stopped`.

### v0.3

- **`range.RangeEnergyTable`**: generic tabulated range-energy relation
  with log-log interpolation in both directions.
- **`emulsion` submodule** with NIST **PSTAR/ASTAR** tables for
  Photographic Emulsion (matno=215) bundled as CSV under
  `energy_loss/data/nist/`. Each CSV's header records source URL,
  density, mean excitation energy, NIST composition and *fetch date*,
  so the data version is always traceable.
- Built-in particles extended: `alpha`, `deuteron`, `triton`, `helion`
  (masses from `scipy.constants`).

### v0.3.1

- Top-level **`energy_from_range`** / **`range_from_energy`** /
  **`compare_energy_from_range`** / **`load_table`** as named by
  `scope.md`, with explicit `model="auto"` / `"nist_pstar"` /
  `"nist_astar"` selection backed by a particle-material-model
  registry (`energy_loss.models`).
- `model="auto"` resolves to `nist_pstar` for protons and `nist_astar`
  for alphas; unknown particle / material / model triples raise a
  clear error rather than silently falling back.
- Fetcher script renamed `scripts/fetch_nist_emulsion.py` →
  **`tools/fetch_nist_tables.py`**, generalised to drive multiple
  `(program, matno, particle)` jobs out of one composition fetch.
- `from_nist_csv` parser tightened (no chained `.replace().isdigit()`
  heuristics).

### v0.4

- **Geant4 backend** as the third model leg suggested by `scope.md`.
  A small C++ CLI under `geant4/` (built separately with CMake) runs
  `G4EmStandardPhysics_option4` and writes a CSV in the same format
  the NIST loader already parses. The Python wrapper
  (`energy_loss.backends.geant4_runner`) finds the executable, runs
  it via `subprocess`, and caches the CSV under
  `energy_loss/data/geant4/` keyed by a hash of the request.
- `model="geant4"` resolves through the registry to the concrete
  `geant4_11_4_1` entry; `compare_energy_from_range(..., models=[
  "nist_pstar", "geant4"])` performs the head-to-head comparison.
- New verification plot `three_legs_emulsion.png` shows NIST tables,
  basic Bethe, and Geant4 on the same axes for proton and alpha in
  nuclear emulsion (NIST and Geant4 overlap within ~1% above 1 MeV).
- New built-in particles (`alpha`, `deuteron`, `triton`, `helion`)
  are routed to Geant4 names via the runner.

### v0.5

- **Geant4 ATIMA** as a second physics-list option in the same
  backend. The C++ generator accepts `--physics {option4,atima}`;
  the YAML schema gets a matching `physics_list:` field; the
  registry exposes the result as a separate entry
  `geant4_atima_11_4_1`, with `model="geant4_atima"` as the alias.
- A unit bug in Geant4 11.4.1's `G4AtimaEnergyLossModel`
  (`ComputeDEDXPerVolume` interprets a MeV/cm result as MeV/mm) is
  compensated in the C++ wrapper. The remaining ~15–25% difference
  vs PSTAR / option4 in the emulsion regime is real model physics
  (ATIMA's `sezi_dedx_e` branch and the lack of density-effect
  correction in `Bethek_dedx_e`).
- The 4-model comparison shows up clearly in the new
  `models_emulsion_low_energy.png` figure (NIST + Geant4 option4 +
  Geant4 ATIMA + basic Bethe for proton & alpha in nuclear emulsion;
  ratios at canonical track lengths annotated).

### Still to come

- Setup YAML extension to drive `propagate_config` with a tabulated
  stopping-power model (Geant4-backed transport).
- SRIM table readers; a true external ATIMA backend (GSI library
  rather than the in-Geant4 model).
- Adaptive step control in the integrator.
- `src/energy_loss/` layout move suggested by `scope.md`.
- Density-effect / shell / Barkas / Bloch corrections in the analytic
  formula; effective charge; nuclear stopping.

Caveats:

- Basic Bethe is reliable roughly for `βγ ≳ 0.05` (about 1 MeV/nucleon
  for protons). Below that the function emits a warning. Use a table
  (SRIM / ATIMA / NIST PSTAR) for short-range work.
- The bundled `nuclear_emulsion` entry uses nominal effective values.
  For real emulsion range-energy work, calibrate against your own batch.
- SRIM / ATIMA source code is **not** included; those tables are
  external comparison data.

## Requirements

- Python ≥ 3.11
- Runtime: `numpy`, `scipy`, `pyyaml`, `periodictable`
- Dev: `pytest`, `ruff`

Exact constraints are pinned by lower bound in `pyproject.toml`.

## Install (with venv)

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -e ".[dev]"
```

## Quick check

```bash
pytest      # runs the test suite
ruff check  # lint
```

## Usage

### Direct call

```python
from energy_loss import get_material, list_materials
from energy_loss.stopping import bethe_mass_stopping_power

# 100 MeV proton in 9Be [MeV cm^2 / g]
s = bethe_mass_stopping_power("proton", 100.0, "Be")
```

### Setup YAML (recommended)

A single YAML file fully captures a beam + target calculation. For
example, the J-PARC E10 case at `examples/configs/jparc_e10_pim.yaml`:

```yaml
beam:
  particle: pion-
  momentum: 1.2
  momentum_unit: GeV/c

target:
  material: Be
  mass_thickness: 3.5
  mass_thickness_unit: g/cm^2
```

```bash
python examples/spectroscopy_loss.py \
       examples/configs/jparc_e10_pim.yaml \
       examples/configs/jparc_e10_km.yaml
```

Optional sections let a config be self-contained — see
`examples/configs/alpha_in_emulsion.yaml` for an inline `particles:` /
`materials:` example.

### Adding your own targets / particles

Either via the same setup YAML (inline `particles:` / `materials:`
blocks) or via a separate registry file:

```yaml
# my_targets.yaml
materials:
  cu_window:
    element: Cu             # Z, A, density auto-resolved from periodictable
    mean_excitation_energy_ev: 322.0

  my_9be:
    element: Be
    isotope: 9              # 9Be(Enriched) — uses the isotopic mass
    mean_excitation_energy_ev: 63.7
```

```python
from energy_loss import load_materials_from_yaml, get_material
load_materials_from_yaml("my_targets.yaml")
get_material("cu_window")
```

For compounds and mixtures, specify `z_over_a` and `density_g_per_cm3`
directly. `mean_excitation_energy_ev` is always required.

### Multi-layer stacks

```yaml
# examples/configs/jparc_e10_stack.yaml
beam:
  particle: pion-
  momentum: 1.2
  momentum_unit: GeV/c
target:
  layers:
    - material: kapton          # entrance window
      thickness: 50.0
      thickness_unit: um
    - material: air             # gap between window and target
      thickness: 5.0
      thickness_unit: cm
    - material: Be              # 9Be target
      mass_thickness: 3.5
      mass_thickness_unit: g/cm^2
```

`propagate_config(cfg)` integrates the energy loss layer by layer
(RK4 along grammage) and returns per-layer `PropagationResult` info.

### Range ↔ energy (top-level API, v0.3.1)

```python
from energy_loss import energy_from_range, range_from_energy

# Auto model selection per scope.md: proton -> NIST PSTAR, alpha -> NIST ASTAR.
t = energy_from_range("proton", "nuclear_emulsion", range_um=35.2)   # MeV
r = range_from_energy("alpha", "nuclear_emulsion", energy_mev=5.486) # um
```

```python
from energy_loss import compare_energy_from_range, list_models

print(list_models(particle="proton", material="nuclear_emulsion"))
# ['nist_pstar']

compare_energy_from_range(
    "proton", "nuclear_emulsion",
    range_um=35.2,
    models=["nist_pstar"],     # add "srim" / "geant4_*" once registered
)
# {'unit': 'MeV', 'nist_pstar': 1.829}
```

Other inputs / output units are first-class:

```python
energy_from_range("proton", "nuclear_emulsion",
                  range_value=10e-3, range_unit="mm",
                  energy_unit="keV")
```

Bundled tables live at `energy_loss/data/nist/*.csv` with full
provenance headers (source URL, density, I, composition, fetch date).
Regenerate / extend them via

```bash
python tools/fetch_nist_tables.py
```

For experiment-specific emulsions, load your own calibration CSV via
`RangeEnergyTable.from_nist_csv(path)` or `from_arrays(T, R)`.

The v0.3 helpers `energy_from_emulsion_range` /
`get_emulsion_range_energy` are kept as thin shims over the same
registry for back-compat.

### Geant4 backend (v0.4)

```bash
# Build the C++ generator once (requires Geant4 >= 11.4 + cmake + Qt).
cmake -S geant4 -B geant4/build \
      -DGeant4_DIR=$HOME/software/geant4/11.4.1/lib/Geant4-11.4.1
cmake --build geant4/build -j

# Generate a table from a YAML job.
python examples/generate_geant4_table.py examples/configs/geant4_table.yaml
```

```python
from energy_loss import compare_energy_from_range

compare_energy_from_range(
    "proton", "nuclear_emulsion",
    range_um=35.2,
    models=["nist_pstar", "geant4"],
)
# {'unit': 'MeV', 'nist_pstar': 1.829, 'geant4_11_4_1': 1.827}
```

The Geant4 backend is **optional** at install time; everything in the
Python package still runs without it. When the executable is missing,
the registry just doesn't advertise the `geant4_*` model and the
3-leg verification plot quietly drops back to 2 legs.

### Verification plots

```bash
python examples/plot_verification.py
```

writes five figures under `examples/figures/`:

- `stopping_power_curves.png` — Bethe curve on 9Be for p / π⁻ / K⁻.
- `jparc_e10_marker.png` — J-PARC E10 working points overlaid.
- `transport_vs_linear.png` — RK4-integrated vs single-point linear
  ΔE through 9Be vs grammage.
- `pstar_vs_bethe_emulsion.png` — NIST PSTAR/ASTAR tables compared
  with this library's basic Bethe CSDA for protons and alphas in
  nuclear emulsion. They agree at high energy and diverge at low
  energy, which is exactly where PSTAR/ASTAR (incorporating shell
  and Barkas corrections) should be used instead of plain Bethe.
- `three_legs_emulsion.png` *(v0.4, requires Geant4 executable)* —
  NIST, Bethe and Geant4 on the same axes. NIST and Geant4 overlap
  within ~1% above 1 MeV; Bethe-only diverges at low energy as
  expected. This is the `scope.md` 3-leg comparison realised.

## Notes on numerical values

No physical constant is hand-typed where a scientific library provides
it: `scipy.constants` is the source for `mₑc²`, the classical electron
radius, Avogadro's number, and the proton/muon masses; `periodictable`
is the source for elemental Z, A, and density. The Bethe prefactor
`K = 4π N_A r_e² mₑc²` is computed from those, not stored. Empirical
mean excitation energies (no library exposes them) and PDG masses for
π / K are the only literal numbers in the package.
