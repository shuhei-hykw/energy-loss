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

## v0.1 scope

- Heavy-charged-particle Bethe formula (no density / shell / Barkas /
  Bloch corrections yet).
- YAML-driven material registry (see `energy_loss/data/materials.yaml`).
- Pure elements: Z, A and standard density come from `periodictable`;
  isotope-specific entries (e.g. `9Be`) are supported via `isotope: <N>`.
  Only the empirical mean excitation energy `I` lives in YAML.
- Particles: proton, π±, K±, μ±, e±. Masses for proton/muon/electron
  come from `scipy.constants` (CODATA); π and K from PDG 2024.
- **Single-file setup YAML** (`load_config(path)`) capturing beam
  (particle + kinetic energy *or* momentum) and target (material +
  optional linear thickness *or* grammage). Optional inline
  `particles:` / `materials:` blocks register custom entries.
- Stopping power is evaluated at the beam's *initial* kinetic energy;
  no step-wise integration of `dT/dx = -S(T(x))` yet (planned for
  v0.2, in the `transport` submodule).
- `transport` / `range` / `emulsion` submodules are **not** in v0.1.

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

### Verification plot

```bash
python examples/plot_verification.py
```

writes `examples/figures/{stopping_power_curves,jparc_e10_marker}.png`
showing the Bethe curve on 9Be for proton / π⁻ / K⁻ and the J-PARC E10
working points overlaid.

## Notes on numerical values

No physical constant is hand-typed where a scientific library provides
it: `scipy.constants` is the source for `mₑc²`, the classical electron
radius, Avogadro's number, and the proton/muon masses; `periodictable`
is the source for elemental Z, A, and density. The Bethe prefactor
`K = 4π N_A r_e² mₑc²` is computed from those, not stored. Empirical
mean excitation energies (no library exposes them) and PDG masses for
π / K are the only literal numbers in the package.
