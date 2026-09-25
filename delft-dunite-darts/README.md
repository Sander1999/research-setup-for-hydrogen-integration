# Delft–dunite hydrogen model with Open-DARTS

A standalone research model for comparing hydrogen generation, retention and deep transport under initially stationary water and an imposed natural hydraulic gradient. It couples finite reactive iron, water consumption, an empirical expansion/damage law and native three-component, two-phase DARTS flow.

The sedimentary sequence is informed by Delft geology. The ideal dunite and its connecting corridor are **hypotheses**, not observed features beneath Delft. The model supports a proposal combining established research areas. Kinetics and damage are scenario assumptions; the code is not a calibrated resource estimate or a resolved fracture model.

## Start here

Use an existing compatible **Python 3.11 / Open-DARTS 2.0.0** environment. See [installation and native build provenance](docs/installation.md) before installing dependencies. A package version alone does not establish native API compatibility.

On the original Mac, the launcher finds the existing shared `darts-py311` environment through your home directory. On another machine, set its interpreter explicitly:

```bash
export DARTS_PYTHON="/absolute/path/to/compatible/environment/bin/python"
bash run_model.sh doctor
bash run_model.sh check
bash run_model.sh pair --preset pilot --output results/pilot
```

The pilot is a one-year, coarse setup check. It must not be compared with the ten-year reference results. The scientific configuration and convergence requirements are unchanged in the standard and detailed runs.

```bash
# Standard ten-year closed / natural-flow comparison
bash run_model.sh pair
# Full eight-case control suite; replaces matching standard cases
bash run_model.sh suite
# Paired timestep and grid checks; can take substantially longer
bash run_model.sh validate
# Water-renewal mechanism audit, after the pair or suite
bash run_model.sh mechanisms
# Deeper host changed from 1 to 100 mD, all other source inputs retained
bash run_model.sh sensitivity
```

Python-only invocation works on any supported platform without Bash:

```text
/absolute/path/to/python run_model.py pair --preset standard
```

On Windows, use the environment's `python.exe`. The native model was verified locally on macOS arm64; other platforms need their own native tests.

## Change the calculation

The full default input is [config.json](hydrogen_research/darts_source/config.json). Supply a JSON override without editing the defaults:

```bash
bash run_model.sh run --config examples/closed.json --output results/custom
bash run_model.sh pair --preset detailed --output results/detailed
```

Relative input/output paths refer to the directory where you issue the command; omitted outputs go under this repository's `results/`. `pair` overrides only the boundary mode between its two cases. `suite` additionally applies the named control overrides. A `layers` override replaces the complete layer list. Do not mix results from different configurations or presets in one output directory. Matching case files are regenerated when rerun.

The standard mesh has 3,648 cells after formation boundaries are honoured. `detailed` uses 20 × 16 horizontal cells, at least 24 vertical cells and a maximum five-day step. Runtime is determined by convergence and resolution; no automatic reduction of fidelity is applied. `validate` refines the saved pair, or creates the standard pair if none exists. Read `verification.json`: command completion is **not** a claim that every sensitivity target passed.

## 3D results and ParaView

Each plotted run creates `history.png`, `source_rock_3d.png` and `transport_section.png`. ParaView is optional and uses its own `pvpython` runtime:

```bash
bash run_model.sh export
export PVPYTHON="/absolute/path/to/ParaView/bin/pvpython"
bash run_model.sh render
```

For a pilot, add `--case results/pilot/natural_flow --output results/paraview_pilot --pilot-label` to `export`, then use `render --output results/paraview_pilot`. Open the generated `source_model.pvd` in ParaView, or load `hydrogen_research.pvsm`. Rebuild the state after moving the output folder. Export verifies every numeric array against native CSVs; rendering reloads the saved state in a separate process. See [ParaView details](paraview/README.md).

## Outputs and checks

Each case saves complete parameters and hashes, status, accepted-step history, initial/final spatial fields, intermediate snapshots and a summary. `summary.json` separates generation, retained dissolved/free gas and **deep-boundary** export. Export is not surface arrival. The mass ledger accounts for hydration/redox water, salt and hydrogen. Failed calculations retain `status: failed` rather than a completed summary.

[Native equations and assumptions](hydrogen_research/darts_source/README.md) describe the finite source, reaction splitting, damage and TPFA flow. [Geology](docs/delft_geological_scenario.md), [source-rock evidence](docs/source_rock_science.md) and [fluid properties](docs/fluid_sources.md) document the scientific basis. Input evidence tables retain their source references.

The compact `reference/` files record earlier standard calculations and refinement checks. The reference generation difference is approximately +0.337%. Free-gas inventory **did not meet** the 10% grid-sensitivity target (15.7–20.1% discrepancy), so gas accumulation and plume shape need further refinement. None of these results validates the hypothetical geology or reaction rate.

[Package verification](docs/package_verification.md) records what was rerun for this standalone copy. Its native core and input file preserve the original SHA-256 fingerprint.

## Repository contents

- `hydrogen_research/darts_source/`: native model, configuration, controls and numerical tests.
- `run_model.py`, `run_model.sh`: portable launchers.
- `paraview/`: export, verification and rendering scripts.
- `examples/`, `docs/`, `data/literature/`: small required inputs and explanations.
- `reference/`: compact reference summaries; full results are generated locally.
- `requirements.txt`, `environment-reference.json`, `build/`: dependency and native-build records.

The manuscript, personal paths, installed engines, environment folders and generated large results are excluded. `.gitignore` keeps new simulation outputs out of Git. This is a local repository-ready folder; no GitHub repository has been created or uploaded. See [provenance and licensing](NOTICE.md) before distributing under a chosen project licence.
