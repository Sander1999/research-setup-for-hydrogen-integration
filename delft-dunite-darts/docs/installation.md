# Python and native DARTS setup

## Reuse a compatible environment first

The tested runtime is Python 3.11.16 on macOS arm64, Open-DARTS 2.0.0 from commit `bb040aa025bfcafff20578654ad0b0d819d463a0`, with CPU/OpenMP engines. The exact checked package versions and engine hash are in `environment-reference.json`. `requirements.txt` includes direct numerical dependencies, including IAPWS. It is not a cross-platform lock of every transitive library.

Select a separate existing DARTS environment; do not install into system Python or Conda base and do not upgrade a shared environment without checking its other consumers. Nothing in this repository creates an environment or downloads an engine automatically. Set `DARTS_PYTHON` to that interpreter, then run:

```bash
"$DARTS_PYTHON" -m pip check
bash run_model.sh doctor
bash run_model.sh check
```

`doctor` checks the actual imports required by this model and reports the native engine hash. `check` solves small native cases, tests conservation, rollback, geometry and permeability updates, and runs the fluid benchmarks. Neither replaces the paired grid/timestep check.

## Set up another machine

If no existing environment is compatible, create a centrally stored Python 3.11 environment using your existing environment manager. For example, with an already installed Conda distribution on macOS/Linux:

```bash
conda create --prefix "$HOME/.local/share/python-envs/darts-py311" python=3.11 pip
export DARTS_PYTHON="$HOME/.local/share/python-envs/darts-py311/bin/python"
"$DARTS_PYTHON" -m pip install -r requirements.txt
"$DARTS_PYTHON" -m pip check
```

Only use that create command if the destination does not already exist. Record the shared environment and this repository as a consumer in your local registry. On Windows, choose a suitable central prefix and its `python.exe`.

[Upstream installation instructions](https://open-darts.readthedocs.io/en/latest/getting_started/installation.html) describe the pip route. Wheel availability and API compatibility vary; this project does not claim that every wheel labelled 2.0.0 is equivalent to the tested source revision. If installation or the required imports fail, build the matching engine from source instead of substituting the unrelated time-series package named `darts`.

## Source build provenance

Clone [Open-DARTS](https://gitlab.com/open-darts/open-darts) into a durable location outside the project and cloud sync, check out the recorded revision, and initialise its submodules:

```bash
git clone https://gitlab.com/open-darts/open-darts.git /path/to/native-sources/open-darts
git -C /path/to/native-sources/open-darts checkout bb040aa025bfcafff20578654ad0b0d819d463a0
git -C /path/to/native-sources/open-darts submodule update --init --recursive
```

Follow the build instructions for that checkout and the [upstream build guide](https://gitlab.com/open-darts/open-darts/-/wikis/Build-instructions), using the selected Python interpreter. The original macOS build used GCC 15, HYPRE/SuperLU, CPU/OpenMP (`OPENDARTS_CONFIG=MT`), and disabled PHREEQC, AMGX and cuDSS. DARTS-flash 0.14.0 and OPM grids were built as separate native dependencies. No equilibrium/mineral thermodynamics backend is required by this reduced model.

`build/macos-packaging.patch` records the two local changes to Open-DARTS: loader-relative macOS RPATHs and inclusion of `.dylib` libraries in package data. Apply it only to a clean matching macOS checkout, before compiling:

```bash
git -C /path/to/native-sources/open-darts apply --check /path/to/this-repository/build/macos-packaging.patch
git -C /path/to/native-sources/open-darts apply /path/to/this-repository/build/macos-packaging.patch
```

It does not alter reservoir physics. This folder contains neither native binaries nor a complete native installer; clean native installation on another machine has not been tested. After building, install this repository's remaining pinned requirements into the selected environment, check dependencies, and run the native tests and paired reference calculation. Compare against `reference/standard_pair.json` before interpreting changes. Do not copy compiled engines between operating systems or Python ABIs.

The repository can run directly without installing its own package. For library imports from elsewhere, optional installation is `"$DARTS_PYTHON" -m pip install --no-deps /path/to/this-repository`, after dependencies have been checked. The main scripts and ParaView workflow are intended to run from a clone.
