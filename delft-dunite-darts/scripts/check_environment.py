"""Check the actual native API and configuration, not just package metadata."""
import hashlib
import importlib.metadata as metadata
import json
import platform
from pathlib import Path
import sys
import darts.engines
from darts.engines import value_vector, well_control_iface
from darts.nonlinear_solvers import NewtonSolver
from hydrogen_research.darts_source.model import Model
from hydrogen_research.darts_source.properties import load_config, solubility_x


def main():
    versions = {name: metadata.version(name) for name in ['open-darts', 'open-darts-flash', 'numpy', 'scipy', 'pandas', 'matplotlib', 'h5py', 'iapws']}
    if versions['open-darts'] != '2.0.0':
        raise RuntimeError('This model targets Open-DARTS 2.0.0; see docs/installation.md for the exact source revision.')
    config = load_config()
    engine = Path(darts.engines.__file__)
    report = {'python': platform.python_version(), 'platform': platform.platform(),
              'packages': versions, 'engine_sha256': hashlib.sha256(engine.read_bytes()).hexdigest(),
              'default_solubility_salt_free_x': float(solubility_x(config['producer_bhp_bar'], config['temperature_k'], config['salinity_molal'])),
              'status': 'imports and configuration passed; run check and pair for native solve verification'}
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
