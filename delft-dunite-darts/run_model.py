#!/usr/bin/env python3
"""Portable launcher for the native DARTS model and its output workflow."""
from __future__ import annotations
import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
MODULE = 'hydrogen_research.darts_source'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['doctor', 'check', 'run', 'pair', 'suite', 'validate', 'mechanisms', 'sensitivity', 'export', 'render'])
    parser.add_argument('--preset', choices=['pilot', 'standard', 'detailed'], default='standard')
    parser.add_argument('--config', type=Path, help='JSON overrides; relative paths are relative to the calling directory')
    parser.add_argument('--output', type=Path, help='Output directory; see README for command defaults')
    parser.add_argument('--reference', type=Path, help='Completed baseline pair, for sensitivity')
    parser.add_argument('--case', type=Path, help='Completed native case, for VTK export')
    parser.add_argument('--no-plots', action='store_true')
    parser.add_argument('--pilot-label', action='store_true', help='Label VTK and ParaView output as pilot only')
    parser.add_argument('--pvpython', default=os.environ.get('PVPYTHON'), help='ParaView pvpython executable (or PVPYTHON)')
    args = parser.parse_args()
    if sys.version_info[:2] != (3, 11):
        parser.error('Use the tested Python 3.11 DARTS environment. See docs/installation.md.')
    for key in ('config', 'output', 'reference', 'case'):
        value = getattr(args, key)
        if value is not None:
            setattr(args, key, value.resolve())
    if (args.config or args.preset != 'standard' or args.no_plots) and args.command not in {'run', 'pair', 'suite'}:
        parser.error('--config, --preset and --no-plots apply only to run, pair or suite')
    default = 'results/paraview' if args.command in {'export', 'render'} else ('results/flow_sensitivity' if args.command == 'sensitivity' else 'results/darts')
    output = args.output or ROOT / default
    # Keep generated caches outside the repository and clean them on exit.
    with tempfile.TemporaryDirectory(prefix='delft-dunite-darts-') as temporary:
        env = os.environ.copy()
        env.update(PYTHONDONTWRITEBYTECODE='1', MPLCONFIGDIR=str(Path(temporary) / 'matplotlib'), OMP_NUM_THREADS='2')
        def execute(arguments, interpreter=None):
            subprocess.run([interpreter or sys.executable, '-B', *map(str, arguments)], cwd=ROOT, env=env, check=True)
        if args.command == 'doctor':
            execute(['-m', 'scripts.check_environment'])
        elif args.command == 'check':
            execute(['-m', 'unittest', MODULE + '.validation.test_model', '-v'])
            execute(['-c', 'from hydrogen_research.darts_source.validation.properties_reference import benchmarks; r=benchmarks(); print("Fluid benchmark rows:", len(r["rows"]))'])
        elif args.command in {'run', 'pair', 'suite'}:
            output.mkdir(parents=True, exist_ok=True)
            common = ['-m', MODULE + '.run', '--output', output, '--preset', args.preset]
            if args.no_plots:
                common.append('--no-plots')
            if args.command == 'pair':
                overrides = json.loads(args.config.read_text()) if args.config else {}
                summaries = {}
                for mode in ['closed', 'natural_flow']:
                    cfg = Path(temporary) / (mode + '.json')
                    cfg.write_text(json.dumps(overrides | {'boundary_mode': mode}))
                    execute(common + ['--config', cfg])
                    summaries[mode] = json.loads((output / mode / 'summary.json').read_text())
                (output / 'suite_summary.json').write_text(json.dumps(summaries, indent=2) + '\n')
                a = summaries['closed']['generated_H2_kg']
                comparison = {'scope': 'Conditional paired model calculation, not field validation.',
                              'generation_change_percent': 100 * (summaries['natural_flow']['generated_H2_kg'] / a - 1) if a else None}
                (output / 'pair_comparison.json').write_text(json.dumps(comparison, indent=2) + '\n')
                print(json.dumps(comparison, indent=2))
            else:
                if args.config:
                    common += ['--config', args.config]
                if args.command == 'suite':
                    common.append('--suite')
                execute(common)
        elif args.command == 'validate':
            execute(['-m', MODULE + '.run', '--validate', '--output', output])
        elif args.command == 'mechanisms':
            execute(['-m', MODULE + '.validation.mechanisms', '--output', output])
        elif args.command == 'sensitivity':
            execute(['-m', MODULE + '.flow_sensitivity', '--reference', args.reference or ROOT / 'results/darts', '--output', output])
        elif args.command == 'export':
            command = ['paraview/export_vtk.py', '--case', args.case or ROOT / 'results/darts/natural_flow', '--output', output]
            if args.pilot_label:
                command.append('--pilot')
            execute(command)
            execute(['paraview/verify_export.py', '--output', output])
        elif args.command == 'render':
            if not args.pvpython:
                parser.error('Set PVPYTHON or pass --pvpython /path/to/ParaView/bin/pvpython')
            execute(['paraview/build_view.py', '--output', output], args.pvpython)
            execute(['paraview/build_view.py', '--output', output, '--verify-state'], args.pvpython)


if __name__ == '__main__':
    main()
