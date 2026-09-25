#!/usr/bin/env bash
set -euo pipefail
MODEL_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
DARTS_PYTHON="${DARTS_PYTHON:-${HOME}/.local/share/python-envs/darts-py311/bin/python}"
[[ -x "$DARTS_PYTHON" ]] || { echo 'Set DARTS_PYTHON to an existing compatible Python 3.11 DARTS interpreter; see docs/installation.md.' >&2; exit 2; }
exec "$DARTS_PYTHON" -B "$MODEL_ROOT/run_model.py" "$@"
