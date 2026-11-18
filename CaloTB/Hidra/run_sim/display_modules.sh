#!/usr/bin/env bash
# Actual interlocking module outlines only; source the HiDRa setup first.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
exec python3 "$script_dir/display.py" --module-borders "$@"
