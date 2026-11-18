#!/usr/bin/env bash
# Source the HiDRa/Key4hep environment before using this portable launcher.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
exec python3 "$script_dir/display.py" "$@"
