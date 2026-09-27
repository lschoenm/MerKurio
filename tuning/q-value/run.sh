#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    cargo run --release --manifest-path "$SCRIPT_DIR/Cargo.toml" -- "$1"
    exit 0
fi

cargo run --release --manifest-path "$SCRIPT_DIR/Cargo.toml" -- "$@"
python3 "$SCRIPT_DIR/analyze.py"

echo "Results written to $SCRIPT_DIR/results"
