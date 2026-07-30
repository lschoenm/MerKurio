#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BINARY="$SCRIPT_DIR/target/release/merkurio-algorithm-selection-tuning"

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    python3 "$SCRIPT_DIR/refine.py" --binary "$BINARY" "$1"
    exit 0
fi

cargo build --release --manifest-path "$SCRIPT_DIR/Cargo.toml"
python3 "$SCRIPT_DIR/refine.py" --binary "$BINARY" "$@"
