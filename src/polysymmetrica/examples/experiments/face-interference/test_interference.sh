#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
cd "$ROOT_DIR"

OPENSCAD_BIN="${OPENSCAD_BIN:-openscad}"
OUT_DIR="${OUT_DIR:-.tmp}"
IMG_SIZE="${IMG_SIZE:-1800,1000}"
SCAD="src/polysymmetrica/examples/experiments/face-interference/test_interference.scad"

mkdir -p "$OUT_DIR"

"$OPENSCAD_BIN" -o "$OUT_DIR/test_interference_face2.png" \
    --imgsize="$IMG_SIZE" --viewall --autocenter \
    "$SCAD"

"$OPENSCAD_BIN" -D 'FACE_IDX=0' \
    -o "$OUT_DIR/test_interference_face0.png" \
    --imgsize="$IMG_SIZE" --viewall --autocenter \
    "$SCAD"

"$OPENSCAD_BIN" -D 'SHOW_CUTTERS=true' \
    -o "$OUT_DIR/test_interference_cutters_face2.png" \
    --imgsize="$IMG_SIZE" --viewall --autocenter \
    "$SCAD"

echo "Wrote:"
echo "  $OUT_DIR/test_interference_face2.png"
echo "  $OUT_DIR/test_interference_face0.png"
echo "  $OUT_DIR/test_interference_cutters_face2.png"
