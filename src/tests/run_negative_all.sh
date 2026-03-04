#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NEG_DIR="${ROOT_DIR}/negative"
OUT_STL="/tmp/ps-neg.stl"
OUT_LOG="/tmp/ps-neg.out"

shopt -s nullglob
files=("${NEG_DIR}"/*.scad)

if [ "${#files[@]}" -eq 0 ]; then
    echo "No negative tests found in ${NEG_DIR}"
    exit 1
fi

failures=0

for f in "${files[@]}"; do
    echo "== $(basename "$f") =="
    if openscad -o "${OUT_STL}" "$f" >"${OUT_LOG}" 2>&1; then
        echo "FAIL: expected failure but command succeeded"
        failures=$((failures + 1))
        continue
    fi

    if rg -q "ERROR: Assertion" "${OUT_LOG}"; then
        echo "PASS: assertion failure as expected"
    else
        echo "FAIL: did not find assertion failure in output"
        echo "--- output ---"
        cat "${OUT_LOG}"
        echo "--------------"
        failures=$((failures + 1))
    fi
done

if [ "${failures}" -ne 0 ]; then
    echo "Negative tests failed: ${failures}"
    exit 1
fi

echo "All negative tests passed."

