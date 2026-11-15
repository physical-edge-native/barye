#!/usr/bin/env bash
# Portable regression check used by CI (Linux / macOS / Windows Git Bash).
# Runs without any hardware and without a GPU.
#
#   1. unit tests of the workspace, including the tests that guard the pinned Dyne crates
#      (border-exchange format and deterministic random numbers, compared bit for bit)
#
# Later milestones add checks here (e.g. a split solution must match the single-patch solution
# bit for bit). Every check must run on all three operating systems.
#
# Usage (from the repository root):  bash scripts/ci-check.sh
set -euo pipefail

fail=0

echo "== [1] unit tests =="
if cargo test --release --workspace --locked; then
  echo "unit tests       OK"
else
  echo "unit tests       FAIL"; fail=1
fi

if [ "$fail" -ne 0 ]; then echo "CI CHECK FAILED"; exit 1; fi
echo "ALL CHECKS PASSED"
