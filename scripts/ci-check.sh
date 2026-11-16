#!/usr/bin/env bash
# Portable regression check used by CI (Linux / macOS / Windows Git Bash).
# Runs without any hardware and without a GPU.
#
#   1. unit tests of the workspace, including the tests that guard the pinned Dyne crates
#      (border-exchange format and deterministic random numbers, compared bit for bit), and the
#      2D and 3D wave tests: split solutions equal the single-patch solution bit for bit, and the
#      native 2D solver reproduces the field hashes of Dyne's WebAssembly wave solver
#   2. the 2D wave equation on a small grid, split 1x1, 2x2, 4x4, 8x8 and 3x5 (uneven) with one
#      thread per patch: every split must print "PASS (bit-exact)", and the field hash must be the
#      one Dyne's own wave scenario produced (18e44331)
#   3. the native solver against Dyne's WebAssembly wave solver (wave2d-actor), built here from the
#      pinned Dyne commit (the .wasm file is not kept in the repository): same field bits and same
#      border-message bytes. Needs the Rust target wasm32-unknown-unknown; skipped with a message
#      if it is missing, except in CI (CI=true), where a missing target is a failure.
#   4. the 3D wave equation on a small grid (25x23x21, 100 steps), split 1x1x1, 2x2x2, 3x2x1,
#      4x4x2 and 1x3x5 (uneven) with one thread per patch: every split must print
#      "PASS (bit-exact)", and the field hash must be the pinned one (2ca26bb5), the same on every OS
#   5. convergence order of the 2D and 3D wave solvers against analytic standing waves (fixed
#      Courant number 0.5, errors at t = 0.5 on refined grids): every observed order must lie in
#      [1.9, 2.1] (band and reasons in crates/barye-core/src/convergence.rs)
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

echo "== [2] 2D wave: split = single patch (bit-exact) =="
out=$(cargo run --release --locked --quiet --bin barye-wave2d -- check --grid 41x41 --steps 400 \
        --splits 1x1,2x2,4x4,8x8,3x5) || true
echo "$out"
passes=$(printf '%s\n' "$out" | grep -c "PASS (bit-exact)" || true)
if [ "$passes" -eq 5 ] && printf '%s\n' "$out" | grep -q "single patch      crc32 18e44331"; then
  echo "wave2d split     OK"
else
  echo "wave2d split     FAIL"; fail=1
fi

echo "== [3] 2D wave: native solver = Dyne's WebAssembly solver =="
if rustup target list --installed 2>/dev/null | grep -q '^wasm32-unknown-unknown'; then
  # The pinned Dyne checkout that cargo already fetched for dyne-wire also holds wave2d-actor.
  wire_manifest=$(cargo metadata --format-version 1 --locked \
                  | grep -o '"manifest_path":"[^"]*dyne-wire[^"]*Cargo.toml"' | sed -n 1p \
                  | sed -e 's/^"manifest_path":"//' -e 's/"$//' -e 's#\\\\#/#g')
  actor_manifest="$(dirname "$(dirname "$wire_manifest")")/wave2d-actor/Cargo.toml"
  wasm="$(pwd)/target/dyne-actors/wasm32-unknown-unknown/release/wave2d_actor.wasm"
  if cargo build --release --locked --target wasm32-unknown-unknown \
       --manifest-path "$actor_manifest" --target-dir target/dyne-actors \
     && BARYE_WAVE2D_WASM="$wasm" cargo test --release --locked --test wasm_actor; then
    echo "wave2d wasm      OK"
  else
    echo "wave2d wasm      FAIL"; fail=1
  fi
elif [ "${CI:-}" = "true" ]; then
  echo "wave2d wasm      FAIL (target wasm32-unknown-unknown missing in CI)"; fail=1
else
  echo "wave2d wasm      SKIPPED: install the target with 'rustup target add wasm32-unknown-unknown'"
fi

echo "== [4] 3D wave: split = single patch (bit-exact) =="
out=$(cargo run --release --locked --quiet --bin barye-wave3d -- check --grid 25x23x21 --steps 100         --splits 1x1x1,2x2x2,3x2x1,4x4x2,1x3x5) || true
echo "$out"
passes=$(printf '%s\n' "$out" | grep -c "PASS (bit-exact)" || true)
if [ "$passes" -eq 5 ] && printf '%s\n' "$out" | grep -q "single patch         crc32 2ca26bb5"; then
  echo "wave3d split     OK"
else
  echo "wave3d split     FAIL"; fail=1
fi

echo "== [5] 2D and 3D wave: convergence order against analytic solutions =="
for dim in 2d 3d; do
  out=$(cargo run --release --locked --quiet --bin barye-wave$dim -- convergence) || true
  echo "$out"
  if printf '%s\n' "$out" | grep -q "^PASS (order within \[1.9, 2.1\])"; then
    echo "wave$dim order    OK"
  else
    echo "wave$dim order    FAIL"; fail=1
  fi
done

if [ "$fail" -ne 0 ]; then echo "CI CHECK FAILED"; exit 1; fi
echo "ALL CHECKS PASSED"
