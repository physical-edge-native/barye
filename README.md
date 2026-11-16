# Barye

**3D field physics on CPUs, built from many cooperating local solvers.**

The name (pronounced /ˈbæri/, "BAR-ee") comes from the barye, the CGS unit of pressure (1 barye = 1 dyne/cm²) — a sibling of
[Dyne](https://github.com/physical-edge-native/dyne), named after the CGS unit of force.

Part of the [Physical-Edge Native Architecture](https://github.com/physical-edge-native) project.

> **Status: early development (2026-11).** The 2D and 3D wave equations run split into
> patches, one thread per patch, and match the single-patch solution bit for bit. The design
> below is a draft and open questions are tracked as issues.

## The question this repository answers

Large physical simulations usually run as one big model on GPUs. This repository
asks how far the opposite approach can go: split the domain into many small
patches, give each patch a small solver (classical or learned), and let patches
cooperate only by exchanging values at their borders — the same scheme that
[Dyne](https://github.com/physical-edge-native/dyne) runs on microcontrollers,
here on ordinary CPU cores and PCs.

Target problems are those where the toolchain matters, i.e. not the plain heat
equation (which has closed-form treatments):

- wave equation (vibration fields), with von Neumann stability analysis
- fluid flow
- mixing classical and learned solvers, including problems where only partial
  stability guarantees are available

## Relationship to Dyne

| | Dyne | Barye |
|---|---|---|
| Runs on | microcontrollers (and a simulator of them) | CPU cores, PCs, clusters |
| Main concern | deployment to real devices, networks, faults | numerical methods and scale |
| Shared | interface of a local solver, border-exchange format, stability conditions | |

## Positioning

Barye complements robot simulators such as Isaac Lab and MuJoCo rather than replacing
them: it aims at field physics (waves, heat, coarse airflow) on CPUs, bit-exact results
across platforms, and mixed classical and learned patches. Rigid-body contact,
reinforcement-learning throughput and rendering are not goals. **The target is 3D physics;
2D problems are development milestones only.** See the design record
[docs/design/2026-11-15-positioning.md](docs/design/2026-11-15-positioning.md).

Barye is written in Rust and uses Dyne's crates (border-exchange format, WebAssembly
solver interface, deterministic simulator parts) as a dependency pinned to one commit,
instead of copying them.

## Building and checking

Requires a stable Rust toolchain; no GPU or hardware is needed.

```sh
cargo test                  # unit tests
bash scripts/ci-check.sh    # the same checks as CI (use Git Bash on Windows)
```

`scripts/ci-check.sh` also compares the native wave solver with Dyne's WebAssembly wave solver,
which it builds from the pinned Dyne commit; this needs `rustup target add wasm32-unknown-unknown`
(the check is skipped with a message otherwise).

## First milestone: the 2D wave equation split into patches

The wave equation `u_tt = c² (u_xx + u_yy)` is solved with the leapfrog scheme (second order in
space and time, stable for Courant numbers `c Δt / h ≤ 1/√2`). The grid is split into blocks;
each block is advanced by its own solver on its own thread, and blocks exchange only the values
on their borders each step, using Dyne's message format. The solver performs the same
floating-point operations as Dyne's WebAssembly wave solver, so both give identical bits (why a
native solver: [docs/design/2026-11-16-native-patch-solver.md](docs/design/2026-11-16-native-patch-solver.md)).

```sh
# split 1x1 ... 8x8 and 3x5 must equal the single-patch solve bit for bit
cargo run --release --bin barye-wave2d -- check --grid 41x41 --steps 400 --splits 1x1,2x2,4x4,8x8,3x5
# strong scaling (fixed grid) and weak scaling (fixed points per patch) over 1..N cores
cargo run --release --bin barye-wave2d -- scaling --mode strong --grid 2048x2048 --cores 1,2,4,8
cargo run --release --bin barye-wave2d -- scaling --mode weak --patch 256x256 --cores 1,2,4,8
```

The decomposition, the face numbering and the border messages are written for any number of
dimensions, so that a 3D solver only adds a new patch solver.

## The 3D wave equation split into patches

`u_tt = c² (u_xx + u_yy + u_zz)` is solved with the same leapfrog scheme and a seven-point
Laplacian (stable for Courant numbers `c Δt / h ≤ 1/√3`), with fixed-value, zero-gradient or
absorbing outer faces. Patches exchange only the values on their six faces, each face as one
list of values (a 2D slab) in the same message format; the seven-point stencil needs no edge or
corner values. The order of every floating-point operation is fixed and documented in
[`wave3d.rs`](crates/barye-core/src/wave3d.rs), so that another implementation (for example a
future WebAssembly solver for Dyne's devices) can be checked against it bit for bit. A 3D problem
that does not vary along z gives the 2D solver's bits in every z slice.

```sh
# even and uneven splits must equal the single-patch solve bit for bit
cargo run --release --bin barye-wave3d -- check --grid 25x23x21 --steps 100 --splits 1x1x1,2x2x2,3x2x1,4x4x2,1x3x5
# strong and weak scaling over 1..N cores
cargo run --release --bin barye-wave3d -- scaling --mode strong --grid 192x192x192 --cores 1,2,4,8
cargo run --release --bin barye-wave3d -- scaling --mode weak --patch 64x64x64 --cores 1,2,4,8
```

## Wider borders: one exchange every k steps

With `--halo k`, each patch keeps `k` layers of its neighbours' values and exchanges them only once
every `k` steps, advancing the overlapping layers itself in between (the correct part shrinks by one
layer per step). Because the stencil's reach grows diagonally, patches then also exchange with the
patches that touch them only at an edge or a corner. Every point still performs the same operations,
so the result equals the single-patch result bit for bit for every `k`; a `k` larger than a patch
along a split axis is rejected. `k = 1` (the default) is the one-layer, every-step exchange of Dyne's
WebAssembly solver, message for message. Details are in [`halo.rs`](crates/barye-core/src/halo.rs).

```sh
cargo run --release --bin barye-wave2d -- check --grid 41x41 --splits 2x2,4x4,3x5 --halo 1,2,3,4
cargo run --release --bin barye-wave3d -- scaling --mode strong --grid 192x192x192 --cores 1,2,4,8 --halo 4
```

Wider borders cut the number of exchanges by a factor of `k`; they do not reduce the memory traffic
of a step (the overlap is even computed twice). In first measurements on a laptop with cores of
different speeds (preliminary, on battery), `k = 2` and `k = 4` gave no measurable speed-up in 2D and
were slower in 3D.

## Accuracy: convergence order against analytic solutions

Bit-exactness shows that splitting does not change the result, not that the result is accurate.
Both solvers are therefore compared with analytic standing waves (`sin(πx)·sin(πy)·cos(ωt)` and its
3D analogue, fixed value 0) at a fixed Courant number of 0.5 on refined grids. The maximum and L2
errors at a fixed time fall as `h²`: every observed order must lie in `[1.9, 2.1]` (the reasons for
this band are given in [`convergence.rs`](crates/barye-core/src/convergence.rs)); the measured orders
are 2.000–2.003 in 2D and 2.001–2.004 in 3D.

```sh
cargo run --release --bin barye-wave2d -- convergence
cargo run --release --bin barye-wave3d -- convergence
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for how to propose changes.

## Open design questions

- Parallel runtime beyond one PC (MPI, WebAssembly on PCs)? The first milestone uses Rust
  threads on one PC.
- How to keep scaling when the update is limited by memory bandwidth? Wider borders (above) cut
  the waiting between patches; computing several steps per cache-sized tile (temporal blocking)
  is not implemented.
- A 3D interface for WebAssembly solvers shared with Dyne (the native 3D solver exists).

## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or https://www.apache.org/licenses/LICENSE-2.0)
- MIT license ([LICENSE-MIT](LICENSE-MIT) or https://opensource.org/licenses/MIT)

at your option.

### Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.
