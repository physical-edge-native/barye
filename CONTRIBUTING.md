# Contributing to Barye

Thank you for your interest. Barye is a research prototype maintained by a
single author, so response times may vary (usually within a week).

## Reporting problems and asking questions

Open a [GitHub issue](https://github.com/physical-edge-native/barye/issues).
Please include your OS, Rust version (`rustc --version`), the exact command,
and its output.

## Proposing changes

1. Open an issue first for anything larger than a small fix, so we can agree on the approach.
2. Fork, create a branch, and make your change.
3. Run `bash scripts/ci-check.sh`. All checks must pass: numerical results are
   compared **bit for bit**, so any change in output is treated as a regression
   unless it is intended and explained in the pull request.
4. Add one changelog fragment `changelog.d/<issue>.<section>.md` (see
   [changelog.d/README.md](changelog.d/README.md)); do not edit `CHANGELOG.md` directly.
5. Open a pull request describing what changed and why, with `Closes #<issue>`.

Issues and pull requests may be written in English or Japanese.

Barye depends on crates of [Dyne](https://github.com/physical-edge-native/dyne) pinned to
one commit (see the root `Cargo.toml`). Changes to the shared parts (border-exchange format,
WebAssembly solver interface) belong in Dyne; moving the pin is a separate pull request.

## Scope

In scope: solvers for field physics (e.g. waves, heat, airflow) split into patches that
exchange only border values, mixing classical and learned solvers, stability conditions at
patch borders, and tools to check results against reference solvers. The target is 3D physics;
2D problems are development milestones. Rigid-body contact, articulated robots and rendering are
out of scope (see [the positioning design record](docs/design/2026-11-15-positioning.md)).
Larger directions are discussed across the
[Physical-Edge Native Architecture](https://github.com/physical-edge-native) organization.

## Code of conduct

This project follows the organization-wide
[Code of Conduct](https://github.com/physical-edge-native/.github/blob/main/CODE_OF_CONDUCT.md).

## License of contributions

Barye is licensed under either of MIT or Apache-2.0, at your option (see [README](README.md#license)).
Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion in this work
by you, as defined in the Apache-2.0 license, shall be dual licensed as above, without any additional terms or conditions.
