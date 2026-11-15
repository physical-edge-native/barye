# Barye

**Physical simulation without GPUs, built from many cooperating local solvers.**

The name comes from the barye, the CGS unit of pressure (1 barye = 1 dyne/cm²) — a sibling of
[Dyne](https://github.com/physical-edge-native/dyne), named after the CGS unit of force.

Part of the [Physical-Edge Native Architecture](https://github.com/physical-edge-native) project.

> **Status: just started (2026-11).** Nothing usable yet. The design below is a draft
> and open questions are tracked as issues.

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

## Open design questions

- Reuse Dyne's crates (wire format, WebAssembly solver interface) as a dependency,
  or keep an independent core and share only the specification?
- Language and parallel runtime (Rust threads, MPI, WebAssembly on PCs)?
- Which first problem and benchmark (e.g. 2D wave equation vs. a reference solver)?

## License

TODO: license not yet chosen.
