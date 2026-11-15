//! Core of Barye: physical simulation from many cooperating local solvers.
//!
//! The domain is split into patches; each patch has its own small solver (classical or learned),
//! and patches cooperate only by exchanging values at their borders. Barye builds on the crates of
//! [Dyne](https://github.com/physical-edge-native/dyne), which runs the same scheme on
//! microcontrollers, so that both projects share one definition of:
//!
//! - the border-exchange messages (wire format version 2), re-exported as [`wire`];
//! - the deterministic simulator parts (random numbers, discrete-event virtual time, virtual
//!   network, loading of WebAssembly solvers), re-exported as [`sim`].
//!
//! Nothing physical is implemented yet; this crate currently fixes the dependency on Dyne and the
//! tests that guard it.

/// Dyne's border-exchange format (wire format version 2) and mailbox conventions.
pub use dyne_wire as wire;

/// Dyne's deterministic simulator parts (random numbers, virtual time, virtual network,
/// WebAssembly solver loading).
pub use dyne_core as sim;
