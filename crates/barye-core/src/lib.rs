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
//! What is implemented:
//!
//! - [`decomp`]: splitting a grid of any dimension into blocks, with face and interface numbering;
//! - [`patch`]: the patch-solver interface and the runners (single patch; one thread per patch,
//!   exchanging only border messages with neighbouring patches);
//! - [`wave2d`]: a native 2D wave-equation patch solver that gives the same bits as Dyne's
//!   WebAssembly wave solver;
//! - [`wasm2d`]: running Dyne's 2D WebAssembly patch solvers as patches;
//! - [`wave3d`]: a native 3D wave-equation patch solver (seven-point Laplacian) with a documented,
//!   fixed order of floating-point operations;
//! - [`convergence`]: the convergence-order check of the wave solvers against analytic standing waves.

pub mod convergence;
pub mod decomp;
pub mod patch;
pub mod wasm2d;
pub mod wave2d;
pub mod wave3d;

/// Dyne's border-exchange format (wire format version 2) and mailbox conventions.
pub use dyne_wire as wire;

/// Dyne's deterministic simulator parts (random numbers, virtual time, virtual network,
/// WebAssembly solver loading).
pub use dyne_core as sim;

/// CRC-32 (IEEE 802.3) of the little-endian bits of a field: the field hash that Dyne's simulators
/// print, so that results can be compared across the two projects.
pub fn field_crc32(field: &[f64]) -> u32 {
    let mut c: u32 = 0xFFFF_FFFF;
    for v in field {
        for b in v.to_bits().to_le_bytes() {
            c ^= b as u32;
            for _ in 0..8 {
                let m = (c & 1).wrapping_neg();
                c = (c >> 1) ^ (0xEDB8_8320 & m);
            }
        }
    }
    !c
}

#[cfg(test)]
mod tests {
    #[test]
    fn crc32_check_value() {
        // The standard check value: CRC-32 of the ASCII bytes "123456789" is 0xCBF43926. Here the
        // input must be whole f64 values, so use "12345678" (one value), whose CRC-32 is 0x9AE0DAAF.
        let v = f64::from_bits(u64::from_le_bytes(*b"12345678"));
        assert_eq!(super::field_crc32(&[v]), 0x9AE0_DAAF);
    }
}
