//! Barye's native 2D wave solver against Dyne's WebAssembly wave solver (`wave2d-actor`), run with
//! the `wasmi` interpreter: same field bits, as one patch and split across threads, and the same
//! border-message bytes at every step.
//!
//! The WebAssembly file is a build product and is not kept in the repository. `scripts/ci-check.sh`
//! builds it from the pinned Dyne commit and passes its path in `BARYE_WAVE2D_WASM`; without that
//! variable this test reports that it was skipped and passes.

use barye_core::decomp::Decomposition;
use barye_core::patch::{run_single, run_threaded, PatchSolver};
use barye_core::wasm2d::Wasm2dPatch;
use barye_core::wave2d::{Initial, Outer, Wave2dPatch, Wave2dProblem, MIN_POINTS};

fn load() -> Option<(wasmi::Engine, wasmi::Module)> {
    let Ok(path) = std::env::var("BARYE_WAVE2D_WASM") else {
        eprintln!(
            "SKIPPED: BARYE_WAVE2D_WASM is not set (run `bash scripts/ci-check.sh`, which builds Dyne's \
             wave2d-actor from the pinned commit and sets it)"
        );
        return None;
    };
    let bytes = std::fs::read(&path).unwrap_or_else(|e| panic!("cannot read {path}: {e}"));
    let engine = wasmi::Engine::default();
    let module = wasmi::Module::new(&engine, &bytes[..]).expect("valid WebAssembly module");
    Some((engine, module))
}

fn problems() -> Vec<(&'static str, Wave2dProblem)> {
    let gauss = Initial::Gauss { x0: 0.3, y0: 0.6, sigma: 0.06, amp: 1.0 };
    let p = |initial, outer| Wave2dProblem { points: [45, 31], h: 1.0 / 44.0, c: 1.0, dt: 0.5 / 44.0, initial, outer };
    vec![
        ("fixed", p(gauss, [Outer::Fixed(0.0); 4])),
        ("zero gradient", p(Initial::CosMode { m: 2.0, n: 1.0, amp: 0.5 }, [Outer::ZeroGradient; 4])),
        ("absorbing", p(gauss, [Outer::Absorbing; 4])),
        ("mixed", p(gauss, [Outer::Absorbing, Outer::Fixed(0.5), Outer::Fixed(0.5), Outer::Absorbing])),
    ]
}

fn bits(v: &[f64]) -> Vec<u64> {
    v.iter().map(|x| x.to_bits()).collect()
}

#[test]
fn native_and_webassembly_solvers_agree_bit_for_bit() {
    let Some((engine, module)) = load() else { return };
    let steps = 200;
    for (name, p) in problems() {
        // One patch.
        let whole = Decomposition::whole(p.points, MIN_POINTS).unwrap();
        let native = run_single(&mut Wave2dPatch::new(&p, &whole.blocks[0]).unwrap(), steps).unwrap();
        let mut w = Wasm2dPatch::new_wave(&engine, &module, &p, &whole.blocks[0]).unwrap();
        let wasm = run_single(&mut w, steps).unwrap();
        assert_eq!(bits(&native), bits(&wasm), "{name}: single patch");

        // Split 3×2 (uneven), one WebAssembly instance per thread.
        let d = Decomposition::new(p.points, [3, 2], MIN_POINTS).unwrap();
        let split = run_threaded(&d, steps, |b| Wasm2dPatch::new_wave(&engine, &module, &p, b)).unwrap();
        assert_eq!(bits(&native), bits(&split.field), "{name}: WebAssembly split 3x2");
    }
}

/// Both solvers write byte-identical border messages at every step (driven in lockstep on one
/// thread, with each patch's messages delivered from the native solvers).
#[test]
fn border_messages_are_byte_identical() {
    let Some((engine, module)) = load() else { return };
    let (_, p) = problems().remove(3);
    let d = Decomposition::new(p.points, [2, 2], MIN_POINTS).unwrap();
    let mut native: Vec<Wave2dPatch> = d.blocks.iter().map(|b| Wave2dPatch::new(&p, b).unwrap()).collect();
    let mut wasm: Vec<Wasm2dPatch> =
        d.blocks.iter().map(|b| Wasm2dPatch::new_wave(&engine, &module, &p, b).unwrap()).collect();
    for step in 0..100 {
        let mut outs = Vec::new();
        for (k, (n, w)) in native.iter_mut().zip(wasm.iter_mut()).enumerate() {
            let (mut a, mut b) = (Vec::new(), Vec::new());
            n.outbox(&mut a).unwrap();
            w.outbox(&mut b).unwrap();
            assert_eq!(a, b, "step {step}, block {k}: outbox bytes differ");
            outs.push(a);
        }
        for (k, b) in d.blocks.iter().enumerate() {
            for i in b.faces.iter().flatten() {
                // Deliver the neighbour's message that crosses interface i.
                let msgs = barye_core::wire::MailboxReader::new(&outs[i.neighbour]);
                for m in msgs {
                    let (h, _, raw) = m.unwrap();
                    if h.src == i.id {
                        let mut framed = (raw.len() as u32).to_le_bytes().to_vec();
                        framed.extend_from_slice(raw);
                        native[k].inbox(&framed).unwrap();
                        wasm[k].inbox(&framed).unwrap();
                    }
                }
            }
        }
        for (n, w) in native.iter_mut().zip(wasm.iter_mut()) {
            n.step().unwrap();
            w.step().unwrap();
        }
    }
}
