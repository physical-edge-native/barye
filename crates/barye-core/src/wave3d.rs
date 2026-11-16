//! The 3D wave equation on one patch: a native Rust solver, the 3D extension of [`crate::wave2d`].
//!
//! Dyne has no 3D WebAssembly solver yet. The order of every floating-point operation is fixed and
//! written down below, so that a future 3D WebAssembly solver can be built to give the same bits and
//! be compared with this one, as the 2D solvers are. A split run gives the same bits as the
//! single-patch run because every point performs the same operations in the same order whatever
//! patch it belongs to (Rust does not fuse multiplications and additions unless asked to).
//!
//! # Method
//!
//! `u_tt = c² (u_xx + u_yy + u_zz)` on a cubic grid (spacing `h` on all axes), seven-point Laplacian
//! and the three-level central difference in time (the leapfrog scheme), `ν = c Δt / h`. The scheme is
//! stable for `ν ≤ 1/√3` on an unbounded grid (von Neumann analysis); the solver does not check this,
//! so that sweeps across the limit stay possible ([`Wave3dProblem::courant`]).
//!
//! # Fixed operation order (the contract for bit-exact comparison)
//!
//! With `u` the current level, `W, E` its neighbours along x (low, high), `S, N` along y and
//! `D, U` along z, every point computes, in this order and with these parentheses:
//!
//! ```text
//! lap = ((W − 2.0·u + E) + (S − 2.0·u + N)) + (D − 2.0·u + U)
//! first step (n = 0):  u^1     = u + half_nu2 · lap
//! later steps:         u^{n+1} = 2.0·u − u^{n−1} + nu2 · lap
//! ```
//!
//! where each bracket `(A − 2.0·u + B)` is `(A − (2.0·u)) + B`, and the coefficients are computed once
//! at set-up: `nu = c · dt / h`, `nu2 = nu · nu`, `half_nu2 = 0.5 · nu2`, `kmur = (nu − 1.0) / (nu + 1.0)`.
//! The x and y terms are exactly those of the 2D solver, so a 3D problem that does not vary along z
//! (with zero-gradient z faces) reproduces the 2D solver bit for bit in every z slice (tested).
//!
//! The first step starts from zero initial velocity. Initial values are computed with the pure-Rust
//! `libm` (identical bits on every platform), from `x = i·h`, `y = j·h`, `z = k·h` with the global
//! point indices `i, j, k` (see [`Initial`] for each expression and its order).
//!
//! # Outer boundaries (per face)
//!
//! - [`Outer::Fixed`]: the face points are reset to a value after every step (Dirichlet condition);
//!   the ghost points of the face hold the same value.
//! - [`Outer::ZeroGradient`]: the ghost point mirrors the first interior point (Neumann condition).
//! - [`Outer::Absorbing`]: first-order Engquist–Majda condition discretised with Mur's method:
//!   `u_0^{n+1} = u_1^n + kmur · (u_1^{n+1} − u_0^n)`. The face is first updated as a zero-gradient
//!   face, then its points are overwritten.
//!
//! One step is: fill the ghost points of the outer faces (faces 0 to 5 in order), update every point
//! (the formula above), apply the absorbing faces in face order 0 to 5, make the new level current,
//! then reset the fixed faces in face order 0 to 5. Where faces meet (edges and corners of the box),
//! the later face wins. This is the 2D order with the two z faces appended.
//!
//! # Faces and border messages
//!
//! Faces are numbered `2·axis + side` (x low 0, x high 1, y low 2, y high 3, z low 4, z high 5). Only
//! the six faces are exchanged (the seven-point stencil uses no edge or corner neighbour): each step,
//! a patch sends the 2D slab of values on each face that touches a neighbour as one IFACE_BAND
//! message (one `f64` per point). The slab is flattened in the block's point order: the lower of the
//! two in-face axes runs fastest (face x: y fastest, then z; face y: x, then z; face z: x, then y).

use crate::decomp::Block;
use crate::patch::{PatchError, PatchSolver};
use dyne_wire as wire;

pub use crate::wave2d::Outer;

/// Smallest number of points per axis of a patch.
pub const MIN_POINTS: usize = 3;

const PI: f64 = core::f64::consts::PI;

/// Initial displacement (the initial velocity is always zero). `x`, `y`, `z` are `i·h`, `j·h`, `k·h`
/// for the global point indices; products are evaluated left to right as written.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Initial {
    /// `amp · sin(mπx) · sin(nπy) · sin(lπz)`: a standing wave in the unit cube with fixed value 0.
    SinMode { m: f64, n: f64, l: f64, amp: f64 },
    /// `amp · exp(−(dx² + dy² + dz²) / (2σ²))` with `dx = x − x0` etc. (sum left to right): a
    /// Gaussian bump.
    Gauss { x0: f64, y0: f64, z0: f64, sigma: f64, amp: f64 },
    /// `amp · cos(mπx) · cos(nπy) · cos(lπz)`: a standing wave in the unit cube with zero gradient.
    CosMode { m: f64, n: f64, l: f64, amp: f64 },
}

/// A 3D wave problem on a grid of `points[0] × points[1] × points[2]` points with spacing `h`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Wave3dProblem {
    /// Number of grid points along x, y and z.
    pub points: [usize; 3],
    /// Grid spacing (the same on all axes).
    pub h: f64,
    /// Wave speed.
    pub c: f64,
    /// Time step.
    pub dt: f64,
    /// Initial displacement.
    pub initial: Initial,
    /// Outer boundary per face `2·axis + side`: x low, x high, y low, y high, z low, z high.
    pub outer: [Outer; 6],
}

impl Wave3dProblem {
    /// Courant number `ν = c Δt / h` (stable for `ν ≤ 1/√3` on an unbounded grid).
    pub fn courant(&self) -> f64 {
        self.c * self.dt / self.h
    }

    /// Initial displacement at global point `(gi, gj, gk)`.
    fn initial_at(&self, gi: usize, gj: usize, gk: usize) -> f64 {
        let x = (gi as f64) * self.h;
        let y = (gj as f64) * self.h;
        let z = (gk as f64) * self.h;
        match self.initial {
            Initial::Gauss { x0, y0, z0, sigma, amp } => {
                let (dx, dy, dz) = (x - x0, y - y0, z - z0);
                amp * libm::exp(-(dx * dx + dy * dy + dz * dz) / (2.0 * sigma * sigma))
            }
            Initial::CosMode { m, n, l, amp } => {
                amp * libm::cos(m * PI * x) * libm::cos(n * PI * y) * libm::cos(l * PI * z)
            }
            Initial::SinMode { m, n, l, amp } => {
                amp * libm::sin(m * PI * x) * libm::sin(n * PI * y) * libm::sin(l * PI * z)
            }
        }
    }
}

/// What one face of a patch is.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Face {
    Interface(u32),
    Fixed(f64),
    ZeroGradient,
    Absorbing,
}

/// The native 3D wave solver for one patch.
pub struct Wave3dPatch {
    n: [usize; 3],
    nu2: f64,
    half_nu2: f64,
    kmur: f64,
    face: [Face; 6],
    /// Current level `u^n` with one layer of ghost points: `(nx+2)·(ny+2)·(nz+2)` values, x fastest.
    u: Vec<f64>,
    /// Previous level `u^{n−1}` (overwritten point by point with `u^{n+1}` during a step).
    p: Vec<f64>,
    got: [bool; 6],
    step: u64,
    /// Scratch space for the values of one face as little-endian bytes.
    edge: Vec<u8>,
}

impl Wave3dPatch {
    /// Sets up the solver of block `b` of `problem` at step 0.
    pub fn new(problem: &Wave3dProblem, b: &Block<3>) -> Result<Self, PatchError> {
        let n = b.size;
        if n.iter().any(|&k| k < MIN_POINTS) {
            return Err(PatchError(format!(
                "patch of {}×{}×{} points; at least {MIN_POINTS} per axis needed",
                n[0], n[1], n[2]
            )));
        }
        if !(problem.h > 0.0) || !(problem.c > 0.0) || !(problem.dt > 0.0) {
            return Err(PatchError("h, c and dt must be positive".into()));
        }
        if let Initial::Gauss { sigma, .. } = problem.initial {
            if !(sigma > 0.0) {
                return Err(PatchError("sigma must be positive".into()));
            }
        }
        let nu = problem.c * problem.dt / problem.h;
        let nu2 = nu * nu;
        let mut face = [Face::ZeroGradient; 6];
        for (f, slot) in face.iter_mut().enumerate() {
            *slot = match (b.faces[f], problem.outer[f]) {
                (Some(i), _) => Face::Interface(i.id),
                (None, Outer::Fixed(v)) => Face::Fixed(v),
                (None, Outer::ZeroGradient) => Face::ZeroGradient,
                (None, Outer::Absorbing) => Face::Absorbing,
            };
        }
        let total = (n[0] + 2) * (n[1] + 2) * (n[2] + 2);
        let largest_face = (n[0] * n[1]).max(n[0] * n[2]).max(n[1] * n[2]);
        let mut s = Wave3dPatch {
            n,
            nu2,
            half_nu2: 0.5 * nu2,
            kmur: (nu - 1.0) / (nu + 1.0),
            face,
            u: vec![0.0; total],
            p: vec![0.0; total],
            got: [false; 6],
            step: 0,
            edge: Vec::with_capacity(8 * largest_face),
        };
        let [ix0, iy0, iz0] = b.origin;
        for k in 0..n[2] {
            for j in 0..n[1] {
                for i in 0..n[0] {
                    let idx = s.at(i as isize, j as isize, k as isize);
                    s.u[idx] = problem.initial_at(ix0 + i, iy0 + j, iz0 + k);
                }
            }
        }
        s.apply_fixed();
        Ok(s)
    }

    /// Number of steps done.
    pub fn steps_done(&self) -> u64 {
        self.step
    }

    #[inline]
    fn at(&self, i: isize, j: isize, k: isize) -> usize {
        let (wx, wy) = (self.n[0] + 2, self.n[1] + 2);
        (((k + 1) as usize) * wy + (j + 1) as usize) * wx + (i + 1) as usize
    }

    /// Sizes of the slab of face `f`: (fast in-face axis, slow in-face axis).
    fn face_dims(&self, f: usize) -> (usize, usize) {
        match f / 2 {
            0 => (self.n[1], self.n[2]),
            1 => (self.n[0], self.n[2]),
            _ => (self.n[0], self.n[1]),
        }
    }

    fn face_len(&self, f: usize) -> usize {
        let (a, b) = self.face_dims(f);
        a * b
    }

    /// Point `(a, b)` of the slab of face `f` (`a` along the fast in-face axis), at `depth` points
    /// inward from the face (−1 = ghost, 0 = face point, 1 = first interior point).
    fn face_point(&self, f: usize, a: usize, b: usize, depth: isize) -> usize {
        let (a, b) = (a as isize, b as isize);
        let along = |axis: usize| -> isize {
            if f % 2 == 0 {
                depth
            } else {
                self.n[axis] as isize - 1 - depth
            }
        };
        match f / 2 {
            0 => self.at(along(0), a, b),
            1 => self.at(a, along(1), b),
            _ => self.at(a, b, along(2)),
        }
    }

    /// Calls `g(a, b)` for every point of the slab of face `f`, in slab order.
    fn for_slab(&self, f: usize, mut g: impl FnMut(usize, usize)) {
        let (na, nb) = self.face_dims(f);
        for b in 0..nb {
            for a in 0..na {
                g(a, b);
            }
        }
    }

    fn apply_fixed(&mut self) {
        for f in 0..6 {
            if let Face::Fixed(v) = self.face[f] {
                let (na, nb) = self.face_dims(f);
                for b in 0..nb {
                    for a in 0..na {
                        let idx = self.face_point(f, a, b, 0);
                        self.u[idx] = v;
                    }
                }
            }
        }
    }

    /// Fills the ghost points of the outer faces from the current level (the ghost points of the
    /// interfaces were filled by [`PatchSolver::inbox`]).
    fn fill_outer_ghosts(&mut self) {
        for f in 0..6 {
            let (na, nb) = self.face_dims(f);
            match self.face[f] {
                Face::Interface(_) => {}
                Face::ZeroGradient | Face::Absorbing => {
                    for b in 0..nb {
                        for a in 0..na {
                            let (g, m) = (self.face_point(f, a, b, -1), self.face_point(f, a, b, 1));
                            self.u[g] = self.u[m];
                        }
                    }
                }
                Face::Fixed(v) => {
                    for b in 0..nb {
                        for a in 0..na {
                            let g = self.face_point(f, a, b, -1);
                            self.u[g] = v;
                        }
                    }
                }
            }
        }
    }

    /// Update of all points (rows along x are independent; every point uses the operations of the
    /// module documentation, in that order).
    fn update_interior(&mut self) {
        let [nx, ny, nz] = self.n;
        let w = nx + 2;
        let plane = w * (ny + 2);
        let (nu2, half_nu2) = (self.nu2, self.half_nu2);
        let first = self.step == 0;
        for k in 1..=nz {
            for j in 1..=ny {
                let row = k * plane + j * w + 1;
                let west = &self.u[row - 1..row - 1 + nx];
                let mid = &self.u[row..row + nx];
                let east = &self.u[row + 1..row + 1 + nx];
                let south = &self.u[row - w..row - w + nx];
                let north = &self.u[row + w..row + w + nx];
                let down = &self.u[row - plane..row - plane + nx];
                let up = &self.u[row + plane..row + plane + nx];
                let out = &mut self.p[row..row + nx];
                if first {
                    for i in 0..nx {
                        let u = mid[i];
                        let lap = ((west[i] - 2.0 * u + east[i]) + (south[i] - 2.0 * u + north[i]))
                            + (down[i] - 2.0 * u + up[i]);
                        out[i] = u + half_nu2 * lap;
                    }
                } else {
                    for i in 0..nx {
                        let u = mid[i];
                        let lap = ((west[i] - 2.0 * u + east[i]) + (south[i] - 2.0 * u + north[i]))
                            + (down[i] - 2.0 * u + up[i]);
                        out[i] = 2.0 * u - out[i] + nu2 * lap;
                    }
                }
            }
        }
    }

    /// Mur's absorbing condition on the points of absorbing faces (new level in `p`).
    fn apply_absorbing(&mut self) {
        for f in 0..6 {
            if self.face[f] != Face::Absorbing {
                continue;
            }
            let (na, nb) = self.face_dims(f);
            for b in 0..nb {
                for a in 0..na {
                    let (e, m) = (self.face_point(f, a, b, 0), self.face_point(f, a, b, 1));
                    self.p[e] = self.u[m] + self.kmur * (self.p[m] - self.u[e]);
                }
            }
        }
    }
}

impl PatchSolver for Wave3dPatch {
    fn outbox(&mut self, out: &mut Vec<u8>) -> Result<(), PatchError> {
        for f in 0..6 {
            let Face::Interface(id) = self.face[f] else { continue };
            let n = self.face_len(f);
            let mut edge = std::mem::take(&mut self.edge);
            edge.clear();
            self.for_slab(f, |a, b| {
                let v = self.u[self.face_point(f, a, b, 0)];
                edge.extend_from_slice(&v.to_bits().to_le_bytes());
            });
            let hdr = wire::Header {
                kind: wire::KIND_IFACE_BAND,
                quantity: wire::QTY_DISPLACEMENT,
                src: id,
                step: self.step,
                t_us: 0,
            };
            let payload = wire::Payload::IfaceBand { data: &edge, count: n as u32, per_point: 1 };
            let start = out.len();
            out.resize(start + 4 + wire::HDR_LEN + 8 + 8 * n, 0);
            let mut mb = wire::MailboxWriter::new(&mut out[start..]);
            let pushed = mb.push(&hdr, &payload);
            let used = mb.used();
            self.edge = edge;
            pushed.map_err(|e| PatchError(format!("cannot encode border message: {e:?}")))?;
            out.truncate(start + used);
        }
        Ok(())
    }

    fn inbox(&mut self, framed: &[u8]) -> Result<(), PatchError> {
        for msg in wire::MailboxReader::new(framed) {
            let (hdr, payload, _) = msg.map_err(|e| PatchError(format!("malformed border message: {e:?}")))?;
            let wire::Payload::IfaceBand { data, count, per_point } = payload else {
                return Err(PatchError(format!("unexpected message kind {}", hdr.kind)));
            };
            let f = (0..6)
                .find(|&f| self.face[f] == Face::Interface(hdr.src))
                .ok_or_else(|| PatchError(format!("message from unknown interface {}", hdr.src)))?;
            if hdr.step != self.step || count as usize != self.face_len(f) || per_point < 1 {
                return Err(PatchError(format!(
                    "border message (interface {}, step {}, {} values) does not fit (step {}, {} values)",
                    hdr.src,
                    hdr.step,
                    count,
                    self.step,
                    self.face_len(f)
                )));
            }
            let pp = per_point as usize;
            let (na, nb) = self.face_dims(f);
            for b in 0..nb {
                for a in 0..na {
                    let g = self.face_point(f, a, b, -1);
                    self.u[g] = wire::band_get(data, (b * na + a) * pp);
                }
            }
            self.got[f] = true;
        }
        Ok(())
    }

    fn step(&mut self) -> Result<(), PatchError> {
        for f in 0..6 {
            if matches!(self.face[f], Face::Interface(_)) && !self.got[f] {
                return Err(PatchError(format!("step {}: border values of face {f} missing", self.step)));
            }
        }
        self.fill_outer_ghosts();
        self.update_interior();
        self.apply_absorbing();
        std::mem::swap(&mut self.u, &mut self.p);
        self.apply_fixed();
        self.step += 1;
        self.got = [false; 6];
        Ok(())
    }

    fn field(&mut self) -> Result<Vec<f64>, PatchError> {
        let [nx, ny, nz] = self.n;
        let mut out = Vec::with_capacity(nx * ny * nz);
        for k in 0..nz {
            for j in 0..ny {
                let row = self.at(0, j as isize, k as isize);
                out.extend_from_slice(&self.u[row..row + nx]);
            }
        }
        Ok(out)
    }
}
