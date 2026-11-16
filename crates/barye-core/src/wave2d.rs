//! The 2D wave equation on one patch: a native Rust solver with the same floating-point operations,
//! in the same order, as Dyne's WebAssembly wave solver (`wave2d-actor` in the Dyne repository,
//! described in Dyne's `docs/2026-11-11-wave2d.md`).
//!
//! Because every operation is the same IEEE 754 operation in the same order (Rust does not fuse
//! multiplications and additions unless asked to), this solver gives the same bits as Dyne's solver
//! running in WebAssembly, and a split run gives the same bits as the single-patch run.
//!
//! # Method
//!
//! `u_tt = c² (u_xx + u_yy)` on a square grid (spacing `h` on both axes), five-point Laplacian and the
//! three-level central difference in time (the leapfrog scheme):
//!
//! ```text
//! u^{n+1} = 2 u^n − u^{n−1} + ν² ((u_W − 2u + u_E) + (u_S − 2u + u_N)),   ν = c Δt / h
//! ```
//!
//! The first step starts from zero initial velocity: `u^1 = u^0 + (ν²/2) (same Laplacian)`.
//! The scheme is stable for `ν ≤ 1/√2` on an unbounded grid (von Neumann analysis); the solver does
//! not check this, so that sweeps across the limit stay possible ([`Wave2dProblem::courant`]).
//!
//! Only the current level `u^n` of a neighbour is used, so a patch needs one layer of ghost points and
//! exchanges one row of values per face and step.
//!
//! # Outer boundaries (per face)
//!
//! - [`Outer::Fixed`]: the edge points are reset to a value after every step (Dirichlet condition).
//! - [`Outer::ZeroGradient`]: the ghost point mirrors the first interior point (Neumann condition).
//! - [`Outer::Absorbing`]: first-order Engquist–Majda condition discretised with Mur's method:
//!   `u_0^{n+1} = u_1^n + k (u_1^{n+1} − u_0^n)`, `k = (ν−1)/(ν+1)`. The face is first updated as a
//!   zero-gradient face, then its edge points are overwritten.
//!
//! Absorbing faces are applied in face order west, east, south, north, then fixed faces in the same
//! order (at a corner, the later face wins), exactly as in Dyne's solver.
//!
//! # Towards 3D
//!
//! The patch interface ([`PatchSolver`]), the decomposition and the border messages do not depend on the
//! dimension. A 3D solver is a new type next to this one (seven-point Laplacian, stability limit
//! `ν ≤ 1/√3`, six faces whose values are sent as flat lists); nothing else needs to change.

use crate::decomp::Block;
use crate::patch::{PatchError, PatchSolver};
use dyne_wire as wire;
use dyne_wire::abi::field2d::{BC_DIRICHLET, BC_INTERFACE, BC_ZERO_GRAD};
use dyne_wire::abi::wave2d as keys;

/// Smallest number of points per axis of a patch (as in Dyne's solver).
pub const MIN_POINTS: usize = 3;

const PI: f64 = core::f64::consts::PI;

/// Initial displacement (the initial velocity is always zero). `x` and `y` are `i·h` and `j·h` for
/// the global point indices `i`, `j`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Initial {
    /// `amp · sin(mπx) · sin(nπy)`: a standing wave on the unit square with fixed value 0.
    SinMode { m: f64, n: f64, amp: f64 },
    /// `amp · exp(−((x−x0)² + (y−y0)²) / (2σ²))`: a Gaussian bump.
    Gauss { x0: f64, y0: f64, sigma: f64, amp: f64 },
    /// `amp · cos(mπx) · cos(nπy)`: a standing wave on the unit square with zero gradient.
    CosMode { m: f64, n: f64, amp: f64 },
}

/// Condition on one face of the outer boundary.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Outer {
    /// Fixed value (Dirichlet condition).
    Fixed(f64),
    /// Zero gradient (Neumann condition, free edge).
    ZeroGradient,
    /// First-order absorbing boundary (Engquist–Majda condition, Mur's discretisation).
    Absorbing,
}

/// A 2D wave problem on a grid of `points[0] × points[1]` points with spacing `h`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Wave2dProblem {
    /// Number of grid points along x and y.
    pub points: [usize; 2],
    /// Grid spacing (the same on both axes).
    pub h: f64,
    /// Wave speed.
    pub c: f64,
    /// Time step.
    pub dt: f64,
    /// Initial displacement.
    pub initial: Initial,
    /// Outer boundary per face: west, east, south, north (face `2·axis + side`).
    pub outer: [Outer; 4],
}

impl Wave2dProblem {
    /// Courant number `ν = c Δt / h` (stable for `ν ≤ 1/√2` on an unbounded grid).
    pub fn courant(&self) -> f64 {
        self.c * self.dt / self.h
    }

    /// The configuration calls (`dyne_cfg(key, value)`) that set up Dyne's WebAssembly solver for this
    /// problem (keys from `dyne_wire::abi::wave2d`).
    pub fn actor_cfg(&self) -> Vec<(u32, f64)> {
        let mut cfg = match self.initial {
            Initial::SinMode { m, n, amp } => vec![
                (keys::CFG_IC, keys::IC_SIN_MODE as f64),
                (keys::CFG_MODE_M, m),
                (keys::CFG_MODE_N, n),
                (keys::CFG_AMP, amp),
            ],
            Initial::Gauss { x0, y0, sigma, amp } => vec![
                (keys::CFG_IC, keys::IC_GAUSS as f64),
                (keys::CFG_GX, x0),
                (keys::CFG_GY, y0),
                (keys::CFG_SIGMA, sigma),
                (keys::CFG_AMP, amp),
            ],
            Initial::CosMode { m, n, amp } => vec![
                (keys::CFG_IC, keys::IC_COS_MODE as f64),
                (keys::CFG_MODE_M, m),
                (keys::CFG_MODE_N, n),
                (keys::CFG_AMP, amp),
            ],
        };
        let absorb: u32 = (0..4).filter(|&f| self.outer[f] == Outer::Absorbing).map(|f| 1 << f).sum();
        cfg.push((keys::CFG_ABSORB, absorb as f64));
        cfg
    }

    /// The face settings (`dyne_set_face(face, kind, param)`) of Dyne's WebAssembly solver for one
    /// block: an interface carries its id, a fixed face its value as `f64` bits.
    pub fn actor_faces(&self, b: &Block<2>) -> [(u32, u64); 4] {
        let mut out = [(0, 0); 4];
        for (f, o) in out.iter_mut().enumerate() {
            *o = match (b.faces[f], self.outer[f]) {
                (Some(i), _) => (BC_INTERFACE, i.id as u64),
                (None, Outer::Fixed(v)) => (BC_DIRICHLET, v.to_bits()),
                // An absorbing face is a zero-gradient face with its bit set in CFG_ABSORB.
                (None, Outer::ZeroGradient | Outer::Absorbing) => (BC_ZERO_GRAD, 0),
            };
        }
        out
    }

    /// Initial displacement at global point `(gi, gj)` (the same expression, in the same order, as
    /// Dyne's solver).
    fn initial_at(&self, gi: usize, gj: usize) -> f64 {
        let x = (gi as f64) * self.h;
        let y = (gj as f64) * self.h;
        match self.initial {
            Initial::Gauss { x0, y0, sigma, amp } => {
                let (dx, dy) = (x - x0, y - y0);
                amp * libm::exp(-(dx * dx + dy * dy) / (2.0 * sigma * sigma))
            }
            Initial::CosMode { m, n, amp } => amp * libm::cos(m * PI * x) * libm::cos(n * PI * y),
            Initial::SinMode { m, n, amp } => amp * libm::sin(m * PI * x) * libm::sin(n * PI * y),
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

/// The native 2D wave solver for one patch.
pub struct Wave2dPatch {
    nx: usize,
    ny: usize,
    nu2: f64,
    half_nu2: f64,
    kmur: f64,
    face: [Face; 4],
    /// Current level `u^n` with one ring of ghost points, row length `nx + 2`.
    u: Vec<f64>,
    /// Previous level `u^{n−1}` (overwritten point by point with `u^{n+1}` during a step).
    p: Vec<f64>,
    got: [bool; 4],
    step: u64,
    /// Scratch space for the values of one face as little-endian bytes.
    edge: Vec<u8>,
}

impl Wave2dPatch {
    /// Sets up the solver of block `b` of `problem` at step 0.
    pub fn new(problem: &Wave2dProblem, b: &Block<2>) -> Result<Self, PatchError> {
        let [nx, ny] = b.size;
        if nx < MIN_POINTS || ny < MIN_POINTS {
            return Err(PatchError(format!("patch of {nx}×{ny} points; at least {MIN_POINTS} per axis needed")));
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
        let mut face = [Face::ZeroGradient; 4];
        for (f, slot) in face.iter_mut().enumerate() {
            *slot = match (b.faces[f], problem.outer[f]) {
                (Some(i), _) => Face::Interface(i.id),
                (None, Outer::Fixed(v)) => Face::Fixed(v),
                (None, Outer::ZeroGradient) => Face::ZeroGradient,
                (None, Outer::Absorbing) => Face::Absorbing,
            };
        }
        let w = nx + 2;
        let mut s = Wave2dPatch {
            nx,
            ny,
            nu2,
            half_nu2: 0.5 * nu2,
            kmur: (nu - 1.0) / (nu + 1.0),
            face,
            u: vec![0.0; w * (ny + 2)],
            p: vec![0.0; w * (ny + 2)],
            got: [false; 4],
            step: 0,
            edge: Vec::with_capacity(8 * nx.max(ny)),
        };
        let [ix0, iy0] = b.origin;
        for j in 0..ny {
            for i in 0..nx {
                let k = s.at(i as isize, j as isize);
                s.u[k] = problem.initial_at(ix0 + i, iy0 + j);
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
    fn at(&self, i: isize, j: isize) -> usize {
        ((j + 1) as usize) * (self.nx + 2) + (i + 1) as usize
    }

    fn face_len(&self, f: usize) -> usize {
        if f < 2 {
            self.ny
        } else {
            self.nx
        }
    }

    /// Point `k` along face `f`, at `depth` points inward from the edge (−1 = ghost, 0 = edge,
    /// 1 = first interior point).
    fn face_point(&self, f: usize, k: usize, depth: isize) -> usize {
        let (nx, ny, k) = (self.nx as isize, self.ny as isize, k as isize);
        match f {
            0 => self.at(depth, k),
            1 => self.at(nx - 1 - depth, k),
            2 => self.at(k, depth),
            _ => self.at(k, ny - 1 - depth),
        }
    }

    fn apply_fixed(&mut self) {
        for f in 0..4 {
            if let Face::Fixed(v) = self.face[f] {
                for k in 0..self.face_len(f) {
                    let idx = self.face_point(f, k, 0);
                    self.u[idx] = v;
                }
            }
        }
    }

    /// Fills the ghost points of the outer faces from the current level (the ghost points of the
    /// interfaces were filled by [`PatchSolver::inbox`]).
    fn fill_outer_ghosts(&mut self) {
        for f in 0..4 {
            match self.face[f] {
                Face::Interface(_) => {}
                Face::ZeroGradient | Face::Absorbing => {
                    for k in 0..self.face_len(f) {
                        let (g, m) = (self.face_point(f, k, -1), self.face_point(f, k, 1));
                        self.u[g] = self.u[m];
                    }
                }
                Face::Fixed(v) => {
                    for k in 0..self.face_len(f) {
                        let g = self.face_point(f, k, -1);
                        self.u[g] = v;
                    }
                }
            }
        }
    }

    /// Interior update of all points (rows are independent; each point uses the same operations
    /// in the same order as Dyne's solver).
    fn update_interior(&mut self) {
        let (nx, ny, w) = (self.nx, self.ny, self.nx + 2);
        let (nu2, half_nu2) = (self.nu2, self.half_nu2);
        let first = self.step == 0;
        for j in 1..=ny {
            let south = &self.u[(j - 1) * w + 1..(j - 1) * w + 1 + nx];
            let north = &self.u[(j + 1) * w + 1..(j + 1) * w + 1 + nx];
            let west = &self.u[j * w..j * w + nx];
            let mid = &self.u[j * w + 1..j * w + 1 + nx];
            let east = &self.u[j * w + 2..j * w + 2 + nx];
            let out = &mut self.p[j * w + 1..j * w + 1 + nx];
            if first {
                for i in 0..nx {
                    let u = mid[i];
                    let lap = (west[i] - 2.0 * u + east[i]) + (south[i] - 2.0 * u + north[i]);
                    out[i] = u + half_nu2 * lap;
                }
            } else {
                for i in 0..nx {
                    let u = mid[i];
                    let lap = (west[i] - 2.0 * u + east[i]) + (south[i] - 2.0 * u + north[i]);
                    out[i] = 2.0 * u - out[i] + nu2 * lap;
                }
            }
        }
    }

    /// Mur's absorbing condition on the edge points of absorbing faces (new level in `p`).
    fn apply_absorbing(&mut self) {
        for f in 0..4 {
            if self.face[f] != Face::Absorbing {
                continue;
            }
            for k in 0..self.face_len(f) {
                let (e, m) = (self.face_point(f, k, 0), self.face_point(f, k, 1));
                self.p[e] = self.u[m] + self.kmur * (self.p[m] - self.u[e]);
            }
        }
    }
}

impl PatchSolver for Wave2dPatch {
    fn outbox(&mut self, out: &mut Vec<u8>) -> Result<(), PatchError> {
        for f in 0..4 {
            let Face::Interface(id) = self.face[f] else { continue };
            let n = self.face_len(f);
            self.edge.clear();
            for k in 0..n {
                let v = self.u[self.face_point(f, k, 0)];
                self.edge.extend_from_slice(&v.to_bits().to_le_bytes());
            }
            let hdr = wire::Header {
                kind: wire::KIND_IFACE_BAND,
                quantity: wire::QTY_DISPLACEMENT,
                src: id,
                step: self.step,
                t_us: 0,
            };
            let payload = wire::Payload::IfaceBand { data: &self.edge, count: n as u32, per_point: 1 };
            let start = out.len();
            out.resize(start + 4 + wire::HDR_LEN + 8 + 8 * n, 0);
            let mut mb = wire::MailboxWriter::new(&mut out[start..]);
            mb.push(&hdr, &payload).map_err(|e| PatchError(format!("cannot encode border message: {e:?}")))?;
            let used = mb.used();
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
            let f = (0..4)
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
            for k in 0..count as usize {
                let g = self.face_point(f, k, -1);
                self.u[g] = wire::band_get(data, k * pp);
            }
            self.got[f] = true;
        }
        Ok(())
    }

    fn step(&mut self) -> Result<(), PatchError> {
        for f in 0..4 {
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
        self.got = [false; 4];
        Ok(())
    }

    fn field(&mut self) -> Result<Vec<f64>, PatchError> {
        let mut out = Vec::with_capacity(self.nx * self.ny);
        let w = self.nx + 2;
        for j in 1..=self.ny {
            out.extend_from_slice(&self.u[j * w + 1..j * w + 1 + self.nx]);
        }
        Ok(out)
    }
}
