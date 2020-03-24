//! D2Q9 Lattice Boltzmann hydrodynamics engine.
//!
//! Provides discrete velocity stencils, equilibrium distribution calculations,
//! macroscopic density and velocity moment extraction, and streaming/collision operators.
//!
//! Note: Zero em dashes are used in this file per coding standards.

use serde::{Deserialize, Serialize};

/// Number of discrete lattice velocity directions in D2Q9.
pub const Q: usize = 9;

/// Lattice weights for D2Q9 stencil.
pub const W: [f64; Q] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

/// Discrete lattice velocity x-components.
pub const CX: [i32; Q] = [0, 1, 0, -1, 0, 1, -1, -1, 1];

/// Discrete lattice velocity y-components.
pub const CY: [i32; Q] = [0, 0, 1, 0, -1, 1, 1, -1, -1];

/// Opposite direction index mapping for bounce-back reflections: c[OPP[i]] = -c[i].
pub const OPP: [usize; Q] = [0, 3, 4, 1, 2, 7, 8, 5, 6];

/// Speed of sound squared: cs^2 = 1/3 in lattice units.
pub const CS2: f64 = 1.0 / 3.0;

/// Inverse speed of sound squared: 1 / cs^2 = 3.0.
pub const INV_CS2: f64 = 3.0;

/// Factor 1 / (2 * cs^4) = 4.5.
pub const INV_2CS4: f64 = 4.5;

/// Factor 1 / (2 * cs^2) = 1.5.
pub const INV_2CS2: f64 = 1.5;

/// Macroscopic state at a single lattice site.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MacroscopicState {
    /// Fluid density rho.
    pub rho: f64,
    /// Fluid horizontal velocity ux.
    pub ux: f64,
    /// Fluid vertical velocity uy.
    pub uy: f64,
}

impl MacroscopicState {
    /// Create new macroscopic state.
    pub fn new(rho: f64, ux: f64, uy: f64) -> Self {
        Self { rho, ux, uy }
    }

    /// Calculate thermodynamic pressure: p = cs^2 * rho.
    pub fn pressure(&self) -> f64 {
        CS2 * self.rho
    }

    /// Velocity magnitude: ||u|| = sqrt(ux^2 + uy^2).
    pub fn speed(&self) -> f64 {
        (self.ux * self.ux + self.uy * self.uy).sqrt()
    }

    /// Local Mach number: Ma = speed / cs.
    pub fn mach(&self) -> f64 {
        self.speed() / CS2.sqrt()
    }
}

/// Compute macroscopic density and velocity from 9 discrete populations.
#[inline]
pub fn compute_moments(f: &[f64; Q]) -> MacroscopicState {
    let rho: f64 = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
    let inv_rho = if rho.abs() > 1e-14 { 1.0 / rho } else { 1.0 };

    let jx = (f[1] - f[3]) + (f[5] - f[6] - f[7] + f[8]);
    let jy = (f[2] - f[4]) + (f[5] + f[6] - f[7] - f[8]);

    MacroscopicState {
        rho,
        ux: jx * inv_rho,
        uy: jy * inv_rho,
    }
}

/// Compute second-order polynomial Maxwell-Boltzmann equilibrium distribution f_i^eq.
#[inline]
pub fn compute_equilibrium(rho: f64, ux: f64, uy: f64) -> [f64; Q] {
    let u_sq = ux * ux + uy * uy;
    let u_sq_term = INV_2CS2 * u_sq;

    let mut feq = [0.0; Q];

    // Rest particle (direction 0)
    feq[0] = W[0] * rho * (1.0 - u_sq_term);

    // Coordinate directions 1..4 and diagonal directions 5..8
    for i in 1..Q {
        let ci_dot_u = (CX[i] as f64) * ux + (CY[i] as f64) * uy;
        feq[i] = W[i] * rho * (1.0 + INV_CS2 * ci_dot_u + INV_2CS4 * ci_dot_u * ci_dot_u - u_sq_term);
    }

    feq
}

/// Single Relaxation Time (SRT / BGK) collision operator.
///
/// f_i^* = f_i - (1 / tau) * (f_i - f_i^eq)
#[inline]
pub fn collide_bgk(f: &mut [f64; Q], feq: &[f64; Q], tau: f64) {
    let omega = 1.0 / tau.max(1e-12);
    for i in 0..Q {
        f[i] -= omega * (f[i] - feq[i]);
    }
}

/// Two-Relaxation-Time (TRT) collision operator with magic parameter Lambda.
///
/// Separates distributions into symmetric (even) and anti-symmetric (odd) components:
///   f_i^+ = 0.5 * (f_i + f_ibar)
///   f_i^- = 0.5 * (f_i - f_ibar)
/// Viscous relaxation omega_plus = 1 / tau, odd relaxation omega_minus from magic Lambda = 1/4.
pub fn collide_trt(f: &mut [f64; Q], feq: &[f64; Q], tau_plus: f64, magic_lambda: f64) {
    let omega_plus = 1.0 / tau_plus.max(1e-12);
    let tau_minus = 0.5 + magic_lambda / (tau_plus - 0.5).max(1e-12);
    let omega_minus = 1.0 / tau_minus.max(1e-12);

    let mut f_post = [0.0; Q];
    f_post[0] = f[0] - omega_plus * (f[0] - feq[0]);

    for i in 1..Q {
        let ibar = OPP[i];
        let f_plus = 0.5 * (f[i] + f[ibar]);
        let f_minus = 0.5 * (f[i] - f[ibar]);

        let feq_plus = 0.5 * (feq[i] + feq[ibar]);
        let feq_minus = 0.5 * (feq[i] - feq[ibar]);

        let f_star_plus = f_plus - omega_plus * (f_plus - feq_plus);
        let f_star_minus = f_minus - omega_minus * (f_minus - feq_minus);

        f_post[i] = f_star_plus + f_star_minus;
    }

    *f = f_post;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lattice_constants() {
        let weight_sum: f64 = W.iter().sum();
        assert!((weight_sum - 1.0).abs() < 1e-12);

        for i in 0..Q {
            let opp_idx = OPP[i];
            assert_eq!(CX[opp_idx], -CX[i]);
            assert_eq!(CY[opp_idx], -CY[i]);
        }
    }

    #[test]
    fn test_equilibrium_moments() {
        let rho_in = 1.05;
        let ux_in = 0.04;
        let uy_in = -0.02;

        let feq = compute_equilibrium(rho_in, ux_in, uy_in);
        let mom = compute_moments(&feq);

        assert!((mom.rho - rho_in).abs() < 1e-12);
        assert!((mom.ux - ux_in).abs() < 1e-12);
        assert!((mom.uy - uy_in).abs() < 1e-12);
    }

    #[test]
    fn test_collision_conservation() {
        let rho_in = 1.0;
        let ux_in = 0.03;
        let uy_in = 0.01;
        let mut f = compute_equilibrium(rho_in, ux_in, uy_in);

        // Perturb non-equilibrium
        f[1] += 0.005;
        f[3] += 0.005;
        f[0] -= 0.010;

        let mom_before = compute_moments(&f);
        let feq = compute_equilibrium(mom_before.rho, mom_before.ux, mom_before.uy);

        collide_bgk(&mut f, &feq, 0.55);
        let mom_after = compute_moments(&f);

        assert!((mom_after.rho - mom_before.rho).abs() < 1e-12);
        assert!((mom_after.ux - mom_before.ux).abs() < 1e-12);
        assert!((mom_after.uy - mom_before.uy).abs() < 1e-12);
    }
}
