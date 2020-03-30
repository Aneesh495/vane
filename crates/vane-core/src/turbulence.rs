//! Subgrid scale Large Eddy Simulation (LES) turbulence models.
//!
//! Implements:
//! - Algebraic Smagorinsky model with resolved non-equilibrium strain rate
//! - Wall-Adapting Local Eddy-viscosity (WALE) model
//! - Local eddy viscosity nu_t and effective relaxation time tau_eff
//!
//! Note: Zero em dashes are used in this file per coding standards.

use crate::lattice::{CX, CY, CS2, Q};

/// Algebraic Smagorinsky LES model.
///
/// Computes local eddy viscosity nu_t = (Cs * Delta)^2 * ||S||
#[derive(Debug, Clone, Copy)]
pub struct Smagorinsky {
    /// Smagorinsky constant Cs (typically 0.10 to 0.18).
    pub cs: f64,
    /// Lattice filter width Delta (default 1.0).
    pub delta: f64,
}

impl Default for Smagorinsky {
    fn default() -> Self {
        Self {
            cs: 0.14,
            delta: 1.0,
        }
    }
}

impl Smagorinsky {
    /// Create new Smagorinsky model.
    pub fn new(cs: f64, delta: f64) -> Self {
        Self { cs, delta }
    }

    /// Compute eddy viscosity nu_t from local non-equilibrium stress tensor.
    pub fn compute_eddy_viscosity(&self, f: &[f64; Q], feq: &[f64; Q], rho: f64, tau0: f64) -> f64 {
        // Non-equilibrium momentum flux: Pi_alpha_beta = sum_i (c_i_alpha * c_i_beta * (f_i - f_i^eq))
        let mut pi_xx = 0.0;
        let mut pi_yy = 0.0;
        let mut pi_xy = 0.0;

        for i in 1..Q {
            let f_neq = f[i] - feq[i];
            let cx = CX[i] as f64;
            let cy = CY[i] as f64;

            pi_xx += cx * cx * f_neq;
            pi_yy += cy * cy * f_neq;
            pi_xy += cx * cy * f_neq;
        }

        // Frobenius norm of Pi: ||Pi|| = sqrt(2 * Pi_ij * Pi_ij)
        let pi_norm = (2.0 * (pi_xx * pi_xx + pi_yy * pi_yy + 2.0 * pi_xy * pi_xy)).sqrt();

        // Strain rate norm: ||S|| = ||Pi|| / (2 * rho * cs^2 * tau0)
        let s_norm = pi_norm / (2.0 * rho * CS2 * tau0).max(1e-12);

        // Eddy viscosity: nu_t = (Cs * Delta)^2 * ||S||
        let filter_len = self.cs * self.delta;
        (filter_len * filter_len) * s_norm
    }

    /// Effective relaxation time: tau_eff = tau0 + 3 * nu_t.
    pub fn effective_tau(&self, f: &[f64; Q], feq: &[f64; Q], rho: f64, tau0: f64) -> f64 {
        let nu_t = self.compute_eddy_viscosity(f, feq, rho, tau0);
        tau0 + 3.0 * nu_t
    }
}

/// Wall-Adapting Local Eddy-viscosity (WALE) model.
///
/// Automatically provides y^3 near-wall damping for turbulent boundary layers.
#[derive(Debug, Clone, Copy)]
pub struct Wale {
    /// WALE model constant Cw (typically 0.325).
    pub cw: f64,
    /// Lattice filter width Delta.
    pub delta: f64,
}

impl Default for Wale {
    fn default() -> Self {
        Self {
            cw: 0.325,
            delta: 1.0,
        }
    }
}

impl Wale {
    /// Create new WALE model.
    pub fn new(cw: f64, delta: f64) -> Self {
        Self { cw, delta }
    }

    /// Compute WALE eddy viscosity from strain and traceless deformation tensors.
    pub fn compute_eddy_viscosity(&self, s_norm: f64, sd_norm: f64) -> f64 {
        let s_sq = s_norm * s_norm;
        let sd_sq = sd_norm * sd_norm;

        let sd_cubed = sd_sq * sd_norm;
        let s_five = s_sq * s_sq * s_norm;
        let sd_five = sd_sq * sd_sq * sd_norm;

        let denom = (s_five + sd_five).powf(1.5).max(1e-14);
        let c_scale = (self.cw * self.delta).powi(2);

        c_scale * (sd_cubed / denom)
    }
}
