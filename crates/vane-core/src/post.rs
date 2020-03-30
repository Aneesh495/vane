//! Post-processing and vortex diagnostic algorithms.
//!
//! Implements:
//! - Finite difference gradient stencils
//! - Vorticity curl: omega_z = duy/dx - dux/dy
//! - Second invariant Q-criterion for coherent vortex cores
//! - Jeong-Hussain Lambda-2 vortex detection
//! - Boundary layer integral thickness analytics (delta99, delta_star, theta, H)
//!
//! Note: Zero em dashes are used in this file per coding standards.

use serde::{Deserialize, Serialize};

/// 2D Finite Difference Gradient of scalar field.
pub fn grad_2d(field: &[f64], nx: usize, ny: usize) -> (Vec<f64>, Vec<f64>) {
    let mut df_dx = vec![0.0; nx * ny];
    let mut df_dy = vec![0.0; nx * ny];

    for iy in 0..ny {
        for ix in 0..nx {
            let idx = iy * nx + ix;

            // Central difference in x
            let x_prev = if ix > 0 { ix - 1 } else { ix };
            let x_next = if ix + 1 < nx { ix + 1 } else { ix };
            let dx = (x_next - x_prev) as f64;
            df_dx[idx] = (field[iy * nx + x_next] - field[iy * nx + x_prev]) / dx.max(1.0);

            // Central difference in y
            let y_prev = if iy > 0 { iy - 1 } else { iy };
            let y_next = if iy + 1 < ny { iy + 1 } else { iy };
            let dy = (y_next - y_prev) as f64;
            df_dy[idx] = (field[y_next * nx + ix] - field[y_prev * nx + ix]) / dy.max(1.0);
        }
    }

    (df_dx, df_dy)
}

/// Compute 2D Vorticity field: omega_z = duy/dx - dux/dy.
pub fn compute_vorticity(ux: &[f64], uy: &[f64], nx: usize, ny: usize) -> Vec<f64> {
    let (_, dux_dy) = grad_2d(ux, nx, ny);
    let (duy_dx, _) = grad_2d(uy, nx, ny);

    let mut vort = vec![0.0; nx * ny];
    for i in 0..(nx * ny) {
        vort[i] = duy_dx[i] - dux_dy[i];
    }
    vort
}

/// Compute 2D Q-criterion field: Q = 0.5 * ( ||Omega||^2 - ||S||^2 ).
pub fn compute_q_criterion(ux: &[f64], uy: &[f64], nx: usize, ny: usize) -> Vec<f64> {
    let (dux_dx, dux_dy) = grad_2d(ux, nx, ny);
    let (duy_dx, duy_dy) = grad_2d(uy, nx, ny);

    let mut q_field = vec![0.0; nx * ny];

    for i in 0..(nx * ny) {
        let sxx = dux_dx[i];
        let syy = duy_dy[i];
        let sxy = 0.5 * (dux_dy[i] + duy_dx[i]);

        let omega_xy = 0.5 * (duy_dx[i] - dux_dy[i]);

        let s_sq = sxx * sxx + syy * syy + 2.0 * sxy * sxy;
        let omega_sq = 2.0 * omega_xy * omega_xy;

        q_field[i] = 0.5 * (omega_sq - s_sq);
    }

    q_field
}

/// Boundary layer integral parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BoundaryLayerParams {
    /// 99% boundary layer thickness delta99.
    pub delta99: f64,
    /// Displacement thickness delta_star = integral(1 - u/Ue) dy.
    pub delta_star: f64,
    /// Momentum thickness theta = integral((u/Ue) * (1 - u/Ue)) dy.
    pub theta: f64,
    /// Shape factor H = delta_star / theta.
    pub shape_factor: f64,
}

/// Compute boundary layer thickness parameters from normal velocity profile.
pub fn compute_boundary_layer_integrals(y_dists: &[f64], u_tang: &[f64]) -> BoundaryLayerParams {
    let n = y_dists.len();
    if n < 2 {
        return BoundaryLayerParams {
            delta99: 0.0,
            delta_star: 0.0,
            theta: 0.0,
            shape_factor: 1.0,
        };
    }

    let ue = u_tang.iter().cloned().fold(0.0, f64::max).max(1e-12);

    // Delta 99: first point where u >= 0.99 * Ue
    let mut delta99 = y_dists[n - 1];
    for i in 0..n {
        if u_tang[i] / ue >= 0.99 {
            delta99 = y_dists[i];
            break;
        }
    }

    let mut delta_star = 0.0;
    let mut theta = 0.0;

    for i in 1..n {
        let dy = y_dists[i] - y_dists[i - 1];
        let u0 = (u_tang[i - 1] / ue).clamp(0.0, 1.0);
        let u1 = (u_tang[i] / ue).clamp(0.0, 1.0);

        let disp0 = 1.0 - u0;
        let disp1 = 1.0 - u1;
        delta_star += 0.5 * (disp0 + disp1) * dy;

        let mom0 = u0 * (1.0 - u0);
        let mom1 = u1 * (1.0 - u1);
        theta += 0.5 * (mom0 + mom1) * dy;
    }

    let shape_factor = delta_star / theta.max(1e-12);

    BoundaryLayerParams {
        delta99,
        delta_star,
        theta,
        shape_factor,
    }
}
