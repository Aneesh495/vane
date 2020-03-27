//! Momentum exchange aerodynamic force and moment evaluation engine.
//!
//! Evaluates instantaneous aerodynamic drag, lift, pitching moment, and
//! non-dimensional coefficients (Cd, Cl, Cm) along immersed solid boundaries.
//!
//! Note: Zero em dashes are used in this file per coding standards.

use crate::lattice::{CX, CY, OPP, Q};
use serde::{Deserialize, Serialize};

/// Aerodynamic forces and non-dimensional coefficients.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct AerodynamicForces {
    /// Dimensional drag force Fx (lattice units).
    pub fx: f64,
    /// Dimensional lift force Fy (lattice units).
    pub fy: f64,
    /// Pitching moment about reference center (positive nose-up).
    pub moment: f64,
    /// Drag coefficient Cd = Fx / (0.5 * rho * u^2 * chord).
    pub cd: f64,
    /// Lift coefficient Cl = Fy / (0.5 * rho * u^2 * chord).
    pub cl: f64,
    /// Pitching moment coefficient Cm = M / (0.5 * rho * u^2 * chord^2).
    pub cm: f64,
}

/// Compute aerodynamic forces and moments using momentum exchange algorithm.
pub fn compute_momentum_exchange(
    f_post: &[f64],
    mask: &[bool],
    nx: usize,
    ny: usize,
    rho0: f64,
    u_inf: f64,
    chord: f64,
    x_ref: (f64, f64),
) -> AerodynamicForces {
    let mut fx_total = 0.0;
    let mut fy_total = 0.0;
    let mut m_total = 0.0;

    let (xref_x, xref_y) = x_ref;

    for iy in 0..ny {
        for ix in 0..nx {
            let f_idx = iy * nx + ix;
            // Only examine fluid nodes
            if mask[f_idx] {
                continue;
            }

            // Check neighbors in directions 1..8
            for i in 1..Q {
                let cx = CX[i];
                let cy = CY[i];
                let nx_pos = ix as i32 + cx;
                let ny_pos = iy as i32 + cy;

                // Check bounds
                if nx_pos >= 0 && nx_pos < (nx as i32) && ny_pos >= 0 && ny_pos < (ny as i32) {
                    let s_idx = (ny_pos as usize) * nx + (nx_pos as usize);
                    // If neighbor is solid, this forms a boundary link
                    if mask[s_idx] {
                        let ibar = OPP[i];
                        let fi = f_post[i * (nx * ny) + f_idx];
                        let fibar = f_post[ibar * (nx * ny) + f_idx];

                        // Momentum transfer per link: Delta P = c_i * (f_i + f_ibar)
                        let delta_p = fi + fibar;
                        let df_x = (cx as f64) * delta_p;
                        let df_y = (cy as f64) * delta_p;

                        fx_total += df_x;
                        fy_total += df_y;

                        // Moment contribution: r x Delta F = (rx * df_y) - (ry * df_x)
                        let rx = (ix as f64 + 0.5 * cx as f64) - xref_x;
                        let ry = (iy as f64 + 0.5 * cy as f64) - xref_y;
                        m_total += (rx * df_y) - (ry * df_x);
                    }
                }
            }
        }
    }

    let q_inf = 0.5 * rho0 * u_inf * u_inf;
    let denom_f = (q_inf * chord).max(1e-12);
    let denom_m = (q_inf * chord * chord).max(1e-12);

    AerodynamicForces {
        fx: fx_total,
        fy: fy_total,
        moment: m_total,
        cd: fx_total / denom_f,
        cl: fy_total / denom_f,
        cm: m_total / denom_m,
    }
}
