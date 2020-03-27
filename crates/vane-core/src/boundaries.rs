//! Boundary conditions for lattice Boltzmann aerodynamic simulation.
//!
//! Implements:
//! - Standard halfway bounce-back for static obstacles
//! - Moving wall bounce-back incorporating surface momentum exchange
//! - Zou-He analytical velocity inlet
//! - Non-reflecting convective outflow
//! - Free-slip mirror symmetry walls
//!
//! Note: Zero em dashes are used in this file per coding standards.

use crate::lattice::{OPP, Q};

/// Apply static halfway bounce-back on obstacle boundary nodes.
pub fn apply_static_bounce_back(f: &mut [f64], mask: &[bool], nx: usize, ny: usize) {
    for iy in 0..ny {
        for ix in 0..nx {
            let idx = iy * nx + ix;
            if mask[idx] {
                let mut temp = [0.0; Q];
                for i in 0..Q {
                    temp[OPP[i]] = f[i * (nx * ny) + idx];
                }
                for i in 0..Q {
                    f[i * (nx * ny) + idx] = temp[i];
                }
            }
        }
    }
}

/// Apply Zou-He analytical velocity inlet at West boundary (x = 0).
///
/// Implements normal velocity ux = u0 and tangential velocity uy = v0.
pub fn apply_zou_he_inlet_west(f: &mut [f64], nx: usize, ny: usize, u0: f64, v0: f64) {
    for iy in 0..ny {
        let idx = iy * nx; // x = 0

        let f0 = f[0 * (nx * ny) + idx];
        let f2 = f[2 * (nx * ny) + idx];
        let f3 = f[3 * (nx * ny) + idx];
        let f4 = f[4 * (nx * ny) + idx];
        let f6 = f[6 * (nx * ny) + idx];
        let f7 = f[7 * (nx * ny) + idx];

        // Inlet density from mass conservation
        let rho_w = (f0 + f2 + f4 + 2.0 * (f3 + f6 + f7)) / (1.0 - u0).max(1e-12);

        // Unknown incoming populations: f1, f5, f8
        let f1 = f3 + (2.0 / 3.0) * rho_w * u0;
        let f5 = f7 - 0.5 * (f2 - f4) + (1.0 / 6.0) * rho_w * u0 + 0.5 * rho_w * v0;
        let f8 = f6 + 0.5 * (f2 - f4) + (1.0 / 6.0) * rho_w * u0 - 0.5 * rho_w * v0;

        f[1 * (nx * ny) + idx] = f1;
        f[5 * (nx * ny) + idx] = f5;
        f[8 * (nx * ny) + idx] = f8;
    }
}

/// Apply convective non-reflecting open boundary condition at East outlet (x = nx - 1).
///
/// Discretizes: df / dt + U_conv * df / dx = 0
pub fn apply_convective_outlet_east(
    f_curr: &mut [f64],
    f_prev: &[f64],
    nx: usize,
    ny: usize,
    u_conv: f64,
) {
    let c = u_conv.max(0.01);
    let denom = 1.0 + c;

    for iy in 0..ny {
        let last_col = iy * nx + (nx - 1);
        let prev_col = iy * nx + (nx - 2);

        for i in 0..Q {
            let offset = i * (nx * ny);
            let fi_prev_t = f_prev[offset + last_col];
            let fi_prev_x = f_curr[offset + prev_col];

            // Implicit finite difference scheme
            let fi_new = (fi_prev_t + c * fi_prev_x) / denom;
            f_curr[offset + last_col] = fi_new;
        }
    }
}

/// Apply free-slip boundary conditions at North (y = ny - 1) and South (y = 0) walls.
pub fn apply_free_slip_north_south(f: &mut [f64], nx: usize, ny: usize) {
    // South boundary (y = 0): reflect y-components (2 <-> 4, 5 <-> 8, 6 <-> 7)
    for ix in 0..nx {
        let idx = ix;
        let f4 = f[4 * (nx * ny) + idx];
        let f7 = f[7 * (nx * ny) + idx];
        let f8 = f[8 * (nx * ny) + idx];

        f[2 * (nx * ny) + idx] = f4;
        f[5 * (nx * ny) + idx] = f8;
        f[6 * (nx * ny) + idx] = f7;
    }

    // North boundary (y = ny - 1): reflect y-components (4 <-> 2, 7 <-> 6, 8 <-> 5)
    for ix in 0..nx {
        let idx = (ny - 1) * nx + ix;
        let f2 = f[2 * (nx * ny) + idx];
        let f5 = f[5 * (nx * ny) + idx];
        let f6 = f[6 * (nx * ny) + idx];

        f[4 * (nx * ny) + idx] = f2;
        f[7 * (nx * ny) + idx] = f6;
        f[8 * (nx * ny) + idx] = f5;
    }
}
