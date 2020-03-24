//! Parametric aerodynamic airfoil geometry generator and discrete grid rasterizer.
//!
//! Implements:
//! - NACA 4-digit analytical profile coordinates (camber line, thickness distribution)
//! - Joukowsky conformal transformation
//! - Subgrid and discrete binary lattice obstacle rasterization
//!
//! Note: Zero em dashes are used in this file per coding standards.

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Coordinates of an aerodynamic profile.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AirfoilCoordinates {
    /// Upper suction surface x-coordinates.
    pub x_upper: Vec<f64>,
    /// Upper suction surface y-coordinates.
    pub y_upper: Vec<f64>,
    /// Lower pressure surface x-coordinates.
    pub x_lower: Vec<f64>,
    /// Lower pressure surface y-coordinates.
    pub y_lower: Vec<f64>,
    /// Chord length c.
    pub chord: f64,
    /// Angle of attack in degrees.
    pub aoa_deg: f64,
}

/// NACA 4-digit airfoil generator.
///
/// Parameters:
/// - m: maximum camber (first digit / 100)
/// - p: position of maximum camber (second digit / 10)
/// - t: maximum thickness ratio (last two digits / 100, e.g. 0.12 for NACA 0012)
#[derive(Debug, Clone, Copy)]
pub struct Naca4 {
    pub m: f64,
    pub p: f64,
    pub t: f64,
    pub chord: f64,
    pub aoa_deg: f64,
}

impl Naca4 {
    /// Create new NACA 4-digit airfoil specification.
    pub fn new(digit1: u32, digit2: u32, digit34: u32, chord: f64, aoa_deg: f64) -> Self {
        Self {
            m: (digit1 as f64) / 100.0,
            p: (digit2 as f64) / 10.0,
            t: (digit34 as f64) / 100.0,
            chord,
            aoa_deg,
        }
    }

    /// Convenience constructor for symmetric NACA 0012.
    pub fn naca0012(chord: f64, aoa_deg: f64) -> Self {
        Self::new(0, 0, 12, chord, aoa_deg)
    }

    /// Convenience constructor for cambered NACA 2412.
    pub fn naca2412(chord: f64, aoa_deg: f64) -> Self {
        Self::new(2, 4, 12, chord, aoa_deg)
    }

    /// Convenience constructor for high-lift NACA 4415.
    pub fn naca4415(chord: f64, aoa_deg: f64) -> Self {
        Self::new(4, 4, 15, chord, aoa_deg)
    }

    /// Generate discretized upper and lower surface coordinates.
    pub fn generate_coordinates(&self, n_points: usize) -> AirfoilCoordinates {
        let n = n_points.max(10);
        let mut x_upper = Vec::with_capacity(n);
        let mut y_upper = Vec::with_capacity(n);
        let mut x_lower = Vec::with_capacity(n);
        let mut y_lower = Vec::with_capacity(n);

        let alpha_rad = self.aoa_deg * PI / 180.0;
        let cos_a = alpha_rad.cos();
        let sin_a = alpha_rad.sin();

        // Quarter-chord moment reference center for rotation
        let x_qc = 0.25 * self.chord;

        for i in 0..n {
            // Cosine spacing clustered at leading and trailing edges
            let beta = (i as f64) * PI / ((n - 1) as f64);
            let xc = 0.5 * (1.0 - beta.cos());
            let x = xc * self.chord;

            // Half-thickness distribution
            let yt = 5.0 * self.t * self.chord
                * (0.2969 * xc.sqrt()
                    - 0.1260 * xc
                    - 0.3516 * xc.powi(2)
                    + 0.2843 * xc.powi(3)
                    - 0.1015 * xc.powi(4));

            // Camber line and gradient
            let (yc, dyc_dx) = if self.p > 1e-6 {
                if xc < self.p {
                    let yc = (self.m / (self.p * self.p)) * (2.0 * self.p * xc - xc * xc) * self.chord;
                    let dyc = (2.0 * self.m / (self.p * self.p)) * (self.p - xc);
                    (yc, dyc)
                } else {
                    let factor = (1.0 - self.p).powi(2);
                    let yc = (self.m / factor) * ((1.0 - 2.0 * self.p) + 2.0 * self.p * xc - xc * xc) * self.chord;
                    let dyc = (2.0 * self.m / factor) * (self.p - xc);
                    (yc, dyc)
                }
            } else {
                (0.0, 0.0)
            };

            let theta = dyc_dx.atan();
            let sin_th = theta.sin();
            let cos_th = theta.cos();

            // Unrotated surface coordinates
            let xu_raw = x - yt * sin_th;
            let yu_raw = yc + yt * cos_th;

            let xl_raw = x + yt * sin_th;
            let yl_raw = yc - yt * cos_th;

            // Rotate about quarter-chord by angle of attack
            let dx_u = xu_raw - x_qc;
            let dy_u = yu_raw;
            let xu_rot = x_qc + (dx_u * cos_a + dy_u * sin_a);
            let yu_rot = -dx_u * sin_a + dy_u * cos_a;

            let dx_l = xl_raw - x_qc;
            let dy_l = yl_raw;
            let xl_rot = x_qc + (dx_l * cos_a + dy_l * sin_a);
            let yl_rot = -dx_l * sin_a + dy_l * cos_a;

            x_upper.push(xu_rot);
            y_upper.push(yu_rot);
            x_lower.push(xl_rot);
            y_lower.push(yl_rot);
        }

        AirfoilCoordinates {
            x_upper,
            y_upper,
            x_lower,
            y_lower,
            chord: self.chord,
            aoa_deg: self.aoa_deg,
        }
    }

    /// Rasterize airfoil geometry into a 2D boolean mask [ny, nx].
    /// True indicates solid obstacle, False indicates fluid.
    pub fn rasterize(&self, nx: usize, ny: usize, origin_x: f64, origin_y: f64) -> Vec<bool> {
        let mut mask = vec![false; ny * nx];

        let alpha_rad = self.aoa_deg * PI / 180.0;
        let cos_a = alpha_rad.cos();
        let sin_a = (-alpha_rad).sin();
        let x_qc = origin_x + 0.25 * self.chord;
        let y_qc = origin_y;

        for iy in 0..ny {
            let y_world = iy as f64;
            for ix in 0..nx {
                let x_world = ix as f64;

                // Unrotate test point relative to quarter-chord
                let dx = x_world - x_qc;
                let dy = y_world - y_qc;
                let x_unrot = x_qc + (dx * cos_a - dy * sin_a);
                let y_unrot = y_qc + (dx * sin_a + dy * cos_a);

                // Relative position along chord [0, c]
                let x_rel = x_unrot - origin_x;
                if x_rel >= 0.0 && x_rel <= self.chord {
                    let xc = x_rel / self.chord;
                    let yt = 5.0 * self.t * self.chord
                        * (0.2969 * xc.sqrt()
                            - 0.1260 * xc
                            - 0.3516 * xc.powi(2)
                            + 0.2843 * xc.powi(3)
                            - 0.1015 * xc.powi(4));

                    let yc = if self.p > 1e-6 {
                        if xc < self.p {
                            (self.m / (self.p * self.p)) * (2.0 * self.p * xc - xc * xc) * self.chord
                        } else {
                            let factor = (1.0 - self.p).powi(2);
                            (self.m / factor) * ((1.0 - 2.0 * self.p) + 2.0 * self.p * xc - xc * xc) * self.chord
                        }
                    } else {
                        0.0
                    };

                    let y_local = y_unrot - origin_y;
                    if (y_local - yc).abs() <= yt {
                        mask[iy * nx + ix] = true;
                    }
                }
            }
        }

        mask
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_naca0012_symmetry() {
        let naca = Naca4::naca0012(100.0, 0.0);
        let coords = naca.generate_coordinates(50);

        for i in 0..50 {
            assert!((coords.x_upper[i] - coords.x_lower[i]).abs() < 1e-10);
            assert!((coords.y_upper[i] + coords.y_lower[i]).abs() < 1e-10);
        }

        // Maximum thickness at ~30% chord should be approx 12 units
        let max_thick = coords.y_upper.iter().cloned().fold(0.0, f64::max);
        assert!((max_thick - 6.0).abs() < 0.2);
    }

    #[test]
    fn test_rasterization_area() {
        let naca = Naca4::naca0012(50.0, 0.0);
        let mask = naca.rasterize(100, 60, 25.0, 30.0);
        let solid_cells = mask.iter().filter(|&&v| v).count();

        // Analytical area of NACA 0012: ~ 0.68 * t * c^2 = 0.68 * 0.12 * 2500 ~ 204 cells
        assert!(solid_cells > 150 && solid_cells < 260);
    }
}
