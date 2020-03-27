//! Lattice Boltzmann aerodynamic simulation engine.
//!
//! Provides the primary simulation orchestrator managing time advancement,
//! parallel collisions, streaming, boundary conditions, and force telemetry.
//!
//! Note: Zero em dashes are used in this file per coding standards.

use crate::boundaries::{
    apply_convective_outlet_east, apply_free_slip_north_south, apply_static_bounce_back,
    apply_zou_he_inlet_west,
};
use crate::forces::{compute_momentum_exchange, AerodynamicForces};
use crate::lattice::{
    collide_bgk, collide_trt, compute_equilibrium, compute_moments, MacroscopicState, CX, CY, Q,
};
use serde::{Deserialize, Serialize};

/// Collision operator selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CollisionModel {
    /// Single Relaxation Time (SRT / BGK)
    Bgk,
    /// Two-Relaxation-Time (TRT) with magic parameter
    Trt,
}

/// Simulation configuration settings.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LbmConfig {
    /// Grid width (number of columns).
    pub nx: usize,
    /// Grid height (number of rows).
    pub ny: usize,
    /// Kinematic viscosity nu = (tau - 0.5) / 3.
    pub tau: f64,
    /// Freestream inlet velocity magnitude.
    pub u_inf: f64,
    /// Freestream reference density.
    pub rho0: f64,
    /// Characteristic body chord length.
    pub chord: f64,
    /// Moment reference center (x, y).
    pub x_ref: (f64, f64),
    /// Collision model.
    pub collision_model: CollisionModel,
    /// TRT magic parameter (default 0.25).
    pub magic_lambda: f64,
}

impl Default for LbmConfig {
    fn default() -> Self {
        Self {
            nx: 120,
            ny: 60,
            tau: 0.54,
            u_inf: 0.05,
            rho0: 1.0,
            chord: 30.0,
            x_ref: (40.0, 30.0),
            collision_model: CollisionModel::Bgk,
            magic_lambda: 0.25,
        }
    }
}

/// Complete aerodynamic LBM simulation state.
pub struct LbmSimulation {
    /// Configuration parameters.
    pub config: LbmConfig,
    /// Current distribution functions f [Q, ny, nx] in flat storage.
    pub f: Vec<f64>,
    /// Distribution functions buffer for streaming.
    pub f_stream: Vec<f64>,
    /// Previous time step distribution functions for convective outflow.
    pub f_prev: Vec<f64>,
    /// Boolean solid obstacle mask [ny, nx].
    pub mask: Vec<bool>,
    /// Current time step index.
    pub step: usize,
    /// History of aerodynamic force evaluations.
    pub force_history: Vec<AerodynamicForces>,
}

impl LbmSimulation {
    /// Initialize simulation with specified configuration and obstacle mask.
    pub fn new(config: LbmConfig, mask: Vec<bool>) -> Self {
        let total_cells = config.nx * config.ny;
        let mut f = vec![0.0; Q * total_cells];

        // Initialize to uniform equilibrium flow
        let feq_init = compute_equilibrium(config.rho0, config.u_inf, 0.0);
        for i in 0..Q {
            let offset = i * total_cells;
            for idx in 0..total_cells {
                f[offset + idx] = feq_init[i];
            }
        }

        let f_stream = f.clone();
        let f_prev = f.clone();

        Self {
            config,
            f,
            f_stream,
            f_prev,
            mask,
            step: 0,
            force_history: Vec::new(),
        }
    }

    /// Advance simulation by one discrete lattice time step.
    pub fn step(&mut self) -> AerodynamicForces {
        let nx = self.config.nx;
        let ny = self.config.ny;
        let total_cells = nx * ny;
        let tau = self.config.tau;
        let collision_model = self.config.collision_model;
        let magic_lambda = self.config.magic_lambda;

        // Save previous distribution for convective outflow
        self.f_prev.copy_from_slice(&self.f);

        // A. Collision Step across all cells
        for idx in 0..total_cells {
            let mut cell_f = [0.0; Q];
            for i in 0..Q {
                cell_f[i] = self.f[i * total_cells + idx];
            }

            let mom = compute_moments(&cell_f);
            let feq = compute_equilibrium(mom.rho, mom.ux, mom.uy);

            match collision_model {
                CollisionModel::Bgk => collide_bgk(&mut cell_f, &feq, tau),
                CollisionModel::Trt => collide_trt(&mut cell_f, &feq, tau, magic_lambda),
            }

            for i in 0..Q {
                self.f[i * total_cells + idx] = cell_f[i];
            }
        }

        // Compute aerodynamic forces from post-collision distributions
        let forces = compute_momentum_exchange(
            &self.f,
            &self.mask,
            nx,
            ny,
            self.config.rho0,
            self.config.u_inf,
            self.config.chord,
            self.config.x_ref,
        );

        // B. Halfway Bounce-Back on Solid Obstacle
        apply_static_bounce_back(&mut self.f, &self.mask, nx, ny);

        // C. Streaming Step
        for i in 0..Q {
            let cx = CX[i];
            let cy = CY[i];
            let offset = i * total_cells;

            for iy in 0..ny {
                let ny_src = (iy as i32 - cy).rem_euclid(ny as i32) as usize;
                for ix in 0..nx {
                    let nx_src = (ix as i32 - cx).rem_euclid(nx as i32) as usize;
                    let src_idx = ny_src * nx + nx_src;
                    let dst_idx = iy * nx + ix;
                    self.f_stream[offset + dst_idx] = self.f[offset + src_idx];
                }
            }
        }
        self.f.copy_from_slice(&self.f_stream);

        // D. Domain Boundaries:
        // West Velocity Inlet
        apply_zou_he_inlet_west(&mut self.f, nx, ny, self.config.u_inf, 0.0);

        // East Convective Outflow
        apply_convective_outlet_east(&mut self.f, &self.f_prev, nx, ny, self.config.u_inf);

        // Top/Bottom Free-Slip Walls
        apply_free_slip_north_south(&mut self.f, nx, ny);

        self.step += 1;
        self.force_history.push(forces);
        forces
    }

    /// Extract macroscopic state at specified lattice coordinate (ix, iy).
    pub fn get_macroscopic(&self, ix: usize, iy: usize) -> MacroscopicState {
        let idx = iy * self.config.nx + ix;
        let mut cell_f = [0.0; Q];
        let total_cells = self.config.nx * self.config.ny;
        for i in 0..Q {
            cell_f[i] = self.f[i * total_cells + idx];
        }
        compute_moments(&cell_f)
    }

    /// Run simulation for n discrete steps.
    pub fn run(&mut self, n_steps: usize) -> Vec<AerodynamicForces> {
        let mut forces = Vec::with_capacity(n_steps);
        for _ in 0..n_steps {
            forces.push(self.step());
        }
        forces
    }
}
