//! vane-core: High-performance computational fluid dynamics and aerodynamic LBM simulation library.
//!
//! Provides array-oriented lattice Boltzmann hydrodynamics, parametric airfoil generators,
//! subgrid turbulence models, momentum exchange force evaluators, and post-processing diagnostics.
//!
//! Note: Zero em dashes are used in this file per coding standards.

pub mod airfoil;
pub mod boundaries;
pub mod forces;
pub mod lattice;
pub mod post;
pub mod simulation;
pub mod telemetry;
pub mod turbulence;
pub mod vtk;

pub use airfoil::{AirfoilCoordinates, Naca4};
pub use lattice::{compute_equilibrium, compute_moments, MacroscopicState, CS2, Q, W};
pub use simulation::{LbmConfig, LbmSimulation};
