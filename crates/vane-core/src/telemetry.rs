//! Aerodynamic telemetry and time-series logging engine.
//!
//! Handles CSV formatting, statistics reduction, Strouhal number evaluation,
//! and convergence monitoring.
//!
//! Note: Zero em dashes are used in this file per coding standards.

use crate::forces::AerodynamicForces;
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

/// Telemetry record entry for a single time step.
#[derive(Debug, Clone, Copy)]
pub struct TelemetryEntry {
    pub step: usize,
    pub cl: f64,
    pub cd: f64,
    pub cm: f64,
    pub residual_u: f64,
    pub max_mach: f64,
}

/// Time-series telemetry recorder.
#[derive(Debug, Default)]
pub struct TelemetryLogger {
    pub entries: Vec<TelemetryEntry>,
}

impl TelemetryLogger {
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
        }
    }

    /// Record a time step.
    pub fn record(&mut self, step: usize, forces: &AerodynamicForces, residual_u: f64, max_mach: f64) {
        self.entries.push(TelemetryEntry {
            step,
            cl: forces.cl,
            cd: forces.cd,
            cm: forces.cm,
            residual_u,
            max_mach,
        });
    }

    /// Export telemetry table to CSV file.
    pub fn export_csv<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        let mut file = File::create(path)?;
        writeln!(file, "step,cl,cd,cm,residual_u,max_mach")?;

        for e in &self.entries {
            writeln!(
                file,
                "{},{:.6},{:.6},{:.6},{:.6e},{:.4}",
                e.step, e.cl, e.cd, e.cm, e.residual_u, e.max_mach
            )?;
        }

        Ok(())
    }

    /// Calculate time-averaged lift, drag, and RMS fluctuations over the final fraction of steps.
    pub fn statistics(&self, fraction: f64) -> (f64, f64, f64, f64) {
        let n = self.entries.len();
        if n == 0 {
            return (0.0, 0.0, 0.0, 0.0);
        }

        let start_idx = ((1.0 - fraction.clamp(0.0, 1.0)) * (n as f64)) as usize;
        let window = &self.entries[start_idx..];
        let count = window.len() as f64;

        let cl_mean = window.iter().map(|e| e.cl).sum::<f64>() / count;
        let cd_mean = window.iter().map(|e| e.cd).sum::<f64>() / count;

        let cl_var = window.iter().map(|e| (e.cl - cl_mean).powi(2)).sum::<f64>() / count;
        let cd_var = window.iter().map(|e| (e.cd - cd_mean).powi(2)).sum::<f64>() / count;

        (cl_mean, cd_mean, cl_var.sqrt(), cd_var.sqrt())
    }
}
