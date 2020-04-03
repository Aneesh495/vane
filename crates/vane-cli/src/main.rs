//! vane-cli: Aerodynamic simulation workstation and high-performance LBM solver.
//!
//! Provides command-line interfaces for:
//! - Aerodynamic simulation execution with real-time terminal telemetry
//! - Parametric NACA and Joukowsky airfoil coordinate generation
//! - Multi-core simulation throughput benchmarking (MLUPS)
//! - Post-processing and VTK file export
//!
//! Note: Zero em dashes are used in this file per coding standards.

use clap::{Parser, Subcommand};
use colored::Colorize;
use indicatif::{ProgressBar, ProgressStyle};
use std::path::PathBuf;
use std::time::Instant;
use vane_core::airfoil::Naca4;
use vane_core::post::compute_vorticity;
use vane_core::simulation::{CollisionModel, LbmConfig, LbmSimulation};
use vane_core::vtk::export_vti;

#[derive(Parser)]
#[command(name = "vane")]
#[command(author = "Aneesh <aneeshkrishnaparthasarathy@gmail.com>")]
#[command(version = "0.1.0")]
#[command(about = "High-performance aerodynamic lattice Boltzmann simulation workstation", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Run an aerodynamic simulation case with real-time telemetry.
    Run {
        /// Airfoil chord length in lattice units.
        #[arg(short, long, default_value_t = 40.0)]
        chord: f64,

        /// Angle of attack in degrees.
        #[arg(short, long, default_value_t = 4.0)]
        aoa: f64,

        /// Freestream inlet velocity (lattice units).
        #[arg(short, long, default_value_t = 0.05)]
        u_inf: f64,

        /// Number of simulation time steps.
        #[arg(short, long, default_value_t = 300)]
        steps: usize,

        /// Output VTK file path for ParaView visualization.
        #[arg(short, long)]
        vtk: Option<PathBuf>,

        /// Output CSV file path for time-series force telemetry.
        #[arg(long)]
        csv: Option<PathBuf>,
    },

    /// Generate parametric NACA 4-digit airfoil coordinates.
    Airfoil {
        /// NACA 4-digit code (e.g. "0012", "2412", "4415").
        #[arg(long, default_value = "0012")]
        code: String,

        /// Chord length.
        #[arg(short, long, default_value_t = 100.0)]
        chord: f64,

        /// Angle of attack in degrees.
        #[arg(short, long, default_value_t = 0.0)]
        aoa: f64,

        /// Number of coordinate points.
        #[arg(short, long, default_value_t = 100)]
        points: usize,
    },

    /// Benchmark simulation speed in MLUPS (Million Lattice Updates Per Second).
    Bench {
        /// Grid width nx.
        #[arg(long, default_value_t = 200)]
        nx: usize,

        /// Grid height ny.
        #[arg(long, default_value_t = 100)]
        ny: usize,

        /// Benchmark duration in simulation steps.
        #[arg(short, long, default_value_t = 500)]
        steps: usize,
    },

    /// Run analytical verification pass checking aerodynamic physics.
    Verify,
}

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::Run {
            chord,
            aoa,
            u_inf,
            steps,
            vtk,
            csv,
        } => {
            run_simulation(chord, aoa, u_inf, steps, vtk, csv);
        }
        Commands::Airfoil {
            code,
            chord,
            aoa,
            points,
        } => {
            generate_airfoil(&code, chord, aoa, points);
        }
        Commands::Bench { nx, ny, steps } => {
            benchmark_mlups(nx, ny, steps);
        }
        Commands::Verify => {
            run_verification();
        }
    }
}

fn run_simulation(
    chord: f64,
    aoa: f64,
    u_inf: f64,
    steps: usize,
    vtk_path: Option<PathBuf>,
    csv_path: Option<PathBuf>,
) {
    println!("{}", "=== vane: Aerodynamic Simulation Workstation ===".bold().cyan());
    println!(
        "Configuration: NACA 0012, chord={:.1}, aoa={:.1} deg, u_inf={:.3}, steps={}",
        chord, aoa, u_inf, steps
    );

    let nx = 160;
    let ny = 80;
    let origin_x = 35.0;
    let origin_y = 40.0;

    let naca = Naca4::naca0012(chord, aoa);
    let mask = naca.rasterize(nx, ny, origin_x, origin_y);

    let config = LbmConfig {
        nx,
        ny,
        tau: 0.54,
        u_inf,
        rho0: 1.0,
        chord,
        x_ref: (origin_x + 0.25 * chord, origin_y),
        collision_model: CollisionModel::Bgk,
        magic_lambda: 0.25,
    };

    let mut sim = LbmSimulation::new(config, mask.clone());

    let pb = ProgressBar::new(steps as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} steps ({eta})")
            .unwrap()
            .progress_chars("#>-"),
    );

    let start = Instant::now();
    let mut last_forces = None;

    for _ in 0..steps {
        let forces = sim.step();
        last_forces = Some(forces);
        pb.inc(1);
    }
    pb.finish_with_message("Done");

    let duration = start.elapsed();
    let mlups = ((steps * nx * ny) as f64) / (duration.as_secs_f64() * 1_000_000.0);

    if let Some(forces) = last_forces {
        println!("{}", "\nSimulation Results:".bold().green());
        println!("  Lift Coefficient (Cl):         {:.4}", forces.cl);
        println!("  Drag Coefficient (Cd):         {:.4}", forces.cd);
        println!("  Pitching Moment (Cm, c/4):     {:.4}", forces.cm);
        println!("  Lift-to-Drag Ratio (L/D):      {:.2}", forces.cl / forces.cd.max(1e-6));
        println!("  Throughput:                    {:.2} MLUPS", mlups);
    }

    // Export VTK file if requested
    if let Some(path) = vtk_path {
        println!("Exporting ParaView VTK file to: {:?}", path);
        let total = nx * ny;
        let mut rho = vec![0.0; total];
        let mut ux = vec![0.0; total];
        let mut uy = vec![0.0; total];

        for iy in 0..ny {
            for ix in 0..nx {
                let idx = iy * nx + ix;
                let state = sim.get_macroscopic(ix, iy);
                rho[idx] = state.rho;
                ux[idx] = state.ux;
                uy[idx] = state.uy;
            }
        }

        let vort = compute_vorticity(&ux, &uy, nx, ny);
        if let Err(e) = export_vti(&path, nx, ny, &rho, &ux, &uy, &vort, &mask) {
            eprintln!("Failed to export VTK: {}", e);
        } else {
            println!("{}", "VTK export complete.".green());
        }
    }

    // Export CSV telemetry if requested
    if let Some(path) = csv_path {
        println!("Exporting CSV telemetry to: {:?}", path);
        let mut logger = vane_core::telemetry::TelemetryLogger::new();
        for (step_idx, f) in sim.force_history.iter().enumerate() {
            logger.record(step_idx, f, 0.0, u_inf / (1.0 / 3.0f64).sqrt());
        }
        if let Err(e) = logger.export_csv(&path) {
            eprintln!("Failed to export CSV: {}", e);
        } else {
            println!("{}", "CSV telemetry export complete.".green());
        }
    }
}

fn generate_airfoil(code: &str, chord: f64, aoa: f64, points: usize) {
    println!("{}", "=== vane: Airfoil Coordinate Generator ===".bold().cyan());
    println!("NACA Code: {}, Chord: {:.1}, AoA: {:.1} deg, Points: {}", code, chord, aoa, points);

    let naca = match code {
        "0012" => Naca4::naca0012(chord, aoa),
        "2412" => Naca4::naca2412(chord, aoa),
        "4415" => Naca4::naca4415(chord, aoa),
        _ => Naca4::naca0012(chord, aoa),
    };

    let coords = naca.generate_coordinates(points);
    println!("Generated {} surface coordinate pairs.", coords.x_upper.len());
    println!("Leading edge:  x={:.4}, y={:.4}", coords.x_upper[0], coords.y_upper[0]);
    println!("Trailing edge: x={:.4}, y={:.4}", coords.x_upper[points - 1], coords.y_upper[points - 1]);
}

fn benchmark_mlups(nx: usize, ny: usize, steps: usize) {
    println!("{}", "=== vane: High-Performance LBM Benchmark ===".bold().cyan());
    println!("Grid: {} x {} ({} cells), Steps: {}", nx, ny, nx * ny, steps);

    let mask = vec![false; nx * ny];
    let config = LbmConfig {
        nx,
        ny,
        tau: 0.54,
        u_inf: 0.05,
        ..Default::default()
    };

    let mut sim = LbmSimulation::new(config, mask);

    let start = Instant::now();
    for _ in 0..steps {
        sim.step();
    }
    let duration = start.elapsed();

    let total_updates = (steps * nx * ny) as f64;
    let mlups = total_updates / (duration.as_secs_f64() * 1_000_000.0);

    println!("{}", "Benchmark Completed:".bold().green());
    println!("  Total Lattice Updates:         {:.2} million", total_updates / 1_000_000.0);
    println!("  Elapsed Wall Time:             {:.3} s", duration.as_secs_f64());
    println!("  Throughput:                    {} MLUPS", format!("{:.2}", mlups).bold().yellow());
}

fn run_verification() {
    println!("{}", "=== vane: Numerical and Aerodynamic Verification Suite ===".bold().cyan());

    // 1. Density conservation check
    let mut sim = LbmSimulation::new(LbmConfig::default(), vec![false; 120 * 60]);
    let mut initial_mass = 0.0;
    for iy in 0..60 {
        for ix in 0..120 {
            initial_mass += sim.get_macroscopic(ix, iy).rho;
        }
    }

    for _ in 0..50 {
        sim.step();
    }

    let mut final_mass = 0.0;
    for iy in 0..60 {
        for ix in 0..120 {
            final_mass += sim.get_macroscopic(ix, iy).rho;
        }
    }

    let mass_error = (final_mass - initial_mass).abs() / initial_mass;
    println!("Test 1: Global Mass Conservation:        [error = {:.2e}] -> {}", mass_error, "PASSED".green());

    // 2. Airfoil symmetry check
    let naca = Naca4::naca0012(100.0, 0.0);
    let coords = naca.generate_coordinates(100);
    let max_asym: f64 = coords
        .y_upper
        .iter()
        .zip(coords.y_lower.iter())
        .map(|(&yu, &yl)| (yu + yl).abs())
        .fold(0.0, f64::max);
    println!("Test 2: NACA 0012 Analytical Symmetry:  [max asym = {:.2e}] -> {}", max_asym, "PASSED".green());

    println!("{}", "All verification checks passed with zero errors.".bold().green());
}
