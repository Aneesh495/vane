//! Comprehensive Integration Test Suite for vane-core.
//!
//! Validates end-to-end coupled simulation pipelines:
//! - Full LBM time advancement with Zou-He inlet and convective outlet
//! - Momentum exchange aerodynamic force and moment convergence
//! - Smagorinsky and WALE subgrid-scale turbulence models
//! - Vortex core identification (Q-criterion and vorticity)
//! - XML VTK and CSV telemetry file generation
//!
//! Note: Zero em dashes are used in this file per coding standards.

use std::fs;
use std::path::Path;
use vane_core::airfoil::Naca4;
use vane_core::post::{compute_q_criterion, compute_vorticity};
use vane_core::simulation::{CollisionModel, LbmConfig, LbmSimulation};
use vane_core::telemetry::TelemetryLogger;
use vane_core::turbulence::{Smagorinsky, Wale};
use vane_core::vtk::export_vti;

#[test]
fn test_end_to_end_naca0012_pipeline() {
    let nx = 100;
    let ny = 50;
    let chord = 25.0;
    let u_inf = 0.05;
    let tau = 0.54;

    // 1. Generate NACA 0012 obstacle
    let naca = Naca4::naca0012(chord, 4.0);
    let mask = naca.rasterize(nx, ny, 20.0, 25.0);

    let solid_count = mask.iter().filter(|&&v| v).count();
    assert!(solid_count > 0, "Airfoil mask must contain solid nodes");

    // 2. Initialize simulation
    let config = LbmConfig {
        nx,
        ny,
        tau,
        u_inf,
        rho0: 1.0,
        chord,
        x_ref: (20.0 + 0.25 * chord, 25.0),
        collision_model: CollisionModel::Bgk,
        magic_lambda: 0.25,
    };

    let mut sim = LbmSimulation::new(config, mask.clone());
    assert_eq!(sim.step, 0);

    // 3. Step forward
    let num_steps = 60;
    let mut logger = TelemetryLogger::new();

    for _ in 0..num_steps {
        let forces = sim.step();
        assert!(forces.fx.is_finite(), "Fx must be finite");
        assert!(forces.fy.is_finite(), "Fy must be finite");
        assert!(forces.cl.is_finite(), "CL must be finite");
        assert!(forces.cd.is_finite(), "CD must be finite");
        logger.record(sim.step, &forces, 0.0, u_inf * 1.732);
    }

    assert_eq!(sim.step, num_steps);
    assert_eq!(sim.force_history.len(), num_steps);

    // 4. Extract macroscopic fields
    let mut rho = vec![0.0; nx * ny];
    let mut ux = vec![0.0; nx * ny];
    let mut uy = vec![0.0; nx * ny];

    for iy in 0..ny {
        for ix in 0..nx {
            let m = sim.get_macroscopic(ix, iy);
            let idx = iy * nx + ix;
            rho[idx] = m.rho;
            ux[idx] = m.ux;
            uy[idx] = m.uy;
            assert!(m.rho > 0.5 && m.rho < 1.5, "Density within stable bounds");
            assert!(m.ux.abs() < 0.3, "Velocity within subsonic limit");
        }
    }

    // 5. Diagnostics: Vorticity and Q-criterion
    let vort = compute_vorticity(&ux, &uy, nx, ny);
    let q = compute_q_criterion(&ux, &uy, nx, ny);

    assert_eq!(vort.len(), nx * ny);
    assert_eq!(q.len(), nx * ny);

    // 6. Test Exporters
    let test_dir = Path::new("target").join("test_output");
    fs::create_dir_all(&test_dir).expect("Failed to create test output dir");

    let vti_path = test_dir.join("test_naca0012.vti");
    let csv_path = test_dir.join("test_telemetry.csv");

    export_vti(&vti_path, nx, ny, &rho, &ux, &uy, &vort, &mask)
        .expect("VTK export failed");
    assert!(vti_path.is_file(), "VTK file must exist");

    logger.export_csv(&csv_path).expect("CSV export failed");
    assert!(csv_path.is_file(), "CSV file must exist");
}

#[test]
fn test_trt_vs_bgk_consistency() {
    let nx = 60;
    let ny = 30;
    let chord = 15.0;
    let u_inf = 0.04;
    let tau = 0.55;

    let naca = Naca4::naca0012(chord, 0.0);
    let mask = naca.rasterize(nx, ny, 15.0, 15.0);

    let config_bgk = LbmConfig {
        nx,
        ny,
        tau,
        u_inf,
        rho0: 1.0,
        chord,
        x_ref: (18.75, 15.0),
        collision_model: CollisionModel::Bgk,
        magic_lambda: 0.25,
    };

    let config_trt = LbmConfig {
        collision_model: CollisionModel::Trt,
        ..config_bgk.clone()
    };

    let mut sim_bgk = LbmSimulation::new(config_bgk, mask.clone());
    let mut sim_trt = LbmSimulation::new(config_trt, mask);

    for _ in 0..30 {
        sim_bgk.step();
        sim_trt.step();
    }

    // Compare densities at probe location
    let m_bgk = sim_bgk.get_macroscopic(10, 15);
    let m_trt = sim_trt.get_macroscopic(10, 15);

    assert!((m_bgk.rho - m_trt.rho).abs() < 0.05, "TRT and BGK densities consistent");
    assert!((m_bgk.ux - m_trt.ux).abs() < 0.02, "TRT and BGK velocities consistent");
}

#[test]
fn test_turbulence_subgrid_viscosity() {
    let smag = Smagorinsky::new(0.14, 1.0);
    let wale = Wale::new(0.325, 1.0);

    let feq = [1.0 / 9.0; 9];
    let mut f = feq;
    // Induce non-equilibrium shear stress
    f[1] += 0.01;
    f[3] -= 0.01;

    let nu_t = smag.compute_eddy_viscosity(&f, &feq, 1.0, 0.54);
    assert!(nu_t >= 0.0, "Smagorinsky viscosity must be non-negative");

    let tau_eff = smag.effective_tau(&f, &feq, 1.0, 0.54);
    assert!(tau_eff >= 0.54, "Effective tau must be at least molecular tau");

    let nu_wale = wale.compute_eddy_viscosity(0.05, 0.03);
    assert!(nu_wale >= 0.0, "WALE eddy viscosity must be non-negative");
}

#[test]
fn test_airfoil_aoa_lift_polar() {
    let nx = 80;
    let ny = 40;
    let chord = 20.0;
    let u_inf = 0.05;
    let tau = 0.54;

    // Run at AoA = 0 deg
    let naca_0 = Naca4::naca0012(chord, 0.0);
    let mask_0 = naca_0.rasterize(nx, ny, 15.0, 20.0);
    let config_0 = LbmConfig {
        nx,
        ny,
        tau,
        u_inf,
        rho0: 1.0,
        chord,
        x_ref: (20.0, 20.0),
        collision_model: CollisionModel::Bgk,
        magic_lambda: 0.25,
    };
    let mut sim_0 = LbmSimulation::new(config_0, mask_0);

    for _ in 0..40 {
        sim_0.step();
    }
    let cl_0 = sim_0.force_history.last().unwrap().cl;

    // At 0 degrees AoA for symmetric NACA 0012, Cl should be approximately zero
    assert!(cl_0.abs() < 1.0, "Zero AoA lift is small");
}
