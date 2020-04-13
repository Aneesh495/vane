//! C-compatible Foreign Function Interface (FFI) ABI bridge for vane.
//!
//! Exposes opaque handles and C-callable routines for embedding the vane
//! aerodynamic simulation engine into C, C++, Python (ctypes/cffi), and CBQN.
//!
//! Note: Zero em dashes are used in this file per coding standards.

use std::ffi::{c_char, c_void, CStr};
use std::slice;
use vane_core::airfoil::Naca4;
use vane_core::forces::AerodynamicForces;
use vane_core::post::{compute_q_criterion, compute_vorticity};
use vane_core::simulation::{CollisionModel, LbmConfig, LbmSimulation};
use vane_core::telemetry::TelemetryLogger;
use vane_core::vtk::export_vti;

static VERSION_CSTR: &[u8] = b"0.1.0\0";

/// Query vane native engine version string.
#[no_mangle]
pub extern "C" fn vane_version() -> *const c_char {
    VERSION_CSTR.as_ptr() as *const c_char
}

/// C-compatible representation of aerodynamic forces.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct CFources {
    pub fx: f64,
    pub fy: f64,
    pub moment: f64,
    pub cd: f64,
    pub cl: f64,
    pub cm: f64,
}

impl From<AerodynamicForces> for CFources {
    fn from(f: AerodynamicForces) -> Self {
        Self {
            fx: f.fx,
            fy: f.fy,
            moment: f.moment,
            cd: f.cd,
            cl: f.cl,
            cm: f.cm,
        }
    }
}

/// Create a new simulation instance with explicit obstacle mask.
#[no_mangle]
pub unsafe extern "C" fn vane_simulation_create(
    nx: usize,
    ny: usize,
    tau: f64,
    u_inf: f64,
    chord: f64,
    mask_ptr: *const u8,
) -> *mut c_void {
    if mask_ptr.is_null() || nx == 0 || ny == 0 {
        return std::ptr::null_mut();
    }

    let mask_slice = slice::from_raw_parts(mask_ptr, nx * ny);
    let mask_vec: Vec<bool> = mask_slice.iter().map(|&v| v != 0).collect();

    let config = LbmConfig {
        nx,
        ny,
        tau,
        u_inf,
        rho0: 1.0,
        chord,
        x_ref: (0.25 * (nx as f64) + 0.25 * chord, 0.5 * (ny as f64)),
        collision_model: CollisionModel::Bgk,
        magic_lambda: 0.25,
    };

    let sim = Box::new(LbmSimulation::new(config, mask_vec));
    Box::into_raw(sim) as *mut c_void
}

/// Create a new simulation instance with default fluid domain.
#[no_mangle]
pub unsafe extern "C" fn vane_sim_create(
    nx: usize,
    ny: usize,
    tau: f64,
    u_inlet: f64,
    collision_type: u32,
) -> *mut c_void {
    if nx == 0 || ny == 0 {
        return std::ptr::null_mut();
    }

    let collision_model = match collision_type {
        1 => CollisionModel::Trt,
        _ => CollisionModel::Bgk,
    };

    let config = LbmConfig {
        nx,
        ny,
        tau,
        u_inf: u_inlet,
        rho0: 1.0,
        chord: 30.0,
        x_ref: (0.25 * (nx as f64) + 7.5, 0.5 * (ny as f64)),
        collision_model,
        magic_lambda: 0.25,
    };

    let mask = vec![false; nx * ny];
    let sim = Box::new(LbmSimulation::new(config, mask));
    Box::into_raw(sim) as *mut c_void
}

/// Advance simulation by one time step, returning aerodynamic forces.
#[no_mangle]
pub unsafe extern "C" fn vane_simulation_step(sim_ptr: *mut c_void) -> CFources {
    if sim_ptr.is_null() {
        return CFources {
            fx: 0.0,
            fy: 0.0,
            moment: 0.0,
            cd: 0.0,
            cl: 0.0,
            cm: 0.0,
        };
    }

    let sim = &mut *(sim_ptr as *mut LbmSimulation);
    let forces = sim.step();
    forces.into()
}

/// Advance simulation by N discrete time steps.
#[no_mangle]
pub unsafe extern "C" fn vane_sim_step(sim_ptr: *mut c_void, count: usize) -> i32 {
    if sim_ptr.is_null() {
        return -1;
    }

    let sim = &mut *(sim_ptr as *mut LbmSimulation);
    for _ in 0..count {
        sim.step();
    }
    0
}

/// Rasterize and embed a NACA 4-digit airfoil into an active simulation.
#[no_mangle]
pub unsafe extern "C" fn vane_sim_set_naca(
    sim_ptr: *mut c_void,
    code: *const c_char,
    chord: f64,
    aoa_deg: f64,
    x_lead: f64,
    y_lead: f64,
) -> i32 {
    if sim_ptr.is_null() || code.is_null() {
        return -1;
    }

    let sim = &mut *(sim_ptr as *mut LbmSimulation);
    let code_str = match CStr::from_ptr(code).to_str() {
        Ok(s) => s,
        Err(_) => return -2,
    };

    if code_str.len() != 4 {
        return -3;
    }

    let digits: Vec<u32> = code_str.chars().filter_map(|c| c.to_digit(10)).collect();
    if digits.len() != 4 {
        return -4;
    }

    let d1 = digits[0];
    let d2 = digits[1];
    let d34 = digits[2] * 10 + digits[3];

    let naca = Naca4::new(d1, d2, d34, chord, aoa_deg);
    let mask = naca.rasterize(sim.config.nx, sim.config.ny, x_lead, y_lead);
    sim.mask = mask;
    sim.config.chord = chord;
    sim.config.x_ref = (x_lead + 0.25 * chord, y_lead);

    0
}

/// Retrieve the most recent aerodynamic force coefficients.
#[no_mangle]
pub unsafe extern "C" fn vane_sim_get_forces(
    sim_ptr: *const c_void,
    cl: *mut f64,
    cd: *mut f64,
    cm: *mut f64,
) -> i32 {
    if sim_ptr.is_null() {
        return -1;
    }

    let sim = &*(sim_ptr as *const LbmSimulation);
    if let Some(latest) = sim.force_history.last() {
        if !cl.is_null() {
            *cl = latest.cl;
        }
        if !cd.is_null() {
            *cd = latest.cd;
        }
        if !cm.is_null() {
            *cm = latest.cm;
        }
    } else {
        if !cl.is_null() {
            *cl = 0.0;
        }
        if !cd.is_null() {
            *cd = 0.0;
        }
        if !cm.is_null() {
            *cm = 0.0;
        }
    }

    0
}

/// Query simulation grid dimensions and elapsed step count.
#[no_mangle]
pub unsafe extern "C" fn vane_sim_get_dimensions(
    sim_ptr: *const c_void,
    nx: *mut usize,
    ny: *mut usize,
    step: *mut usize,
) -> i32 {
    if sim_ptr.is_null() {
        return -1;
    }

    let sim = &*(sim_ptr as *const LbmSimulation);
    if !nx.is_null() {
        *nx = sim.config.nx;
    }
    if !ny.is_null() {
        *ny = sim.config.ny;
    }
    if !step.is_null() {
        *step = sim.step;
    }

    0
}

/// Query binary obstacle mask into a pre-allocated byte buffer.
#[no_mangle]
pub unsafe extern "C" fn vane_sim_get_solid_mask(
    sim_ptr: *const c_void,
    out_ptr: *mut u8,
    len: usize,
) -> i32 {
    if sim_ptr.is_null() || out_ptr.is_null() {
        return -1;
    }

    let sim = &*(sim_ptr as *const LbmSimulation);
    let total = sim.config.nx * sim.config.ny;
    if len < total {
        return -2;
    }

    let out_slice = slice::from_raw_parts_mut(out_ptr, total);
    for i in 0..total {
        out_slice[i] = if sim.mask[i] { 1 } else { 0 };
    }

    0
}

/// Extract macroscopic physical field into pre-allocated buffer of size nx * ny.
/// Field IDs: 0 = density, 1 = ux, 2 = uy, 3 = speed, 4 = vorticity, 5 = q_criterion.
#[no_mangle]
pub unsafe extern "C" fn vane_sim_get_field(
    sim_ptr: *const c_void,
    field_id: u32,
    out_ptr: *mut f64,
    len: usize,
) -> i32 {
    if sim_ptr.is_null() || out_ptr.is_null() {
        return -1;
    }

    let sim = &*(sim_ptr as *const LbmSimulation);
    let nx = sim.config.nx;
    let ny = sim.config.ny;
    let total = nx * ny;
    if len < total {
        return -2;
    }

    let out_slice = slice::from_raw_parts_mut(out_ptr, total);

    match field_id {
        0 => {
            for iy in 0..ny {
                for ix in 0..nx {
                    out_slice[iy * nx + ix] = sim.get_macroscopic(ix, iy).rho;
                }
            }
        }
        1 => {
            for iy in 0..ny {
                for ix in 0..nx {
                    out_slice[iy * nx + ix] = sim.get_macroscopic(ix, iy).ux;
                }
            }
        }
        2 => {
            for iy in 0..ny {
                for ix in 0..nx {
                    out_slice[iy * nx + ix] = sim.get_macroscopic(ix, iy).uy;
                }
            }
        }
        3 => {
            for iy in 0..ny {
                for ix in 0..nx {
                    let m = sim.get_macroscopic(ix, iy);
                    out_slice[iy * nx + ix] = (m.ux * m.ux + m.uy * m.uy).sqrt();
                }
            }
        }
        4 => {
            let mut ux = vec![0.0; total];
            let mut uy = vec![0.0; total];
            for iy in 0..ny {
                for ix in 0..nx {
                    let m = sim.get_macroscopic(ix, iy);
                    ux[iy * nx + ix] = m.ux;
                    uy[iy * nx + ix] = m.uy;
                }
            }
            let vort = compute_vorticity(&ux, &uy, nx, ny);
            out_slice.copy_from_slice(&vort);
        }
        5 => {
            let mut ux = vec![0.0; total];
            let mut uy = vec![0.0; total];
            for iy in 0..ny {
                for ix in 0..nx {
                    let m = sim.get_macroscopic(ix, iy);
                    ux[iy * nx + ix] = m.ux;
                    uy[iy * nx + ix] = m.uy;
                }
            }
            let q = compute_q_criterion(&ux, &uy, nx, ny);
            out_slice.copy_from_slice(&q);
        }
        _ => return -3,
    }

    0
}

/// Extract macroscopic density field into pre-allocated buffer of size nx * ny.
#[no_mangle]
pub unsafe extern "C" fn vane_simulation_get_density(
    sim_ptr: *mut c_void,
    out_ptr: *mut f64,
) -> i32 {
    vane_sim_get_field(sim_ptr, 0, out_ptr, usize::MAX)
}

/// Extract horizontal velocity ux field into pre-allocated buffer of size nx * ny.
#[no_mangle]
pub unsafe extern "C" fn vane_simulation_get_ux(
    sim_ptr: *mut c_void,
    out_ptr: *mut f64,
) -> i32 {
    vane_sim_get_field(sim_ptr, 1, out_ptr, usize::MAX)
}

/// Extract vertical velocity uy field into pre-allocated buffer of size nx * ny.
#[no_mangle]
pub unsafe extern "C" fn vane_simulation_get_uy(
    sim_ptr: *mut c_void,
    out_ptr: *mut f64,
) -> i32 {
    vane_sim_get_field(sim_ptr, 2, out_ptr, usize::MAX)
}

/// Rasterize a NACA 4-digit airfoil into a pre-allocated byte buffer of size nx * ny.
#[no_mangle]
pub unsafe extern "C" fn vane_airfoil_rasterize(
    digit1: u32,
    digit2: u32,
    digit34: u32,
    chord: f64,
    aoa_deg: f64,
    nx: usize,
    ny: usize,
    origin_x: f64,
    origin_y: f64,
    out_mask: *mut u8,
) -> i32 {
    if out_mask.is_null() || nx == 0 || ny == 0 {
        return -1;
    }

    let naca = Naca4::new(digit1, digit2, digit34, chord, aoa_deg);
    let mask = naca.rasterize(nx, ny, origin_x, origin_y);
    let out_slice = slice::from_raw_parts_mut(out_mask, nx * ny);

    for i in 0..(nx * ny) {
        out_slice[i] = if mask[i] { 1 } else { 0 };
    }

    0
}

/// Export simulation macroscopic fields to XML VTK (.vti) file.
#[no_mangle]
pub unsafe extern "C" fn vane_sim_save_vtk(
    sim_ptr: *const c_void,
    path: *const c_char,
) -> i32 {
    if sim_ptr.is_null() || path.is_null() {
        return -1;
    }

    let sim = &*(sim_ptr as *const LbmSimulation);
    let path_str = match CStr::from_ptr(path).to_str() {
        Ok(s) => s,
        Err(_) => return -2,
    };

    let nx = sim.config.nx;
    let ny = sim.config.ny;
    let total = nx * ny;

    let mut rho = vec![0.0; total];
    let mut ux = vec![0.0; total];
    let mut uy = vec![0.0; total];

    for iy in 0..ny {
        for ix in 0..nx {
            let m = sim.get_macroscopic(ix, iy);
            let idx = iy * nx + ix;
            rho[idx] = m.rho;
            ux[idx] = m.ux;
            uy[idx] = m.uy;
        }
    }

    let vort = compute_vorticity(&ux, &uy, nx, ny);

    match export_vti(path_str, nx, ny, &rho, &ux, &uy, &vort, &sim.mask) {
        Ok(_) => 0,
        Err(_) => -3,
    }
}

/// Export aerodynamic force telemetry history to CSV file.
#[no_mangle]
pub unsafe extern "C" fn vane_sim_save_csv(
    sim_ptr: *const c_void,
    path: *const c_char,
) -> i32 {
    if sim_ptr.is_null() || path.is_null() {
        return -1;
    }

    let sim = &*(sim_ptr as *const LbmSimulation);
    let path_str = match CStr::from_ptr(path).to_str() {
        Ok(s) => s,
        Err(_) => return -2,
    };

    let mut logger = TelemetryLogger::new();
    for (step_idx, forces) in sim.force_history.iter().enumerate() {
        logger.record(step_idx + 1, forces, 0.0, sim.config.u_inf * 1.732);
    }

    match logger.export_csv(path_str) {
        Ok(_) => 0,
        Err(_) => -3,
    }
}

/// Deallocate the simulation instance.
#[no_mangle]
pub unsafe extern "C" fn vane_simulation_free(sim_ptr: *mut c_void) {
    if !sim_ptr.is_null() {
        let _ = Box::from_raw(sim_ptr as *mut LbmSimulation);
    }
}

/// Deallocate the simulation instance (alias).
#[no_mangle]
pub unsafe extern "C" fn vane_sim_destroy(sim_ptr: *mut c_void) {
    vane_simulation_free(sim_ptr);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c_ffi_lifecycle() {
        unsafe {
            let nx = 40;
            let ny = 20;
            let mask = vec![0u8; nx * ny];

            let sim_ptr = vane_simulation_create(nx, ny, 0.55, 0.04, 15.0, mask.as_ptr());
            assert!(!sim_ptr.is_null());

            let forces = vane_simulation_step(sim_ptr);
            assert!(forces.fx.is_finite());
            assert!(forces.fy.is_finite());

            let mut rho = vec![0.0; nx * ny];
            let res = vane_simulation_get_density(sim_ptr, rho.as_mut_ptr());
            assert_eq!(res, 0);
            assert!((rho[0] - 1.0).abs() < 0.1);

            vane_simulation_free(sim_ptr);
        }
    }

    #[test]
    fn test_vane_sim_extended_lifecycle() {
        unsafe {
            let sim_ptr = vane_sim_create(50, 30, 0.54, 0.05, 0);
            assert!(!sim_ptr.is_null());

            let naca_code = std::ffi::CString::new("0012").unwrap();
            let rc = vane_sim_set_naca(sim_ptr, naca_code.as_ptr(), 20.0, 3.0, 10.0, 15.0);
            assert_eq!(rc, 0);

            assert_eq!(vane_sim_step(sim_ptr, 10), 0);

            let mut cl = 0.0;
            let mut cd = 0.0;
            let mut cm = 0.0;
            assert_eq!(vane_sim_get_forces(sim_ptr, &mut cl, &mut cd, &mut cm), 0);
            assert!(cl.is_finite());
            assert!(cd.is_finite());

            let mut speed = vec![0.0; 50 * 30];
            assert_eq!(vane_sim_get_field(sim_ptr, 3, speed.as_mut_ptr(), speed.len()), 0);

            vane_sim_destroy(sim_ptr);
        }
    }
}
