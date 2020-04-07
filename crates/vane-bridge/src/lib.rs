//! C-compatible Foreign Function Interface (FFI) ABI bridge for vane.
//!
//! Exposes opaque handles and C-callable routines for embedding the vane
//! aerodynamic simulation engine into C, C++, Python (ctypes/cffi), and CBQN.
//!
//! Note: Zero em dashes are used in this file per coding standards.

use std::ffi::c_void;
use std::slice;
use vane_core::airfoil::Naca4;
use vane_core::forces::AerodynamicForces;
use vane_core::simulation::{CollisionModel, LbmConfig, LbmSimulation};

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

/// Create a new simulation instance.
///
/// Returns an opaque pointer to the heap-allocated simulation.
/// Ownership must be returned to `vane_simulation_free` when finished.
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

/// Advance simulation by one time step.
///
/// Returns the instantaneous aerodynamic forces.
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

/// Extract macroscopic density field into pre-allocated buffer of size nx * ny.
#[no_mangle]
pub unsafe extern "C" fn vane_simulation_get_density(
    sim_ptr: *mut c_void,
    out_ptr: *mut f64,
) -> i32 {
    if sim_ptr.is_null() || out_ptr.is_null() {
        return -1;
    }

    let sim = &*(sim_ptr as *mut LbmSimulation);
    let nx = sim.config.nx;
    let ny = sim.config.ny;
    let out_slice = slice::from_raw_parts_mut(out_ptr, nx * ny);

    for iy in 0..ny {
        for ix in 0..nx {
            let idx = iy * nx + ix;
            out_slice[idx] = sim.get_macroscopic(ix, iy).rho;
        }
    }

    0
}

/// Extract horizontal velocity ux field into pre-allocated buffer of size nx * ny.
#[no_mangle]
pub unsafe extern "C" fn vane_simulation_get_ux(
    sim_ptr: *mut c_void,
    out_ptr: *mut f64,
) -> i32 {
    if sim_ptr.is_null() || out_ptr.is_null() {
        return -1;
    }

    let sim = &*(sim_ptr as *mut LbmSimulation);
    let nx = sim.config.nx;
    let ny = sim.config.ny;
    let out_slice = slice::from_raw_parts_mut(out_ptr, nx * ny);

    for iy in 0..ny {
        for ix in 0..nx {
            let idx = iy * nx + ix;
            out_slice[idx] = sim.get_macroscopic(ix, iy).ux;
        }
    }

    0
}

/// Extract vertical velocity uy field into pre-allocated buffer of size nx * ny.
#[no_mangle]
pub unsafe extern "C" fn vane_simulation_get_uy(
    sim_ptr: *mut c_void,
    out_ptr: *mut f64,
) -> i32 {
    if sim_ptr.is_null() || out_ptr.is_null() {
        return -1;
    }

    let sim = &*(sim_ptr as *mut LbmSimulation);
    let nx = sim.config.nx;
    let ny = sim.config.ny;
    let out_slice = slice::from_raw_parts_mut(out_ptr, nx * ny);

    for iy in 0..ny {
        for ix in 0..nx {
            let idx = iy * nx + ix;
            out_slice[idx] = sim.get_macroscopic(ix, iy).uy;
        }
    }

    0
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

/// Deallocate the simulation instance.
#[no_mangle]
pub unsafe extern "C" fn vane_simulation_free(sim_ptr: *mut c_void) {
    if !sim_ptr.is_null() {
        let _ = Box::from_raw(sim_ptr as *mut LbmSimulation);
    }
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
}
