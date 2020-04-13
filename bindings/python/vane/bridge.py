"""
Low-level ctypes bindings loader for vane-bridge native shared library.
"""

import ctypes
import os
import sys
from pathlib import Path


def find_library_path() -> Path:
    """
    Search for libvane_bridge shared library across build artifacts and system paths.
    """
    system = sys.platform
    if system == "darwin":
        lib_name = "libvane_bridge.dylib"
    elif system.startswith("win"):
        lib_name = "vane_bridge.dll"
    else:
        lib_name = "libvane_bridge.so"

    search_roots = [
        # Relative to current file: in repository target directories
        Path(__file__).resolve().parent.parent.parent.parent / "target" / "release",
        Path(__file__).resolve().parent.parent.parent.parent / "target" / "debug",
        # Package local directory
        Path(__file__).resolve().parent,
        # Environment override
        Path(os.environ.get("VANE_LIB_DIR", ".")),
    ]

    for search_dir in search_roots:
        candidate = search_dir / lib_name
        if candidate.is_file():
            return candidate

    raise FileNotFoundError(
        f"Could not locate {lib_name}. Please build vane-bridge with "
        "`cargo build --release -p vane-bridge` first."
    )


def load_vane_bridge(lib_path: str | Path | None = None) -> ctypes.CDLL:
    """
    Load and initialize ctypes interface for vane-bridge shared library.
    """
    path = Path(lib_path) if lib_path else find_library_path()
    lib = ctypes.CDLL(str(path))

    # Opaque pointer handle
    handle_p = ctypes.c_void_p

    # vane_version
    lib.vane_version.restype = ctypes.c_char_p
    lib.vane_version.argtypes = []

    # vane_sim_create
    lib.vane_sim_create.restype = handle_p
    lib.vane_sim_create.argtypes = [
        ctypes.c_size_t,  # nx
        ctypes.c_size_t,  # ny
        ctypes.c_double,  # tau
        ctypes.c_double,  # u_inlet
        ctypes.c_uint32,  # collision_type
    ]

    # vane_sim_destroy
    lib.vane_sim_destroy.restype = None
    lib.vane_sim_destroy.argtypes = [handle_p]

    # vane_sim_set_naca
    lib.vane_sim_set_naca.restype = ctypes.c_int32
    lib.vane_sim_set_naca.argtypes = [
        handle_p,
        ctypes.c_char_p,  # code
        ctypes.c_double,  # chord
        ctypes.c_double,  # aoa_deg
        ctypes.c_double,  # x_lead
        ctypes.c_double,  # y_lead
    ]

    # vane_sim_step
    lib.vane_sim_step.restype = ctypes.c_int32
    lib.vane_sim_step.argtypes = [handle_p, ctypes.c_size_t]

    # vane_sim_get_forces
    lib.vane_sim_get_forces.restype = ctypes.c_int32
    lib.vane_sim_get_forces.argtypes = [
        handle_p,
        ctypes.POINTER(ctypes.c_double),  # cl
        ctypes.POINTER(ctypes.c_double),  # cd
        ctypes.POINTER(ctypes.c_double),  # cm
    ]

    # vane_sim_get_dimensions
    lib.vane_sim_get_dimensions.restype = ctypes.c_int32
    lib.vane_sim_get_dimensions.argtypes = [
        handle_p,
        ctypes.POINTER(ctypes.c_size_t),  # nx
        ctypes.POINTER(ctypes.c_size_t),  # ny
        ctypes.POINTER(ctypes.c_size_t),  # step
    ]

    # vane_sim_get_field
    lib.vane_sim_get_field.restype = ctypes.c_int32
    lib.vane_sim_get_field.argtypes = [
        handle_p,
        ctypes.c_uint32,                  # field_id
        ctypes.POINTER(ctypes.c_double),  # out_buffer
        ctypes.c_size_t,                  # buffer_len
    ]

    # vane_sim_get_solid_mask
    lib.vane_sim_get_solid_mask.restype = ctypes.c_int32
    lib.vane_sim_get_solid_mask.argtypes = [
        handle_p,
        ctypes.POINTER(ctypes.c_uint8),   # out_buffer
        ctypes.c_size_t,                  # buffer_len
    ]

    # vane_sim_save_vtk
    lib.vane_sim_save_vtk.restype = ctypes.c_int32
    lib.vane_sim_save_vtk.argtypes = [handle_p, ctypes.c_char_p]

    # vane_sim_save_csv
    lib.vane_sim_save_csv.restype = ctypes.c_int32
    lib.vane_sim_save_csv.argtypes = [handle_p, ctypes.c_char_p]

    return lib
