"""
High-level object-oriented Python interface for vane CFD simulations.
"""

import ctypes
from pathlib import Path
from typing import Dict, Optional, Tuple, Union

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

from .bridge import load_vane_bridge


FIELD_IDS = {
    "density": 0,
    "ux": 1,
    "uy": 2,
    "speed": 3,
    "vorticity": 4,
    "q_criterion": 5,
}


class VaneSimulation:
    """
    Object-oriented wrapper around the native vane CFD lattice Boltzmann engine.

    Supports collision models:
    - 0: Single-Relaxation-Time (BGK)
    - 1: Two-Relaxation-Time (TRT)
    - 2: Multiple-Relaxation-Time (MRT)
    """

    def __init__(
        self,
        nx: int,
        ny: int,
        tau: float = 0.54,
        u_inlet: float = 0.05,
        collision_type: int = 0,
        lib_path: Optional[Union[str, Path]] = None,
    ):
        """
        Initialize simulation domain and allocate lattice distributions.

        Parameters
        ----------
        nx : int
            Number of lattice nodes along x-axis (streamwise).
        ny : int
            Number of lattice nodes along y-axis (transverse).
        tau : float
            LBM relaxation time (controls kinematic viscosity: nu = (tau - 0.5) / 3).
        u_inlet : float
            Inlet free-stream velocity in lattice units (dx/dt).
        collision_type : int
            Collision operator: 0 = BGK, 1 = TRT, 2 = MRT.
        """
        self._lib = load_vane_bridge(lib_path)
        self.nx = int(nx)
        self.ny = int(ny)
        self.tau = float(tau)
        self.u_inlet = float(u_inlet)
        self.collision_type = int(collision_type)

        self._handle = self._lib.vane_sim_create(
            self.nx,
            self.ny,
            self.tau,
            self.u_inlet,
            self.collision_type,
        )
        if not self._handle:
            raise RuntimeError("Failed to allocate vane simulation handle via FFI.")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        """Free native simulation handle."""
        if hasattr(self, "_handle") and self._handle:
            self._lib.vane_sim_destroy(self._handle)
            self._handle = None

    def __del__(self):
        self.close()

    def _check_alive(self):
        if not self._handle:
            raise RuntimeError("VaneSimulation instance has been destroyed.")

    @property
    def version(self) -> str:
        """Query native vane engine version string."""
        return self._lib.vane_version().decode("utf-8")

    @property
    def current_step(self) -> int:
        """Query current simulation time step count."""
        self._check_alive()
        nx_val = ctypes.c_size_t(0)
        ny_val = ctypes.c_size_t(0)
        step_val = ctypes.c_size_t(0)
        rc = self._lib.vane_sim_get_dimensions(
            self._handle,
            ctypes.byref(nx_val),
            ctypes.byref(ny_val),
            ctypes.byref(step_val),
        )
        if rc != 0:
            raise RuntimeError(f"Error querying dimensions: code {rc}")
        return int(step_val.value)

    def set_naca_airfoil(
        self,
        code: str = "0012",
        chord: float = 40.0,
        aoa_deg: float = 0.0,
        x_lead: Optional[float] = None,
        y_lead: Optional[float] = None,
    ) -> "VaneSimulation":
        """
        Rasterize and embed a 4-digit NACA airfoil into the simulation domain.

        Parameters
        ----------
        code : str
            Four-digit NACA designation, e.g. "0012", "2412", "4415".
        chord : float
            Airfoil chord length in lattice units.
        aoa_deg : float
            Angle of attack in degrees (positive pitches leading edge up).
        x_lead : float, optional
            Leading edge x-coordinate (defaults to 0.25 * nx).
        y_lead : float, optional
            Leading edge y-coordinate (defaults to 0.50 * ny).
        """
        self._check_alive()
        if x_lead is None:
            x_lead = self.nx * 0.25
        if y_lead is None:
            y_lead = self.ny * 0.50

        c_code = code.encode("utf-8")
        rc = self._lib.vane_sim_set_naca(
            self._handle,
            c_code,
            float(chord),
            float(aoa_deg),
            float(x_lead),
            float(y_lead),
        )
        if rc != 0:
            raise ValueError(f"Failed to set NACA airfoil '{code}': error code {rc}")
        return self

    def step(self, count: int = 1) -> "VaneSimulation":
        """
        Advance simulation forward by a given number of lattice time steps.
        """
        self._check_alive()
        if count <= 0:
            return self
        rc = self._lib.vane_sim_step(self._handle, int(count))
        if rc != 0:
            raise RuntimeError(f"Simulation step failed: error code {rc}")
        return self

    def get_forces(self) -> Tuple[float, float, float]:
        """
        Compute aerodynamic lift, drag, and pitch moment coefficients.

        Returns
        -------
        cl : float
            Lift coefficient normalized by free-stream dynamic pressure and chord.
        cd : float
            Drag coefficient.
        cm : float
            Quarter-chord aerodynamic pitching moment coefficient.
        """
        self._check_alive()
        cl = ctypes.c_double(0.0)
        cd = ctypes.c_double(0.0)
        cm = ctypes.c_double(0.0)
        rc = self._lib.vane_sim_get_forces(
            self._handle,
            ctypes.byref(cl),
            ctypes.byref(cd),
            ctypes.byref(cm),
        )
        if rc != 0:
            raise RuntimeError(f"Failed to compute aerodynamic forces: code {rc}")
        return float(cl.value), float(cd.value), float(cm.value)

    def get_solid_mask(self):
        """
        Retrieve 2D binary obstacle mask (1 = solid obstacle node, 0 = fluid node).
        """
        self._check_alive()
        total_cells = self.nx * self.ny
        c_buf = (ctypes.c_uint8 * total_cells)()
        rc = self._lib.vane_sim_get_solid_mask(self._handle, c_buf, total_cells)
        if rc != 0:
            raise RuntimeError(f"Failed to read solid mask: code {rc}")

        if HAS_NUMPY:
            arr = np.ctypeslib.as_array(c_buf).copy()
            return arr.reshape((self.ny, self.nx))
        else:
            flat = list(c_buf)
            return [flat[y * self.nx : (y + 1) * self.nx] for y in range(self.ny)]

    def get_field(self, name: str):
        """
        Extract 2D scalar/vector component field from the simulation domain.

        Available field names:
        - "density"
        - "ux"
        - "uy"
        - "speed"
        - "vorticity"
        - "q_criterion"
        """
        self._check_alive()
        name_lower = name.lower()
        if name_lower not in FIELD_IDS:
            valid_names = ", ".join(FIELD_IDS.keys())
            raise ValueError(f"Unknown field '{name}'. Valid fields: {valid_names}")

        field_id = FIELD_IDS[name_lower]
        total_cells = self.nx * self.ny
        c_buf = (ctypes.c_double * total_cells)()
        rc = self._lib.vane_sim_get_field(self._handle, field_id, c_buf, total_cells)
        if rc != 0:
            raise RuntimeError(f"Failed to extract field '{name}': code {rc}")

        if HAS_NUMPY:
            arr = np.ctypeslib.as_array(c_buf).copy()
            return arr.reshape((self.ny, self.nx))
        else:
            flat = list(c_buf)
            return [flat[y * self.nx : (y + 1) * self.nx] for y in range(self.ny)]

    def get_all_fields(self) -> Dict[str, any]:
        """Retrieve dictionary containing all available physical fields."""
        return {name: self.get_field(name) for name in FIELD_IDS.keys()}

    def save_vtk(self, filepath: Union[str, Path]) -> "VaneSimulation":
        """
        Export domain macroscopic fields to XML VTK Image Data (.vti).
        """
        self._check_alive()
        c_path = str(filepath).encode("utf-8")
        rc = self._lib.vane_sim_save_vtk(self._handle, c_path)
        if rc != 0:
            raise IOError(f"Failed to write VTK file '{filepath}': code {rc}")
        return self

    def save_csv(self, filepath: Union[str, Path]) -> "VaneSimulation":
        """
        Export aerodynamic force telemetry history to CSV file.
        """
        self._check_alive()
        c_path = str(filepath).encode("utf-8")
        rc = self._lib.vane_sim_save_csv(self._handle, c_path)
        if rc != 0:
            raise IOError(f"Failed to write CSV file '{filepath}': code {rc}")
        return self

    def __repr__(self) -> str:
        status = "active" if self._handle else "destroyed"
        return (
            f"<VaneSimulation {self.nx}x{self.ny} "
            f"step={self.current_step if self._handle else 0} "
            f"status={status}>"
        )
