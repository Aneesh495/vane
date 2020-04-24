# System Architecture and Memory Subsystem

vane is structured as a multi-tier computational fluid dynamics workstation designed for high-resolution aerodynamic flows, boundary layer analytics, and vortex dynamics. The system combines array-native functional programming in BQN with cache-optimized systems engineering in Rust, exposed through high-level Python bindings and a command-line interface.

```
+-------------------------------------------------------------------------+
|                                vane-cli                                 |
|         Interactive Terminal Workstation & Real-Time Telemetry          |
+--------------------+--------------------------------+-------------------+
                     |                                |
                     v                                v
+--------------------+--------+              +--------+-------------------+
|        bindings/python      |              |          bqn/cases         |
|   ctypes FFI + Visualizers  |              | Numerical Benchmark Suites |
+--------------------+--------+              +--------+-------------------+
                     |                                |
                     v                                v
+--------------------+--------+              +--------+-------------------+
|         vane-bridge         |              |         bqn/lbm            |
|   C-ABI Dynamic Shared Lib  |              | Array Physics Engine (BQN) |
+--------------------+--------+              +----------------------------+
                     |
                     v
+--------------------+----------------------------------------------------+
|                               vane-core                                 |
|  - D2Q9 / D3Q19 Lattice Storage (Flat SoA)                              |
|  - BGK / TRT / MRT Collision Operators                                  |
|  - Zou-He & Convective Outflow Boundaries                               |
|  - Momentum Exchange Aerodynamic Force Integration                      |
|  - Smagorinsky & WALE Subgrid Scale LES                                 |
|  - XML VTK (.vti) & CSV Telemetry Serializers                           |
+-------------------------------------------------------------------------+
```

---

## 1. Multi-Tier Subsystem Layout

1. **`bqn/` (Array Physics Engine)**
   Pure functional implementation of discrete lattice Boltzmann stencils, multi-relaxation-time collision matrices, Joukowsky conformal maps, aeroacoustic directivities, and 2-DOF aeroelastic flutter integrators. Serves as the mathematical reference formulation where algorithms are expressed in rank-polymorphic array primitives.

2. **`crates/vane-core` (High-Performance Engine)**
   Memory-aligned Rust implementation with cache-friendly streaming loops, multi-threaded collision kernels powered by Rayon, SIMD vectorization, and algebraic subgrid turbulence models.

3. **`crates/vane-cli` (Terminal Workstation)**
   Interactive terminal simulation environment featuring animated Unicode progress bars, live aerodynamic coefficient readouts (CL, CD, CM, L/D), automated NACA geometry rasterization, and automated verification suites.

4. **`crates/vane-bridge` (C-ABI Dynamic Library)**
   Shared dynamic library (`libvane_bridge.dylib` / `.so`) providing an `extern "C"` foreign function interface. Encapsulates simulation lifecycle management, in-place macroscopic field extraction, and VTK/CSV export for external runtimes.

5. **`bindings/python/` (Python Interface & Visualization)**
   NumPy-integrated ctypes wrapper offering an object-oriented Python client (`VaneSimulation`), terminal ASCII field renderer (`AsciiFlowRenderer`), and raster graphics export (`FlowPlotter`).

---

## 2. Memory Topology and Layout

### Structure of Arrays (SoA) Storage
Lattice distribution functions are stored in a contiguous 1D vector representing a 3D tensor of shape `[Q, Ny, Nx]`, where `Q = 9` for D2Q9:

$$\text{Index}(i, y, x) = i \cdot (N_y \cdot N_x) + y \cdot N_x + x$$

This layout guarantees that:
- For a fixed discrete velocity direction $i$, all spatial nodes $(y, x)$ are laid out sequentially in memory with unit stride (`stride = 1`).
- During collision, cell-local populations $f_0, \dots, f_8$ are gathered into CPU registers or L1 cache lines without stride contention.
- During streaming, reading from neighbor offset $(y - c_{i,y}, x - c_{i,x})$ corresponds to a constant memory offset, maximizing hardware prefetcher efficiency.

### Ping-Pong Double Buffering
To prevent race conditions during streaming, vane employs two distribution buffers:
- `f`: Active post-collision population field.
- `f_stream`: Target post-streaming population field.

After streaming completes across all cells, the pointers are swapped or `f` is updated in-place via cache-line block copies.

---

## 3. Boundary Condition Pipelines

```
       +-------------------------------------------------------------+
       |                  Top: Free-Slip Wall (dvy/dy = 0)           |
       +-------------------------------------------------------------+
       |                                                             |
Inlet  | --->              [Airfoil Solid Obstacle]             ---> | Outlet
Zou-He | --->              Halfway Bounce-Back (No-Slip)        ---> | Convective
u = u0 | --->              Momentum Exchange Force Integration  ---> | Outflow
       |                                                             |
       +-------------------------------------------------------------+
       |                 Bottom: Free-Slip Wall (dvy/dy = 0)         |
       +-------------------------------------------------------------+
```

1. **West Boundary (Zou-He Velocity Inlet):**
   Specifies constant velocity $\mathbf{u} = (u_{\text{inf}}, 0)$. Missing populations $f_1, f_5, f_8$ entering from the west boundary are computed from density $\rho_{\text{in}}$ and known populations using momentum conservation constraints.

2. **East Boundary (Convective Outflow):**
   Implements first-order convective boundary conditions:
   $$\frac{\partial f_i}{\partial t} + U_{\text{conv}} \frac{\partial f_i}{\partial x} = 0$$
   This prevents non-physical acoustic reflections from bouncing back into the computational domain.

3. **North and South Boundaries (Free-Slip Symmetry Walls):**
   Vertical velocity normal to the boundary is set to zero while horizontal tangential velocity has zero gradient, simulating wind tunnel test section boundaries without artificial wall friction.

4. **Solid Obstacle Surface (Halfway Bounce-Back):**
   Nodes located on the obstacle boundary reflect incoming distributions back along their reverse direction $\bar{i}$:
   $$f_{\bar{i}}(\mathbf{x}_f, t + \Delta t) = f_i^*(\mathbf{x}_f, t) - 2 w_i \rho_0 \frac{\mathbf{c}_i \cdot \mathbf{u}_w}{c_s^2}$$
   Yielding no-slip boundary behavior with second-order spatial accuracy.

---

## 4. Foreign Function Interface (FFI) Architecture

The FFI boundary is designed for zero memory leaks and thread safety:
- **Opaque Pointer Handles:** The simulation state is boxed on the Rust heap as `Box<LbmSimulation>` and returned to external runtimes as an opaque pointer `*mut c_void`.
- **Buffer Protocol:** Field extraction methods take caller-allocated pointers and buffer lengths, preventing internal reallocation and enabling zero-copy wrapping via NumPy's `as_array`.
- **Exception Boundary:** All C-ABI functions contain panic-safe barriers returning integer status codes ($0$ for success, negative for invalid parameters or bounds errors).
