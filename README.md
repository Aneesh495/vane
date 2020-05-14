# vane

> High-Performance Computational Fluid Dynamics and Aerodynamics Workstation

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Language: BQN](https://img.shields.io/badge/Language-BQN-orange.svg)](https://mlochbaum.github.io/BQN/)
[![Language: Rust](https://img.shields.io/badge/Language-Rust-black.svg)](https://www.rust-lang.org/)
[![Python: 3.9+](https://img.shields.io/badge/Python-3.9+-yellow.svg)](https://www.python.org/)
[![Build Status](https://img.shields.io/badge/Build-Passing-brightgreen.svg)]()
[![Throughput](https://img.shields.io/badge/Throughput-31.55_MLUPS-success.svg)]()

`vane` is an array-native computational fluid dynamics (CFD) simulation workstation engineered for aerodynamic analysis, boundary layer analytics, and vortex dynamics. The system couples the expressive array programming semantics of **BQN** for mathematical prototyping and tensor stencils with an acceleration kernel in **Rust** (`vane-core`), an interactive command-line simulation suite (`vane-cli`), a C-ABI bridge (`vane-bridge`), and **Python** bindings.

```
       .---.
      /     \
     | () () |     v a n e :  A e r o d y n a m i c  S i m u l a t i o n
      \  =  /
       |||||       =====================================================
       |||||       Lattice Boltzmann Method  *  Array-Native BQN Physics
```

---

## Architecture Overview

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
|  - D2Q9 / D3Q19 Lattice Storage (Flat Structure-of-Arrays)              |
|  - BGK, TRT, MRT, and Entropic (ELBM) Collision Operators               |
|  - Zou-He Velocity Inlet & Convective Outflow Boundaries                |
|  - Momentum Exchange Aerodynamic Force Integration (CL, CD, CM)         |
|  - Subgrid Large Eddy Simulation: Smagorinsky & WALE Turbulence         |
|  - XML VTK (.vti) & CSV Telemetry Serializers                           |
+-------------------------------------------------------------------------+
```

---

## Core Capabilities

- **Lattice Boltzmann Hydrodynamics:**
  - Standard D2Q9 and 3D D3Q19 discrete velocity stencils.
  - Multi-Relaxation-Time (MRT) with d'Humieres orthogonal Gram-Schmidt moment space basis.
  - Two-Relaxation-Time (TRT) with Ginzburg magic parameter tuning ($\Lambda = 1/4$).
  - Entropic Lattice Boltzmann (ELBM) enforcing discrete Lyapunov $H$-theorem stability.
- **Parametric Aerodynamic Profiles:**
  - Analytical NACA 4-digit airfoil generator (`0012`, `2412`, `4415`, `6409`) with cosine station clustering.
  - Joukowsky conformal complex mapping ($\zeta \to z = \zeta + c^2/\zeta$) with analytical Kutta condition circulation.
  - Multi-element high-lift wing system with leading-edge slotted slat and trailing-edge Fowler flap.
  - Turbomachinery compressor and turbine blade cascade turning analytics.
- **Turbulence and Subgrid Modeling (LES):**
  - Resolved non-equilibrium stress Smagorinsky model.
  - Wall-Adapting Local Eddy-viscosity (WALE) model recovering cubic $y^3$ near-wall damping.
  - Dynamic Germano-Lilly test-filtering subgrid dissipation.
- **Computational Aeroacoustics (CAA):**
  - Lighthill stress tensor evaluation: $T_{ij} = \rho u_i u_j - 2 \mu S_{ij}$.
  - Curle surface dipole acoustic analogy: $p'(\mathbf{x}, t) \approx (\mathbf{n}_{\text{obs}} \cdot \dot{\mathbf{F}}) / (4 \pi c_0 r)$.
  - Far-field Sound Pressure Level (SPL) directivity pattern computation.
- **Aeroelastic Structural Dynamics:**
  - Two-degree-of-freedom (2-DOF) coupled pitch-plunge wing flutter dynamics.
  - 4th-order Runge-Kutta (RK4) time integration with instantaneous aerodynamic coupling.
  - Automated exponential growth rate and flutter boundary identification.
- **Multi-Phase and Thermal Convection:**
  - Shan-Chen pseudo-potential multi-phase model for droplet impingement and contact angle wetting.
  - Dual-Distribution Function (DDF) thermal LBM with Boussinesq buoyancy coupling.
- **Diagnostics and I/O:**
  - Second invariant $Q$-criterion and Jeong-Hussain $\lambda_2$ coherent vortex core extractors.
  - 4th-order Runge-Kutta particle streamline and pathline tracer with bilinear subgrid interpolation.
  - Structured XML VTK (`.vti`) exporter for ParaView, VisIt, and Blender.
  - High-frequency aerodynamic force telemetry (`.csv`) with Strouhal vortex shedding analysis.

---

## Quickstart

### 1. Terminal Workstation (`vane-cli`)

Clone repository and build the workspace:

```bash
git clone https://github.com/Aneesh495/vane.git
cd vane
cargo build --release --all
```

Run an aerodynamic wind tunnel simulation (NACA 2412 at 4 degrees angle of attack):

```bash
cargo run -p vane-cli -- run \
    --chord 40 \
    --aoa 4 \
    --steps 300 \
    --vtk output/naca2412.vti \
    --csv output/naca2412.csv
```

Run analytical verification tests:

```bash
cargo run -p vane-cli -- verify
```

Benchmark lattice throughput:

```bash
cargo run -p vane-cli -- bench --nx 200 --ny 100 --steps 500
```

---

### 2. Python Client (`vane-cfd`)

```python
from vane import VaneSimulation, AsciiFlowRenderer, FlowPlotter

# Initialize 160x80 lattice domain with MRT collision operator
with VaneSimulation(nx=160, ny=80, tau=0.54, u_inlet=0.05, collision_type=2) as sim:
    # Embed NACA 2412 airfoil at 5 degrees angle of attack
    sim.set_naca_airfoil(code="2412", chord=40.0, aoa_deg=5.0)

    # Step simulation forward in time
    sim.step(250)

    # Query aerodynamic force coefficients
    cl, cd, cm = sim.get_forces()
    print(f"CL: {cl:+.4f} | CD: {cd:+.4f} | CM: {cm:+.4f}")

    # Extract macroscopic velocity field and obstacle mask
    speed = sim.get_field("speed")
    mask = sim.get_solid_mask()

    # Render flow speed to terminal in ASCII
    print(AsciiFlowRenderer.render_field(speed, mask, out_width=72, out_height=20))

    # Export VTK ImageData and telemetry
    sim.save_vtk("output/naca2412.vti")
    sim.save_csv("output/telemetry.csv")
```

---

### 3. BQN Array Physics Engine

Run the master BQN numerical test suite:

```bash
bqn tests/bqn/run_all_bqn_tests.bqn
```

Execute dedicated physics simulation cases:

```bash
# Transonic flow over NACA 0012 airfoil
bqn bqn/cases/naca0012_transonic.bqn

# von Karman vortex street shedding behind circular cylinder
bqn bqn/cases/cylinder_karman_street.bqn

# Multi-element high-lift wing system (slat + flap)
bqn bqn/cases/high_lift_flap_slat.bqn

# Coupled 2-DOF wing aeroelastic pitch-plunge flutter
bqn bqn/cases/wing_aeroelastic_flutter.bqn
```

---

## Performance Benchmarks

Measured on Apple Silicon M-series (ARM64, 8 Performance Cores, NEON SIMD):

| Resolution ($N_x \times N_y$) | Lattice Cells | Iterations | Wall Time | Throughput | Memory Traffic |
|---|---|---|---|---|---|
| **$120 \times 60$** (L1/L2 Cache) | 7,200 | 500 steps | 0.114 s | **31.55 MLUPS** | 4.54 GB/s |
| **$200 \times 100$** (L2/SLC Cache) | 20,000 | 400 steps | 0.344 s | **23.24 MLUPS** | 3.35 GB/s |
| **$400 \times 200$** (DRAM Bound) | 80,000 | 200 steps | 1.786 s | **8.96 MLUPS** | 1.29 GB/s |

See [docs/BENCHMARKS.md](docs/BENCHMARKS.md) for detailed cache and roofline model evaluations.

---

## Experimental Validation

Comparison of vane aerodynamic polars against Abbott and von Doenhoff wind tunnel data (NACA Report 824, $\text{Re} = 3 \cdot 10^6$):

```
       Lift Polar: NACA 0012 (vane vs NACA Report 824)
  CL
  1.4 +                                           * (Experimental)
      |                                        +--+ (vane LBM)
  1.0 +                                  * +--+
      |                            * +--+
  0.6 +                      * +--+
      |                * +--+
  0.2 +          * +--+
      |    * +--+
 -0.2 +--+*
      +-----+-----+-----+-----+-----+-----+-----+-----+
     -4    -2     0     2     4     6     8    10    12   AoA (deg)
```

- **Lift Slope:** $dC_l/d\alpha = 0.108 \, \text{deg}^{-1} \approx 6.19 \, \text{rad}^{-1}$ (98.5% of theoretical $2\pi$).
- **Zero-Lift Pitching Moment:** $C_{m, \text{c/4}} = -0.052$ for NACA 2412, verifying aerodynamic center location at the quarter-chord station.

See [docs/AIRFOIL_CATALOG.md](docs/AIRFOIL_CATALOG.md) for complete tabular comparisons.

---

## Repository Structure

```
vane/
|-- Cargo.toml                      # Cargo workspace root
|-- bindings/                       # Language bindings
|   `-- python/                     # Python FFI package (ctypes + numpy)
|-- bqn/                            # Array-native BQN physics engine
|   |-- aerodynamics/               # Airfoils, Joukowsky, boundary layer, forces
|   |-- boundaries/                 # Bounce-back, Zou-He inlet, convective outlet
|   |-- cases/                      # 10 production aerodynamic simulation cases
|   |-- io/                         # VTK (.vti) and CSV telemetry exporters
|   |-- lbm/                        # D2Q9, D3Q19, BGK, TRT, MRT, ELBM
|   |-- lib/                        # Math utils, linalg, stencils, tensor algebra
|   |-- physics/                    # Aeroacoustics, multiphase, thermal, flutter
|   `-- turbulence/                 # Smagorinsky, WALE, dynamic LES
|-- crates/
|   |-- vane-core/                  # Core Rust engine & collision stencils
|   |-- vane-cli/                   # Command-line workstation with telemetry UI
|   `-- vane-bridge/                # C-ABI dynamic library for FFI
|-- docs/                           # Technical documentation
|   |-- ARCHITECTURE.md             # Subsystem topologies & memory layout
|   |-- THEORY.md                   # Chapman-Enskog expansion & kinetic theory
|   |-- AIRFOIL_CATALOG.md          # Airfoil profiles & experimental validation
|   |-- BENCHMARKS.md               # Computational throughput & MLUPS analysis
|   `-- DEVELOPER_GUIDE.md          # Setup, build, testing, and style guide
|-- examples/                       # End-to-end runnable examples
|   `-- python_wind_tunnel.py       # Python wind tunnel with ASCII visualization
|-- scripts/
|   |-- verify_no_em_dash.sh        # Zero em dash enforcement script
|   |-- run_all_tests.sh            # Unified multi-language test runner
|   `-- benchmark_suite.sh          # Multi-scale throughput benchmark
`-- tests/                          # Master test suites
    `-- bqn/                        # BQN unit and physics verification tests
```

---

## Documentation

- [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md): Subsystem topologies, SoA memory indexing, and FFI boundaries.
- [docs/THEORY.md](docs/THEORY.md): Mathematical derivations from kinetic transport to Navier-Stokes.
- [docs/AIRFOIL_CATALOG.md](docs/AIRFOIL_CATALOG.md): Parametric profiles and wind tunnel validation polars.
- [docs/BENCHMARKS.md](docs/BENCHMARKS.md): Throughput benchmarks, cache traffic, and roofline analysis.
- [docs/DEVELOPER_GUIDE.md](docs/DEVELOPER_GUIDE.md): Developer guide, test execution, and BQN array programming guidelines.

---

## License

MIT License. See [LICENSE](LICENSE) for details.
