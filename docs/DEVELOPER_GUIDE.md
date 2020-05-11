# Developer and Contributor Guide

Welcome to the development guide for vane. This document covers toolchain setup, compilation instructions, test execution, benchmark procedures, and BQN array programming style guidelines.

---

## 1. Prerequisites and Toolchains

Ensure the following compilers and interpreters are available in your path:

1. **Rust Toolchain:**
   - Minimum Supported Rust Version (MSRV): `1.75.0` or newer.
   - Install via rustup: `curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh`
   - Target components: `cargo`, `rustc`, `rust-std`.

2. **CBQN / BQN Interpreter:**
   - `cbqn` (native C-based BQN) or `bqn` on PATH.
   - Homebrew (macOS): `brew install bqn` or build from source via `git clone https://github.com/dzaima/CBQN`.

3. **Python (Bindings & Post-Processing):**
   - Python `3.9+` with standard ctypes.
   - Recommended packages: `numpy`, `pillow`, `matplotlib` (`pip install numpy pillow matplotlib`).

---

## 2. Workspace Build Instructions

Compile all crates in the workspace in optimized release mode:

```bash
# Build Rust core, CLI workstation, and C-ABI FFI bridge
cargo build --release --all
```

This compiles:
- `target/release/vane-cli`: Standalone terminal simulation workstation executable.
- `target/release/libvane_bridge.dylib` (or `.so` on Linux): Shared C-ABI library for Python and external bindings.
- `crates/vane-core`: Core simulation library and numerical algorithms.

---

## 3. Test Suites and Verification

### A. Run Master Test Suite
To execute the comprehensive verification pass across all languages in a single command:

```bash
./scripts/run_all_tests.sh
```

This script executes:
1. Zero em dash verification across all tracked files.
2. Rust unit and integration tests (`cargo test --all`).
3. Master BQN array physics test runner (`bqn tests/bqn/run_all_bqn_tests.bqn`).
4. Python FFI wind tunnel simulation (`python3 examples/python_wind_tunnel.py`).
5. Terminal CLI smoke test and aerodynamic verification pass.

### B. Run BQN Physics Tests Individually
```bash
# Mathematical utilities and linear algebra
bqn tests/bqn/test_math.bqn

# Lattice Boltzmann dynamics (D2Q9, MRT, TRT, D3Q19, ELBM)
bqn tests/bqn/test_lbm.bqn

# Aerodynamics, boundary layer, and momentum exchange
bqn tests/bqn/test_aerodynamics.bqn

# Multi-physics (aeroacoustics, multiphase, flutter)
bqn tests/bqn/test_physics.bqn

# Post-processing (Q-criterion, Lambda-2, streamlines)
bqn tests/bqn/test_post.bqn
```

### C. Run Throughput Benchmarks
```bash
./scripts/benchmark_suite.sh
```

---

## 4. Running Simulation Cases

### A. Terminal Interactive CLI
```bash
# Simulate NACA 2412 at 4 degrees AoA for 300 steps with VTK and CSV output
cargo run -p vane-cli -- run --chord 40 --aoa 4 --steps 300 --vtk output/naca2412.vti --csv output/naca2412.csv

# Analytical verification pass
cargo run -p vane-cli -- verify

# Profile geometric station generation
cargo run -p vane-cli -- airfoil --code 0012 --chord 50 --aoa 6
```

### B. Python Wind Tunnel Script
```bash
python3 examples/python_wind_tunnel.py
```

Outputs:
- Instantaneous aerodynamic lift, drag, and moment readouts.
- Terminal ASCII flow speed visualization.
- XML VTK ImageData (`.vti`) file for visualization in ParaView.
- Time-series force telemetry (`.csv`).
- Color raster plot (`.png`).

---

## 5. Coding Standards and Style Guidelines

### Absolute Zero Em Dash Rule
To preserve pristine typography across all operating systems, terminals, and compilers:
- **Never use em dashes (`—`) or en dashes (`–`) anywhere in comments, docstrings, documentation, markdown, or code strings.**
- Use hyphens (`-`), colons (`:`), or parentheses for separation.
- Always run `./scripts/verify_no_em_dash.sh` before committing changes.

### BQN Array Programming Conventions
1. **Case-Sensitive Identifier Roles:**
   - Uppercase identifiers (`Sin`, `Moments`, `StepRK4`) denote functions and operations.
   - Lowercase identifiers (`rho`, `ux`, `chord`, `tau`) denote subjects (arrays, scalars, namespaces).
2. **Rank-Polymorphic Expressions:**
   - Prefer array-level operations over scalar indexing loops.
   - Use `+´` (sum reduction) and `⌈´` (max reduction) on flat arrays (`⥊ array`).
   - Use table outer product (`⌜`) and rank cells (`˘`) for tensor contractions.
3. **Branch-Free Choices:**
   - Avoid branch divergence where possible.
   - Use branch-free selection `condition ⊑ ⟨false_val, true_val⟩` for smooth numerical evaluation.
