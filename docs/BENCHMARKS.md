# Computational Performance and Throughput Benchmarks

This report benchmarks the computational throughput, memory bandwidth utilization, and scaling characteristics of the vane lattice Boltzmann aerodynamic simulation engine.

---

## 1. Throughput Benchmarks across Grid Resolutions

Simulations were executed across three standard grid resolutions on an Apple Silicon M-series unified memory architecture (ARM64, 8 performance cores, NEON SIMD). Throughput is quantified in **MLUPS** (Million Lattice Updates Per Second):

$$\text{MLUPS} = \frac{N_{\text{cells}} \cdot N_{\text{steps}}}{\Delta t_{\text{wall}} \cdot 10^6}$$

| Domain Grid ($N_x \times N_y$) | Total Lattice Cells | Time Steps | Wall Time | Throughput |
|---|---|---|---|---|
| **$120 \times 60$ (L1/L2 Cache Resident)** | 7,200 | 500 | 0.114 s | **31.55 MLUPS** |
| **$200 \times 100$ (L2/SLC Cache Resident)** | 20,000 | 400 | 0.344 s | **23.24 MLUPS** |
| **$400 \times 200$ (DRAM Bandwidth Bound)** | 80,000 | 200 | 1.786 s | **8.96 MLUPS** |

---

## 2. Memory Traffic and Roofline Analysis

In Lattice Boltzmann CFD, each lattice node requires updating $Q = 9$ double-precision distribution functions:
- **Read Phase (Streaming + Moments):** 9 values $\times 8$ bytes $= 72$ bytes/cell.
- **Write Phase (Post-Collision Store):** 9 values $\times 8$ bytes $= 72$ bytes/cell.
- **Minimum Memory Bandwidth:** $144$ bytes per lattice site update.

### Effective Memory Bandwidth
At 31.55 MLUPS in cache:
$$\text{Bandwidth}_{\text{eff}} = 31.55 \cdot 10^6 \, \frac{\text{cells}}{\text{s}} \times 144 \, \frac{\text{bytes}}{\text{cell}} \approx 4.54 \, \text{GB/s}$$

For the $400 \times 200$ grid:
$$\text{Bandwidth}_{\text{DRAM}} = 8.96 \cdot 10^6 \times 144 \approx 1.29 \, \text{GB/s}$$

### Operational Intensity
A standard D2Q9 BGK collision step requires approximately 38 floating-point operations per cell:
- Equilibrium calculation: $\sim 24$ FLOPs
- Relaxation update: $9 \times 2 = 18$ FLOPs
- Operational intensity:
  $$I = \frac{38 \, \text{FLOPs}}{144 \, \text{bytes}} \approx 0.264 \, \frac{\text{FLOP}}{\text{byte}}$$

Because $I < 1.0 \, \text{FLOP/byte}$, LBM streaming and collision is memory-bandwidth bound rather than compute bound. The Structure-of-Arrays (SoA) layout adopted in vane maintains sequential unit-stride access patterns during streaming, allowing hardware memory prefetchers to achieve near-optimal saturation of memory channels.

---

## 3. Memory Footprint per Grid Scale

Contiguous buffer allocations for D2Q9 double-precision float arrays (`f`, `f_stream`, `f_prev`) plus macroscopic fields and obstacle masks:

| Grid Size | Cells | Distributions ($3 \times 9 \times 8$B) | Macroscopic ($3 \times 8$B) | Mask & Aux (2B) | Total Memory |
|---|---|---|---|---|---|
| $120 \times 60$ | 7,200 | 1.55 MB | 0.17 MB | 0.01 MB | **1.73 MB** |
| $200 \times 100$ | 20,000 | 4.32 MB | 0.48 MB | 0.04 MB | **4.84 MB** |
| $400 \times 200$ | 80,000 | 17.28 MB | 1.92 MB | 0.16 MB | **19.36 MB** |
| $1000 \times 500$ | 500,000 | 108.0 MB | 12.0 MB | 1.00 MB | **121.0 MB** |

The compact memory footprint allows domain grids up to $200 \times 100$ to fit entirely within L2 and System Level Caches (SLC), explaining the sustained $>20$ MLUPS throughput.
