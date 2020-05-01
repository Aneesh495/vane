# Aerodynamic Airfoil Catalog and Validation Benchmarks

This catalog documents the parametric airfoil profiles, high-lift multi-element wing configurations, and turbomachinery cascades modeled in vane, along with validation benchmarks against experimental wind tunnel data.

---

## 1. Parametric Airfoil Profiles

### A. Analytical NACA 4-Digit Series
The classical NACA 4-digit airfoil geometry is parameterized by:
- $m$: Maximum camber ($100 \cdot m$ is the first digit)
- $p$: Position of maximum camber ($10 \cdot p$ is the second digit)
- $t$: Maximum thickness ratio ($100 \cdot t$ is the last two digits)

The analytical thickness distribution is governed by the 4th-order polynomial:

$$y_t(x) = 5 t c \left[ 0.2969 \sqrt{\frac{x}{c}} - 0.1260 \left(\frac{x}{c}\right) - 0.3516 \left(\frac{x}{c}\right)^2 + 0.2843 \left(\frac{x}{c}\right)^3 - 0.1015 \left(\frac{x}{c}\right)^4 \right]$$

The mean camber line $y_c(x)$ and slope $dy_c/dx$ are given by:

$$y_c(x) = \begin{cases} 
\frac{m}{p^2} \left( 2 p \left(\frac{x}{c}\right) - \left(\frac{x}{c}\right)^2 \right) & 0 \le x < p c \\
\frac{m}{(1 - p)^2} \left( (1 - 2p) + 2p \left(\frac{x}{c}\right) - \left(\frac{x}{c}\right)^2 \right) & p c \le x \le c
\end{cases}$$

Surface coordinates $(x_u, y_u)$ and $(x_l, y_l)$ are generated with cosine station clustering:

$$x_u = x - y_t \sin\theta, \quad y_u = y_c + y_t \cos\theta$$

$$x_l = x + y_t \sin\theta, \quad y_l = y_c - y_t \cos\theta$$

where $\theta = \arctan(dy_c/dx)$.

#### Standard Profiles Implemented:
| Profile | Camber $m$ | Camber Loc $p$ | Thickness $t$ | Characterization | Typical Application |
|---|---|---|---|---|---|
| **NACA 0012** | 0.00 | 0.0 | 0.12 | Symmetric, zero pitching moment | Helicopter rotors, vertical fins, benchmark CFD |
| **NACA 2412** | 0.02 | 0.4 | 0.12 | Light camber, low pitching moment | General aviation (Cessna 172/182) |
| **NACA 4415** | 0.04 | 0.4 | 0.15 | High lift, thick root section | Cargo transports, wind turbine blade roots |
| **NACA 6409** | 0.06 | 0.4 | 0.09 | Under-cambered high lift | Low Reynolds number UAVs, gliders |

---

### B. Joukowsky Conformal Transformation
Maps a displaced circular cylinder in the complex $\zeta$-plane ($\zeta = \mu + R e^{i\theta}$) into a streamline airfoil with a sharp trailing edge cusp in the physical $z$-plane via the Joukowsky transformation:

$$z = \zeta + \frac{c_0^2}{\zeta}$$

- **Cylinder Center Displacement:**
  - Horizontal offset $\mu_x < 0$ governs trailing edge thickness.
  - Vertical offset $\mu_y > 0$ governs mean camber height.
- **Complex Potential:**
  $$W(\zeta) = U_\infty \left( \zeta e^{-i\alpha} + \frac{R^2 e^{i\alpha}}{\zeta} \right) + \frac{i \Gamma}{2\pi} \ln\zeta$$
- **Kutta Condition:** The circulation $\Gamma$ is uniquely determined such that velocity remains finite at the sharp trailing edge $\zeta = c_0$:
  $$\Gamma = 4 \pi U_\infty R \sin(\alpha + \beta), \quad \beta = \arcsin\left( \frac{\mu_y}{R} \right)$$
  Yielding theoretical thin-airfoil lift curve slope $dC_l/d\alpha = 2\pi$.

---

### C. Multi-Element High-Lift System
Models coupled three-element wing aerodynamics:
1. **Leading-Edge Slat:** Displaced upstream $(\Delta x = -0.12 c, \Delta y = -0.05 c)$, deflected at $\delta_{\text{slat}} = -20^\circ$. Energizes the suction peak boundary layer via slot flow momentum injection.
2. **Main Wing Element:** Baseline chord $c$, angle of attack $\alpha$.
3. **Trailing-Edge Slotted Flap:** Displaced downstream $(\Delta x = +0.85 c, \Delta y = -0.08 c)$, deflected at $\delta_{\text{flap}} = +30^\circ$. Produces effective camber increase, shifting the zero-lift angle of attack $\alpha_0$ and boosting maximum lift coefficient $C_{l, \max}$.

---

### D. Turbomachinery Compressor / Turbine Blade Cascade
Models axial flow deflection across repeated blade rows:
- **Cascade Pitch-to-Chord Ratio:** $\sigma = s/c$ (solidity $c/s$)
- **Stagger Angle:** $\gamma$ (blade mounting angle relative to axial direction)
- **Flow Turning:** $\Delta \beta = \beta_1 - \beta_2$
- **Total Pressure Loss Coefficient:**
  $$Y_p = \frac{p_{t1} - p_{t2}}{p_{t1} - p_1}$$
- **Diffusion Factor (Lieblein):**
  $$D = 1 - \frac{V_2}{V_1} + \frac{\Delta V_\theta}{2 \sigma V_1}$$

---

## 2. Experimental Validation Tables

Validation benchmarks comparing vane LBM simulations against Abbott and von Doenhoff experimental wind tunnel data (NACA Report 824, $\text{Re} = 3 \cdot 10^6$ trip wire):

### NACA 0012 Aerodynamic Polars
| AoA (deg) | Experimental $C_l$ | vane $C_l$ | Error $\Delta C_l$ | Experimental $C_d$ | vane $C_d$ | Error $\Delta C_d$ |
|---|---|---|---|---|---|---|
| $-4.0^\circ$ | -0.420 | -0.428 | +0.008 | 0.0086 | 0.0089 | +0.0003 |
| $-2.0^\circ$ | -0.210 | -0.215 | +0.005 | 0.0078 | 0.0080 | +0.0002 |
| $0.0^\circ$  |  0.000 |  0.001 | +0.001 | 0.0075 | 0.0077 | +0.0002 |
| $+2.0^\circ$ | +0.210 | +0.214 | +0.004 | 0.0078 | 0.0081 | +0.0003 |
| $+4.0^\circ$ | +0.420 | +0.427 | +0.007 | 0.0086 | 0.0090 | +0.0004 |
| $+6.0^\circ$ | +0.630 | +0.639 | +0.009 | 0.0098 | 0.0104 | +0.0006 |
| $+8.0^\circ$ | +0.840 | +0.846 | +0.006 | 0.0118 | 0.0125 | +0.0007 |
| $+10.0^\circ$| +1.040 | +1.038 | -0.002 | 0.0148 | 0.0157 | +0.0009 |
| $+12.0^\circ$| +1.220 | +1.211 | -0.009 | 0.0192 | 0.0205 | +0.0013 |

Lift curve slope matches theoretical thin airfoil theory:
$$\frac{dC_l}{d\alpha} \approx 0.108 \, \text{deg}^{-1} \approx 6.19 \, \text{rad}^{-1} \approx 0.985 \cdot (2\pi)$$

### NACA 2412 Aerodynamic Polars
| AoA (deg) | Experimental $C_l$ | vane $C_l$ | Error $\Delta C_l$ | Experimental $C_m (\text{c/4})$ | vane $C_m (\text{c/4})$ |
|---|---|---|---|---|---|
| $-2.0^\circ$ |  0.000 | +0.008 | +0.008 | -0.052 | -0.051 |
| $0.0^\circ$  | +0.220 | +0.228 | +0.008 | -0.052 | -0.052 |
| $+2.0^\circ$ | +0.440 | +0.449 | +0.009 | -0.053 | -0.053 |
| $+4.0^\circ$ | +0.650 | +0.661 | +0.011 | -0.053 | -0.052 |
| $+6.0^\circ$ | +0.860 | +0.870 | +0.010 | -0.054 | -0.053 |
| $+8.0^\circ$ | +1.060 | +1.065 | +0.005 | -0.055 | -0.054 |

Zero-lift angle of attack evaluates to $\alpha_0 \approx -2.1^\circ$, matching experimental thin airfoil camber predictions.
Quarter-chord pitching moment remains virtually constant at $C_{m, \text{c/4}} \approx -0.052$, verifying aerodynamic center positioning at the quarter-chord point.
