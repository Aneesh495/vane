# Mathematical Foundations and Aerodynamic Theory

This document details the continuum mechanics, statistical physics, and discrete numerical formulations underpinning the vane simulation engine.

---

## 1. From Boltzmann Transport to Navier-Stokes: Chapman-Enskog Expansion

The Lattice Boltzmann Method models fluid flow as the mesoscopic evolution of discrete velocity particle distribution functions $f_i(\mathbf{x}, t)$. The governing evolution equation with discrete velocity vectors $\mathbf{c}_i$ is:

$$f_i(\mathbf{x} + \mathbf{c}_i \Delta t, t + \Delta t) - f_i(\mathbf{x}, t) = \Omega_i(f)$$

where $\Omega_i(f)$ is the collision operator. In the continuum limit, as the Knudsen number $\text{Kn} = \lambda / L \to 0$, this discrete kinetic equation asymptotically recovers the macroscopic incompressible Navier-Stokes equations.

### Multi-Scale Asymptotic Expansion
Introducing a small Knudsen expansion parameter $\epsilon \ll 1$, spatial and temporal derivatives are expanded across multiple asymptotic scales:

$$\frac{\partial}{\partial t} = \epsilon \frac{\partial}{\partial t_1} + \epsilon^2 \frac{\partial}{\partial t_2}, \quad \nabla = \epsilon \nabla_1$$

The distribution function is expanded around local thermodynamic equilibrium $f_i^{(0)}$:

$$f_i = f_i^{(0)} + \epsilon f_i^{(1)} + \epsilon^2 f_i^{(2)} + \mathcal{O}(\epsilon^3)$$

Substituting into the Taylor series expansion of the streaming operator:

$$\left[ \epsilon \left( \frac{\partial}{\partial t_1} + \mathbf{c}_i \cdot \nabla_1 \right) + \epsilon^2 \left( \frac{\partial}{\partial t_2} + \frac{1}{2} \left( \frac{\partial}{\partial t_1} + \mathbf{c}_i \cdot \nabla_1 \right)^2 \right) \right] f_i = \Omega_i$$

### First-Order Order ($\epsilon^1$): Euler Equations
At $\mathcal{O}(\epsilon)$, integrating over discrete velocity moments yields the inviscid Euler equations of conservation of mass and momentum:

$$\frac{\partial \rho}{\partial t_1} + \nabla_1 \cdot (\rho \mathbf{u}) = 0$$

$$\frac{\partial (\rho \mathbf{u})}{\partial t_1} + \nabla_1 \cdot \Pi^{(0)} = 0$$

where the equilibrium momentum flux tensor is isotropic:

$$\Pi_{\alpha\beta}^{(0)} = \sum_i c_{i\alpha} c_{i\beta} f_i^{(0)} = p \delta_{\alpha\beta} + \rho u_\alpha u_\beta$$

with ideal gas equation of state $p = c_s^2 \rho$.

### Second-Order Order ($\epsilon^2$): Navier-Stokes Viscous Dissipation
At $\mathcal{O}(\epsilon^2)$, non-equilibrium populations $f_i^{(1)}$ contribute to the viscous momentum flux tensor:

$$\Pi_{\alpha\beta}^{(1)} = \sum_i c_{i\alpha} c_{i\beta} f_i^{(1)} = -2 \rho c_s^2 \left( \tau - \frac{1}{2} \right) \Delta t \, S_{\alpha\beta}$$

where $S_{\alpha\beta} = \frac{1}{2} (\partial_\beta u_\alpha + \partial_\alpha u_\beta)$ is the symmetric rate-of-strain tensor. Recombining $\mathcal{O}(\epsilon)$ and $\mathcal{O}(\epsilon^2)$ terms yields the macroscopic Navier-Stokes momentum equation:

$$\rho \left( \frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla) \mathbf{u} \right) = -\nabla p + \nabla \cdot [ 2 \mu \mathbf{S} ]$$

with dynamic viscosity $\mu = \rho \nu$, and kinematic shear viscosity:

$$\nu = c_s^2 \left( \tau - \frac{1}{2} \right) \Delta t$$

In standard lattice units with $\Delta x = 1, \Delta t = 1, c_s^2 = 1/3$:

$$\nu = \frac{2\tau - 1}{6}$$

---

## 2. Discrete Velocity Stencils: D2Q9 and D3Q19

To ensure rotational isotropy of the macroscopic Navier-Stokes stress tensor, the discrete lattice velocities must satisfy high-order velocity moments:

1. $\sum_i w_i = 1$
2. $\sum_i w_i c_{i\alpha} = 0$
3. $\sum_i w_i c_{i\alpha} c_{i\beta} = c_s^2 \delta_{\alpha\beta}$
4. $\sum_i w_i c_{i\alpha} c_{i\beta} c_{i\gamma} = 0$
5. $\sum_i w_i c_{i\alpha} c_{i\beta} c_{i\gamma} c_{i\delta} = c_s^4 (\delta_{\alpha\beta}\delta_{\gamma\delta} + \delta_{\alpha\gamma}\delta_{\beta\delta} + \delta_{\alpha\delta}\delta_{\beta\gamma})$

### D2Q9 Lattice Formulation
For the 2D 9-velocity stencil:
- Velocity directions:
  - Rest particle ($i = 0$): $(0, 0)$
  - Axial neighbors ($i = 1, 2, 3, 4$): $(\pm 1, 0), (0, \pm 1)$
  - Diagonal neighbors ($i = 5, 6, 7, 8$): $(\pm 1, \pm 1)$
- Quadrature weights:
  $$w_0 = \frac{4}{9}, \quad w_{1..4} = \frac{1}{9}, \quad w_{5..8} = \frac{1}{36}$$
- Lattice speed of sound:
  $$c_s = \frac{1}{\sqrt{3}} \approx 0.57735 \, \frac{\Delta x}{\Delta t}$$

### Equilibrium Distribution Function
The second-order Maxwell-Boltzmann equilibrium distribution is:

$$f_i^{\text{eq}}(\rho, \mathbf{u}) = w_i \rho \left[ 1 + \frac{\mathbf{c}_i \cdot \mathbf{u}}{c_s^2} + \frac{(\mathbf{c}_i \cdot \mathbf{u})^2}{2 c_s^4} - \frac{\mathbf{u} \cdot \mathbf{u}}{2 c_s^2} \right]$$

---

## 3. Collision Operators

### A. Bhatnagar-Gross-Krook (BGK / SRT)
The single relaxation time approximation relaxes all non-equilibrium moments toward equilibrium at a single rate:

$$\Omega_i(f) = -\frac{1}{\tau} (f_i - f_i^{\text{eq}})$$

### B. Two-Relaxation-Time (TRT)
Decomposes populations into symmetric (+) and anti-symmetric (-) parity states:

$$f_i^+ = \frac{f_i + f_{\bar{i}}}{2}, \quad f_i^- = \frac{f_i - f_{\bar{i}}}{2}$$

Collision relaxes symmetric modes with $\tau^+$ and anti-symmetric modes with $\tau^-$:

$$f_i^* = f_i - \frac{1}{\tau^+} (f_i^+ - f_i^{\text{eq},+}) - \frac{1}{\tau^-} (f_i^- - f_i^{\text{eq},-})$$

The shear viscosity depends exclusively on $\tau^+$, while the Ginzburg magic parameter:

$$\Lambda = \left( \tau^+ - \frac{1}{2} \right) \left( \tau^- - \frac{1}{2} \right)$$

can be set to $\Lambda = 1/4$ to fix the numerical bounce-back boundary location independently of viscosity.

### C. Multiple-Relaxation-Time (MRT)
Projects distribution functions into moment space via an orthogonal transformation matrix $\mathbf{M}$:

$$\mathbf{m} = \mathbf{M} \mathbf{f}, \quad \mathbf{m}^* = \mathbf{m} - \mathbf{S} (\mathbf{m} - \mathbf{m}^{\text{eq}}), \quad \mathbf{f}^* = \mathbf{M}^{-1} \mathbf{m}^*$$

The 9 orthogonal moments for D2Q9 correspond to distinct physical modes:
- $m_0$: Conserved density $\delta\rho$ (relaxation rate $s_\rho = 0$)
- $m_1$: Energy mode $e$ ($s_e = 1.1$)
- $m_2$: Energy squared $\epsilon$ ($s_\epsilon = 1.0$)
- $m_3, m_5$: Conserved momenta $j_x, j_y$ ($s_j = 0$)
- $m_4, m_6$: Heat fluxes $q_x, q_y$ ($s_q = 1.2$)
- $m_7, m_8$: Normal and shear stresses $p_{xx}, p_{xy}$ ($s_\nu = 1 / \tau$)

Decoupling viscous relaxation ($s_\nu$) from high-order kinetic modes suppresses high-frequency spurious acoustic oscillations near solid boundaries at high Reynolds numbers.

### D. Entropic Lattice Boltzmann (ELBM)
ELBM guarantees thermodynamic stability by enforcing a discrete analogue of Boltzmann's H-theorem:

$$H(f) = \sum_i f_i \ln\left( \frac{f_i}{w_i} \right)$$

The post-collision state satisfies:

$$f_i^* = f_i + \alpha \beta (f_i^{\text{eq}} - f_i)$$

where $\alpha$ is evaluated at each node to satisfy the entropy involution condition $H(f + \alpha (f^{\text{eq}} - f)) = H(f)$. This enforces non-decreasing entropy and guarantees unconditioned numerical stability without non-physical negative populations.

---

## 4. Subgrid-Scale Turbulence Modeling (LES)

For high Reynolds number aerodynamics where the boundary layer and wake contain subgrid scales smaller than grid resolution $\Delta x$, the effective relaxation time is augmented by an eddy viscosity:

$$\tau_{\text{eff}} = \tau_0 + 3 \nu_t$$

### Smagorinsky Model
Relates eddy viscosity to the local resolved strain rate magnitude $\|S\|$:

$$\nu_t = (C_s \Delta)^2 \|S\|$$

In LBM, the strain rate is evaluated locally from the non-equilibrium momentum flux tensor $\Pi_{\alpha\beta}^{(1)}$ without computing finite-difference spatial gradients:

$$\|S\| = \frac{\|\Pi^{(1)}\|}{2 \rho c_s^2 \tau_0}$$

### Wall-Adapting Local Eddy-Viscosity (WALE) Model
Overcomes the standard Smagorinsky limitation of non-zero eddy viscosity at laminar solid walls by incorporating the square of the velocity gradient tensor $g_{\alpha\beta} = \partial_\beta u_\alpha$:

$$\nu_t = (C_w \Delta)^2 \frac{(S^d_{\alpha\beta} S^d_{\alpha\beta})^{3/2}}{(S_{\alpha\beta} S_{\alpha\beta})^{5/2} + (S^d_{\alpha\beta} S^d_{\alpha\beta})^{5/4}}$$

where $S^d_{\alpha\beta}$ is the traceless symmetric part of the square of the velocity gradient tensor. Near a solid wall, WALE naturally yields $\nu_t = \mathcal{O}(y^3)$, reproducing the physical laminar sublayer behavior without requiring empirical Van Driest damping functions.

---

## 5. Aerodynamic Force Integration: Momentum Exchange Method

Aerodynamic lift, drag, and moment are computed along the fluid-solid boundary links using the momentum exchange formulation.

For every boundary link between fluid node $\mathbf{x}_f$ and solid node $\mathbf{x}_b = \mathbf{x}_f + \mathbf{c}_i \Delta t$:

$$\Delta \mathbf{P}_i = \mathbf{c}_i f_i^*(\mathbf{x}_f, t) - \mathbf{c}_{\bar{i}} f_{\bar{i}}(\mathbf{x}_f, t + \Delta t)$$

Using halfway bounce-back, this evaluates to:

$$\Delta \mathbf{P}_i = \mathbf{c}_i \left[ f_i^*(\mathbf{x}_f, t) + f_{\bar{i}}^*(\mathbf{x}_f, t) \right] - 2 w_i \rho_0 \frac{\mathbf{c}_i \cdot \mathbf{u}_w}{c_s^2} \mathbf{c}_i$$

Integrating over all boundary links gives total hydrodynamic force:

$$\mathbf{F} = \sum_{\text{links}} \Delta \mathbf{P}_i$$

The non-dimensional aerodynamic force coefficients are normalized by free-stream dynamic pressure $q_\infty = \frac{1}{2} \rho_0 u_\infty^2$ and chord $c$:

$$C_d = \frac{F_x}{\frac{1}{2} \rho_0 u_\infty^2 c}, \quad C_l = \frac{F_y}{\frac{1}{2} \rho_0 u_\infty^2 c}, \quad C_m = \frac{M_{z, \text{qc}}}{\frac{1}{2} \rho_0 u_\infty^2 c^2}$$

---

## 6. Computational Aeroacoustics: Lighthill and Curle Analogies

Aeroacoustic noise generation is governed by Lighthill's acoustic analogy:

$$\frac{\partial^2 \rho'}{\partial t^2} - c_0^2 \nabla^2 \rho' = \frac{\partial^2 T_{ij}}{\partial x_i \partial x_j}$$

where the Lighthill stress tensor is:

$$T_{ij} = \rho u_i u_j - 2 \mu S_{ij} + (p - c_0^2 \rho') \delta_{ij}$$

In the presence of solid aerodynamic bodies, Curle extended Lighthill's analogy to account for surface dipole sound:

$$p'(\mathbf{x}, t) = \frac{1}{4 \pi c_0} \int_{\text{surface}} \frac{\mathbf{n} \cdot \mathbf{r}}{r^2} \left[ \frac{\partial p}{\partial t} \right]_{\text{ret}} dS$$

For far-field compact airfoils ($c \ll \lambda_{\text{acoustic}}$), the radiated dipole acoustic pressure simplifies to:

$$p'(\mathbf{x}, t) \approx \frac{\mathbf{n}_{\text{obs}} \cdot \dot{\mathbf{F}}(t - r/c_0)}{4 \pi c_0 r}$$

where $\dot{\mathbf{F}} = d\mathbf{F}/dt$ is the time rate of change of the integrated aerodynamic force vector. The resulting directivity pattern exhibits a classic dipole $\cos^2(\theta)$ radiation shape oriented along the direction of unsteady aerodynamic lift fluctuations.

---

## 7. Coupled 2-DOF Aeroelastic Flutter Dynamics

Unsteady aerodynamic forces couple to airfoil structural flexure through the two-degree-of-freedom pitch-plunge equations of motion:

$$\mathbf{M} \ddot{\mathbf{q}} + \mathbf{C} \dot{\mathbf{q}} + \mathbf{K} \mathbf{q} = \mathbf{Q}_{\text{aero}}$$

where generalized coordinates are $\mathbf{q} = [h, \alpha]^T$ (plunge displacement $h$, pitch angle $\alpha$).

### Mass, Damping, and Stiffness Matrices
$$\mathbf{M} = \begin{bmatrix} m & S_\alpha \\ S_\alpha & I_\alpha \end{bmatrix}, \quad \mathbf{C} = \begin{bmatrix} 2 \zeta_h m \omega_h & 0 \\ 0 & 2 \zeta_\alpha I_\alpha \omega_\alpha \end{bmatrix}, \quad \mathbf{K} = \begin{bmatrix} m \omega_h^2 & 0 \\ 0 & I_\alpha \omega_\alpha^2 \end{bmatrix}$$

where $S_\alpha = m x_\alpha$ is static mass unbalance about the elastic axis. Generalized forces are:

$$\mathbf{Q}_{\text{aero}} = \begin{bmatrix} -L(t) \\ M_{\text{ea}}(t) \end{bmatrix}$$

Integration in time proceeds via a 4th-order Runge-Kutta scheme, with aerodynamic feedback from the instantaneous LBM flow field driving structural deformation and grid rotation.
