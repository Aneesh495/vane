# Dual-Distribution Thermal Convection and Aerodynamic Heat Transfer Engine
#
# Simulates coupled hydrodynamic and thermodynamic energy transport using the
# Dual-Distribution Function (DDF) thermal Lattice Boltzmann Method with
# Boussinesq buoyancy coupling.
#
# References:
# He, X., Chen, S., and Doolen, G. D., "A novel thermal model for the lattice
# Boltzmann method in incompressible limit", J. Comp. Phys. 146, 1998.
# Peng, Y., Shu, C., and Chew, Y. T., "Simplified thermal lattice Boltzmann model
# for incompressible thermal flows", Phys. Rev. E 68, 2003.
#
# Note: Zero em dashes are used in this file per coding standards.

math ← •Import "../lib/math_utils.bqn"
tensor ← •Import "../lib/tensor.bqn"
d2q9 ← •Import "../lbm/d2q9.bqn"

# Compute macroscopic temperature T = sum(g_i)
# Arguments:
#   g: list of 9 thermal distribution matrices [g0, ..., g8]
# Returns:
#   T: [Ny, Nx] temperature field
Temperature ⇐ {
  g ← 𝕩
  +´ g
}

# Thermal equilibrium distribution functions g_i^eq(T, ux, uy)
# Arguments:
#   T: [Ny, Nx] temperature field
#   ux: [Ny, Nx] horizontal velocity
#   uy: [Ny, Nx] vertical velocity
# Returns:
#   geq: list of 9 thermal equilibrium matrices
ThermalEquilibrium ⇐ {
  ⟨temperature, ux, uy⟩ ← 𝕩
  cxVec ← d2q9.cx
  cyVec ← d2q9.cy
  w ← d2q9.w

  uSq ← (ux × ux) + (uy × uy)

  geq ← {
    i ← 𝕩
    wi ← i ⊑ w
    cx ← 0.0 + (i ⊑ cxVec)
    cy ← 0.0 + (i ⊑ cyVec)

    ciDotU ← (cx × ux) + (cy × uy)
    # g_i^eq = T * w_i * [ 1 + 3*(c.u) + 4.5*(c.u)^2 - 1.5*u^2 ]
    poly ← 1.0 + ((3.0 × ciDotU) + ((4.5 × ciDotU × ciDotU) - (1.5 × uSq)))
    (temperature × wi) × poly
  }¨ (↕9)

  geq
}

# Thermal BGK collision operator
# Arguments:
#   g: list of 9 thermal distribution matrices
#   geq: list of 9 thermal equilibrium matrices
#   tauT: thermal relaxation time (related to thermal diffusivity alpha = (tauT - 0.5)/3)
# Returns:
#   gPost: list of 9 post-collision thermal distribution matrices
CollideThermalBGK ⇐ {
  ⟨g, geq, tauT⟩ ← 𝕩
  invTauT ← 1.0 ÷ (tauT ⌈ 1e¯12)
  g - (invTauT × (g - geq))
}

# Boussinesq buoyancy force field F_buoyant = [0, rho * g_beta * (T - T0)]
# Arguments:
#   rho: [Ny, Nx] fluid density field
#   temperature: [Ny, Nx] temperature field
#   t0: reference reference temperature
#   gBeta: gravitational acceleration times thermal expansion coefficient g * beta
# Returns:
#   ⟨fBuoyantX, fBuoyantY⟩
BoussinesqForce ⇐ {
  ⟨rho, temperature, t0, gBeta⟩ ← 𝕩
  ⟨ny, nx⟩ ← ≢temperature

  fBuoyantX ← ny‿nx ⥊ 0.0
  fBuoyantY ← (rho × gBeta) × (temperature - t0)
  ⟨fBuoyantX, fBuoyantY⟩
}

# Isothermal Dirichlet wall boundary condition (Anti-Bounce-Back)
# For solid boundary node with specified temperature Tw:
# g_ibar(x, t+dt) = -g_i^*(x, t) + 2 * w_i * Tw
# Arguments:
#   gPost: list of 9 post-collision thermal distributions
#   solidMask: [Ny, Nx] boolean mask (1 = solid wall, 0 = fluid)
#   tWall: wall temperature
# Returns:
#   gBC: list of 9 distributions with isothermal wall condition applied
IsothermalWallBC ⇐ {
  ⟨gPost, solidMask, tWall⟩ ← 𝕩
  cxVec ← d2q9.cx
  cyVec ← d2q9.cy
  opp ← d2q9.opp
  w ← d2q9.w

  fluidMask ← 1 - solidMask

  gBC ← {
    i ← 𝕩
    giPost ← i ⊑ gPost
    cx ← i ⊑ cxVec
    cy ← i ⊑ cyVec
    ibar ← i ⊑ opp
    gibarPost ← ibar ⊑ gPost
    wi ← i ⊑ w

    solidNeighbor ← ((-cy) ⌽ ((-cx) ⌽˘ solidMask))
    linkMask ← fluidMask × solidNeighbor

    # Anti-bounce-back: g_ibar = -g_i^* + 2 * w_i * Tw
    antiBB ← (0.0 - gibarPost) + ((2.0 × wi) × tWall)
    (linkMask × antiBB) + ((1 - linkMask) × giPost)
  }¨ (↕9)

  gBC
}

# Local and average Nusselt numbers along heated surface
# Nu(x) = (H / deltaT) * |dT/dy|_wall
# Arguments:
#   temperature: [Ny, Nx] temperature field
#   yWallIdx: row index of heated wall
#   heightH: characteristic domain height
#   deltaT: characteristic temperature difference (Th - Tc)
# Returns:
#   ⟨nuLocal, nuAvg⟩
ComputeNusselt ⇐ {
  ⟨temperature, yWallIdx, heightH, deltaT⟩ ← 𝕩
  ⟨ny, nx⟩ ← ≢temperature

  yFluidIdx ← (ny - 1) ⌊ (yWallIdx + 1)
  tWallRow ← yWallIdx ⊏ temperature
  tFluidRow ← yFluidIdx ⊏ temperature

  # Temperature gradient normal to wall: |dT/dy| = |T_fluid - T_wall| / dy
  # dy = 1 in lattice units
  gradT ← | (tFluidRow - tWallRow)
  prefactor ← heightH ÷ (deltaT ⌈ 1e¯12)
  nuLocal ← prefactor × gradT
  nuAvg ← (+´ nuLocal) ÷ (nx ⌈ 1)

  ⟨nuLocal, nuAvg⟩
}
