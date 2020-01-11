# Shan-Chen Multi-Phase and Surface Wetting Lattice Boltzmann Model
#
# Simulates liquid-gas phase separation, surface tension, droplet impingement,
# contact angles, and aerodynamic surface wetting phenomena via the Shan-Chen
# pseudo-potential interaction method.
#
# Reference:
# Shan, X. and Chen, H., "Lattice Boltzmann model for simulating flows with multiple
# phases and components", Phys. Rev. E 47, 1993.
# Shan, X. and Chen, H., "Simulation of nonideal gases and liquid-gas phase transitions",
# Phys. Rev. E 49, 1994.
#
# Note: Zero em dashes are used in this file per coding standards.

math ← •Import "../lib/math_utils.bqn"
tensor ← •Import "../lib/tensor.bqn"
d2q9 ← •Import "../lbm/d2q9.bqn"

# Pseudo-potential function psi(rho) = psi0 * (1 - exp(-rho / rho0))
# Default: psi0 = 1.0, rho0 = 1.0
Psi ⇐ {
  rho ← 𝕩
  1.0 - (⋆ (0.0 - rho))
}

# Intermolecular interaction force F_int(x)
# F_int(x) = -G * psi(x) * sum_i [ w_i * psi(x + c_i * dt) * c_i ]
# Arguments:
#   rho: [Ny, Nx] density field
#   gParam: interaction strength parameter G (G < -4 for phase separation)
#
# Returns list:
#   ⟨fIntX, fIntY⟩
InteractionForce ⇐ {
  ⟨rho, gParam⟩ ← 𝕩
  ⟨ny, nx⟩ ← ≢rho

  cxVec ← d2q9.cx
  cyVec ← d2q9.cy
  w ← d2q9.w

  psiField ← Psi rho

  # Directional neighbor accumulation
  acc ← +´ {
    i ← 𝕩
    cx ← i ⊑ cxVec
    cy ← i ⊑ cyVec
    wi ← i ⊑ w

    # Shift psi field to sample neighbor psi(x + ci)
    # cy shift along rows, cx shift along cols
    psiNeighbor ← ((-cy) ⌽ ((-cx) ⌽˘ psiField))
    weightNeighbor ← wi × psiNeighbor

    ⟨(0.0 + cx) × weightNeighbor, (0.0 + cy) × weightNeighbor⟩
  }¨ (1 + ↕8)

  # F_int = -G * psi * acc
  prefactor ← (0.0 - gParam) × psiField
  fIntX ← prefactor × (0 ⊑ acc)
  fIntY ← prefactor × (1 ⊑ acc)

  ⟨fIntX, fIntY⟩
}

# Fluid-solid wetting surface force F_wall(x)
# Adjusts contact angle and droplet wettability against solid walls
# F_wall(x) = -G_wall * psi(x) * sum_i [ w_i * solidMask(x + c_i) * c_i ]
# Arguments:
#   rho: [Ny, Nx] density field
#   solidMask: [Ny, Nx] solid mask (1 = solid, 0 = fluid)
#   gWall: wetting parameter (positive for hydrophobic, negative for hydrophilic)
#
# Returns list:
#   ⟨fWallX, fWallY⟩
SurfaceWettingForce ⇐ {
  ⟨rho, solidMask, gWall⟩ ← 𝕩
  cxVec ← d2q9.cx
  cyVec ← d2q9.cy
  w ← d2q9.w

  psiField ← Psi rho

  acc ← +´ {
    i ← 𝕩
    cx ← i ⊑ cxVec
    cy ← i ⊑ cyVec
    wi ← i ⊑ w

    solidNeighbor ← ((-cy) ⌽ ((-cx) ⌽˘ solidMask))
    weightNeighbor ← wi × solidNeighbor

    ⟨(0.0 + cx) × weightNeighbor, (0.0 + cy) × weightNeighbor⟩
  }¨ (1 + ↕8)

  prefactor ← (0.0 - gWall) × psiField
  fWallX ← prefactor × (0 ⊑ acc)
  fWallY ← prefactor × (1 ⊑ acc)

  ⟨fWallX, fWallY⟩
}

# Equilibrium populations with Shan-Chen velocity shift
# Arguments:
#   rho: [Ny, Nx] density field
#   ux: [Ny, Nx] bare fluid x-velocity
#   uy: [Ny, Nx] bare fluid y-velocity
#   fxTotal: [Ny, Nx] total force x-component
#   fyTotal: [Ny, Nx] total force y-component
#   tau: relaxation time
#
# Returns list:
#   feq: 9 equilibrium population matrices
EquilibriumWithForce ⇐ {
  ⟨rho, ux, uy, fxTotal, fyTotal, tau⟩ ← 𝕩

  denomRho ← rho ⌈ 1e¯12
  # Shifted velocity: u_eq = u + (tau * F) / rho
  uxEq ← ux + ((tau × fxTotal) ÷ denomRho)
  uyEq ← uy + ((tau × fyTotal) ÷ denomRho)

  feq ← d2q9.Equilibrium ⟨rho, uxEq, uyEq⟩
  feq
}

# Real physical fluid velocity u_phys = u + F / (2 * rho)
# Arguments:
#   ux: [Ny, Nx] bare fluid x-velocity
#   uy: [Ny, Nx] bare fluid y-velocity
#   fxTotal: [Ny, Nx] total force x-component
#   fyTotal: [Ny, Nx] total force y-component
#   rho: [Ny, Nx] density field
PhysicalVelocity ⇐ {
  ⟨ux, uy, fxTotal, fyTotal, rho⟩ ← 𝕩
  denomTwoRho ← (2.0 × rho) ⌈ 1e¯12
  uPhysX ← ux + (fxTotal ÷ denomTwoRho)
  uPhysY ← uy + (fyTotal ÷ denomTwoRho)
  ⟨uPhysX, uPhysY⟩
}

# Phase Interface Identification and Droplet Radius
# Finds boundary interface where density equals average of liquid and vapor densities
# Arguments:
#   rho: [Ny, Nx] density field
#   rhoLiquid: reference liquid density (e.g. 2.0)
#   rhoVapor: reference vapor density (e.g. 0.1)
#
# Returns list:
#   ⟨interfaceMask, equivalentRadius⟩
PhaseInterface ⇐ {
  ⟨rho, rhoLiquid, rhoVapor⟩ ← 𝕩
  rhoMid ← 0.5 × (rhoLiquid + rhoVapor)
  liquidMask ← rho ≥ rhoMid

  # Droplet volume/area = total number of liquid cells
  area ← +´ ⥊ liquidMask
  pi ← 3.14159265358979323846
  # Equivalent 2D circle radius R = sqrt(Area / pi)
  equivRadius ← ((0.0 ⌈ area) ÷ pi) ⋆ 0.5

  ⟨liquidMask, equivRadius⟩
}
