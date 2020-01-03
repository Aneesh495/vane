# Boundary Layer Velocity Profiles and Integral Thickness Analytics
#
# Analyzes viscous boundary layers over aerodynamic bodies, computing:
#   - 99% boundary layer thickness delta99
#   - Displacement thickness delta_star = integral(1 - u/Ue) dy
#   - Momentum thickness theta = integral(u/Ue * (1 - u/Ue)) dy
#   - Shape factor H = delta_star / theta (Blasius ~ 2.59, Turbulent ~ 1.4)
#   - Wall shear stress tau_w and skin friction coefficient Cf
#
# Reference:
# Schlichting, H. and Gersten, K., "Boundary-Layer Theory", 9th ed., Springer.
#
# Note: Zero em dashes are used in this file per coding standards.

math ← •Import "../lib/math_utils.bqn"
tensor ← •Import "../lib/tensor.bqn"

# Compute boundary layer integral thicknesses from discrete wall-normal velocity profile
# Arguments:
#   yNorm: [N] wall-normal distances (lattice units, starting at 0 at wall)
#   uTang: [N] tangential velocity values (starting at 0 at wall)
#   nu: kinematic viscosity
#   rho0: freestream reference density
#   uInf: freestream reference velocity magnitude
#
# Returns list:
#   delta99: 99% boundary layer thickness
#   deltaStar: displacement thickness
#   theta: momentum thickness
#   shapeFactor: H = deltaStar / theta
#   tauWall: wall shear stress
#   cf: skin friction coefficient
IntegralParameters ⇐ {
  ⟨yNorm, uTang, nu, rho0, uInf⟩ ← 𝕩
  nPoints ← ≠yNorm

  # Freestream edge velocity Ue is max tangential velocity in outer profile
  ueRaw ← ⌈´ uTang
  ue ← ueRaw ⌈ 1e¯12

  # Normalized velocity u / Ue clamped to [0, 1]
  uOverUe ← 0.0 ⌈ (1.0 ⌊ (uTang ÷ ue))

  # 1. Delta 99: find first height where u >= 0.99 * Ue
  above99 ← (uOverUe ≥ 0.99)
  hasAbove99 ← ∨´ above99
  firstIdx ← hasAbove99 ⊑ ⟨nPoints - 1, 0 ⊑ ((above99 = 1) ⊸/ (↕nPoints))⟩
  delta99 ← firstIdx ⊑ yNorm

  # 2. Integrands:
  # Displacement integrand: 1 - u/Ue
  dispIntegrand ← 1.0 - uOverUe
  # Momentum integrand: (u/Ue) * (1 - u/Ue)
  momIntegrand ← uOverUe × (1.0 - uOverUe)

  # Trapezoidal integration across profile
  dy ← (1 ↓ yNorm) - (¯1 ↓ yNorm)
  dispMid ← 0.5 × ((1 ↓ dispIntegrand) + (¯1 ↓ dispIntegrand))
  momMid ← 0.5 × ((1 ↓ momIntegrand) + (¯1 ↓ momIntegrand))

  deltaStar ← +´ (dispMid × dy)
  theta ← +´ (momMid × dy)

  # Shape factor H = deltaStar / theta
  denomTheta ← theta ⌈ 1e¯12
  shapeFactor ← deltaStar ÷ denomTheta

  # 3. Wall shear stress tau_w = mu * (du/dy)_wall
  # Finite difference at wall using first two points
  dy0 ← ((1 ⊑ yNorm) - (0 ⊑ yNorm)) ⌈ 1e¯12
  du0 ← (1 ⊑ uTang) - (0 ⊑ uTang)
  dudyWall ← du0 ÷ dy0
  mu ← rho0 × nu
  tauWall ← mu × dudyWall

  # Skin friction coefficient Cf = tau_w / (0.5 * rho0 * uInf^2)
  qInf ← (0.5 × rho0) × (uInf × uInf)
  denomQ ← qInf ⌈ 1e¯12
  cf ← tauWall ÷ denomQ

  ⟨delta99, deltaStar, theta, shapeFactor, tauWall, cf⟩
}

# Extract upper surface wall-normal velocity profile at specified column colIdx
# Arguments:
#   ux: [Ny, Nx] horizontal velocity
#   uy: [Ny, Nx] vertical velocity
#   solidMask: [Ny, Nx] solid mask
#   colIdx: column index to sample
#   profileHeight: number of lattice cells normal to wall to sample
#
# Returns list:
#   ⟨yDistances, uTangential⟩
ExtractUpperProfile ⇐ {
  ⟨ux, uy, solidMask, colIdx, profileHeight⟩ ← 𝕩
  ⟨ny, nx⟩ ← ≢solidMask

  colMask ← colIdx ⊏˘ solidMask
  solidY ← (colMask = 1) ⊸/ (↕ny)
  nSolid ← ≠solidY

  # If solid exists in column, find highest solid cell (upper surface)
  hasSolid ← nSolid > 0
  safeSolidY ← hasSolid ⊑ ⟨⟨0⟩, solidY⟩
  yWall ← hasSolid ⊑ ⟨0, ⌈´ safeSolidY⟩

  # Sample upwards from yWall
  samples ← {
    step ← 𝕩
    ySample ← (ny - 1) ⌊ (yWall + step)
    coord ← ySample‿colIdx
    uVal ← coord ⊑ ux
    # Normal distance
    yDist ← 0.0 + step
    ⟨yDist, (step = 0) ⊑ ⟨uVal, 0.0⟩⟩
  }¨ (↕profileHeight)

  yDistances ← { 0 ⊑ 𝕩 }¨ samples
  uTangential ← { 1 ⊑ 𝕩 }¨ samples

  ⟨yDistances, uTangential⟩
}

# Sweep chordwise stations and return boundary layer parameters along upper surface
# Arguments:
#   ux: [Ny, Nx] horizontal velocity
#   uy: [Ny, Nx] vertical velocity
#   solidMask: [Ny, Nx] solid mask
#   xLead: chord leading edge x-coordinate
#   chord: reference chord length
#   nStations: number of chordwise stations
#   profileHeight: sampling height
#   nu: kinematic viscosity
#   rho0: reference density
#   uInf: freestream velocity
ChordwiseSweep ⇐ {
  ⟨ux, uy, solidMask, xLead, chord, nStations, profileHeight, nu, rho0, uInf⟩ ← 𝕩
  ⟨ny, nx⟩ ← ≢solidMask

  stepX ← (chord ÷ nStations) ⌈ 1.0

  sweepData ← {
    i ← 𝕩
    colX ← (nx - 2) ⌊ (0 ⌈ (⌊ (xLead + (i × stepX))))
    prof ← ExtractUpperProfile ⟨ux, uy, solidMask, colX, profileHeight⟩
    ⟨yDists, uTangs⟩ ← prof
    params ← IntegralParameters ⟨yDists, uTangs, nu, rho0, uInf⟩
    xFrac ← (((colX - xLead)) ÷ (chord ⌈ 1.0))
    ⟨xFrac⟩ ∾ params
  }¨ (↕nStations)

  sweepData
}
