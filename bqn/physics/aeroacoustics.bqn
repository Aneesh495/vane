# Computational Aeroacoustics (CAA) Lighthill and Curle Analogy Engine
#
# Computes acoustic source terms and far-field radiated sound from aerodynamic flow:
#   - Lighthill stress tensor T_ij = rho * u_i * u_j - 2 * mu * S_ij
#   - Curle surface dipole acoustic radiation to far-field observer locations
#   - Sound Pressure Level (SPL) spectrum in decibels (dB)
#   - Far-field acoustic directivity pattern SPL(theta)
#
# References:
# Lighthill, M. J., "On sound generated aerodynamically. I. General theory",
# Proc. R. Soc. Lond. A 211, 1952.
# Curle, N., "The influence of solid boundaries upon aerodynamic sound",
# Proc. R. Soc. Lond. A 231, 1955.
#
# Note: Zero em dashes are used in this file per coding standards.

math ← •Import "../lib/math_utils.bqn"
stencils ← •Import "../lib/stencils.bqn"

# Sound speed in lattice units: cs = 1 / sqrt(3)
c0 ⇐ 1.0 ÷ (3.0 ⋆ 0.5)

# Reference acoustic pressure threshold: p_ref (lattice units)
pRef ⇐ 1e¯5

# Compute Lighthill Stress Tensor T_ij
# Arguments:
#   rho: [Ny, Nx] density field
#   ux: [Ny, Nx] x-velocity
#   uy: [Ny, Nx] y-velocity
#   nu: kinematic viscosity
#
# Returns list:
#   ⟨Txx, Txy, Tyy⟩ components of Lighthill tensor
LighthillTensor ⇐ {
  ⟨rho, ux, uy, nu⟩ ← 𝕩

  # Compute strain rate tensor components with grid spacing dx=1, dy=1
  ⟨sxx, syy, sxy⟩ ← 1.0‿1.0 stencils.StrainRate2D ⟨ux, uy⟩

  # Viscous stress: tau_ij = 2 * mu * S_ij = 2 * rho * nu * S_ij
  mu ← rho × nu
  twoMu ← 2.0 × mu
  tauXx ← twoMu × sxx
  tauXy ← twoMu × sxy
  tauYy ← twoMu × syy

  # Reynolds stress: rho * u_i * u_j
  reyXx ← rho × (ux × ux)
  reyXy ← rho × (ux × uy)
  reyYy ← rho × (uy × uy)

  # Lighthill tensor: T_ij = rho * u_i * u_j - tau_ij
  txx ← reyXx - tauXx
  txy ← reyXy - tauXy
  tyy ← reyYy - tauYy

  ⟨txx, txy, tyy⟩
}

# Curle Far-Field Dipole Radiated Acoustic Pressure
# Computes acoustic sound pressure at observer position xObs from time series of forces
# Arguments:
#   timeHistoryF: [N, 2] list of force vectors ⟨Fx(t), Fy(t)⟩ over time
#   dt: time step
#   xObs: [2] observer position relative to acoustic source center [x_obs, y_obs]
#
# Returns list:
#   pPrime: [N-1] radiated acoustic pressure time series at observer
CurleDipoleFarField ⇐ {
  ⟨timeHistoryF, dt, xObs⟩ ← 𝕩

  rDist ← (((0 ⊑ xObs) ⋆ 2) + ((1 ⊑ xObs) ⋆ 2)) ⋆ 0.5
  denomR ← rDist ⌈ 1e¯12
  # Unit observer direction vector: [nx, ny]
  nObsX ← (0 ⊑ xObs) ÷ denomR
  nObsY ← (1 ⊑ xObs) ÷ denomR

  # Time derivative of aerodynamic force: dF/dt via central difference
  nSteps ← ≠timeHistoryF
  safeSteps ← (nSteps > 2) ⊑ ⟨2, nSteps⟩

  fxList ← { 0 ⊑ 𝕩 }¨ timeHistoryF
  fyList ← { 1 ⊑ 𝕩 }¨ timeHistoryF

  # First derivative dFx/dt and dFy/dt
  dfxDt ← ((1 ↓ fxList) - (¯1 ↓ fxList)) ÷ dt
  dfyDt ← ((1 ↓ fyList) - (¯1 ↓ fyList)) ÷ dt

  # Dipole source projection: (n_obs . dF/dt)
  dipoleProj ← (nObsX × dfxDt) + (nObsY × dfyDt)

  # Far-field 2D Curle prefactor: 1 / (4 * pi * c0 * r)
  pi ← 3.14159265358979323846
  prefactor ← 1.0 ÷ ((4.0 × pi) × (c0 × denomR))

  # Radiated acoustic pressure fluctuation p'(t)
  pPrime ← prefactor × dipoleProj
  pPrime
}

# Sound Pressure Level (SPL) Spectrum and Overall SPL (OASPL)
# Arguments:
#   pSignal: [N] acoustic pressure signal p'(t)
#
# Returns list:
#   ⟨oaspl, rmsPressure⟩
#   oaspl: Overall Sound Pressure Level in decibels (dB re pRef)
#   rmsPressure: root-mean-square acoustic pressure
AcousticSPL ⇐ {
  pSignal ← 𝕩
  n ← ≠pSignal

  # Mean and root mean square pressure
  pMean ← (+´ pSignal) ÷ (n ⌈ 1)
  pFluc ← pSignal - pMean
  pSqMean ← (+´ (pFluc × pFluc)) ÷ (n ⌈ 1)
  pRMS ← pSqMean ⋆ 0.5

  # OASPL = 20 * log10(pRMS / pRef) = 10 * log10(pSqMean / pRef^2)
  # In BQN, log10(x) = Ln(x) / Ln(10)
  ratio ← (pRMS ÷ pRef) ⌈ 1e¯12
  lnRatio ← ⋆⁼ ratio
  ln10 ← ⋆⁼ 10.0
  log10Ratio ← lnRatio ÷ ln10
  oaspl ← 20.0 × log10Ratio

  ⟨oaspl, pRMS⟩
}

# Far-field Acoustic Directivity Pattern SPL(theta)
# Evaluates radiated acoustic pressure amplitude across 360 degrees
# Arguments:
#   timeHistoryF: [N, 2] list of forces ⟨Fx(t), Fy(t)⟩
#   dt: time step
#   radius: observer radial distance R
#   nAngles: number of angular observer stations (e.g. 36)
#
# Returns list:
#   ⟨anglesRad, splTheta⟩
DirectivityPattern ⇐ {
  ⟨timeHistoryF, dt, radius, nAngles⟩ ← 𝕩
  pi ← 3.14159265358979323846
  twoPi ← 2.0 × pi

  angles ← { (𝕩 × twoPi) ÷ nAngles }¨ (↕nAngles)

  splPattern ← {
    theta ← 𝕩
    xObs ← ⟨radius × (math.Cos theta), radius × (math.Sin theta)⟩
    pPrime ← CurleDipoleFarField ⟨timeHistoryF, dt, xObs⟩
    splRes ← AcousticSPL pPrime
    0 ⊑ splRes
  }¨ angles

  ⟨angles, splPattern⟩
}
