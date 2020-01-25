# Jeong-Hussain Lambda-2 Vortex Core Identification Engine
#
# Identifies coherent vortex cores based on the second eigenvalue of the tensor:
#   M = S^2 + Omega^2
# A vortex core is present where lambda_2 < 0.
#
# Reference:
# Jeong, J. and Hussain, F., "On the identification of a vortex",
# Journal of Fluid Mechanics, Vol. 285, pp. 69-94, 1995.
#
# Note: Zero em dashes are used in this file per coding standards.

math ← •Import "../lib/math_utils.bqn"
stencils ← •Import "../lib/stencils.bqn"
linalg ← •Import "../lib/linalg.bqn"

# Compute 2D Lambda-2 scalar field from velocity components
# Arguments:
#   ux: [Ny, Nx] horizontal velocity
#   uy: [Ny, Nx] vertical velocity
#
# Returns:
#   lambda2: [Ny, Nx] second eigenvalue field (vortex cores where lambda2 < 0)
Lambda2_2D ⇐ {
  ⟨ux, uy⟩ ← 𝕩

  # S_ij: [sxx, syy, sxy]
  ⟨sxx, syy, sxy⟩ ← 1.0‿1.0 stencils.StrainRate2D ⟨ux, uy⟩

  # Vorticity curl omega_z = duy/dx - dux/dy
  # Omega_xy = 0.5 * omega_z
  vortZ ← 1.0‿1.0 stencils.Curl2D ⟨ux, uy⟩
  omegaXy ← 0.5 × vortZ
  omegaSq ← omegaXy × omegaXy

  # Matrix M = S^2 + Omega^2
  # M_xx = sxx^2 + sxy^2 - omega_xy^2
  # M_yy = sxy^2 + syy^2 - omega_xy^2
  # M_xy = sxy * (sxx + syy)
  mxx ← ((sxx × sxx) + (sxy × sxy)) - omegaSq
  myy ← ((sxy × sxy) + (syy × syy)) - omegaSq
  mxy ← sxy × (sxx + syy)

  # Analytical eigenvalues of 2x2 symmetric matrix M:
  # lambda = (mxx + myy)/2 +/- sqrt( ((mxx - myy)/2)^2 + mxy^2 )
  traceHalf ← 0.5 × (mxx + myy)
  diffHalf ← 0.5 × (mxx - myy)
  disc ← ((diffHalf × diffHalf) + (mxy × mxy)) ⋆ 0.5

  # Second eigenvalue lambda_2 is the smaller eigenvalue:
  lambda2 ← traceHalf - disc
  lambda2
}

# Binary Vortex Core Identification Mask for Lambda-2 criterion
# Arguments:
#   lambda2Field: [Ny, Nx] lambda-2 field
#   threshold: negative threshold (e.g. -1e-4) or 0.0
#
# Returns:
#   mask: [Ny, Nx] boolean mask (1 = vortex core, 0 = background)
VortexCoreMask ⇐ {
  ⟨lambda2Field, threshold⟩ ← 𝕩
  lambda2Field < threshold
}

# Summary Statistics of Identified Vortex Structures
# Arguments:
#   lambda2Field: [Ny, Nx] lambda-2 field
#   threshold: negative threshold for vortex detection
#
# Returns list:
#   ⟨vortexArea, minLambda2, meanCoreIntensity⟩
VortexStatistics ⇐ {
  ⟨lambda2Field, threshold⟩ ← 𝕩
  mask ← VortexCoreMask ⟨lambda2Field, threshold⟩

  area ← +´ ⥊ mask
  minVal ← ⌊´ ⥊ lambda2Field

  # Mean intensity inside core
  coreSum ← +´ ⥊ (mask × lambda2Field)
  denomArea ← area ⌈ 1
  meanIntensity ← coreSum ÷ (0.0 + denomArea)

  ⟨area, minVal, meanIntensity⟩
}
