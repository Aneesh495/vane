# Q-Criterion Coherent Vortex Core Identification and Diagnostics Engine
#
# Identifies coherent turbulent vortices and separated shear layers using the
# second invariant of the velocity gradient tensor:
#   Q = 0.5 * ( ||Omega||^2 - ||S||^2 )
#
# Reference:
# Hunt, J. C. R., Wray, A. A., and Moin, P., "Eddies, streams, and convergence zones
# in turbulent flows", Center for Turbulence Research Report CTR-S88, 1988.
#
# Note: Zero em dashes are used in this file per coding standards.

math ← •Import "../lib/math_utils.bqn"
stencils ← •Import "../lib/stencils.bqn"
tensor ← •Import "../lib/tensor.bqn"

# Compute 2D Q-Criterion field from velocity components
# Arguments:
#   ux: [Ny, Nx] horizontal velocity field
#   uy: [Ny, Nx] vertical velocity field
#
# Returns:
#   qField: [Ny, Nx] second invariant Q
QCriterion2D ⇐ {
  ⟨ux, uy⟩ ← 𝕩

  # Rate-of-strain tensor components: [sxx, syy, sxy]
  ⟨sxx, syy, sxy⟩ ← 1.0‿1.0 stencils.StrainRate2D ⟨ux, uy⟩

  # Vorticity curl scalar omega_z = duy/dx - dux/dy
  vortZ ← 1.0‿1.0 stencils.Curl2D ⟨ux, uy⟩

  # In 2D: ||Omega||^2 = 0.5 * omega_z^2
  omegaSqNorm ← 0.5 × (vortZ × vortZ)

  # ||S||^2 = sxx^2 + syy^2 + 2 * sxy^2
  sSqNorm ← ((sxx × sxx) + (syy × syy)) + (2.0 × (sxy × sxy))

  # Q = 0.5 * ( ||Omega||^2 - ||S||^2 )
  qField ← 0.5 × (omegaSqNorm - sSqNorm)
  qField
}

# Normalized Q-Criterion Q* in range [-0.5, 0.5]
# Q* = Q / ( ||S||^2 + ||Omega||^2 + eps )
# Useful for thresholding across widely varying local shear magnitudes
NormalizedQ2D ⇐ {
  ⟨ux, uy⟩ ← 𝕩

  ⟨sxx, syy, sxy⟩ ← 1.0‿1.0 stencils.StrainRate2D ⟨ux, uy⟩
  vortZ ← 1.0‿1.0 stencils.Curl2D ⟨ux, uy⟩

  omegaSqNorm ← 0.5 × (vortZ × vortZ)
  sSqNorm ← ((sxx × sxx) + (syy × syy)) + (2.0 × (sxy × sxy))

  qNum ← 0.5 × (omegaSqNorm - sSqNorm)
  denom ← (sSqNorm + omegaSqNorm) ⌈ 1e¯12

  qNorm ← qNum ÷ denom
  qNorm
}

# Binary Vortex Core Identification Mask
# Arguments:
#   qField: [Ny, Nx] Q-criterion field
#   threshold: positive threshold for vortex core boundary (e.g. 1e-4)
#
# Returns:
#   mask: [Ny, Nx] boolean mask (1 = vortex core, 0 = background)
VortexCoreMask ⇐ {
  ⟨qField, threshold⟩ ← 𝕩
  qField > threshold
}

# Compute Vortex Core Area and Circulation Integral
# Arguments:
#   qField: [Ny, Nx] Q field
#   ux: [Ny, Nx] horizontal velocity
#   uy: [Ny, Nx] vertical velocity
#   threshold: Q threshold
#
# Returns list:
#   ⟨totalVortexArea, circulation⟩
#   circulation = integral(omega_z) over vortex cores
VortexCirculation ⇐ {
  ⟨qField, ux, uy, threshold⟩ ← 𝕩
  mask ← VortexCoreMask ⟨qField, threshold⟩
  vortZ ← 1.0‿1.0 stencils.Curl2D ⟨ux, uy⟩

  area ← +´ ⥊ mask
  circ ← +´ ⥊ (mask × vortZ)

  ⟨area, circ⟩
}
