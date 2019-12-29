# Surface Pressure Distribution and Cp Chordwise Integration Engine
#
# Extracts surface pressure profiles along upper (suction) and lower (pressure)
# airfoil surfaces, computes non-dimensional pressure coefficient Cp(x),
# and integrates normal force coefficient Cn and pitching moment Cm,p.
#
# Definition:
#   Cp(x) = (p(x) - p_inf) / (0.5 * rho_inf * u_inf^2)
#   p = cs2 * rho = (1/3) * rho in lattice Boltzmann hydrodynamics
#
# Note: Zero em dashes are used in this file per coding standards.

math ← •Import "../lib/math_utils.bqn"
tensor ← •Import "../lib/tensor.bqn"
d2q9 ← •Import "../lbm/d2q9.bqn"

# Extract surface pressure coefficient Cp distribution along chord
# Arguments:
#   rho: [Ny, Nx] fluid density field
#   solidMask: [Ny, Nx] solid obstacle mask (1 = solid, 0 = fluid)
#   xLead: chord leading edge x-coordinate (lattice units)
#   chord: reference chord length
#   rho0: freestream reference density
#   uInf: freestream reference velocity magnitude
#
# Returns list:
#   xNorm: normalized chordwise stations x / c (range [0, 1])
#   cpUpper: Cp values along suction surface
#   cpLower: Cp values along pressure surface
#   pUpper: dimensional pressure along suction surface
#   pLower: dimensional pressure along pressure surface
#   cn: integrated normal force coefficient
#   cmQuarter: integrated pitching moment about c/4
ExtractCp ⇐ {
  ⟨rho, solidMask, xLead, chord, rho0, uInf⟩ ← 𝕩
  ⟨ny, nx⟩ ← ≢solidMask

  cs2 ← 1.0 ÷ 3.0
  pInf ← cs2 × rho0
  qInf ← (0.5 × rho0) × (uInf × uInf)
  denom ← qInf ⌈ 1e¯12

  xTrail ← xLead + chord
  xStart ← 0 ⌈ (⌊ xLead)
  xEnd ← (nx - 1) ⌊ (⌈ xTrail)
  nSamples ← 1 ⌈ ((xEnd - xStart) + 1)

  # Sample columns across airfoil chord
  stationData ← {
    colIdx ← xStart + 𝕩
    colMask ← colIdx ⊏˘ solidMask

    solidY ← (colMask = 1) ⊸/ (↕ny)
    nSolid ← ≠solidY
    isSolid ← nSolid > 0

    # Guard against empty solid indices with default fallback
    safeSolidY ← isSolid ⊑ ⟨⟨0⟩, solidY⟩
    yMinSolid ← ⌊´ safeSolidY
    yMaxSolid ← ⌈´ safeSolidY

    yUpper ← (ny - 1) ⌊ (yMaxSolid + 1)
    yLower ← 0 ⌈ (yMinSolid - 1)

    # Densities at adjacent fluid layers
    rhoUpper ← (yUpper‿colIdx) ⊑ rho
    rhoLower ← (yLower‿colIdx) ⊑ rho

    pUpperLoc ← cs2 × rhoUpper
    pLowerLoc ← cs2 × rhoLower

    cpUp ← (pUpperLoc - pInf) ÷ denom
    cpLow ← (pLowerLoc - pInf) ÷ denom
    xFrac ← (((colIdx - xLead)) ÷ (chord ⌈ 1.0))

    # If not solid, zero out Cp and return base pressure
    ⟨
      xFrac,
      isSolid ⊑ ⟨0.0, cpUp⟩,
      isSolid ⊑ ⟨0.0, cpLow⟩,
      isSolid ⊑ ⟨pInf, pUpperLoc⟩,
      isSolid ⊑ ⟨pInf, pLowerLoc⟩,
      isSolid
    ⟩
  }¨ (↕nSamples)

  # Filter stations where solid was present
  validMask ← { 5 ⊑ 𝕩 }¨ stationData
  validIndices ← (validMask = 1) ⊸/ (↕nSamples)
  nValid ← ≠validIndices

  # Extract filtered vectors
  hasValid ← nValid > 1
  safeIndices ← hasValid ⊑ ⟨⟨0, 0⟩, validIndices⟩

  xNorm ← { 0 ⊑ (𝕩 ⊑ stationData) }¨ safeIndices
  cpUpper ← { 1 ⊑ (𝕩 ⊑ stationData) }¨ safeIndices
  cpLower ← { 2 ⊑ (𝕩 ⊑ stationData) }¨ safeIndices
  pUpper ← { 3 ⊑ (𝕩 ⊑ stationData) }¨ safeIndices
  pLower ← { 4 ⊑ (𝕩 ⊑ stationData) }¨ safeIndices

  # Trapezoidal integration for normal force and pitching moment
  dx ← (1 ↓ xNorm) - (¯1 ↓ xNorm)
  deltaCpMid ← 0.5 × (((1 ↓ cpLower) - (1 ↓ cpUpper)) + ((¯1 ↓ cpLower) - (¯1 ↓ cpUpper)))
  cnRaw ← +´ (deltaCpMid × dx)

  xMid ← 0.5 × ((1 ↓ xNorm) + (¯1 ↓ xNorm))
  armMid ← xMid - 0.25
  momentCpMid ← 0.5 × (((1 ↓ cpUpper) - (1 ↓ cpLower)) + ((¯1 ↓ cpUpper) - (¯1 ↓ cpLower)))
  cmRaw ← +´ ((momentCpMid × armMid) × dx)

  cn ← hasValid ⊑ ⟨0.0, cnRaw⟩
  cmQuarter ← hasValid ⊑ ⟨0.0, cmRaw⟩

  ⟨xNorm, cpUpper, cpLower, pUpper, pLower, cn, cmQuarter⟩
}
