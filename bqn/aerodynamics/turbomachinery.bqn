# vane/bqn/aerodynamics/turbomachinery.bqn
# Axial turbine and compressor blade cascade parametric geometry and periodic rasterizer.
# Note: Zero em dashes are used in this file per coding standards.

t ← •Import "../lib/tensor.bqn"
m ← •Import "../lib/math_utils.bqn"
Linspace ← t.Linspace
pi ← m.pi

# Generate Turbomachinery Blade Surface Coordinates
# chord: true chord length
# staggerDeg: stagger angle gamma (setting angle of chord line)
# camberAngleDeg: total flow deflection / camber turning angle theta = beta1 - beta2
# maxThick: maximum profile thickness ratio t/c
# nPoints: discretization stations
GenerateTurbomachineryBlade ⇐ {
  # 𝕩: ⟨chord, staggerDeg, camberAngleDeg, maxThick, nPoints⟩
  ⟨chord, staggerDeg, camberAngleDeg, maxThick, nPoints⟩ ← 𝕩
  
  staggerRad ← staggerDeg × (pi ÷ 180.0)
  camberRad ← camberAngleDeg × (pi ÷ 180.0)
  
  # Circular arc camber line radius R and center
  # theta = 2 * asin(c / (2*R)) -> R = c / (2 * sin(theta/2))
  halfCamber ← 0.5 × camberRad
  sinHalfCamber ← •math.Sin halfCamber
  safeSin ← 0.001 ⌈ sinHalfCamber
  radiusArc ← chord ÷ (2.0 × safeSin)
  
  # Coordinate parameter s along chord
  sNorm ← Linspace ⟨0.0, 1.0, nPoints⟩
  
  # Camber line height: yc(x) = sqrt(R^2 - (x - c/2)^2) - (R*cos(theta/2))
  cosHalfCamber ← •math.Cos halfCamber
  hBaseline ← radiusArc × cosHalfCamber
  
  yc ← {
    s ← 𝕩
    xFromCenter ← (s - 0.5) × chord
    termUnderRoot ← 0.0 ⌈ ((radiusArc × radiusArc) - (xFromCenter × xFromCenter))
    (√termUnderRoot) - hBaseline
  }¨ sNorm
  
  # Standard C4 thickness distribution for turbomachinery
  yt ← {
    s ← 𝕩
    maxThick × chord × (2.969 × (√s) - (1.260 × s) - (3.516 × s × s) + (2.843 × s × s × s) - (1.015 × s × s × s × s))
  }¨ sNorm
  
  # Upper and lower suction / pressure surfaces
  xUpper ← (sNorm × chord)
  yUpper ← yc + yt
  
  xLower ← (sNorm × chord)
  yLower ← yc - yt
  
  # Rotate by stagger angle gamma
  cosStagger ← •math.Cos staggerRad
  sinStagger ← •math.Sin staggerRad
  
  xuRot ← (xUpper × cosStagger) - (yUpper × sinStagger)
  yuRot ← (xUpper × sinStagger) + (yUpper × cosStagger)
  
  xlRot ← (xLower × cosStagger) - (yLower × sinStagger)
  ylRot ← (xLower × sinStagger) + (yLower × cosStagger)
  
  ⟨xuRot, yuRot, xlRot, ylRot⟩
}

# Rasterize Periodic Turbine Blade Cascade Row on Domain ny x nx
# Pitch spacing: pitchSpacing (vertical distance s between adjacent blades)
RasterizeBladeCascade2D ⇐ {
  # 𝕨: ⟨ny, nx, originX, chord, pitchSpacing, staggerDeg, camberAngleDeg, maxThick⟩
  # 𝕩: dummy
  ⟨ny, nx, originX, chord, pitchSpacing, staggerDeg, camberAngleDeg, maxThick⟩ ← 𝕨
  
  staggerRad ← staggerDeg × (pi ÷ 180.0)
  cosStagger ← •math.Cos (-staggerRad)
  sinStagger ← •math.Sin (-staggerRad)
  
  camberRad ← camberAngleDeg × (pi ÷ 180.0)
  halfCamber ← 0.5 × camberRad
  radiusArc ← chord ÷ (2.0 × (0.001 ⌈ •math.Sin halfCamber))
  hBaseline ← radiusArc × •math.Cos halfCamber
  
  # Vertical center of first blade
  yCenter0 ← 0.5 × pitchSpacing
  
  gridRows ← {
    iy ← 𝕩
    {
      ix ← 𝕩
      # Map iy into periodic cascade pitch interval [0, pitchSpacing)
      # dy to nearest blade center:
      yRel ← (iy - yCenter0) - (pitchSpacing × ⌊ ((iy - yCenter0) ÷ pitchSpacing))
      dyFromBladeCenter ← yRel - (0.5 × pitchSpacing)
      dxFromBladeLead ← ix - originX
      
      # Un-rotate to blade coordinates
      xBlade ← (dxFromBladeLead × cosStagger) - (dyFromBladeCenter × sinStagger)
      yBlade ← (dxFromBladeLead × sinStagger) + (dyFromBladeCenter × cosStagger)
      
      sNorm ← xBlade ÷ chord
      inChord ← (sNorm ≥ 0.0) ∧ (sNorm ≤ 1.0)
      
      safeS ← 0.0 ⌈ (1.0 ⌊ sNorm)
      xFromCenter ← (safeS - 0.5) × chord
      termUnderRoot ← 0.0 ⌈ ((radiusArc × radiusArc) - (xFromCenter × xFromCenter))
      ycLoc ← (√termUnderRoot) - hBaseline
      ytLoc ← maxThick × chord × (2.969 × (√safeS) - (1.260 × safeS) - (3.516 × safeS × safeS) + (2.843 × safeS × safeS × safeS) - (1.015 × safeS × safeS × safeS × safeS))
      
      isSolid ← inChord ∧ ((| (yBlade - ycLoc)) ≤ ytLoc)
      isSolid ⊑ ⟨0.0, 1.0⟩
    }¨ ↕nx
  }¨ ↕ny
  
  > gridRows
}
