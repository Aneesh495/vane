# vane/bqn/aerodynamics/airfoil_naca.bqn
# Analytical NACA 4-digit airfoil geometry generation, angle-of-attack rotation, and grid rasterization.
# Note: Zero em dashes are used in this file per coding standards.

t ← •Import "../lib/tensor.bqn"
m ← •Import "../lib/math_utils.bqn"
Linspace ← t.Linspace
pi ← m.pi

# Generate 2D Surface Coordinates for NACA 4-Digit Profile
# m: max camber (e.g. 0.02 for NACA 2412)
# p: camber position (e.g. 0.4 for NACA 2412)
# th: thickness (e.g. 0.12 for NACA 2412)
# chord: chord length
# aoaDeg: angle of attack in degrees
# nPoints: number of surface chordwise discretization stations
GenerateNACA4 ⇐ {
  # 𝕩: ⟨camberM, camberP, thickT, chord, aoaDeg, nPoints⟩
  ⟨camberM, camberP, thickT, chord, aoaDeg, nPoints⟩ ← 𝕩
  
  # Cosine clustering of chordwise coordinates for fine resolution at leading edge
  beta ← Linspace ⟨0.0, pi, nPoints⟩
  xcNorm ← 0.5 × (1.0 - •math.Cos beta)
  xc ← chord × xcNorm
  
  # Half-thickness distribution:
  yt ← {
    x ← 𝕩
    5.0 × thickT × chord × ((0.2969 × √x) - (0.1260 × x) - (0.3516 × x × x) + (0.2843 × x × x × x) - (0.1015 × x × x × x × x))
  }¨ xcNorm
  
  # Camber line yc and camber slope dyc/dx
  pSafe ← 0.001 ⌈ camberP
  ycAndTheta ← {
    x ← 𝕩
    isBeforeP ← x ≤ pSafe
    
    yc1 ← (camberM ÷ (pSafe × pSafe)) × ((2.0 × pSafe × x) - (x × x))
    dy1 ← (2.0 × camberM ÷ (pSafe × pSafe)) × (pSafe - x)
    
    oneMinusP ← 1.0 - pSafe
    yc2 ← (camberM ÷ (oneMinusP × oneMinusP)) × ((1.0 - (2.0 × pSafe)) + (2.0 × pSafe × x) - (x × x))
    dy2 ← (2.0 × camberM ÷ (oneMinusP × oneMinusP)) × (pSafe - x)
    
    ycVal ← (camberM = 0.0) ⊑ ⟨isBeforeP ⊑ ⟨yc2, yc1⟩, 0.0⟩
    dyVal ← (camberM = 0.0) ⊑ ⟨isBeforeP ⊑ ⟨dy2, dy1⟩, 0.0⟩
    
    thVal ← •math.Atan dyVal
    ⟨chord × ycVal, thVal⟩
  }¨ xcNorm
  
  yc ← { 0 ⊑ 𝕩 }¨ ycAndTheta
  thetaAngle ← { 1 ⊑ 𝕩 }¨ ycAndTheta
  
  # Upper and lower surface unrotated coordinates
  sinTh ← •math.Sin¨ thetaAngle
  cosTh ← •math.Cos¨ thetaAngle
  
  xu ← xc - (yt × sinTh)
  yu ← yc + (yt × cosTh)
  
  xl ← xc + (yt × sinTh)
  yl ← yc - (yt × cosTh)
  
  # Rotate by angle of attack alpha (aoaDeg)
  alphaRad ← aoaDeg × (pi ÷ 180.0)
  cosA ← •math.Cos alphaRad
  sinA ← •math.Sin alphaRad
  
  xqc ← 0.25 × chord
  
  xuRot ← xqc + (((xu - xqc) × cosA) + (yu × sinA))
  yuRot ← (-(xu - xqc) × sinA) + (yu × cosA)
  
  xlRot ← xqc + (((xl - xqc) × cosA) + (yl × sinA))
  ylRot ← (-(xl - xqc) × sinA) + (yl × cosA)
  
  ⟨xuRot, yuRot, xlRot, ylRot⟩
}

# Rasterize Airfoil Geometry into Boolean Discrete Lattice Obstacle Matrix
RasterizeAirfoil2D ⇐ {
  # 𝕨: ⟨ny, nx, originX, originY, chord, aoaDeg, camberM, camberP, thickT⟩
  ⟨ny, nx, originX, originY, chord, aoaDeg, camberM, camberP, thickT⟩ ← 𝕨
  
  # Generate high-resolution airfoil coordinates
  ⟨xu, yu, xl, yl⟩ ← GenerateNACA4 ⟨camberM, camberP, thickT, chord, aoaDeg, 200⟩
  
  xuShift ← originX + xu
  yuShift ← originY + yu
  xlShift ← originX + xl
  ylShift ← originY + yl
  
  minX ← ⌊´ (xuShift ∾ xlShift)
  maxX ← ⌈´ (xuShift ∾ xlShift)
  minY ← ⌊´ (yuShift ∾ ylShift)
  maxY ← ⌈´ (yuShift ∾ ylShift)
  
  alphaRad ← aoaDeg × (pi ÷ 180.0)
  cosA ← •math.Cos (-alphaRad)
  sinA ← •math.Sin (-alphaRad)
  xqc ← originX + (0.25 × chord)
  yqc ← originY
  
  gridRows ← {
    iy ← 𝕩
    {
      ix ← 𝕩
      dx ← ix - xqc
      dy ← iy - yqc
      xUnrot ← xqc + ((dx × cosA) - (dy × sinA))
      yUnrot ← yqc + ((dx × sinA) + (dy × cosA))
      xNorm ← (xUnrot - originX) ÷ chord
      
      inBox ← (ix ≥ minX) ∧ ((ix ≤ maxX) ∧ ((iy ≥ minY) ∧ ((iy ≤ maxY) ∧ ((xNorm ≥ 0.0) ∧ (xNorm ≤ 1.0)))))
      
      safeX ← 0.0 ⌈ (1.0 ⌊ xNorm)
      ytLoc ← 5.0 × thickT × chord × ((0.2969 × √safeX) - (0.1260 × safeX) - (0.3516 × safeX × safeX) + (0.2843 × safeX × safeX × safeX) - (0.1015 × safeX × safeX × safeX × safeX))
      
      pSafe ← 0.001 ⌈ camberP
      isBeforeP ← safeX ≤ pSafe
      yc1 ← (camberM ÷ (pSafe × pSafe)) × ((2.0 × pSafe × safeX) - (safeX × safeX))
      oneMinusP ← 1.0 - pSafe
      yc2 ← (camberM ÷ (oneMinusP × oneMinusP)) × ((1.0 - (2.0 × pSafe)) + (2.0 × pSafe × safeX) - (safeX × safeX))
      ycLoc ← originY + (chord × ((camberM = 0.0) ⊑ ⟨isBeforeP ⊑ ⟨yc2, yc1⟩, 0.0⟩))
      
      isSolid ← inBox ∧ ((| (yUnrot - ycLoc)) ≤ ytLoc)
      isSolid ⊑ ⟨0.0, 1.0⟩
    }¨ ↕nx
  }¨ ↕ny
  
  > gridRows
}
