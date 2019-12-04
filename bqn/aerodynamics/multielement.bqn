# vane/bqn/aerodynamics/multielement.bqn
# High-lift multi-element wing configurations: main wing with leading-edge slat and slotted flap.
# Note: Zero em dashes are used in this file per coding standards.

af ← •Import "airfoil_naca.bqn"
RasterizeAirfoil2D ← af.RasterizeAirfoil2D

# Multi-Element High-Lift Wing Rasterizer
# Generates composite discrete lattice mask representing:
# 1. Main Wing: chord c_main, origin (x0, y0), angle aoa
# 2. Leading-Edge Slat: chord c_slat, deflection delta_slat, slot gap, overlap
# 3. Trailing-Edge Flap: chord c_flap, deflection delta_flap, slot gap, overlap
RasterizeMultiElement2D ⇐ {
  # 𝕨: ⟨ny, nx, originX, originY, mainChord, aoaDeg, slatDeflectDeg, flapDeflectDeg⟩
  # 𝕩: dummy
  ⟨ny, nx, originX, originY, mainChord, aoaDeg, slatDeflectDeg, flapDeflectDeg⟩ ← 𝕨
  
  # 1. Main Wing Element (NACA 2412 baseline)
  mainMask ← ⟨ny, nx, originX, originY, mainChord, aoaDeg, 0.02, 0.4, 0.12⟩ RasterizeAirfoil2D 0
  
  # 2. Leading-Edge Slat
  # Slat chord is ~15% of main chord
  slatChord ← 0.15 × mainChord
  # Positioned ahead and below leading edge of main wing
  slatOriginX ← originX - (0.12 × mainChord)
  slatOriginY ← originY - (0.04 × mainChord)
  slatEffectiveAoa ← aoaDeg + slatDeflectDeg
  slatMask ← ⟨ny, nx, slatOriginX, slatOriginY, slatChord, slatEffectiveAoa, 0.04, 0.2, 0.10⟩ RasterizeAirfoil2D 0
  
  # 3. Trailing-Edge Slotted Flap
  # Flap chord is ~30% of main chord
  flapChord ← 0.30 × mainChord
  # Positioned near trailing edge of main wing (x = originX + 0.90 * mainChord)
  flapOriginX ← originX + (0.88 × mainChord)
  flapOriginY ← originY - (0.05 × mainChord)
  flapEffectiveAoa ← aoaDeg + flapDeflectDeg
  flapMask ← ⟨ny, nx, flapOriginX, flapOriginY, flapChord, flapEffectiveAoa, 0.03, 0.3, 0.12⟩ RasterizeAirfoil2D 0
  
  # Composite boolean union of all three elements (solid if any element is solid)
  compositeMask ← 1.0 ⌊ (mainMask + (slatMask + flapMask))
  compositeMask
}
