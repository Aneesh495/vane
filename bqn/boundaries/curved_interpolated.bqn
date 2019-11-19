# vane/bqn/boundaries/curved_interpolated.bqn
# Bouzidi-Firdaouss-Lallemand (BFL) interpolated bounce-back for sub-grid curved aerodynamic surfaces.
# Note: Zero em dashes are used in this file per coding standards.

d2q9 ← •Import "../lbm/d2q9.bqn"
t ← •Import "../lib/tensor.bqn"
Roll2D ← t.Roll2D
opp ← d2q9.opp
cx ← d2q9.cx
cy ← d2q9.cy

# Bouzidi Interpolated Boundary Condition
# Given:
# - fPost: post-collision population matrices
# - isBoundaryFluid: boolean mask of fluid nodes adjacent to solid wall
# - qFraction: matrix of link distance fractions q in (0, 1) to solid surface
ApplyBouzidiBoundary ⇐ {
  # 𝕨: ⟨isBoundaryFluid, qFraction⟩
  # 𝕩: populations fPost
  ⟨isBound, qMat⟩ ← 𝕨
  fPost ← 𝕩
  
  fOut ← {
    i ← 𝕩
    iOpp ← i ⊑ opp
    cxi ← i ⊑ cx
    cyi ← i ⊑ cy
    
    fi ← i ⊑ fPost
    fiOpp ← iOpp ⊑ fPost
    
    # Upstream neighbor f_i*(x_f - c_i)
    fiNeighbor ← (cyi‿cxi Roll2D fi)
    
    # Bouzidi scheme branching on q < 0.5 vs q >= 0.5
    # For q < 0.5: f_opp = 2*q*fi + (1 - 2*q)*fiNeighbor
    # For q >= 0.5: f_opp = (1/(2*q))*fi + ((2*q - 1)/(2*q))*fiOpp
    safeQ ← 0.001 ⌈ (0.999 ⌊ qMat)
    
    termLowQ ← (2.0 × safeQ × fi) + ((1.0 - (2.0 × safeQ)) × fiNeighbor)
    
    invTwoQ ← 1.0 ÷ (2.0 × safeQ)
    termHighQ ← (invTwoQ × fi) + ((1.0 - invTwoQ) × fiOpp)
    
    isLowQ ← safeQ < 0.5
    fInterpolated ← (isLowQ × termLowQ) + ((1.0 - isLowQ) × termHighQ)
    
    # Apply only at boundary fluid nodes
    (isBound × fInterpolated) + ((1.0 - isBound) × fi)
  }¨ ↕9
  
  fOut
}
