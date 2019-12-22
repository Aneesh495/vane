# vane/bqn/turbulence/dynamic_les.bqn
# Dynamic Germano-Lilly Subgrid Scale model with local test-filtering.
# Note: Zero em dashes are used in this file per coding standards.

t ← •Import "../lib/tensor.bqn"
st ← •Import "../lib/stencils.bqn"
d2q9 ← •Import "../lbm/d2q9.bqn"
Roll2D ← t.Roll2D
Grad2DPeriodic ← st.Grad2DPeriodic
Moments ← d2q9.Moments

# Discrete 2D Top-Hat Test Filter (Filter width 2*Delta)
# Stencil weights: center = 1/4, 4 cardinal neighbors = 1/8, 4 diagonals = 1/16
TestFilter2D ⇐ {
  # 𝕩: 2D field phi
  phi ← 𝕩
  
  # Center
  wCenter ← 0.25 × phi
  
  # Cardinal neighbors
  wEast ← 0.125 × (0‿¯1 Roll2D phi)
  wWest ← 0.125 × (0‿1 Roll2D phi)
  wNorth ← 0.125 × (¯1‿0 Roll2D phi)
  wSouth ← 0.125 × (1‿0 Roll2D phi)
  
  # Diagonal neighbors
  wNE ← 0.0625 × (¯1‿¯1 Roll2D phi)
  wNW ← 0.0625 × (¯1‿1 Roll2D phi)
  wSE ← 0.0625 × (1‿¯1 Roll2D phi)
  wSW ← 0.0625 × (1‿1 Roll2D phi)
  
  wCenter + (wEast + (wWest + (wNorth + (wSouth + (wNE + (wNW + (wSE + wSW)))))))
}

# Dynamic Calculation of Local Smagorinsky Coefficient Cs^2(x, y)
ComputeDynamicCs2 ⇐ {
  # 𝕨: grid spacing dx
  # 𝕩: populations f
  dx ← 𝕨
  f ← 𝕩
  
  ⟨rho, ux, uy⟩ ← Moments f
  
  # Filtered velocities: u_hat
  uxHat ← TestFilter2D ux
  uyHat ← TestFilter2D uy
  
  # Leonard stress tensor: L_ij = (u_i * u_j)_hat - u_i_hat * u_j_hat
  lXX ← (TestFilter2D (ux × ux)) - (uxHat × uxHat)
  lYY ← (TestFilter2D (uy × uy)) - (uyHat × uyHat)
  lXY ← (TestFilter2D (ux × uy)) - (uxHat × uyHat)
  
  # Grid-scale strain rate: S_ij
  ⟨duxDx, duxDy⟩ ← dx‿dx Grad2DPeriodic ux
  ⟨duyDx, duyDy⟩ ← dx‿dx Grad2DPeriodic uy
  sXX ← duxDx
  sYY ← duyDy
  sXY ← 0.5 × (duxDy + duyDx)
  sMag ← √ ((2.0 × ((sXX × sXX) + (sYY × sYY))) + (4.0 × (sXY × sXY)))
  
  # Test-scale strain rate: S_hat_ij
  ⟨duxHatDx, duxHatDy⟩ ← dx‿dx Grad2DPeriodic uxHat
  ⟨duyHatDx, duyHatDy⟩ ← dx‿dx Grad2DPeriodic uyHat
  sHatXX ← duxHatDx
  sHatYY ← duyHatDy
  sHatXY ← 0.5 × (duyHatDx + duxHatDy)
  sHatMag ← √ ((2.0 × ((sHatXX × sHatXX) + (sHatYY × sHatYY))) + (4.0 × (sHatXY × sHatXY)))
  
  # Germano model tensor: M_ij = 2 * dx^2 * ( (|S|*S_ij)_hat - 4 * |S_hat| * S_hat_ij )
  mXX ← 2.0 × dx × dx × ((TestFilter2D (sMag × sXX)) - (4.0 × (sHatMag × sHatXX)))
  mYY ← 2.0 × dx × dx × ((TestFilter2D (sMag × sYY)) - (4.0 × (sHatMag × sHatYY)))
  mXY ← 2.0 × dx × dx × ((TestFilter2D (sMag × sXY)) - (4.0 × (sHatMag × sHatXY)))
  
  # Lilly contraction: Cs^2 = 0.5 * (L_ij * M_ij) / (M_ij * M_ij)
  lDotM ← (lXX × mXX) + ((lYY × mYY) + (2.0 × (lXY × mXY)))
  mDotM ← (mXX × mXX) + ((mYY × mYY) + (2.0 × (mXY × mXY)))
  
  cs2Raw ← 0.5 × (lDotM ÷ (mDotM + 1e¯12))
  
  # Clip to [0, 0.25] to prevent backscatter instability and unphysical damping
  cs2Clipped ← 0.0 ⌈ (0.05 ⌊ cs2Raw)
  cs2Clipped
}

# Step with Dynamic Smagorinsky LES Model
StepDynamicLES ⇐ {
  # 𝕨: ⟨tau0, dx⟩
  # 𝕩: populations f
  ⟨tau0, dx⟩ ← 𝕨
  f ← 𝕩
  
  cs2Local ← dx ComputeDynamicCs2 f
  csLocal ← √cs2Local
  
  # Compute local tau_eff with dynamic Cs
  smag ← •Import "smagorinsky.bqn"
  tauEff ← ⟨tau0, csLocal⟩ smag.ComputeTauEffective f
  tauEff d2q9.StepBGK f
}
