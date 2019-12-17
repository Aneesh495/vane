# vane/bqn/turbulence/wale.bqn
# Wall-Adapting Local Eddy-viscosity (WALE) model with asymptotic cubic wall scaling.
# Note: Zero em dashes are used in this file per coding standards.

st ← •Import "../lib/stencils.bqn"
d2q9 ← •Import "../lbm/d2q9.bqn"
Grad2DPeriodic ← st.Grad2DPeriodic
Moments ← d2q9.Moments

defaultCw ⇐ 0.325
eps ⇐ 1e¯12

# Compute WALE Turbulent Eddy Viscosity Field
ComputeWaleEddyViscosity ⇐ {
  # 𝕨: ⟨cwVal, dx⟩
  # 𝕩: populations f
  ⟨cwVal, dx⟩ ← 𝕨
  f ← 𝕩
  
  ⟨rho, ux, uy⟩ ← Moments f
  
  # Velocity gradients g_ij = du_i/dx_j
  ⟨duxDx, duxDy⟩ ← dx‿dx Grad2DPeriodic ux
  ⟨duyDx, duyDy⟩ ← dx‿dx Grad2DPeriodic uy
  
  # Symmetric strain rate tensor S_ij
  sxx ← duxDx
  syy ← duyDy
  sxy ← 0.5 × (duxDy + duyDx)
  
  # S_ij * S_ij contraction
  sSq ← (sxx × sxx) + ((syy × syy) + (2.0 × (sxy × sxy)))
  sMag5 ← sSq ⋆ 2.5
  
  # Square of velocity gradient tensor: g2_ij = g_ik * g_kj
  # In 2D:
  # g2_xx = duxDx*duxDx + duxDy*duyDx
  # g2_yy = duyDx*duxDy + duyDy*duyDy
  # g2_xy = duxDx*duxDy + duxDy*duyDy
  # g2_yx = duyDx*duxDx + duyDy*duyDx
  g2XX ← (duxDx × duxDx) + (duxDy × duyDx)
  g2YY ← (duyDx × duxDy) + (duyDy × duyDy)
  g2XY ← (duxDx × duxDy) + (duxDy × duyDy)
  g2YX ← (duyDx × duxDx) + (duyDy × duyDx)
  
  # Traceless symmetric tensor: S^d_ij = 0.5*(g2_ij + g2_ji) - (1/3)*delta_ij*tr(g2)
  trG2 ← g2XX + g2YY
  oneThirdTr ← (1.0 ÷ 3.0) × trG2
  
  sdXX ← g2XX - oneThirdTr
  sdYY ← g2YY - oneThirdTr
  sdXY ← 0.5 × (g2XY + g2YX)
  
  # Sd_ij * Sd_ij contraction
  sdSq ← (sdXX × sdXX) + ((sdYY × sdYY) + (2.0 × (sdXY × sdXY)))
  
  sdMag3 ← sdSq ⋆ 1.5
  sdMag2pt5 ← sdSq ⋆ 1.25
  
  # WALE eddy viscosity: (Cw * Delta)^2 * (Sd_ij * Sd_ij)^(3/2) / ((S_ij*S_ij)^(5/2) + (Sd_ij*Sd_ij)^(5/4))
  lengthScaleSq ← (cwVal × dx) × (cwVal × dx)
  denom ← sMag5 + (sdMag2pt5 + eps)
  
  nuT ← lengthScaleSq × (sdMag3 ÷ denom)
  nuT
}

# Step with WALE LES Turbulence Model
StepWaleLES ⇐ {
  # 𝕨: ⟨tau0, cwVal, dx⟩
  # 𝕩: populations f
  ⟨tau0, cwVal, dx⟩ ← 𝕨
  f ← 𝕩
  
  nuT ← ⟨cwVal, dx⟩ ComputeWaleEddyViscosity f
  tauEff ← tau0 + (3.0 × nuT)
  
  # Apply collision with local tauEff
  tauEff d2q9.StepBGK f
}
