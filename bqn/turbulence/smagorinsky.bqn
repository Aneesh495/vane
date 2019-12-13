# vane/bqn/turbulence/smagorinsky.bqn
# Algebraic Subgrid Scale Smagorinsky turbulence model formulated from non-equilibrium LBM moments.
# Note: Zero em dashes are used in this file per coding standards.

d2q9 ← •Import "../lbm/d2q9.bqn"
StressTensor ← d2q9.StressTensor
Moments ← d2q9.Moments

# Classical Smagorinsky Constant (typically 0.10 to 0.17 for aerodynamic flows)
defaultCs ⇐ 0.14

# Local Effective Relaxation Time via Algebraic Non-Equilibrium Stress Tensor
# Solves quadratic equation: tau_eff = tau0 + 0.5 * (sqrt(tau0^2 + 18*sqrt(2)*Cs^2*|Pi|/rho) - tau0)
ComputeTauEffective ⇐ {
  # 𝕨: ⟨tau0, csVal⟩
  # 𝕩: populations f
  ⟨tau0, csVal⟩ ← 𝕨
  f ← 𝕩
  
  ⟨rho, dummyUx, dummyUy⟩ ← Moments f
  ⟨piXX, piYY, piXY⟩ ← StressTensor f
  
  # Frobenius magnitude of non-equilibrium stress tensor: |Pi| = sqrt(Pi_xx^2 + Pi_yy^2 + 2*Pi_xy^2)
  piSq ← (piXX × piXX) + ((piYY × piYY) + (2.0 × (piXY × piXY)))
  piMag ← √piSq
  
  # Constant factor: 18 * sqrt(2) * Cs^2 = 25.45584412271571 * Cs^2
  constFactor ← 25.45584412271571 × (csVal × csVal)
  
  # Term under square root: tau0^2 + constFactor * |Pi| / rho
  termStress ← constFactor × (piMag ÷ (rho + 1e¯12))
  tau0Sq ← tau0 × tau0
  rootTerm ← √ (tau0Sq + termStress)
  
  tauEff ← tau0 + (0.5 × (rootTerm - tau0))
  tauEff
}

# Turbulent Eddy Viscosity nu_t
# nu_t = (tau_eff - tau0) / 3.0
ComputeEddyViscosity ⇐ {
  # 𝕨: ⟨tau0, csVal⟩
  # 𝕩: populations f
  ⟨tau0, csVal⟩ ← 𝕨
  tauEff ← 𝕨 ComputeTauEffective 𝕩
  nuT ← (tauEff - tau0) ÷ 3.0
  nuT
}

# Smagorinsky LES Time Step: Computes local tau_eff then applies BGK collision and streaming
StepSmagorinskyLES ⇐ {
  # 𝕨: ⟨tau0, csVal⟩
  # 𝕩: populations f
  tauEff ← 𝕨 ComputeTauEffective 𝕩
  tauEff d2q9.StepBGK 𝕩
}
