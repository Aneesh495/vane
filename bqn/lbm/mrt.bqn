# vane/bqn/lbm/mrt.bqn
# Multiple Relaxation Time (MRT) collision model for robust high Reynolds number aerodynamics.
# Note: Zero em dashes are used in this file per coding standards.

d2q9 ← •Import "d2q9.bqn"
la ← •Import "../lib/linalg.bqn"
t ← •Import "../lib/tensor.bqn"

# d'Humieres D2Q9 Orthogonal Transformation Matrix (9x9)
mMat ⇐ 9‿9 ⥊ ⟨
  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
 ¯4.0, ¯1.0, ¯1.0, ¯1.0, ¯1.0,  2.0,  2.0,  2.0,  2.0,
  4.0, ¯2.0, ¯2.0, ¯2.0, ¯2.0,  1.0,  1.0,  1.0,  1.0,
  0.0,  1.0,  0.0, ¯1.0,  0.0,  1.0, ¯1.0, ¯1.0,  1.0,
  0.0, ¯2.0,  0.0,  2.0,  0.0,  1.0, ¯1.0, ¯1.0,  1.0,
  0.0,  0.0,  1.0,  0.0, ¯1.0,  1.0,  1.0, ¯1.0, ¯1.0,
  0.0,  0.0, ¯2.0,  0.0,  2.0,  1.0,  1.0, ¯1.0, ¯1.0,
  0.0,  1.0, ¯1.0,  1.0, ¯1.0,  0.0,  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,  0.0,  0.0,  1.0, ¯1.0,  1.0, ¯1.0
⟩

# Analytical Inverse Transformation Matrix M^-1 = M^T * D^-1
rowNormSq ← { +´ 𝕩 × 𝕩 }¨ <˘ mMat
invDiag ← 1.0 ÷ rowNormSq
mInv ⇐ ⍉ (> (<˘ mMat) × invDiag)

# Equilibrium Moments m_k^eq(rho, ux, uy)
EquilibriumMoments ⇐ {
  # 𝕩: ⟨rho, ux, uy⟩
  ⟨rho, ux, uy⟩ ← 𝕩
  uxSq ← ux × ux
  uySq ← uy × uy
  uSq ← uxSq + uySq
  
  # 9 equilibrium moments
  # m0: density
  m0 ← rho
  # m1: energy e
  m1 ← rho × (¯2.0 + (3.0 × uSq))
  # m2: energy squared epsilon
  m2 ← rho × (1.0 - (3.0 × uSq))
  # m3: momentum jx
  m3 ← rho × ux
  # m4: heat flux qx
  m4 ← - rho × ux
  # m5: momentum jy
  m5 ← rho × uy
  # m6: heat flux qy
  m6 ← - rho × uy
  # m7: normal stress p_xx
  m7 ← rho × (uxSq - uySq)
  # m8: shear stress p_xy
  m8 ← rho × ux × uy
  
  ⟨m0, m1, m2, m3, m4, m5, m6, m7, m8⟩
}

# Forward Transform: populations f (list of 9 matrices) -> moments m (list of 9 matrices)
# m_k = sum_i M_ki * f_i
FToM ⇐ {
  f ← 𝕩
  moments ← {
    k ← 𝕩
    rowK ← k ⊏ mMat
    +´ { (𝕩 ⊑ rowK) × (𝕩 ⊑ f) }¨ ↕9
  }¨ ↕9
  moments
}

# Inverse Transform: moments m (list of 9 matrices) -> populations f (list of 9 matrices)
# f_i = sum_k (M^-1)_ik * m_k
MToF ⇐ {
  m ← 𝕩
  populations ← {
    i ← 𝕩
    rowI ← i ⊏ mInv
    +´ { (𝕩 ⊑ rowI) × (𝕩 ⊑ m) }¨ ↕9
  }¨ ↕9
  populations
}

# Tuned Relaxation Rates Vector S
# w: kinematic viscosity relaxation tau_nu
DefaultRelaxationRates ⇐ {
  tauNu ← 𝕩
  sNu ← 1.0 ÷ tauNu
  # Conserved moments have s = 0
  # Non-hydrodynamic moments are tuned for optimal numerical stability
  s0 ← 0.0
  s1 ← 1.63
  s2 ← 1.14
  s3 ← 0.0
  s4 ← 1.92
  s5 ← 0.0
  s6 ← 1.92
  s7 ← sNu
  s8 ← sNu
  ⟨s0, s1, s2, s3, s4, s5, s6, s7, s8⟩
}

# Multiple Relaxation Time Collision Operator
CollideMRT ⇐ {
  # 𝕨: relaxation rates vector S (or tauNu)
  # 𝕩: populations f
  isScalar ← (≠𝕨) ≤ 1 ⋄ rates ← isScalar ⊑ ⟨𝕨, DefaultRelaxationRates 𝕨⟩
  f ← 𝕩
  
  # Extract macroscopic moments
  mom ← d2q9.Moments f
  meq ← EquilibriumMoments mom
  mRaw ← FToM f
  
  # Relax moments: m* = m - S * (m - meq)
  mPost ← {
    k ← 𝕩
    sk ← k ⊑ rates
    mk ← k ⊑ mRaw
    meqk ← k ⊑ meq
    mk - (sk × (mk - meqk))
  }¨ ↕9
  
  # Project back to population space: f* = M^-1 * m*
  MToF mPost
}

# Full D2Q9 MRT Time Step: Collision followed by Streaming
StepMRT ⇐ {
  # 𝕨: relaxation rates vector or tauNu
  # 𝕩: populations f
  rates ← 𝕨
  f ← 𝕩
  fPost ← rates CollideMRT f
  d2q9.Stream fPost
}
