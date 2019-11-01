# vane/bqn/lbm/d3q19.bqn
# D3Q19 three-dimensional Lattice Boltzmann hydrodynamic solver for 3D aerodynamic wakes.
# Note: Zero em dashes are used in this file per coding standards.

t ← •Import "../lib/tensor.bqn"
Roll3D ← t.Roll3D

# Lattice Weights
w0 ← 1.0 ÷ 3.0
w1 ← 1.0 ÷ 18.0
w2 ← 1.0 ÷ 36.0
w ⇐ ⟨w0, w1, w1, w1, w1, w1, w1, w2, w2, w2, w2, w2, w2, w2, w2, w2, w2, w2, w2⟩

# Discrete Velocity Vectors (19 directions)
cx ⇐ ⟨0,  1, ¯1,  0,  0,  0,  0,  1, ¯1,  1, ¯1,  1, ¯1,  1, ¯1,  0,  0,  0,  0⟩
cy ⇐ ⟨0,  0,  0,  1, ¯1,  0,  0,  1, ¯1, ¯1,  1,  0,  0,  0,  0,  1, ¯1,  1, ¯1⟩
cz ⇐ ⟨0,  0,  0,  0,  0,  1, ¯1,  0,  0,  0,  0,  1, ¯1, ¯1,  1,  1, ¯1, ¯1,  1⟩
opp ⇐ ⟨0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17⟩

cs2 ⇐ 1.0 ÷ 3.0
eps ⇐ 1e¯12

# Macroscopic Moments: Density and 3D Velocity Field
Moments ⇐ {
  # 𝕩: list of 19 population 3D tensors [f0, f1, ..., f18]
  f ← 𝕩
  
  # Density: rho = sum_i f_i
  rho ← +´ f
  
  # Momentum: rho * u_alpha = sum_i c_i_alpha * f_i
  momX ← +´ { (𝕩 ⊑ cx) × (𝕩 ⊑ f) }¨ ↕19
  momY ← +´ { (𝕩 ⊑ cy) × (𝕩 ⊑ f) }¨ ↕19
  momZ ← +´ { (𝕩 ⊑ cz) × (𝕩 ⊑ f) }¨ ↕19
  
  invRho ← 1.0 ÷ (rho + eps)
  ux ← momX × invRho
  uy ← momY × invRho
  uz ← momZ × invRho
  
  ⟨rho, ux, uy, uz⟩
}

# 3D Equilibrium Distribution Function f_i^eq(rho, ux, uy, uz)
Equilibrium ⇐ {
  # 𝕩: ⟨rho, ux, uy, uz⟩
  ⟨rho, ux, uy, uz⟩ ← 𝕩
  uSq ← (ux × ux) + ((uy × uy) + (uz × uz))
  
  feq ← {
    i ← 𝕩
    wi ← i ⊑ w
    cxi ← i ⊑ cx
    cyi ← i ⊑ cy
    czi ← i ⊑ cz
    
    ciDotU ← (cxi × ux) + ((cyi × uy) + (czi × uz))
    wi × rho × (1.0 + (3.0 × ciDotU) + (4.5 × ciDotU × ciDotU) - (1.5 × uSq))
  }¨ ↕19
  
  feq
}

# 3D BGK Collision Operator
CollideBGK ⇐ {
  # 𝕨: relaxation time tau
  # 𝕩: list of 19 population tensors f
  tau ← 𝕨
  f ← 𝕩
  mom ← Moments f
  feq ← Equilibrium mom
  invTau ← 1.0 ÷ tau
  
  fPost ← {
    i ← 𝕩
    fi ← i ⊑ f
    feqi ← i ⊑ feq
    fi - (invTau × (fi - feqi))
  }¨ ↕19
  
  fPost
}

# 3D Streaming Operator
Stream ⇐ {
  # 𝕩: list of 19 post-collision population tensors fPost
  fPost ← 𝕩
  fStreamed ← {
    i ← 𝕩
    cxi ← i ⊑ cx
    cyi ← i ⊑ cy
    czi ← i ⊑ cz
    czi‿cyi‿cxi Roll3D (i ⊑ fPost)
  }¨ ↕19
  
  fStreamed
}

# Full D3Q19 BGK Time Step
StepBGK ⇐ {
  # 𝕨: relaxation time tau
  # 𝕩: list of 19 population tensors f
  tau ← 𝕨
  f ← 𝕩
  fPost ← tau CollideBGK f
  Stream fPost
}

# 3D Non-Equilibrium Viscous Stress Tensor: 6 independent components
StressTensor3D ⇐ {
  # 𝕩: list of 19 population tensors f
  f ← 𝕩
  mom ← Moments f
  feq ← Equilibrium mom
  
  piXX ← +´ { i ← 𝕩 ⋄ (i ⊑ cx) × (i ⊑ cx) × ((i ⊑ f) - (i ⊑ feq)) }¨ ↕19
  piYY ← +´ { i ← 𝕩 ⋄ (i ⊑ cy) × (i ⊑ cy) × ((i ⊑ f) - (i ⊑ feq)) }¨ ↕19
  piZZ ← +´ { i ← 𝕩 ⋄ (i ⊑ cz) × (i ⊑ cz) × ((i ⊑ f) - (i ⊑ feq)) }¨ ↕19
  piXY ← +´ { i ← 𝕩 ⋄ (i ⊑ cx) × (i ⊑ cy) × ((i ⊑ f) - (i ⊑ feq)) }¨ ↕19
  piYZ ← +´ { i ← 𝕩 ⋄ (i ⊑ cy) × (i ⊑ cz) × ((i ⊑ f) - (i ⊑ feq)) }¨ ↕19
  piXZ ← +´ { i ← 𝕩 ⋄ (i ⊑ cx) × (i ⊑ cz) × ((i ⊑ f) - (i ⊑ feq)) }¨ ↕19
  
  ⟨piXX, piYY, piZZ, piXY, piYZ, piXZ⟩
}
