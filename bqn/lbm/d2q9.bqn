# vane/bqn/lbm/d2q9.bqn
# D2Q9 Lattice Boltzmann hydrodynamics engine with BGK collision and streaming.
# Note: Zero em dashes are used in this file per coding standards.

t ← •Import "../lib/tensor.bqn"
Roll2D ← t.Roll2D

# Lattice Parameters
# 9 discrete velocity directions
w ⇐ ⟨4.0÷9.0, 1.0÷9.0, 1.0÷9.0, 1.0÷9.0, 1.0÷9.0, 1.0÷36.0, 1.0÷36.0, 1.0÷36.0, 1.0÷36.0⟩
cx ⇐ ⟨0, 1, 0, ¯1, 0, 1, ¯1, ¯1, 1⟩
cy ⇐ ⟨0, 0, 1, 0, ¯1, 1, 1, ¯1, ¯1⟩
opp ⇐ ⟨0, 3, 4, 1, 2, 7, 8, 5, 6⟩

cs2 ⇐ 1.0 ÷ 3.0
invCs2 ⇐ 3.0
inv2Cs4 ⇐ 4.5
eps ⇐ 1e¯12

# Macroscopic Density and Velocity Moments
Moments ⇐ {
  # 𝕩: list of 9 population matrices [f0, f1, ..., f8]
  f ← 𝕩
  
  # Density: rho = sum(f_i)
  rho ← +´ f
  
  # Momentum: rho * u = sum(c_i * f_i)
  momX ← +´ { (𝕩 ⊑ cx) × (𝕩 ⊑ f) }¨ ↕9
  momY ← +´ { (𝕩 ⊑ cy) × (𝕩 ⊑ f) }¨ ↕9
  
  # Velocity: u = momentum / rho
  invRho ← 1.0 ÷ (rho + eps)
  ux ← momX × invRho
  uy ← momY × invRho
  
  ⟨rho, ux, uy⟩
}

# Equilibrium Distribution Function f_i^eq(rho, ux, uy)
Equilibrium ⇐ {
  # 𝕩: ⟨rho, ux, uy⟩
  ⟨rho, ux, uy⟩ ← 𝕩
  uSq ← (ux × ux) + (uy × uy)
  
  feq ← {
    i ← 𝕩
    wi ← i ⊑ w
    cxi ← i ⊑ cx
    cyi ← i ⊑ cy
    ciDotU ← (cxi × ux) + (cyi × uy)
    wi × rho × (1.0 + (3.0 × ciDotU) + (4.5 × ciDotU × ciDotU) - (1.5 × uSq))
  }¨ ↕9
  
  feq
}

# BGK Single Relaxation Time Collision Operator
CollideBGK ⇐ {
  # 𝕨: relaxation time tau (or effective local tau matrix)
  # 𝕩: list of 9 population matrices f
  tau ← 𝕨
  f ← 𝕩
  mom ← Moments f
  feq ← Equilibrium mom
  invTau ← 1.0 ÷ tau
  
  # f* = f - (f - feq) / tau
  fPost ← {
    i ← 𝕩
    fi ← i ⊑ f
    feqi ← i ⊑ feq
    fi - (invTau × (fi - feqi))
  }¨ ↕9
  
  fPost
}

# Streaming Operator: f_i(x + c_i, t + dt) = f_i*(x, t)
Stream ⇐ {
  # 𝕩: list of 9 post-collision population matrices fPost
  fPost ← 𝕩
  fStreamed ← {
    i ← 𝕩
    cxi ← i ⊑ cx
    cyi ← i ⊑ cy
    cyi‿cxi Roll2D (i ⊑ fPost)
  }¨ ↕9
  
  fStreamed
}

# Full D2Q9 Time Step: Collision followed by Streaming
StepBGK ⇐ {
  # 𝕨: relaxation time tau
  # 𝕩: list of 9 population matrices f
  tau ← 𝕨
  f ← 𝕩
  fPost ← tau CollideBGK f
  Stream fPost
}

# Guo External Body Force Scheme
GuoForce ⇐ {
  # 𝕨: ⟨fx, fy, tau⟩
  # 𝕩: ⟨rho, ux, uy⟩
  ⟨fx, fy, tau⟩ ← 𝕨
  ⟨rho, ux, uy⟩ ← 𝕩
  invTau ← 1.0 ÷ tau
  factor ← 1.0 - (0.5 × invTau)
  
  sourceTerms ← {
    i ← 𝕩
    wi ← i ⊑ w
    cxi ← i ⊑ cx
    cyi ← i ⊑ cy
    
    ciDotU ← (cxi × ux) + (cyi × uy)
    
    # Bracketed term: (c_i - u)/cs^2 + (c_i . u)/cs^4 * c_i
    termX ← (3.0 × (cxi - ux)) + (9.0 × ciDotU × cxi)
    termY ← (3.0 × (cyi - uy)) + (9.0 × ciDotU × cyi)
    
    dotF ← (termX × fx) + (termY × fy)
    wi × factor × dotF
  }¨ ↕9
  
  sourceTerms
}

# Non-Equilibrium Viscous Stress Tensor Pi^(1)
# Pi_alpha_beta = sum_i c_i_alpha * c_i_beta * (f_i - f_i^eq)
StressTensor ⇐ {
  # 𝕩: list of 9 population matrices f
  f ← 𝕩
  mom ← Moments f
  feq ← Equilibrium mom
  
  # Compute Pi_xx, Pi_yy, Pi_xy
  piXX ← +´ {
    i ← 𝕩
    (i ⊑ cx) × (i ⊑ cx) × ((i ⊑ f) - (i ⊑ feq))
  }¨ ↕9
  
  piYY ← +´ {
    i ← 𝕩
    (i ⊑ cy) × (i ⊑ cy) × ((i ⊑ f) - (i ⊑ feq))
  }¨ ↕9
  
  piXY ← +´ {
    i ← 𝕩
    (i ⊑ cx) × (i ⊑ cy) × ((i ⊑ f) - (i ⊑ feq))
  }¨ ↕9
  
  ⟨piXX, piYY, piXY⟩
}
