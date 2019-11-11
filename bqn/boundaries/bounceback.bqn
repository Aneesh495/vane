# vane/bqn/boundaries/bounceback.bqn
# Half-way and full-way bounce-back boundary conditions for complex aerodynamic bodies.
# Note: Zero em dashes are used in this file per coding standards.

d2q9 ← •Import "../lbm/d2q9.bqn"
opp ← d2q9.opp
w ← d2q9.w
cx ← d2q9.cx
cy ← d2q9.cy

# Standard Half-Way Bounce-Back for Static Obstacles
# At solid nodes, populations reflect back along opposite lattice directions:
# f_opp(x_f, t + dt) = f_i*(x_f, t)
ApplyStaticBounceBack ⇐ {
  # 𝕨: boolean obstacle mask (1 for solid, 0 for fluid)
  # 𝕩: populations fPost (post-collision or post-streaming)
  mask ← 𝕨
  fIn ← 𝕩
  
  # Reflect each direction into its opposite
  fReflected ← {
    i ← 𝕩
    iOpp ← i ⊑ opp
    fi ← i ⊑ fIn
    fiOpp ← iOpp ⊑ fIn
    
    # Where mask is 1 (solid), swap with opposite population
    # Where mask is 0 (fluid), keep current fi
    (mask × fiOpp) + ((1.0 - mask) × fi)
  }¨ ↕9
  
  fReflected
}

# Moving Wall Bounce-Back (e.g. flapping wings, rotor blades)
# Incorporates momentum exchange from moving solid boundary:
# f_opp(x_f, t + dt) = f_i*(x_f, t) - 2 * w_i * rho * (c_i . u_w) / cs^2
ApplyMovingBounceBack ⇐ {
  # 𝕨: ⟨mask, rho, uwX, uwY⟩
  # 𝕩: populations fPost
  ⟨mask, rho, uwX, uwY⟩ ← 𝕨
  fIn ← 𝕩
  
  fOut ← {
    i ← 𝕩
    iOpp ← i ⊑ opp
    wi ← i ⊑ w
    cxi ← i ⊑ cx
    cyi ← i ⊑ cy
    
    fi ← i ⊑ fIn
    fiOpp ← iOpp ⊑ fIn
    
    # Momentum transfer term: 2 * w_i * rho * 3 * (c_i . u_w)
    ciDotUw ← (cxi × uwX) + (cyi × uwY)
    deltaMom ← 6.0 × wi × rho × ciDotUw
    
    fSolid ← fiOpp - deltaMom
    (mask × fSolid) + ((1.0 - mask) × fi)
  }¨ ↕9
  
  fOut
}

# Step with Static Obstacle Bounce-Back:
# Collide -> Stream -> Apply Bounce-Back at solid boundaries
StepWithObstacle ⇐ {
  # 𝕨: ⟨tau, mask⟩
  # 𝕩: populations f
  ⟨tau, mask⟩ ← 𝕨
  f ← 𝕩
  
  # Collision
  fPost ← tau d2q9.CollideBGK f
  
  # Stream
  fStreamed ← d2q9.Stream fPost
  
  # Half-way bounce-back at solid nodes
  mask ApplyStaticBounceBack fStreamed
}
