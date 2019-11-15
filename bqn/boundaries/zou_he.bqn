# vane/bqn/boundaries/zou_he.bqn
# Zou-He analytical non-equilibrium velocity inlet and pressure outlet boundary conditions.
# Note: Zero em dashes are used in this file per coding standards.

# Zou-He Velocity Inlet at West Boundary (x = 0)
ApplyZouHeInletWest ⇐ {
  # 𝕨: ⟨u0x, u0y⟩
  # 𝕩: populations f (list of 9 matrices)
  ⟨u0x, u0y⟩ ← 𝕨
  fIn ← 𝕩
  
  f0 ← 0 ⊑ fIn
  f1 ← 1 ⊑ fIn
  f2 ← 2 ⊑ fIn
  f3 ← 3 ⊑ fIn
  f4 ← 4 ⊑ fIn
  f5 ← 5 ⊑ fIn
  f6 ← 6 ⊑ fIn
  f7 ← 7 ⊑ fIn
  f8 ← 8 ⊑ fIn
  
  f0Col ← 0 ⊏˘ f0
  f2Col ← 0 ⊏˘ f2
  f3Col ← 0 ⊏˘ f3
  f4Col ← 0 ⊏˘ f4
  f6Col ← 0 ⊏˘ f6
  f7Col ← 0 ⊏˘ f7
  
  # rho = (f0 + f2 + f4 + 2*(f3 + f6 + f7)) / (1 - u0x)
  rhoIn ← (f0Col + f2Col + f4Col + (2.0 × (f3Col + f6Col + f7Col))) ÷ (1.0 - u0x)
  
  diff24 ← 0.5 × (f2Col - f4Col)
  fluxX ← (1.0 ÷ 6.0) × (rhoIn × u0x)
  fluxY ← 0.5 × (rhoIn × u0y)
  
  # Explicitly parenthesized left-to-right terms
  f1Col ← f3Col + ((2.0 ÷ 3.0) × (rhoIn × u0x))
  f5Col ← ((f7Col - diff24) + fluxX) + fluxY
  f8Col ← ((f6Col + diff24) + fluxX) - fluxY
  
  f1New ← f1 ↩ f1Col ⌾(0⊸⊏˘) f1
  f5New ← f5 ↩ f5Col ⌾(0⊸⊏˘) f5
  f8New ← f8 ↩ f8Col ⌾(0⊸⊏˘) f8
  
  ⟨f0, f1New, f2, f3, f4, f5New, f6, f7, f8New⟩
}

# Zou-He Pressure Outlet at East Boundary (x = nx - 1)
ApplyZouHeOutletEast ⇐ {
  # 𝕨: rhoOut (scalar)
  # 𝕩: populations f (list of 9 matrices)
  rhoOut ← 𝕨
  fIn ← 𝕩
  
  ⟨ny, nx⟩ ← ≢(0 ⊑ fIn)
  lastCol ← nx - 1
  
  f0 ← 0 ⊑ fIn
  f1 ← 1 ⊑ fIn
  f2 ← 2 ⊑ fIn
  f3 ← 3 ⊑ fIn
  f4 ← 4 ⊑ fIn
  f5 ← 5 ⊑ fIn
  f6 ← 6 ⊑ fIn
  f7 ← 7 ⊑ fIn
  f8 ← 8 ⊑ fIn
  
  f0Col ← lastCol ⊏˘ f0
  f1Col ← lastCol ⊏˘ f1
  f2Col ← lastCol ⊏˘ f2
  f4Col ← lastCol ⊏˘ f4
  f5Col ← lastCol ⊏˘ f5
  f8Col ← lastCol ⊏˘ f8
  
  uxOut ← ¯1.0 + ((f0Col + f2Col + f4Col + (2.0 × (f1Col + f5Col + f8Col))) ÷ rhoOut)
  uyOut ← 0.0 × uxOut
  
  diff24 ← 0.5 × (f2Col - f4Col)
  fluxX ← (1.0 ÷ 6.0) × (rhoOut × uxOut)
  fluxY ← 0.5 × (rhoOut × uyOut)
  
  f3Col ← f1Col - ((2.0 ÷ 3.0) × (rhoOut × uxOut))
  f6Col ← ((f8Col - diff24) - fluxX) + fluxY
  f7Col ← ((f5Col + diff24) - fluxX) - fluxY
  
  f3New ← f3 ↩ f3Col ⌾(lastCol⊸⊏˘) f3
  f6New ← f6 ↩ f6Col ⌾(lastCol⊸⊏˘) f6
  f7New ← f7 ↩ f7Col ⌾(lastCol⊸⊏˘) f7
  
  ⟨f0, f1, f2, f3New, f4, f5, f6New, f7New, f8⟩
}
