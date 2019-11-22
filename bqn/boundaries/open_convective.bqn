# vane/bqn/boundaries/open_convective.bqn
# Convective non-reflecting outflow boundary conditions for wake vortex shedding.
# Note: Zero em dashes are used in this file per coding standards.

# Apply Non-Reflecting Convective Boundary at East Outlet (x = nx - 1)
# Equation: df/dt + U_conv * df/dx = 0
ApplyConvectiveOutletEast ⇐ {
  # 𝕨: ⟨uConv, fPrev⟩
  # 𝕩: populations f (current time step)
  ⟨uConv, fPrev⟩ ← 𝕨
  fCurr ← 𝕩
  
  ⟨ny, nx⟩ ← ≢(0 ⊑ fCurr)
  lastCol ← nx - 1
  prevCol ← nx - 2
  
  cFactor ← uConv ⌈ 0.01
  denom ← 1.0 + cFactor
  
  fOut ← {
    i ← 𝕩
    fiCurr ← i ⊑ fCurr
    fiPrev ← i ⊑ fPrev
    
    colPrevT ← lastCol ⊏˘ fiPrev
    colUpstream ← prevCol ⊏˘ fiCurr
    
    colNewLast ← (colPrevT + (cFactor × colUpstream)) ÷ denom
    
    fiCurr ↩ colNewLast ⌾(lastCol⊸⊏˘) fiCurr
  }¨ ↕9
  
  fOut
}

# Free-Slip Boundary Condition at North (y = ny - 1) and South (y = 0) Walls
# Normal velocity is zero, tangential shear is zero (mirror reflection)
ApplyFreeSlipNorthSouth ⇐ {
  # 𝕩: populations f (list of 9 matrices)
  fIn ← 𝕩
  ⟨ny, nx⟩ ← ≢(0 ⊑ fIn)
  lastRow ← ny - 1
  
  f0 ← 0 ⊑ fIn
  f1 ← 1 ⊑ fIn
  f2 ← 2 ⊑ fIn
  f3 ← 3 ⊑ fIn
  f4 ← 4 ⊑ fIn
  f5 ← 5 ⊑ fIn
  f6 ← 6 ⊑ fIn
  f7 ← 7 ⊑ fIn
  f8 ← 8 ⊑ fIn
  
  # South wall (y = 0): incoming populations are f2, f5, f6 reflecting from f4, f8, f7
  # North wall (y = lastRow): incoming populations are f4, f7, f8 reflecting from f2, f6, f5
  f2Row0 ← 0 ⊏ f4
  f5Row0 ← 0 ⊏ f8
  f6Row0 ← 0 ⊏ f7
  
  f4RowLast ← lastRow ⊏ f2
  f7RowLast ← lastRow ⊏ f6
  f8RowLast ← lastRow ⊏ f5
  
  f2New ← f2 ↩ f2Row0 ⌾(0⊸⊏) f2
  f5New ← f5 ↩ f5Row0 ⌾(0⊸⊏) f5
  f6New ← f6 ↩ f6Row0 ⌾(0⊸⊏) f6
  
  f4New ← f4 ↩ f4RowLast ⌾(lastRow⊸⊏) f4
  f7New ← f7 ↩ f7RowLast ⌾(lastRow⊸⊏) f7
  f8New ← f8 ↩ f8RowLast ⌾(lastRow⊸⊏) f8
  
  ⟨f0, f1, f2New, f3, f4New, f5New, f6New, f7New, f8New⟩
}
