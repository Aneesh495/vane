# vane - Lattice Boltzmann Collision and Lattice Physics Test Suite
# Zero em dashes used in comments or strings.

la   ← •Import "../../bqn/lib/linalg.bqn"
d2q9 ← •Import "../../bqn/lbm/d2q9.bqn"
mrt  ← •Import "../../bqn/lbm/mrt.bqn"
trt  ← •Import "../../bqn/lbm/trt.bqn"
d3q19 ← •Import "../../bqn/lbm/d3q19.bqn"
elbm ← •Import "../../bqn/lbm/entropic.bqn"

# Assertion helpers
Assert ← { 𝕨 ! 𝕩 }
AssertNear ← {
  diff ← | 𝕩 - 𝕨
  "Assertion near-equality failed" ! diff < 1e¯4
}

TestD2Q9Lattice ← {
  dummy ← 𝕩
  # 1. Weights sum exactly to 1.0
  sumW ← +´ d2q9.w
  1.0 AssertNear sumW

  # 2. Sound speed squared is 1/3
  (1.0 ÷ 3.0) AssertNear d2q9.cs2

  # 3. Test Equilibrium Moment Recovery
  # At rho = 1.05, ux = 0.08, uy = ¯0.02
  rhoIn ← 1.05
  uxIn  ← 0.08
  uyIn  ← ¯0.02
  # Grid of 2x2 cells
  rhoGrid ← 2‿2 ⥊ rhoIn
  uxGrid  ← 2‿2 ⥊ uxIn
  uyGrid  ← 2‿2 ⥊ uyIn
  feq ← d2q9.Equilibrium ⟨rhoGrid, uxGrid, uyGrid⟩

  # Recover moments from feq
  ⟨rhoRec, uxRec, uyRec⟩ ← d2q9.Moments feq
  rhoIn AssertNear (0‿0 ⊑ rhoRec)
  uxIn  AssertNear (0‿0 ⊑ uxRec)
  uyIn  AssertNear (0‿0 ⊑ uyRec)

  •Out "  PASS: D2Q9 lattice weights and moments validated."
  1
}

TestMRTTransform ← {
  dummy ← 𝕩
  # Transformation matrix M should satisfy M * M^-1 = Identity
  prod ← mrt.mMat la.MatMul mrt.mInv
  ident ← la.Eye 9
  maxDiff ← ⌈´ ⥊ | prod - ident
  "MRT Orthogonal Invertibility" Assert (maxDiff < 1e¯6)

  # Check equilibrium moments at rest (rho=1.0, u=0):
  # m0 = 1.0, m1 = -2.0, m2 = 1.0, all others 0.0
  mEqRest ← mrt.EquilibriumMoments ⟨1.0, 0.0, 0.0⟩
  1.0  AssertNear (0 ⊑ mEqRest)
  ¯2.0 AssertNear (1 ⊑ mEqRest)
  1.0  AssertNear (2 ⊑ mEqRest)
  0.0  AssertNear (3 ⊑ mEqRest)
  0.0  AssertNear (5 ⊑ mEqRest)

  •Out "  PASS: MRT orthogonal transformation validated."
  1
}

TestTRTOperator ← {
  dummy ← 𝕩
  # Magic parameter Lambda = 1/4 with tauPlus = 0.8:
  # (tauPlus - 0.5) * (tauMinus - 0.5) = 0.25 -> (0.3) * (tauMinus - 0.5) = 0.25
  # tauMinus - 0.5 = 0.25 / 0.3 = 5/6 = 0.833333 -> tauMinus = 1.333333
  tauMinus ← 0.25 trt.TauMinusFromLambda 0.8
  (0.5 + (0.25 ÷ 0.3)) AssertNear tauMinus

  # Test Symmetric / Anti-symmetric decomposition:
  # f = fSym + fAnti
  fDummy ← { 2‿2 ⥊ (𝕩 + 1.0) × 0.1 }¨ ↕9
  ⟨fSym, fAnti⟩ ← trt.DecomposeSymAnti fDummy
  fRecon ← fSym + fAnti
  maxDiff ← ⌈´ { ⌈´ ⥊ | (𝕩 ⊑ fDummy) - (𝕩 ⊑ fRecon) }¨ ↕9
  "TRT Decomposition Identity" Assert (maxDiff < 1e¯10)

  •Out "  PASS: TRT symmetry decomposition validated."
  1
}

TestD3Q19Lattice ← {
  dummy ← 𝕩
  # 19 velocity directions
  "19 directions" Assert (19 = ≠ d3q19.w)
  # Weights sum to 1.0
  1.0 AssertNear (+´ d3q19.w)

  •Out "  PASS: D3Q19 3D lattice parameters validated."
  1
}

TestEntropicELBM ← {
  dummy ← 𝕩
  # Calculate discrete H for uniform equilibrium at rest
  feqRest ← { 2‿2 ⥊ (𝕩 ⊑ d2q9.w) }¨ ↕9
  hVal ← elbm.EntropyH feqRest
  # For f = w, f/w = 1, ln(f/w) = 0, so H = sum(w * ln(1)) = 0.0
  0.0 AssertNear (0‿0 ⊑ hVal)

  •Out "  PASS: Entropic ELBM H-theorem properties validated."
  1
}

# Run all test cases in this file
•Out "Running Lattice Boltzmann Physics Tests:"
TestD2Q9Lattice ⟨⟩
TestMRTTransform ⟨⟩
TestTRTOperator ⟨⟩
TestD3Q19Lattice ⟨⟩
TestEntropicELBM ⟨⟩
•Out "All LBM tests passed successfully."
