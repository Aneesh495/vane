# vane - Aerodynamics, Airfoil Geometry, and Boundary Layer Test Suite
# Zero em dashes used in comments or strings.

naca ← •Import "../../bqn/aerodynamics/airfoil_naca.bqn"
bl   ← •Import "../../bqn/aerodynamics/boundary_layer.bqn"
me   ← •Import "../../bqn/aerodynamics/momentum_exchange.bqn"
d2q9 ← •Import "../../bqn/lbm/d2q9.bqn"

# Assertion helpers
Assert ← { 𝕨 ! 𝕩 }
AssertNear ← {
  diff ← | 𝕩 - 𝕨
  "Assertion near-equality failed" ! diff < 0.05
}

TestNACA4Airfoil ← {
  dummy ← 𝕩
  # Generate NACA 0012: chord=100.0, aoa=0.0, 50 points
  coords ← naca.GenerateNACA4 ⟨0.0, 0.0, 0.12, 100.0, 0.0, 51⟩
  # Coordinates format: <xUpper, yUpper, xLower, yLower, camberY, xc>
  ⟨xuRot, yuRot, xlRot, ylRot⟩ ← coords

  # Symmetric airfoil upper and lower thickness should be identical and opposite
  symDiff ← ⌈´ | (yuRot + ylRot)
  "NACA 0012 upper and lower surface symmetry" Assert (symDiff < 1e¯6)

  # Check total thickness near 30% chord (max thickness for NACA 0012 is 12 lu)
  halfThickStation ← 15 ⊑ yuRot
  "NACA 0012 half thickness check" Assert (halfThickStation > 0.04 × 100.0)

  •Out "  PASS: NACA 4-digit airfoil geometry generation validated."
  1
}

TestBoundaryLayerIntegral ← {
  dummy ← 𝕩
  # Create linear boundary layer velocity profile up to delta = 10:
  # y = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
  # u = [0, 0.1, 0.2, ..., 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
  # For linear profile:
  # deltaStar = 10 / 2 = 5.0
  # theta = 10 / 6 = 1.6667
  # Shape factor H = 5.0 / 1.6667 = 3.0
  yNorm ← (↕16) × 1.0
  uTang ← 1.0 ⌊ ((↕16) × 0.1)

  params ← bl.IntegralParameters ⟨yNorm, uTang, 0.01, 1.0, 0.1⟩
  ⟨delta99, deltaStar, theta, shapeFactor, tauWall, cf⟩ ← params

  10.0 AssertNear delta99
  5.0  AssertNear deltaStar
  (10.0 ÷ 6.0) AssertNear theta
  3.0  AssertNear shapeFactor

  •Out "  PASS: Boundary layer integral thickness analytics validated."
  1
}

TestMomentumExchangeSymmetry ← {
  dummy ← 𝕩
  # Create a 20x20 domain with uniform horizontal flow at rest,
  # and a symmetric square obstacle centered at (10, 10)
  ny ← 20
  nx ← 20
  # Set 4x4 square obstacle
  solidMask ← > {
    iy ← 𝕩
    {
      ix ← 𝕩
      ((iy ≥ 8) ∧ (iy ≤ 11) ∧ (ix ≥ 8) ∧ (ix ≤ 11)) ⊑ ⟨0, 1⟩
    }¨ ↕nx
  }¨ ↕ny

  # Equilibrium distribution at u_x = 0.05, u_y = 0.0
  rhoGrid ← ny‿nx ⥊ 1.0
  uxGrid  ← ny‿nx ⥊ 0.05
  uyGrid  ← ny‿nx ⥊ 0.0
  feq     ← d2q9.Equilibrium ⟨rhoGrid, uxGrid, uyGrid⟩

  ⟨fx, fy, moment, cd, cl, cm⟩ ← me.ComputeForces ⟨feq, solidMask, 0.0‿0.0, 1.0, 0.05, 4.0, 10.0‿10.0⟩

  # Under symmetric horizontal flow, lift Cl should be virtually 0
  0.0 AssertNear cl

  •Out "  PASS: Momentum exchange boundary force integration validated."
  1
}

# Run all test cases in this file
•Out "Running Aerodynamics & Boundary Layer Tests:"
TestNACA4Airfoil ⟨⟩
TestBoundaryLayerIntegral ⟨⟩
TestMomentumExchangeSymmetry ⟨⟩
•Out "All Aerodynamics tests passed successfully."
