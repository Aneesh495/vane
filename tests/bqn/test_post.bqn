# vane - Post-Processing and Vortex Core Diagnostics Test Suite
# Tests: Q-criterion, Jeong-Hussain Lambda-2, Streamline Tracer.
# Zero em dashes used in comments or strings.

qc   ← •Import "../../bqn/post/q_criterion.bqn"
l2   ← •Import "../../bqn/post/lambda2.bqn"
str  ← •Import "../../bqn/post/streamlines.bqn"

# Assertion helpers
Assert ← { 𝕨 ! 𝕩 }
AssertNear ← {
  diff ← | 𝕩 - 𝕨
  "Assertion near-equality failed" ! diff < 1e¯3
}

TestQCriterionRigidRotation ← {
  dummy ← 𝕩
  # Construct pure rigid-body vortex on 15x15 grid:
  # u_x = -omega * (y - y0)
  # u_y = +omega * (x - x0)
  # with omega = 0.1, center at (7, 7)
  ny ← 15
  nx ← 15
  omega ← 0.1
  ux ← > { iy ← 𝕩 ⋄ { ix ← 𝕩 ⋄ (0.0 - omega) × (iy - 7.0) }¨ ↕nx }¨ ↕ny
  uy ← > { iy ← 𝕩 ⋄ { ix ← 𝕩 ⋄ omega × (ix - 7.0) }¨ ↕nx }¨ ↕ny

  qField ← qc.QCriterion2D ⟨ux, uy⟩

  # In pure rigid rotation, S = 0 and ||Omega||^2 > 0 -> Q > 0 everywhere in interior
  # Check center node
  qCenter ← 7‿7 ⊑ qField
  "Rigid vortex positive Q-criterion" Assert (qCenter > 0.0)

  •Out "  PASS: Q-criterion vortex core detection validated."
  1
}

TestLambda2RigidRotation ← {
  dummy ← 𝕩
  ny ← 15
  nx ← 15
  omega ← 0.1
  ux ← > { iy ← 𝕩 ⋄ { ix ← 𝕩 ⋄ (0.0 - omega) × (iy - 7.0) }¨ ↕nx }¨ ↕ny
  uy ← > { iy ← 𝕩 ⋄ { ix ← 𝕩 ⋄ omega × (ix - 7.0) }¨ ↕nx }¨ ↕ny

  lam2Field ← l2.Lambda2_2D ⟨ux, uy⟩

  # Vortex core condition: lambda_2 < 0
  lam2Center ← 7‿7 ⊑ lam2Field
  "Rigid vortex negative Lambda-2" Assert (lam2Center < 0.0)

  •Out "  PASS: Jeong-Hussain Lambda-2 vortex detection validated."
  1
}

TestBilinearInterpolation ← {
  dummy ← 𝕩
  # Create simple 10x10 field f(x, y) = 10*y + x
  ny ← 10
  nx ← 10
  field ← > { iy ← 𝕩 ⋄ { ix ← 𝕩 ⋄ (10.0 × iy) + (0.0 + ix) }¨ ↕nx }¨ ↕ny

  # Exact node interpolation at x=3.0, y=4.0 -> 43.0
  valNode ← str.BilinearInterp ⟨field, 3.0‿4.0⟩
  43.0 AssertNear valNode

  # Midway interpolation at x=3.5, y=4.5 -> 10*4.5 + 3.5 = 48.5
  valMid ← str.BilinearInterp ⟨field, 3.5‿4.5⟩
  48.5 AssertNear valMid

  •Out "  PASS: Subgrid bilinear field interpolation validated."
  1
}

# Run all test cases in this file
•Out "Running Post-Processing & Vortex Diagnostics Tests:"
TestQCriterionRigidRotation ⟨⟩
TestLambda2RigidRotation ⟨⟩
TestBilinearInterpolation ⟨⟩
•Out "All Post-Processing tests passed successfully."
