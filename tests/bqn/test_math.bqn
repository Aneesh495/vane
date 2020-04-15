# vane - Mathematical Utilities and Linear Algebra Test Suite
# Zero em dashes used in comments or strings.

math ← •Import "../../bqn/lib/math_utils.bqn"
la   ← •Import "../../bqn/lib/linalg.bqn"
st   ← •Import "../../bqn/lib/stencils.bqn"

# Assertion helpers
Assert ← { 𝕨 ! 𝕩 }
AssertNear ← {
  diff ← | 𝕩 - 𝕨
  "Assertion near-equality failed" ! diff < 1e¯4
}

TestMathUtils ← {
  dummy ← 𝕩
  # Test NormL2 on 3-4 vector
  5.0 AssertNear (math.NormL2 3.0‿4.0)

  # Test NormL1
  7.0 AssertNear (math.NormL1 3.0‿¯4.0)

  # Test Dot Product
  11.0 AssertNear (1.0‿2.0 math.Dot 3.0‿4.0)

  # Test Cross2D: [1, 0] x [0, 1] = 1.0
  1.0 AssertNear (1.0‿0.0 math.Cross2D 0.0‿1.0)

  # Test 1D Linear Interpolation
  15.0 AssertNear (0.0‿1.0‿10.0‿20.0 math.Lerp1D 0.5)

  # Test Error Function Erf(0) = 0
  0.0 AssertNear (math.Erf 0.0)

  # Test Simpson 1D integration of f(x) = 1 over interval [0, 2] -> 2.0
  2.0 AssertNear (0.5 math.Simpson1D 5⥊1.0)

  •Out "  PASS: Mathematical utilities validated."
  1
}

TestLinAlg ← {
  dummy ← 𝕩
  # Test Matrix Determinant
  m ← 2‿2 ⥊ 4.0‿7.0‿2.0‿6.0
  10.0 AssertNear (la.Det2x2 m)

  # Test Matrix Inverse
  invM ← la.Inv2x2 m
  prod ← m la.MatMul invM
  ident ← la.Eye 2
  maxDiff ← ⌈´ ⥊ | prod - ident
  "Matrix Inverse Identity" Assert (maxDiff < 1e¯6)

  # Test Trace: 4 + 6 = 10
  10.0 AssertNear (la.Trace m)

  # Test Frobenius Norm: sqrt(16 + 49 + 4 + 36) = sqrt(105)
  (√105.0) AssertNear (la.NormFrobenius m)

  # Test 2x2 Eigenvalues of symmetric matrix [[2, 1], [1, 2]] -> eigenvalues 3 and 1
  symM ← 2‿2 ⥊ 2.0‿1.0‿1.0‿2.0
  eigs ← la.Eigen2x2Sym symM
  3.0 AssertNear (0 ⊑ eigs)
  1.0 AssertNear (1 ⊑ eigs)

  •Out "  PASS: Linear algebra routines validated."
  1
}

TestStencils ← {
  dummy ← 𝕩
  # Compute gradient of field
  field ← (10‿10 ⥊ (↕100)) × 0.1
  ⟨dfdx, dfdy⟩ ← 1.0‿1.0 st.Grad2DPeriodic field

  # Check shapes
  "Grad2D dfdx shape" Assert ((≢ dfdx) ≡ 10‿10)
  "Grad2D dfdy shape" Assert ((≢ dfdy) ≡ 10‿10)

  •Out "  PASS: Finite difference stencils validated."
  1
}

# Run all test cases in this file
•Out "Running Math & Stencil Tests:"
TestMathUtils ⟨⟩
TestLinAlg ⟨⟩
TestStencils ⟨⟩
•Out "All Math tests passed successfully."
