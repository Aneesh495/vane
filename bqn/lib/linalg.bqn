# vane/bqn/lib/linalg.bqn
# High-performance linear algebra and tensor contraction routines for aerodynamic workflows.
# Note: Zero em dashes are used in this file per coding standards.

eps ← 1e¯12

# Matrix Multiplication
# A (w) is M x K matrix, B (x) is K x N matrix
MatMul ⇐ {
  (<˘ 𝕨) {+´𝕨×𝕩}⌜ (<˘ ⍉𝕩)
}

# Matrix-Vector Multiplication
# A (w) is M x N matrix, v (x) is N-element vector
MatVec ⇐ {
  (<˘ 𝕨) {+´𝕨×𝕩}¨ <𝕩
}

# Matrix Transpose
Transpose ⇐ {
  ⍉𝕩
}

# Identity Matrix of dimension n
Eye ⇐ {
  (↕𝕩) =⌜ (↕𝕩)
}

# Trace of Square Matrix
Trace ⇐ {
  +´ ⥊ ((↕≠𝕩) =⌜ (↕≠𝕩)) × 𝕩
}

# Frobenius Norm of Matrix or Tensor
NormFrobenius ⇐ {
  √+´ ⥊ (𝕩 × 𝕩)
}

# Outer Product: u (w) outer v (x) -> u_i * v_j
Outer ⇐ {
  𝕨 ×⌜ 𝕩
}

# 2x2 Determinant
Det2x2 ⇐ {
  ((0‿0 ⊑ 𝕩) × (1‿1 ⊑ 𝕩)) - ((0‿1 ⊑ 𝕩) × (1‿0 ⊑ 𝕩))
}

# 2x2 Matrix Inverse
Inv2x2 ⇐ {
  d ← Det2x2 𝕩
  (1e¯12 ≥ |d) ◶ ⟨{
    (2‿2 ⥊ ⟨1‿1 ⊑ 𝕩, -0‿1 ⊑ 𝕩, -1‿0 ⊑ 𝕩, 0‿0 ⊑ 𝕩⟩) ÷ d
  }, {0×𝕩}⟩ 𝕩
}

# Solve 2x2 Linear System: A x = b
Solve2x2 ⇐ {
  (Inv2x2 𝕨) MatVec 𝕩
}

# 3x3 Determinant via Triple Scalar Product of Rows
Det3x3 ⇐ {
  ⟨a, b, c⟩ ← 0 ⊏ 𝕩
  ⟨d, e, f⟩ ← 1 ⊏ 𝕩
  ⟨g, h, i⟩ ← 2 ⊏ 𝕩
  (a × ((e × i) - (f × h))) - ((b × ((d × i) - (f × g))) - (c × ((d × h) - (e × g))))
}

# 3x3 Matrix Inverse via Cofactor Matrix
Inv3x3 ⇐ {
  detVal ← Det3x3 𝕩
  (1e¯12 ≥ |detVal) ◶ ⟨{
    ⟨a, b, c⟩ ← 0 ⊏ 𝕩
    ⟨d, e, f⟩ ← 1 ⊏ 𝕩
    ⟨g, h, i⟩ ← 2 ⊏ 𝕩
    
    c00 ← (e × i) - (f × h)
    c01 ← (f × g) - (d × i)
    c02 ← (d × h) - (e × g)
    
    c10 ← (c × h) - (b × i)
    c11 ← (a × i) - (c × g)
    c12 ← (b × g) - (a × h)
    
    c20 ← (b × f) - (c × e)
    c21 ← (c × d) - (a × f)
    c22 ← (a × e) - (b × d)
    
    det ← (a × c00) + ((b × c01) + (c × c02))
    adj ← 3‿3 ⥊ ⟨c00, c10, c20, c01, c11, c21, c02, c12, c22⟩
    adj ÷ det
  }, {0×𝕩}⟩ 𝕩
}

# Solve 3x3 Linear System: A x = b
Solve3x3 ⇐ {
  (Inv3x3 𝕨) MatVec 𝕩
}

# Eigenvalues of Symmetric 2x2 Tensor
Eigen2x2Sym ⇐ {
  a ← 0‿0 ⊑ 𝕩
  b ← 0‿1 ⊑ 𝕩
  c ← 1‿1 ⊑ 𝕩
  tr ← a + c
  gap ← √ (((a - c) × (a - c)) + (4.0 × (b × b)))
  ⟨(tr + gap) ÷ 2.0, (tr - gap) ÷ 2.0⟩
}

# Eigenvalues of Symmetric 3x3 Tensor via Analytical Cardano Method
Eigen3x3Sym ⇐ {
  ⟨a, b, c⟩ ← 0 ⊏ 𝕩
  ⟨d0, d, e⟩ ← 1 ⊏ 𝕩
  ⟨f0, f1, f⟩ ← 2 ⊏ 𝕩
  
  p1 ← (b×b) + ((c×c) + (e×e))
  
  (p1 ≤ eps) ◶ ⟨{
    q1 ← a + (d + f)
    q ← q1 ÷ 3.0
    p2 ← (((a - q) × (a - q)) + (((d - q) × (d - q)) + ((f - q) × (f - q)))) + (2.0 × p1)
    p ← √(p2 ÷ 6.0)
    
    invP ← 1.0 ÷ (p ⌈ eps)
    b00 ← (a - q) × invP
    b11 ← (d - q) × invP
    b22 ← (f - q) × invP
    b01 ← b × invP
    b02 ← c × invP
    b12 ← e × invP
    
    detB ← (b00 × ((b11 × b22) - (b12 × b12))) - ((b01 × ((b01 × b22) - (b12 × b02))) - (b02 × ((b01 × b12) - (b11 × b02))))
    r ← ¯1.0 ⌈ (1.0 ⌊ (detB ÷ 2.0))
    
    phi ← (•math.Acos r) ÷ 3.0
    piThird ← 1.0471975511965976
    
    eig0 ← q + (2.0 × (p × •math.Cos phi))
    eig2 ← q + (2.0 × (p × •math.Cos (phi + (2.0 × piThird))))
    eig1 ← (3.0 × q) - (eig0 + eig2)
    
    ⟨eig0, eig1, eig2⟩
  }, {
    diag ← ⟨a, d, f⟩
    (⍒ diag) ⊏ diag
  }⟩ 𝕩
}
