# vane/bqn/lib/math_utils.bqn
# High-precision mathematical and numerical utilities for aerodynamics and fluid simulation.
# Note: Zero em dashes are used in this file per coding standards.

# Basic Constants
pi ⇐ 3.14159265358979323846
twoPi ⇐ 6.28318530717958647692
halfPi ⇐ 1.57079632679489661923
eps ⇐ 1e¯12

Sin ⇐ •math.Sin
Cos ⇐ •math.Cos

# Vector Norms
NormL2 ⇐ {
  # Computes Euclidean L2 norm of vector x
  √+´𝕩×𝕩
}

NormL1 ⇐ {
  # Computes Manhattan L1 norm: sum(|x_i|)
  +´|𝕩
}

NormLInf ⇐ {
  # Computes Chebyshev L-infinity norm: max(|x_i|)
  ⌈´|𝕩
}

Normalize ⇐ {
  # Normalizes vector to unit length; returns zero vector if magnitude is near zero
  (NormL2 𝕩) ≤ eps ? 0×𝕩 ; 𝕩 ÷ NormL2 𝕩
}

Dot ⇐ {
  # Computes Euclidean dot product between two vectors w and x
  +´𝕨×𝕩
}

Cross3D ⇐ {
  # Computes 3D vector cross product between 3-element vectors w and x
  ⟨w0, w1, w2⟩ ← 𝕨
  ⟨x0, x1, x2⟩ ← 𝕩
  ⟨ (w1×x2) - (w2×x1),
    (w2×x0) - (w0×x2),
    (w0×x1) - (w1×x0) ⟩
}

Cross2D ⇐ {
  # Computes 2D scalar cross product (z-component): w0*x1 - w1*x0
  ⟨w0, w1⟩ ← 𝕨
  ⟨x0, x1⟩ ← 𝕩
  (w0×x1) - (w1×x0)
}

AngleBetween ⇐ {
  # Computes angle in radians between vectors w and x using safe arccos
  ((NormL2 𝕨) × (NormL2 𝕩)) ≤ eps ? 0.0 ; •math.Acos (¯1.0 ⌈ 1.0 ⌊ ((𝕨 Dot 𝕩) ÷ ((NormL2 𝕨) × (NormL2 𝕩))))
}

ProjectVector ⇐ {
  # Projects vector x onto vector w: proj_w(x) = ((x . w) / (w . w)) * w
  (𝕨 Dot 𝕨) ≤ eps ? 0×𝕩 ; ((𝕩 Dot 𝕨) ÷ (𝕨 Dot 𝕨)) × 𝕨
}

# Coordinate Transformations and Rotations
Rot2D ⇐ {
  # Generates 2D rotation matrix for angle theta in radians
  c ← •math.Cos 𝕩
  s ← •math.Sin 𝕩
  2‿2 ⥊ ⟨c, -s, s, c⟩
}

RotateVec2D ⇐ {
  # Rotates 2D vector x by angle theta (w)
  ⟨x0, x1⟩ ← 𝕩
  r ← Rot2D 𝕨
  ⟨ ((0‿0 ⊑ r) × x0) + ((0‿1 ⊑ r) × x1),
    ((1‿0 ⊑ r) × x0) + ((1‿1 ⊑ r) × x1) ⟩
}

Rodrigues3D ⇐ {
  # Applies Rodrigues rotation formula: rotates 3D vector x around unit axis k (w) by angle theta
  ⟨k_vec, theta⟩ ← 𝕨
  k ← Normalize k_vec
  c ← •math.Cos theta
  s ← •math.Sin theta
  vCrossK ← k Cross3D 𝕩
  vDotK ← k Dot 𝕩
  (𝕩 × c) + (vCrossK × s) + (k × (vDotK × (1.0 - c)))
}

# Numerical Integration
Trapezoid1D ⇐ {
  # Integrates array of samples y (x) with uniform spacing dx (w)
  (≠𝕩) < 2 ? 0.0 ; (0.5 × ((0 ⊑ 𝕩) + (((≠𝕩) - 1) ⊑ 𝕩)) + (+´1 ↓ ¯1 ↓ 𝕩)) × 𝕨
}

Simpson1D ⇐ {
  # Integrates array of samples y (x) with uniform spacing dx (w) via Simpson 1/3 rule
  (≠𝕩) < 3 ? (𝕨 Trapezoid1D 𝕩) ; {
    n ← ≠𝕩
    idx ← ↕n
    weights ← {
      i ← 𝕩
      ((i = 0) ∨ (i = n - 1)) ⊑ ⟨((1 = 2 | i) ⊑ ⟨2.0, 4.0⟩), 1.0⟩
    }¨ idx
    (𝕨 ÷ 3.0) × +´weights × 𝕩
  }
}

# Polynomial Evaluation via Horner Scheme
Horner ⇐ {
  # w: coefficients array [c0, c1, ..., cn]
  # x: scalar x value
  coeffsRev ← ⌽𝕨
  {𝕩 × 𝕨 + 𝕩}´ coeffsRev
}

# Special Functions: Error Function (Erf) via Rational Approximation
Erf ⇐ {
  val ← 𝕩
  sign ← (val < 0) ⊑ ⟨1.0, ¯1.0⟩
  ax ← |val
  p ← 0.3275911
  t ← 1.0 ÷ (1.0 + p × ax)
  a1 ← 0.254829592
  a2 ← ¯0.284496736
  a3 ← 1.421413741
  a4 ← ¯1.453152027
  a5 ← 1.061405429
  poly ← t × (a1 + t × (a2 + t × (a3 + t × (a4 + t × a5))))
  res ← 1.0 - (poly × ⋆(¯1.0 × ax × ax))
  sign × res
}

Erfc ⇐ {
  # Complementary error function erfc(x) = 1 - erf(x)
  1.0 - Erf 𝕩
}

# Linear Interpolation (1D)
Lerp1D ⇐ {
  # w: ⟨x0, x1, y0, y1⟩
  # x: target x
  ⟨x0, x1, y0, y1⟩ ← 𝕨
  t ← (𝕩 - x0) ÷ (x1 - x0)
  y0 + t × (y1 - y0)
}

# Bilinear Interpolation (2D)
Bilinear ⇐ {
  # w: ⟨x0, x1, y0, y1, q00, q10, q01, q11⟩
  # x: ⟨x_target, y_target⟩
  ⟨x0, x1, y0, y1, q00, q10, q01, q11⟩ ← 𝕨
  ⟨xt, yt⟩ ← 𝕩
  tx ← (xt - x0) ÷ (x1 - x0)
  ty ← (yt - y0) ÷ (y1 - y0)
  r0 ← q00 + tx × (q10 - q00)
  r1 ← q01 + tx × (q11 - q01)
  r0 + ty × (r1 - r0)
}

# 1D Root Finding: Secant Method
_secantSolve ⇐ {
  # 𝔽: function to find root of
  # 𝕩: ⟨x0, x1, tol, maxIter⟩
  ⟨x0, x1, tol, maxI⟩ ← 𝕩
  f0 ← 𝔽 x0
  f1 ← 𝔽 x1
  _step ← {
    ⟨p0, p1, q0, q1, k⟩ ← 𝕩
    den ← q1 - q0
    pNext ← p1 - (q1 × (p1 - p0)) ÷ den
    ⟨p1, pNext, q1, 𝔽 pNext, k + 1⟩
  }
  _loop ← {
    ⟨p0, p1, q0, q1, k⟩ ← 𝕩
    stop ← ((k ≥ maxI) ∨ ((|q1) ≤ tol)) ∨ ((|q1 - q0) < 1e¯15)
    stop ? p1 ; 𝔽 _loop 𝔽 _step 𝕩
  }
  𝔽 _loop ⟨x0, x1, f0, f1, 0⟩
}

# Numerical Derivative via 5-point Central Difference
Deriv5Pt ⇐ {
  # 𝕨: grid spacing h
  # 𝕩: array of sample values [y_{i-2}, y_{i-1}, y_i, y_{i+1}, y_{i+2}]
  ⟨y0, y1, y2, y3, y4⟩ ← 𝕩
  ((-y4) + (8.0 × y3) - (8.0 × y1) + y0) ÷ (12.0 × 𝕨)
}
