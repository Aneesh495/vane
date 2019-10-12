# vane/bqn/lib/tensor.bqn
# Multidimensional tensor slicing, rolling, padding, and coordinate utilities.
# Note: Zero em dashes are used in this file per coding standards.

# Linspace: n evenly spaced points between a and b inclusive
Linspace ⇐ {
  ⟨a, b, n⟩ ← 𝕩
  (n < 2) ◶ ⟨{ a + (↕n) × (b - a) ÷ (n - 1) }, { ⟨a⟩ }⟩ 𝕩
}

# Logspace: n points spaced logarithmically between 10^a and 10^b
Logspace ⇐ {
  ⟨a, b, n⟩ ← 𝕩
  lin ← Linspace ⟨a, b, n⟩
  10.0 ⋆ lin
}

# 2D Meshgrid: given nx and ny, generates 2D coordinate matrices x and y
# shape of returned matrices is ny x nx
Meshgrid2D ⇐ {
  ⟨nx, ny⟩ ← 𝕩
  xMat ← ny‿nx ⥊ ∾ ny ⥊ < ↕nx
  yMat ← ny‿nx ⥊ ∾ nx ⥊¨ ↕ny
  ⟨xMat, yMat⟩
}

# 2D Bounded Meshgrid: given [x_min, x_max, nx, y_min, y_max, ny]
MeshgridBounds2D ⇐ {
  ⟨xMin, xMax, nx, yMin, yMax, ny⟩ ← 𝕩
  xVec ← Linspace ⟨xMin, xMax, nx⟩
  yVec ← Linspace ⟨yMin, yMax, ny⟩
  xMat ← ny‿nx ⥊ ∾ ny ⥊ < xVec
  yMat ← ny‿nx ⥊ ∾ nx ⥊¨ yVec
  ⟨xMat, yMat⟩
}

# 2D Periodic Roll: shifts 2D array by [dy, dx]
Roll2D ⇐ {
  # w: [dy, dx]
  # x: 2D array
  ⟨dy, dx⟩ ← 𝕨
  (-dy) ⌽ ((-dx) ⌽˘ 𝕩)
}

# 3D Periodic Roll: shifts 3D array by [dz, dy, dx]
Roll3D ⇐ {
  # w: [dz, dy, dx]
  # x: 3D array of shape nz x ny x nx
  ⟨dz, dy, dx⟩ ← 𝕨
  (-dz) ⌽ ((-dy) ⌽˘ ((-dx) ⌽⎉1 𝕩))
}

# Clamp array elements between minVal and maxVal
Clamp ⇐ {
  # w: [minVal, maxVal]
  # x: array
  ⟨minVal, maxVal⟩ ← 𝕨
  minVal ⌈ (maxVal ⌊ 𝕩)
}

# Flatten multidimensional tensor into 1D vector
Flatten ⇐ {
  ⥊𝕩
}

# Reshape 1D vector or array into specified shape
Reshape ⇐ {
  # w: target shape vector
  # x: source array
  𝕨 ⥊ 𝕩
}

# Arithmetic Mean of Array
Mean ⇐ {
  (+´ ⥊𝕩) ÷ (≠ ⥊𝕩)
}

# Maximum element index in 1D array
ArgMax ⇐ {
  0 ⊑ ⍒ 𝕩
}

# Minimum element index in 1D array
ArgMin ⇐ {
  0 ⊑ ⍋ 𝕩
}

# Sum along axis 0
SumAxis0 ⇐ {
  +˝ 𝕩
}

# 2D Zero Pad: pads array with zeros by padY (rows) and padX (cols)
PadZero2D ⇐ {
  # w: [padY, padX]
  # x: 2D array
  ⟨padY, padX⟩ ← 𝕨
  src ← 𝕩
  ⟨ny, nx⟩ ← ≢src
  outNy ← ny + (2 × padY)
  outNx ← nx + (2 × padX)
  
  yOut ← {
    iy ← 𝕩
    isInnerY ← (iy ≥ padY) ∧ (iy < (padY + ny))
    rIdx ← 0 ⌈ ((ny - 1) ⌊ (iy - padY))
    rowIn ← rIdx ⊏ src
    isInnerY ⊑ ⟨outNx ⥊ 0.0, (padX ⥊ 0.0) ∾ rowIn ∾ (padX ⥊ 0.0)⟩
  }¨ ↕outNy
  > yOut
}
