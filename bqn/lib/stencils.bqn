# vane/bqn/lib/stencils.bqn
# Discrete differential operators and finite-difference spatial stencils for aerodynamics.
# Note: Zero em dashes are used in this file per coding standards.

t ← •Import "tensor.bqn"
Roll2D ← t.Roll2D
Roll3D ← t.Roll3D

# 2D Periodic Gradient: d(phi)/dx and d(phi)/dy
Grad2DPeriodic ⇐ {
  # w: [dy, dx]
  # x: scalar field phi of shape ny x nx
  ⟨dy, dx⟩ ← 𝕨
  phi ← 𝕩
  dPhiDx ← ((0‿¯1 Roll2D phi) - (0‿1 Roll2D phi)) ÷ (2.0 × dx)
  dPhiDy ← ((¯1‿0 Roll2D phi) - (1‿0 Roll2D phi)) ÷ (2.0 × dy)
  ⟨dPhiDx, dPhiDy⟩
}

# 2D Bounded Gradient: central differences in interior, 1st-order forward/backward at walls
Grad2DBounded ⇐ {
  # w: [dy, dx]
  # x: scalar field phi of shape ny x nx
  ⟨dy, dx⟩ ← 𝕨
  phi ← 𝕩
  ⟨ny, nx⟩ ← ≢phi
  
  # Interior periodic stencil
  ⟨pDx, pDy⟩ ← 𝕨 Grad2DPeriodic phi
  
  # Fix boundary columns for dx:
  # col 0: (phi[y, 1] - phi[y, 0]) / dx
  # col nx-1: (phi[y, nx-1] - phi[y, nx-2]) / dx
  col0Dx ← (((1 ⌽˘ phi) - phi) ÷ dx)
  colLastDx ← ((phi - (¯1 ⌽˘ phi)) ÷ dx)
  
  # Fix boundary rows for dy:
  # row 0: (phi[1, x] - phi[0, x]) / dy
  # row ny-1: (phi[ny-1, x] - phi[ny-2, x]) / dy
  row0Dy ← (((1 ⌽ phi) - phi) ÷ dy)
  rowLastDy ← ((phi - (¯1 ⌽ phi)) ÷ dy)
  
  # Composite gradient
  # Mask for boundaries
  idxX ← ↕nx
  idxY ← ↕ny
  
  dPhiDx ← {
    iy ← 𝕩
    pRow ← iy ⊏ pDx
    c0Row ← iy ⊏ col0Dx
    cLRow ← iy ⊏ colLastDx
    {
      ix ← 𝕩
      (ix = 0) ⊑ ⟨((ix = (nx - 1)) ⊑ ⟨ix ⊏ pRow, ix ⊏ cLRow⟩), ix ⊏ c0Row⟩
    }¨ idxX
  }¨ idxY
  
  dPhiDy ← {
    iy ← 𝕩
    (iy = 0) ⊑ ⟨((iy = (ny - 1)) ⊑ ⟨iy ⊏ pDy, iy ⊏ rowLastDy⟩), iy ⊏ row0Dy⟩
  }¨ idxY
  
  ⟨> dPhiDx, > dPhiDy⟩
}

# 2D Divergence of Vector Field: div(u) = du_x/dx + du_y/dy
Div2D ⇐ {
  # w: [dy, dx]
  # x: [ux, uy]
  ⟨dy, dx⟩ ← 𝕨
  ⟨ux, uy⟩ ← 𝕩
  ⟨duxDx, dummyY⟩ ← 𝕨 Grad2DPeriodic ux
  ⟨dummyX, duyDy⟩ ← 𝕨 Grad2DPeriodic uy
  duxDx + duyDy
}

# 2D Vorticity / Curl of Velocity Field: omega_z = duy/dx - dux/dy
Curl2D ⇐ {
  # w: [dy, dx]
  # x: [ux, uy]
  ⟨dy, dx⟩ ← 𝕨
  ⟨ux, uy⟩ ← 𝕩
  ⟨duyDx, dummyY2⟩ ← 𝕨 Grad2DPeriodic uy
  ⟨dummyX2, duxDy⟩ ← 𝕨 Grad2DPeriodic ux
  duyDx - duxDy
}

# 2D Discrete 5-Point Laplacian: nabla^2 phi
Laplacian2D ⇐ {
  # w: [dy, dx]
  # x: scalar field phi
  ⟨dy, dx⟩ ← 𝕨
  phi ← 𝕩
  phiE ← 0‿¯1 Roll2D phi
  phiW ← 0‿1 Roll2D phi
  phiN ← ¯1‿0 Roll2D phi
  phiS ← 1‿0 Roll2D phi
  
  d2x ← (phiE + phiW - (2.0 × phi)) ÷ (dx × dx)
  d2y ← (phiN + phiS - (2.0 × phi)) ÷ (dy × dy)
  d2x + d2y
}

# 2D Symmetric Strain Rate Tensor S_ij = 0.5 * (du_i/dx_j + du_j/dx_i)
StrainRate2D ⇐ {
  # w: [dy, dx]
  # x: [ux, uy]
  ⟨dy, dx⟩ ← 𝕨
  ⟨ux, uy⟩ ← 𝕩
  ⟨duxDx, duxDy⟩ ← 𝕨 Grad2DPeriodic ux
  ⟨duyDx, duyDy⟩ ← 𝕨 Grad2DPeriodic uy
  
  sxx ← duxDx
  syy ← duyDy
  sxy ← 0.5 × (duxDy + duyDx)
  
  # Returns [sxx, syy, sxy]
  ⟨sxx, syy, sxy⟩
}

# Magnitude of 2D Strain Rate Tensor: |S| = sqrt(2 * S_ij * S_ij)
StrainRateMag2D ⇐ {
  # w: [dy, dx]
  # x: [ux, uy]
  ⟨sxx, syy, sxy⟩ ← 𝕨 StrainRate2D 𝕩
  #(2 * (sxx^2 + syy^2 + 2*sxy^2))^0.5
  √(2.0 × ((sxx × sxx) + ((syy × syy) + (2.0 × (sxy × sxy)))))
}

# 3D Vorticity / Curl: omega = nabla x u
Curl3D ⇐ {
  # w: [dz, dy, dx]
  # x: [ux, uy, uz]
  ⟨dz, dy, dx⟩ ← 𝕨
  ⟨ux, uy, uz⟩ ← 𝕩
  
  # Derivatives via periodic roll
  duzDy ← ((0‿¯1‿0 Roll3D uz) - (0‿1‿0 Roll3D uz)) ÷ (2.0 × dy)
  duyDz ← ((¯1‿0‿0 Roll3D uy) - (1‿0‿0 Roll3D uy)) ÷ (2.0 × dz)
  
  duxDz ← ((¯1‿0‿0 Roll3D ux) - (1‿0‿0 Roll3D ux)) ÷ (2.0 × dz)
  duzDx ← ((0‿0‿¯1 Roll3D uz) - (0‿0‿1 Roll3D uz)) ÷ (2.0 × dx)
  
  duyDx ← ((0‿0‿¯1 Roll3D uy) - (0‿0‿1 Roll3D uy)) ÷ (2.0 × dx)
  duxDy ← ((0‿¯1‿0 Roll3D ux) - (0‿1‿0 Roll3D ux)) ÷ (2.0 × dy)
  
  wx ← duzDy - duyDz
  wy ← duxDz - duzDx
  wz ← duyDx - duxDy
  
  ⟨wx, wy, wz⟩
}
