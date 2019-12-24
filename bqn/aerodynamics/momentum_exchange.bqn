# Momentum Exchange Aerodynamic Force and Moment Evaluation Engine
#
# Computes instantaneous aerodynamic drag, lift, pitching moment, and non-dimensional
# force coefficients (Cd, Cl, Cm) using the momentum exchange method along fluid-solid
# boundary interfaces.
#
# Reference:
# Ladd, A. J. C., "Numerical simulations of particulate suspensions via a discretized
# Boltzmann equation. Part 1. Theoretical foundation", J. Fluid Mech. 271, 1994.
# Mei, R. et al., "Force evaluations in lattice Boltzmann simulations", J. Comp. Phys. 155, 1999.

math ← •Import "../lib/math_utils.bqn"
tensor ← •Import "../lib/tensor.bqn"
d2q9 ← •Import "../lbm/d2q9.bqn"

# Evaluate aerodynamic forces on solid obstacle from post-collision populations
# Arguments:
#   fPost: [9, Ny, Nx] post-collision distribution functions
#   solidMask: [Ny, Nx] boolean mask (1 = solid, 0 = fluid)
#   uWall: [2] wall velocity (usually [0, 0] for static airfoil)
#   rho0: freestream reference density
#   uInf: freestream reference velocity magnitude
#   chord: reference chord length (lattice units)
#   xRef: [2] moment reference center [x_quarter_chord, y_quarter_chord]
#
# Returns namespace containing:
#   fx, fy: dimensional total forces in lattice units (drag, lift)
#   moment: pitching moment about xRef (positive nose-up)
#   cd: drag coefficient
#   cl: lift coefficient
#   cm: pitching moment coefficient
ComputeForces ⇐ {
  # Arguments unpacked: ⟨fPost, solidMask, uWall, rho0, uInf, chord, xRef⟩
  ⟨fPost, solidMask, uWall, rho0, uInf, chord, xRef⟩ ← 𝕩

  ⟨ny, nx⟩ ← ≢solidMask
  cxVec ← d2q9.cx
  cyVec ← d2q9.cy
  opp ← d2q9.opp
  w ← d2q9.w

  # Initialize accumulated forces and moments
  # Iterate over all 8 non-zero directions (indices 1..8)
  # For each direction i, a boundary link occurs where fluid at (y, x) has solid neighbor at (y+cy, x+cx)
  fluidMask ← 1 - solidMask

  # Loop over directions 1 to 8 to compute link momentum transfers
  acc ← +´ {
    i ← 𝕩
    cx ← i ⊑ cxVec
    cy ← i ⊑ cyVec
    ci ← ⟨cx, cy⟩
    ibar ← i ⊑ opp
    wi ← i ⊑ w

    # Neighbor coordinate shifts
    # Fluid node (y, x) looking toward solid node (y+cy, x+cx)
    # Solid mask shifted so that solidMask at (y+cy, x+cx) aligns with (y, x)
    solidShifted ← (cy ⌽ solidMask)
    solidNeighbor ← (cx ⌽˘ solidShifted)

    # Active link mask: current node is fluid AND neighbor in direction i is solid
    # Mask out boundary wrap-around edges to avoid periodic false links
    linkMask ← fluidMask × solidNeighbor

    # Momentum transfer per link:
    # Delta P = ci * (f_i*(x_f) + f_ibar(x_f)) - 2*wi*rho*(ci . u_wall)/cs2
    fiPost ← i ⊑ fPost
    fibarPost ← ibar ⊑ fPost

    # Wall velocity correction term
    # uWallDotCi = (cx * uWall[0]) + (cy * uWall[1])
    uWallDotCi ← ((cx × (0 ⊑ uWall)) + (cy × (1 ⊑ uWall)))
    wallCorr ← (6.0 × wi) × (rho0 × uWallDotCi)

    # Momentum transfer per node in this direction
    deltaP ← (linkMask × ((fiPost + fibarPost) - wallCorr))

    # Total force contribution from this direction:
    # F = sum(deltaP) * ci
    sumP ← +´ ⥊ deltaP
    fDir ← (sumP × ci)

    # Pitching moment contribution:
    # Moment = sum( (x_b - x_ref) x Delta F )
    # Approximate link boundary location as midpoint: x_b = x + 0.5*cx, y_b = y + 0.5*cy
    # 2D cross product: (r_x * F_y) - (r_y * F_x)
    # Generate coordinates grids
    ⟨xGrid, yGrid⟩ ← tensor.Meshgrid2D ⟨nx, ny⟩
    rx ← (((xGrid + (0.5 × cx)) - (0 ⊑ xRef)))
    ry ← (((yGrid + (0.5 × cy)) - (1 ⊑ xRef)))

    momentMap ← (linkMask × deltaP) × (((rx × cy) - (ry × cx)))
    sumM ← +´ ⥊ momentMap

    # Return vector: [Fx, Fy, Moment]
    ⟨0 ⊑ fDir, 1 ⊑ fDir, sumM⟩
  }¨ (1 + ↕8)

  fxTotal ← 0 ⊑ acc
  fyTotal ← 1 ⊑ acc
  mTotal ← 2 ⊑ acc

  # Dynamic pressure: qInf = 0.5 * rho0 * uInf^2
  qInf ← (0.5 × rho0) × (uInf × uInf)
  denomForce ← (qInf × chord) ⌈ 1e¯12
  denomMoment ← ((qInf × chord) × chord) ⌈ 1e¯12

  cd ← (fxTotal ÷ denomForce)
  cl ← (fyTotal ÷ denomForce)
  cm ← (mTotal ÷ denomMoment)

  ⟨fxTotal, fyTotal, mTotal, cd, cl, cm⟩
}

# Sectional chordwise force distribution dFx/dx, dFy/dx
# Bins forces along chordwise sections for aerodynamic loading analysis
SectionalDistribution ⇐ {
  ⟨fPost, solidMask, xLead, chord, nBins⟩ ← 𝕩
  ⟨ny, nx⟩ ← ≢solidMask
  cxVec ← d2q9.cx
  cyVec ← d2q9.cy
  opp ← d2q9.opp
  fluidMask ← 1 - solidMask

  binWidth ← (chord ÷ nBins)
  ⟨xGrid, dummyY⟩ ← tensor.Meshgrid2D ⟨nx, ny⟩

  # For each bin, calculate localized forces
  bins ← {
    binIdx ← 𝕩
    xMin ← (xLead + (binIdx × binWidth))
    xMax ← (xMin + binWidth)
    inBin ← ((xGrid ≥ xMin) × (xGrid < xMax))

    accBin ← +´ {
      i ← 𝕩
      cx ← i ⊑ cxVec
      cy ← i ⊑ cyVec
      ci ← ⟨cx, cy⟩
      ibar ← i ⊑ opp
      solidNeighbor ← (cx ⌽˘ (cy ⌽ solidMask))
      linkMask ← (fluidMask × solidNeighbor) × inBin

      fiPost ← i ⊑ fPost
      fibarPost ← ibar ⊑ fPost
      deltaP ← (linkMask × (fiPost + fibarPost))
      sumP ← +´ ⥊ deltaP
      (sumP × ci)
    }¨ (1 + ↕8)

    ⟨0 ⊑ accBin, 1 ⊑ accBin⟩
  }¨ (↕nBins)

  bins
}
