# Particle Streamline and Pathline Integration Engine
#
# Computes instantaneous streamlines and unsteady particle pathlines across
# aerodynamic flow fields using 4th-order Runge-Kutta (RK4) integration with
# subgrid bilinear velocity interpolation.
#
# Equations of motion:
#   dx / dt = u(x, y)
#   dy / dt = v(x, y)
#
# Note: Zero em dashes are used in this file per coding standards.

math ← •Import "../lib/math_utils.bqn"
tensor ← •Import "../lib/tensor.bqn"

# Bilinear interpolation of 2D scalar field at continuous coordinate (x, y)
# Arguments:
#   field: [Ny, Nx] 2D matrix
#   coord: [2] continuous coordinate ⟨x, y⟩
#
# Returns:
#   interpolated scalar value
BilinearInterp ⇐ {
  ⟨field, coord⟩ ← 𝕩
  ⟨ny, nx⟩ ← ≢field
  xCont ← 0 ⊑ coord
  yCont ← 1 ⊑ coord

  # Clamp coordinates inside domain
  xClamped ← 0.0 ⌈ ((nx - 1.001) ⌊ xCont)
  yClamped ← 0.0 ⌈ ((ny - 1.001) ⌊ yCont)

  x0 ← ⌊ xClamped
  y0 ← ⌊ yClamped
  x1 ← (nx - 1) ⌊ (x0 + 1)
  y1 ← (ny - 1) ⌊ (y0 + 1)

  fx ← xClamped - (0.0 + x0)
  fy ← yClamped - (0.0 + y0)

  # Four corner values
  v00 ← (y0‿x0) ⊑ field
  v01 ← (y0‿x1) ⊑ field
  v10 ← (y1‿x0) ⊑ field
  v11 ← (y1‿x1) ⊑ field

  # Bilinear blend
  top ← ((1.0 - fx) × v00) + (fx × v01)
  bot ← ((1.0 - fx) × v10) + (fx × v11)
  val ← ((1.0 - fy) × top) + (fy × bot)
  val
}

# Interpolate 2D velocity vector [ux, uy] at continuous coordinate
# Arguments:
#   ux: [Ny, Nx] horizontal velocity
#   uy: [Ny, Nx] vertical velocity
#   coord: [2] continuous coordinate ⟨x, y⟩
InterpolateVelocity ⇐ {
  ⟨ux, uy, coord⟩ ← 𝕩
  uVal ← BilinearInterp ⟨ux, coord⟩
  vVal ← BilinearInterp ⟨uy, coord⟩
  ⟨uVal, vVal⟩
}

# Single RK4 step for particle position
# Arguments:
#   coord: [2] current particle position ⟨x, y⟩
#   ux: [Ny, Nx] velocity field
#   uy: [Ny, Nx] velocity field
#   dt: integration step size (can be normalized step ds / ||u|| for streamlines)
#
# Returns:
#   nextCoord: [2] next position ⟨x_next, y_next⟩
StepRK4 ⇐ {
  ⟨coord, ux, uy, dt⟩ ← 𝕩

  k1 ← InterpolateVelocity ⟨ux, uy, coord⟩

  coord2 ← coord + ((0.5 × dt) × k1)
  k2 ← InterpolateVelocity ⟨ux, uy, coord2⟩

  coord3 ← coord + ((0.5 × dt) × k2)
  k3 ← InterpolateVelocity ⟨ux, uy, coord3⟩

  coord4 ← coord + (dt × k3)
  k4 ← InterpolateVelocity ⟨ux, uy, coord4⟩

  stepDelta ← (dt ÷ 6.0) × (k1 + ((2.0 × k2) + ((2.0 × k3) + k4)))
  coord + stepDelta
}

# Trace a single streamline from seed point for nSteps or until domain exit
# Arguments:
#   seed: [2] starting coordinate ⟨x0, y0⟩
#   ux: [Ny, Nx] velocity
#   uy: [Ny, Nx] velocity
#   dt: integration step
#   maxSteps: maximum step count
#
# Returns list:
#   trajectory: list of coordinates [ ⟨x0, y0⟩, ⟨x1, y1⟩, ... ]
TraceStreamline ⇐ {
  ⟨seed, ux, uy, dt, maxSteps⟩ ← 𝕩
  ⟨ny, nx⟩ ← ≢ux

  StepFn ← {
    state ← 𝕩
    curr ← 0 ⊑ state
    traj ← 1 ⊑ state
    isAlive ← 2 ⊑ state

    # Check domain bounds
    x ← 0 ⊑ curr
    y ← 1 ⊑ curr
    inBounds ← ((x ≥ 0) ∧ (x < (nx - 1))) ∧ ((y ≥ 0) ∧ (y < (ny - 1)))
    stillAlive ← isAlive ∧ inBounds

    nextPos ← stillAlive ⊑ ⟨curr, (StepRK4 ⟨curr, ux, uy, dt⟩)⟩
    nextTraj ← stillAlive ⊑ ⟨traj, (traj ∾ ⟨nextPos⟩)⟩

    ⟨nextPos, nextTraj, stillAlive⟩
  }

  init ← ⟨seed, ⟨seed⟩, 1⟩
  res ← StepFn⍟maxSteps init
  1 ⊑ res
}

# Generate rake of seeds and trace all streamlines
# Arguments:
#   seedList: list of seed coordinates [ ⟨x0, y0⟩, ⟨x1, y1⟩, ... ]
#   ux: [Ny, Nx] velocity
#   uy: [Ny, Nx] velocity
#   dt: step size
#   maxSteps: maximum steps per line
#
# Returns list of streamline trajectories
TraceRake ⇐ {
  ⟨seedList, ux, uy, dt, maxSteps⟩ ← 𝕩
  { TraceStreamline ⟨𝕩, ux, uy, dt, maxSteps⟩ }¨ seedList
}
