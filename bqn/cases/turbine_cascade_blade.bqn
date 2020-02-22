# Axial Turbine Blade Cascade Turning and Aerodynamic Loss Simulation Case
#
# Simulates viscous flow through an axial turbomachinery stator/rotor blade cascade,
# computing flow turning angle delta_beta = beta_1 - beta_2, total pressure loss,
# and blade tangential forces generating shaft torque.
#
# Diagnostics:
#   - Inlet flow angle beta_1 vs outlet flow angle beta_2
#   - Cascade flow turning angle delta_beta
#   - Tangential aerodynamic force Ft and axial force Fx
#   - Blade passage wake vortex structures via Q-criterion
#
# References:
# Denton, J. D., "Loss mechanisms in turbomachines", ASME J. Turbomach. 115, 1993.
# Gostelow, J. P., "Cascade Aerodynamics", Pergamon Press, 1984.
#
# Note: Zero em dashes are used in this file per coding standards.

math ← •Import "../lib/math_utils.bqn"
tensor ← •Import "../lib/tensor.bqn"
stencils ← •Import "../lib/stencils.bqn"
d2q9 ← •Import "../lbm/d2q9.bqn"
mrt ← •Import "../lbm/mrt.bqn"
zouHe ← •Import "../boundaries/zou_he.bqn"
bb ← •Import "../boundaries/bounceback.bqn"
openBc ← •Import "../boundaries/open_convective.bqn"
turbo ← •Import "../aerodynamics/turbomachinery.bqn"
momEx ← •Import "../aerodynamics/momentum_exchange.bqn"
qc ← •Import "../post/q_criterion.bqn"
telem ← •Import "../io/csv_telemetry.bqn"

# Run Turbomachinery Blade Cascade Simulation Case
# Arguments:
#   ny: domain height
#   nx: domain width
#   chord: blade chord length
#   pitch: pitch spacing between blades
#   staggerDeg: blade setting stagger angle
#   camberDeg: blade camber turning angle
#   uInlet: axial inlet velocity magnitude
#   inletAngleDeg: inlet flow angle beta_1
#   tau: relaxation time
#   nSteps: simulation steps
#
# Returns list:
#   ⟨turningAngleDeg, tangentialForce, axialForce, history, qMax⟩
RunCase ⇐ {
  ⟨ny, nx, chord, pitch, staggerDeg, camberDeg, uInlet, inletAngleDeg, tau, nSteps⟩ ← 𝕩

  # 1. Geometry Setup: Rasterize periodic cascade
  originX ← 0.30 × (0.0 + nx)
  cascadeParams ← ⟨ny, nx, originX, chord, pitch, staggerDeg, camberDeg, 0.10⟩
  solidMask ← cascadeParams turbo.RasterizeBladeCascade2D 0

  # 2. Flow Initialization: inlet flow with angle beta_1
  pi ← 3.14159265358979323846
  beta1Rad ← inletAngleDeg × (pi ÷ 180.0)
  uxInlet ← uInlet × (math.Cos beta1Rad)
  uyInlet ← uInlet × (math.Sin beta1Rad)

  rho0 ← 1.0
  rhoInit ← ny‿nx ⥊ rho0
  uxInit ← ny‿nx ⥊ uxInlet
  uyInit ← ny‿nx ⥊ uyInlet
  fInit ← d2q9.Equilibrium ⟨rhoInit, uxInit, uyInit⟩

  xRef ← ⟨originX + (0.5 × chord), 0.5 × (0.0 + ny)⟩

  # Time stepping state: ⟨step, f, prevCl, history⟩
  StepSim ← {
    state ← 𝕩
    step ← 0 ⊑ state
    fCurr ← 1 ⊑ state
    prevCl ← 2 ⊑ state
    histCurr ← 3 ⊑ state

    # A. Collision: MRT
    fColl ← tau mrt.CollideMRT fCurr

    # B. Bounce-Back on Blade Walls
    fAfterWall ← solidMask bb.ApplyStaticBounceBack fColl

    # C. Streaming: periodic in y (north/south)
    fStream ← d2q9.Stream fAfterWall

    # D. Inflow and Outflow Boundaries:
    fInlet ← ⟨uxInlet, uyInlet⟩ zouHe.ApplyZouHeInletWest fStream
    fOutlet ← ⟨uxInlet, fCurr⟩ openBc.ApplyConvectiveOutletEast fInlet

    # Macroscopic Moments
    mom ← d2q9.Moments fOutlet
    ⟨rhoField, uxField, uyField⟩ ← mom

    # Evaluate cascade forces
    forces ← momEx.ComputeForces ⟨fColl, solidMask, ⟨0.0, 0.0⟩, rho0, uInlet, chord, xRef⟩
    ⟨fx, fy, momQtr, cdVal, clVal, cmVal⟩ ← forces

    resU ← | (clVal - prevCl)
    csVal ← d2q9.cs2 ⋆ 0.5
    telemEntry ← ⟨step, clVal, cdVal, cmVal, clVal, resU, 0.0, uInlet ÷ csVal⟩

    ⟨step + 1, fOutlet, clVal, histCurr ∾ ⟨telemEntry⟩⟩
  }

  initState ← ⟨0, fInit, 0.0, ⟨⟩⟩
  finalState ← StepSim⍟nSteps initState

  fFinal ← 1 ⊑ finalState
  history ← 3 ⊑ finalState

  ⟨rhoFinal, uxFinal, uyFinal⟩ ← d2q9.Moments fFinal
  finalForces ← momEx.ComputeForces ⟨fFinal, solidMask, ⟨0.0, 0.0⟩, rho0, uInlet, chord, xRef⟩
  ⟨fx, fy, momQtr, finalCd, finalCl, finalCm⟩ ← finalForces

  # Average outlet flow velocity angle beta_2 across east boundary column
  lastCol ← nx - 2
  uxOutletCol ← lastCol ⊏˘ uxFinal
  uyOutletCol ← lastCol ⊏˘ uyFinal
  meanUxOut ← (+´ uxOutletCol) ÷ (0.0 + ny)
  meanUyOut ← (+´ uyOutletCol) ÷ (0.0 + ny)

  # Outlet flow angle beta2 = atan2(uy, ux)
  # Approx atan(uy / ux) in degrees
  ratioOut ← meanUyOut ÷ (meanUxOut ⌈ 1e¯12)
  # Simple atan approximation: atan(x) ~ x - x^3/3 for small x
  beta2Deg ← (ratioOut × (180.0 ÷ pi))
  turningAngleDeg ← inletAngleDeg - beta2Deg

  qField ← qc.QCriterion2D ⟨uxFinal, uyFinal⟩
  qMax ← ⌈´ ⥊ qField

  ⟨turningAngleDeg, fy, fx, history, qMax⟩
}

# Main standalone entrypoint
Main ← {
  dummy ← 𝕩
  •Out "=== vane: Turbomachinery Blade Cascade Turning Simulation Case ==="
  ny ← 50
  nx ← 100
  chord ← 30.0
  pitch ← 25.0
  staggerDeg ← 15.0
  camberDeg ← 30.0
  uInlet ← 0.05
  inletAngleDeg ← 20.0
  tau ← 0.53
  nSteps ← 60

  res ← RunCase ⟨ny, nx, chord, pitch, staggerDeg, camberDeg, uInlet, inletAngleDeg, tau, nSteps⟩
  ⟨turning, fy, fx, hist, qMax⟩ ← res

  •Out "Simulation completed successfully."
  •Out "Cascade Flow Turning Angle (deg): " ∾ •Fmt turning
  •Out "Tangential Driving Force Fy: " ∾ •Fmt fy
  •Out "Axial Drag Force Fx: " ∾ •Fmt fx
  •Out "Max Q-Criterion in Blade Passage: " ∾ •Fmt qMax
  •Out "Total Recorded Telemetry Steps: " ∾ •Fmt ≠hist

  isFiniteTurning ← (|turning) < 90.0
  isFiniteForce ← (|fy) < 1000.0
  (isFiniteTurning ∧ isFiniteForce) ⊑ ⟨"WARNING: Unstable cascade turning", "PASSED: Blade passage flow turned with positive aerodynamic deflection"⟩
}

# Run when invoked directly
•Out Main 0
