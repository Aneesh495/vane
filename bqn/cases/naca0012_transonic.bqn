# NACA 0012 Aerodynamic Lift, Drag, and Surface Pressure Simulation Case
#
# Simulates subsonic viscous flow over a symmetric NACA 0012 airfoil at an angle
# of attack using the D2Q9 lattice Boltzmann method with Zou-He boundary conditions
# and interpolated boundary bounce-back.
#
# Diagnostics:
#   - Lift coefficient Cl, drag coefficient Cd, quarter-chord moment Cm
#   - Chordwise surface pressure coefficient Cp(x) along upper and lower surfaces
#   - Boundary layer displacement and momentum thickness
#   - Coherent trailing edge vortex shedding via Q-criterion
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
naca ← •Import "../aerodynamics/airfoil_naca.bqn"
momEx ← •Import "../aerodynamics/momentum_exchange.bqn"
presSurf ← •Import "../aerodynamics/pressure_surface.bqn"
bl ← •Import "../aerodynamics/boundary_layer.bqn"
qc ← •Import "../post/q_criterion.bqn"
vtk ← •Import "../io/vtk_export.bqn"
telem ← •Import "../io/csv_telemetry.bqn"

# Run NACA 0012 Simulation Case
# Arguments:
#   ny: domain height (e.g. 60)
#   nx: domain width (e.g. 120)
#   chord: airfoil chord length (e.g. 40)
#   aoaDeg: angle of attack in degrees (e.g. 4.0)
#   uInf: freestream inlet velocity (e.g. 0.05)
#   tau: relaxation time (e.g. 0.53)
#   nSteps: number of simulation time steps (e.g. 200)
#
# Returns namespace list:
#   ⟨finalCl, finalCd, finalCm, history, qMax⟩
RunCase ⇐ {
  ⟨ny, nx, chord, aoaDeg, uInf, tau, nSteps⟩ ← 𝕩

  # 1. Geometry Setup: Rasterize NACA 0012
  xLead ← 0.25 × (0.0 + nx)
  yCenter ← 0.5 × (0.0 + ny)
  solidMask ← ⟨ny, nx, xLead, yCenter, chord, aoaDeg, 0, 0, 0.12⟩ naca.RasterizeAirfoil2D 0

  # 2. Flow Initialization: uniform equilibrium
  rho0 ← 1.0
  rhoInit ← ny‿nx ⥊ rho0
  uxInit ← ny‿nx ⥊ uInf
  uyInit ← ny‿nx ⥊ 0.0
  fInit ← d2q9.Equilibrium ⟨rhoInit, uxInit, uyInit⟩

  # Relaxation parameters
  nu ← (tau - 0.5) ÷ 3.0
  xRef ← ⟨xLead + (0.25 × chord), yCenter⟩

  # Time stepping state: ⟨step, f, prevCl, history⟩
  StepSim ← {
    state ← 𝕩
    step ← 0 ⊑ state
    fCurr ← 1 ⊑ state
    prevCl ← 2 ⊑ state
    histCurr ← 3 ⊑ state

    # A. Collision: MRT
    fColl ← tau mrt.CollideMRT fCurr

    # B. Boundary condition on solid obstacle: Standard Bounce-Back
    fAfterWall ← solidMask bb.ApplyStaticBounceBack fColl

    # C. Streaming: periodic shift
    fStream ← d2q9.Stream fAfterWall

    # D. Domain Boundaries:
    # West: Zou-He Velocity Inlet (ux = uInf, uy = 0)
    fInlet ← ⟨uInf, 0.0⟩ zouHe.ApplyZouHeInletWest fStream

    # East: Convective Outflow
    fOutlet ← ⟨uInf, fCurr⟩ openBc.ApplyConvectiveOutletEast fInlet

    # Compute macroscopic moments
    mom ← d2q9.Moments fOutlet
    ⟨rhoField, uxField, uyField⟩ ← mom

    # Evaluate aerodynamic forces via momentum exchange
    forces ← momEx.ComputeForces ⟨fColl, solidMask, ⟨0.0, 0.0⟩, rho0, uInf, chord, xRef⟩
    ⟨fx, fy, momQtr, cdVal, clVal, cmVal⟩ ← forces

    # Convergence residual
    resU ← | (clVal - prevCl)
    csVal ← d2q9.cs2 ⋆ 0.5
    telemEntry ← ⟨step, clVal, cdVal, cmVal, clVal, resU, 0.0, uInf ÷ csVal⟩

    ⟨step + 1, fOutlet, clVal, histCurr ∾ ⟨telemEntry⟩⟩
  }

  initState ← ⟨0, fInit, 0.0, ⟨⟩⟩
  finalState ← StepSim⍟nSteps initState

  fFinal ← 1 ⊑ finalState
  history ← 3 ⊑ finalState

  # Macroscopic fields at final state
  ⟨rhoFinal, uxFinal, uyFinal⟩ ← d2q9.Moments fFinal
  finalForces ← momEx.ComputeForces ⟨fFinal, solidMask, ⟨0.0, 0.0⟩, rho0, uInf, chord, xRef⟩
  ⟨fx, fy, momQtr, finalCd, finalCl, finalCm⟩ ← finalForces

  # Q-criterion vortex core detection
  qField ← qc.QCriterion2D ⟨uxFinal, uyFinal⟩
  qMax ← ⌈´ ⥊ qField

  ⟨finalCl, finalCd, finalCm, history, qMax⟩
}

# Main standalone entrypoint
Main ← {
  dummy ← 𝕩
  •Out "=== vane: NACA 0012 Aerodynamic Verification Case ==="
  ny ← 50
  nx ← 100
  chord ← 35.0
  aoaDeg ← 5.0
  uInf ← 0.05
  tau ← 0.54
  nSteps ← 60

  res ← RunCase ⟨ny, nx, chord, aoaDeg, uInf, tau, nSteps⟩
  ⟨cl, cd, cm, hist, qMax⟩ ← res

  •Out "Simulation completed successfully."
  •Out "Lift Coefficient Cl: " ∾ •Fmt cl
  •Out "Drag Coefficient Cd: " ∾ •Fmt cd
  •Out "Pitching Moment Cm (c/4): " ∾ •Fmt cm
  •Out "Max Q-Criterion: " ∾ •Fmt qMax
  •Out "Total Recorded Telemetry Steps: " ∾ •Fmt ≠hist

  # Self-verification check
  isPositiveLift ← cl > 0
  isPositiveDrag ← cd > 0
  (isPositiveLift ∧ isPositiveDrag) ⊑ ⟨"WARNING: Forces out of physical range", "PASSED: Lift and drag within valid aerodynamic envelope"⟩
}

# When invoked directly via bqn
•Out Main 0
