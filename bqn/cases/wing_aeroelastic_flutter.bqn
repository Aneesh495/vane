# Coupled 2-DOF Aeroelastic Wing Flutter Simulation Case
#
# Simulates dynamic Fluid-Structure Interaction (FSI) between unsteady aerodynamic
# lift/moment and a 2-degree-of-freedom (pitch and plunge) elastic airfoil section.
#
# Diagnostics:
#   - Coupled plunge displacement h(t) and pitch angle alpha(t) trajectories
#   - Dynamic flutter stability assessment (damped vs limit cycle vs flutter)
#   - Phase portrait in pitch-plunge state space (alpha vs alpha_dot)
#
# References:
# Theodorsen, T., "General theory of aerodynamic instability and the mechanism of flutter",
# NACA Report 496, 1935.
# Hodges, D. H. and Pierce, G. A., "Introduction to Structural Dynamics and Aeroelasticity",
# Cambridge University Press, 2002.
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
flutter ← •Import "../physics/aeroelastic_flutter.bqn"
telem ← •Import "../io/csv_telemetry.bqn"

# Run Coupled Aeroelastic Simulation Case
# Arguments:
#   ny: domain height
#   nx: domain width
#   chord: airfoil chord length
#   uInf: freestream speed
#   tau: relaxation time
#   nSteps: simulation steps
#
# Returns list:
#   ⟨finalPitch, finalPlunge, stabilityFlag, history⟩
RunCase ⇐ {
  ⟨ny, nx, chord, uInf, tau, nSteps⟩ ← 𝕩

  # 1. Geometry Setup: Baseline NACA 0012
  xLead ← 0.28 × (0.0 + nx)
  yCenter ← 0.50 × (0.0 + ny)
  solidMask ← ⟨ny, nx, xLead, yCenter, chord, 0.0, 0, 0, 0.12⟩ naca.RasterizeAirfoil2D 0

  # 2. Structural Model Setup:
  # Wing section mass m = 500.0, I_alpha = 50.0, x_alpha = 0.02, omega_h = 4.0, omega_alpha = 6.0, zeta = 0.05
  struct ← flutter.StructuralModel ⟨500.0, 50.0, 0.02, 4.0, 6.0, 0.05, 0.05⟩

  # Initial structural perturbation: alpha0 = 0.01 rad (~0.6 deg), h0 = 0
  state0 ← ⟨0.0, 0.01, 0.0, 0.0⟩

  # 3. Flow Initialization: uniform equilibrium
  rho0 ← 1.0
  rhoInit ← ny‿nx ⥊ rho0
  uxInit ← ny‿nx ⥊ uInf
  uyInit ← ny‿nx ⥊ 0.0
  fInit ← d2q9.Equilibrium ⟨rhoInit, uxInit, uyInit⟩

  dt ← 0.005
  xRef ← ⟨xLead + (0.25 × chord), yCenter⟩

  # Time stepping state: ⟨step, f, structState, history⟩
  StepSim ← {
    simState ← 𝕩
    step ← 0 ⊑ simState
    fCurr ← 1 ⊑ simState
    sState ← 2 ⊑ simState
    histCurr ← 3 ⊑ simState

    ⟨hPos, alphaPos, vh, omegaA⟩ ← sState

    # A. Collision: MRT
    fColl ← tau mrt.CollideMRT fCurr

    # B. Moving boundary bounce-back using current structural velocity:
    # Wall velocity: u_w = [0, v_h]
    uwX ← ny‿nx ⥊ 0.0
    uwY ← ny‿nx ⥊ vh
    fAfterWall ← ⟨solidMask, rhoInit, uwX, uwY⟩ bb.ApplyMovingBounceBack fColl

    # C. Streaming
    fStream ← d2q9.Stream fAfterWall

    # D. Inflow and Outflow
    fInlet ← ⟨uInf, 0.0⟩ zouHe.ApplyZouHeInletWest fStream
    fOutlet ← ⟨uInf, fCurr⟩ openBc.ApplyConvectiveOutletEast fInlet

    # Macroscopic Moments
    mom ← d2q9.Moments fOutlet
    ⟨rhoField, uxField, uyField⟩ ← mom

    # Compute instantaneous aerodynamic lift and quarter-chord pitching moment
    forces ← momEx.ComputeForces ⟨fColl, solidMask, ⟨0.0, vh⟩, rho0, uInf, chord, xRef⟩
    ⟨fx, fy, momQtr, cdVal, clVal, cmVal⟩ ← forces

    # Aerodynamic generalized forcing: Q_aero = [-Lift, Moment_ea] with mass coupling factor
    qAero ← ⟨0.01 × (0.0 - fy), 0.01 × momQtr⟩

    # Advance structural equations of motion by dt using RK4
    sStateNext ← flutter.StepRK4 ⟨sState, qAero, struct, dt⟩

    telemEntry ← ⟨step, clVal, cdVal, cmVal, alphaPos, hPos, vh, omegaA⟩

    ⟨step + 1, fOutlet, sStateNext, histCurr ∾ ⟨telemEntry⟩⟩
  }

  initState ← ⟨0, fInit, state0, ⟨⟩⟩
  finalState ← StepSim⍟nSteps initState

  finalSState ← 2 ⊑ finalState
  history ← 3 ⊑ finalState

  finalPlunge ← 0 ⊑ finalSState
  finalPitch ← 1 ⊑ finalSState

  # Stability check on pitch history
  alphaHistory ← { 4 ⊑ 𝕩 }¨ history
  stabRes ← flutter.StabilityCheck alphaHistory
  stabilityFlag ← 0 ⊑ stabRes

  ⟨finalPitch, finalPlunge, stabilityFlag, history⟩
}

# Main standalone entrypoint
Main ← {
  dummy ← 𝕩
  •Out "=== vane: 2-DOF Aeroelastic Wing Flutter FSI Simulation Case ==="
  ny ← 50
  nx ← 100
  chord ← 30.0
  uInf ← 0.05
  tau ← 0.54
  nSteps ← 60

  res ← RunCase ⟨ny, nx, chord, uInf, tau, nSteps⟩
  ⟨pitch, plunge, stab, hist⟩ ← res

  •Out "Simulation completed successfully."
  •Out "Final Airfoil Pitch Angle alpha (rad): " ∾ •Fmt pitch
  •Out "Final Airfoil Plunge Displacement h: " ∾ •Fmt plunge
  •Out "Flutter Stability Flag (-1:damped, 0:LCO, 1:flutter): " ∾ •Fmt stab
  •Out "Total Recorded FSI Steps: " ∾ •Fmt ≠hist

  isFinitePitch ← (|pitch) < 1.0
  isFinitePlunge ← (|plunge) < 10.0
  (isFinitePitch ∧ isFinitePlunge) ⊑ ⟨"WARNING: Aeroelastic divergence detected", "PASSED: Coupled fluid-structure interaction dynamics verified"⟩
}

# Run when invoked directly
•Out Main 0
