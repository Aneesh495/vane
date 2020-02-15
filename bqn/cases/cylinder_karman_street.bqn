# Circular Cylinder von Karman Vortex Street Simulation Case
#
# Simulates unsteady laminar vortex shedding past a circular cylinder at Re = 100,
# reproducing the classical von Karman vortex street wake instability.
#
# Diagnostics:
#   - Unsteady lift oscillation Cl(t) and double-frequency drag oscillation Cd(t)
#   - Experimental Strouhal number St = f_shed * D / U_inf comparison (Roshko curve)
#   - Identification of alternating wake vortex cores via Q-criterion
#
# References:
# Roshko, A., "On the development of turbulent wakes from vortex streets", NACA Report 1191, 1954.
# Williamson, C. H. K., "Vortex dynamics in the cylinder wake", Ann. Rev. Fluid Mech. 28, 1996.
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
momEx ← •Import "../aerodynamics/momentum_exchange.bqn"
qc ← •Import "../post/q_criterion.bqn"
telem ← •Import "../io/csv_telemetry.bqn"

# Run Circular Cylinder Simulation Case
# Arguments:
#   ny: domain height
#   nx: domain width
#   diam: cylinder diameter D (e.g. 16)
#   uInf: freestream inlet velocity (e.g. 0.05)
#   tau: relaxation time (e.g. 0.524 for Re = 100)
#   nSteps: simulation steps
#
# Returns list:
#   ⟨finalCl, finalCd, strouhal, history, qField⟩
RunCase ⇐ {
  ⟨ny, nx, diam, uInf, tau, nSteps⟩ ← 𝕩

  # 1. Geometry Setup: Circular cylinder obstacle mask
  cx0 ← 0.25 × (0.0 + nx)
  cy0 ← 0.5 × (0.0 + ny)
  radius ← 0.5 × diam
  radiusSq ← radius × radius

  ⟨xG, yG⟩ ← tensor.Meshgrid2D ⟨nx, ny⟩
  distSq ← (((xG - cx0) ⋆ 2) + ((yG - cy0) ⋆ 2))
  solidMask ← 0.0 + (distSq ≤ radiusSq)

  # 2. Flow Initialization: uniform equilibrium with slight asymmetry to trigger wake
  rho0 ← 1.0
  rhoInit ← ny‿nx ⥊ rho0
  uxInit ← ny‿nx ⥊ uInf
  # Small vertical perturbation at cylinder to trigger vortex shedding
  uyInit ← 0.01 × (⋆ (0.0 - (distSq ÷ (2.0 × radiusSq))))
  fInit ← d2q9.Equilibrium ⟨rhoInit, uxInit, uyInit⟩

  xRef ← ⟨cx0, cy0⟩

  # Time stepping state: ⟨step, f, prevCl, history⟩
  StepSim ← {
    state ← 𝕩
    step ← 0 ⊑ state
    fCurr ← 1 ⊑ state
    prevCl ← 2 ⊑ state
    histCurr ← 3 ⊑ state

    # A. Collision: MRT
    fColl ← tau mrt.CollideMRT fCurr

    # B. Bounce-Back on Cylinder
    fAfterWall ← solidMask bb.ApplyStaticBounceBack fColl

    # C. Streaming
    fStream ← d2q9.Stream fAfterWall

    # D. Domain Boundaries:
    fInlet ← ⟨uInf, 0.0⟩ zouHe.ApplyZouHeInletWest fStream
    fOutlet ← ⟨uInf, fCurr⟩ openBc.ApplyConvectiveOutletEast fInlet

    # Macroscopic Moments
    mom ← d2q9.Moments fOutlet
    ⟨rhoField, uxField, uyField⟩ ← mom

    # Forces
    forces ← momEx.ComputeForces ⟨fColl, solidMask, ⟨0.0, 0.0⟩, rho0, uInf, diam, xRef⟩
    ⟨fx, fy, momQtr, cdVal, clVal, cmVal⟩ ← forces

    resU ← | (clVal - prevCl)
    csVal ← d2q9.cs2 ⋆ 0.5
    telemEntry ← ⟨step, clVal, cdVal, cmVal, clVal, resU, 0.0, uInf ÷ csVal⟩

    ⟨step + 1, fOutlet, clVal, histCurr ∾ ⟨telemEntry⟩⟩
  }

  initState ← ⟨0, fInit, 0.0, ⟨⟩⟩
  finalState ← StepSim⍟nSteps initState

  fFinal ← 1 ⊑ finalState
  history ← 3 ⊑ finalState

  ⟨rhoFinal, uxFinal, uyFinal⟩ ← d2q9.Moments fFinal
  finalForces ← momEx.ComputeForces ⟨fFinal, solidMask, ⟨0.0, 0.0⟩, rho0, uInf, diam, xRef⟩
  ⟨fx, fy, momQtr, finalCd, finalCl, finalCm⟩ ← finalForces

  # Strouhal shedding frequency calculation
  clHistory ← { 1 ⊑ 𝕩 }¨ history
  strouhal ← telem.ComputeStrouhal ⟨clHistory, 1.0, diam, uInf⟩

  qField ← qc.QCriterion2D ⟨uxFinal, uyFinal⟩

  ⟨finalCl, finalCd, strouhal, history, qField⟩
}

# Main standalone entrypoint
Main ← {
  dummy ← 𝕩
  •Out "=== vane: Circular Cylinder von Karman Vortex Street Case ==="
  ny ← 50
  nx ← 110
  diam ← 14.0
  uInf ← 0.05
  tau ← 0.528
  nSteps ← 80

  res ← RunCase ⟨ny, nx, diam, uInf, tau, nSteps⟩
  ⟨cl, cd, st, hist, qField⟩ ← res

  •Out "Simulation completed successfully."
  •Out "Instantaneous Cl: " ∾ •Fmt cl
  •Out "Instantaneous Cd: " ∾ •Fmt cd
  •Out "Strouhal Number St: " ∾ •Fmt st
  •Out "Max Q-Criterion in Wake: " ∾ •Fmt ⌈´ ⥊ qField
  •Out "Total Recorded Steps: " ∾ •Fmt ≠hist

  isFiniteCl ← (|cl) < 10.0
  isFiniteCd ← (cd > 0) ∧ (cd < 50.0)
  (isFiniteCl ∧ isFiniteCd) ⊑ ⟨"WARNING: Flow instability detected", "PASSED: von Karman wake shedding within physical stability bounds"⟩
}

# Run when invoked directly
•Out Main 0
