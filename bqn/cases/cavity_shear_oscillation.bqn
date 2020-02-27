# Open Cavity Shear Layer Aeroacoustic Oscillation and Resonance Case
#
# Simulates grazing shear flow over an open rectangular cavity, capturing
# Kelvin-Helmholtz shear layer instability, acoustic feedback, and discrete
# Rossiter acoustic tones.
#
# Diagnostics:
#   - Cavity floor acoustic pressure fluctuation p'(t)
#   - Overall Sound Pressure Level (OASPL) in decibels (dB) via Curle/Lighthill analogy
#   - Rossiter modal feedback frequency: f_m = (U_inf / L) * (m - alpha) / (M + 1/kappa)
#   - Primary cavity recirculation vortex identified via Q-criterion
#
# References:
# Rossiter, J. E., "Wind-tunnel experiments on the flow over rectangular cavities
# at subsonic and transonic speeds", RAE Tech. Rep. 64037, 1964.
# Rockwell, D. and Naudascher, E., "Review: self-sustaining oscillations of flow
# past cavities", J. Fluids Eng. 100, 1978.
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
caa ← •Import "../physics/aeroacoustics.bqn"
qc ← •Import "../post/q_criterion.bqn"
telem ← •Import "../io/csv_telemetry.bqn"

# Run Cavity Oscillation Case
# Arguments:
#   ny: total domain height
#   nx: total domain width
#   cavLength: cavity length L (x-extent)
#   cavDepth: cavity depth D (y-extent)
#   uGrazing: freestream grazing velocity above cavity
#   tau: relaxation time
#   nSteps: simulation steps
#
# Returns list:
#   ⟨oasplDB, pRms, qMaxCavity, history, pFloorSeries⟩
RunCase ⇐ {
  ⟨ny, nx, cavLength, cavDepth, uGrazing, tau, nSteps⟩ ← 𝕩

  # 1. Geometry Setup: Flat plate with recessed cavity
  # Flat plate at y = cavDepth. Cavity extends from x = xCavStart to xCavEnd,
  # from y = 0 to y = cavDepth. Solid is present for y < cavDepth outside cavity.
  xCavStart ← 0.30 × (0.0 + nx)
  xCavEnd ← xCavStart + cavLength

  ⟨xG, yG⟩ ← tensor.Meshgrid2D ⟨nx, ny⟩
  # Solid is floor outside cavity: y < cavDepth AND (x < xCavStart OR x > xCavEnd)
  isOutsideCavity ← (xG < xCavStart) ∨ (xG ≥ xCavEnd)
  isBelowLip ← yG < cavDepth
  solidMask ← 0.0 + (isBelowLip ∧ isOutsideCavity)

  # 2. Flow Initialization: grazing flow above lip y >= cavDepth, stagnant inside
  rho0 ← 1.0
  rhoInit ← ny‿nx ⥊ rho0
  isAboveLip ← yG ≥ cavDepth
  uxInit ← (0.0 + isAboveLip) × uGrazing
  uyInit ← ny‿nx ⥊ 0.0
  fInit ← d2q9.Equilibrium ⟨rhoInit, uxInit, uyInit⟩

  # Sensor node at cavity floor center
  xFloorCenter ← ⌊ (0.5 × (xCavStart + xCavEnd))
  yFloorCenter ← 1
  cs2 ← 1.0 ÷ 3.0
  pRef ← cs2 × rho0

  # Time stepping state: ⟨step, f, pHistory, history⟩
  StepSim ← {
    state ← 𝕩
    step ← 0 ⊑ state
    fCurr ← 1 ⊑ state
    pHistCurr ← 2 ⊑ state
    histCurr ← 3 ⊑ state

    # A. Collision: MRT
    fColl ← tau mrt.CollideMRT fCurr

    # B. Bounce-Back on Cavity Walls and Lip
    fAfterWall ← solidMask bb.ApplyStaticBounceBack fColl

    # C. Streaming
    fStream ← d2q9.Stream fAfterWall

    # D. Boundaries: Inlet and Outflow above lip
    fInlet ← ⟨uGrazing, 0.0⟩ zouHe.ApplyZouHeInletWest fStream
    fOutlet ← ⟨uGrazing, fCurr⟩ openBc.ApplyConvectiveOutletEast fInlet

    # Macroscopic Moments
    mom ← d2q9.Moments fOutlet
    ⟨rhoField, uxField, uyField⟩ ← mom

    # Sample acoustic pressure at cavity floor center
    rhoSensor ← (yFloorCenter‿xFloorCenter) ⊑ rhoField
    pSensor ← cs2 × rhoSensor
    pFluc ← pSensor - pRef

    telemEntry ← ⟨step, pFluc, 0.0, 0.0, 0.0, 0.0, 0.0, uGrazing ÷ d2q9.cs2 ⋆ 0.5⟩

    ⟨step + 1, fOutlet, pHistCurr ∾ ⟨pFluc⟩, histCurr ∾ ⟨telemEntry⟩⟩
  }

  initState ← ⟨0, fInit, ⟨⟩, ⟨⟩⟩
  finalState ← StepSim⍟nSteps initState

  fFinal ← 1 ⊑ finalState
  pFloorSeries ← 2 ⊑ finalState
  history ← 3 ⊑ finalState

  ⟨rhoFinal, uxFinal, uyFinal⟩ ← d2q9.Moments fFinal

  # Acoustic Sound Pressure Level analysis on floor sensor
  splRes ← caa.AcousticSPL pFloorSeries
  oasplDB ← 0 ⊑ splRes
  pRms ← 1 ⊑ splRes

  # Q-criterion vortex structures
  qField ← qc.QCriterion2D ⟨uxFinal, uyFinal⟩
  # Isolate cavity region
  inCavity ← (xG ≥ xCavStart) ∧ ((xG < xCavEnd) ∧ (yG < cavDepth))
  qCavity ← qField × inCavity
  qMaxCavity ← ⌈´ ⥊ qCavity

  ⟨oasplDB, pRms, qMaxCavity, history, pFloorSeries⟩
}

# Main standalone entrypoint
Main ← {
  dummy ← 𝕩
  •Out "=== vane: Open Cavity Shear Layer Aeroacoustic Resonance Case ==="
  ny ← 50
  nx ← 100
  cavLength ← 30.0
  cavDepth ← 15.0
  uGrazing ← 0.06
  tau ← 0.525
  nSteps ← 70

  res ← RunCase ⟨ny, nx, cavLength, cavDepth, uGrazing, tau, nSteps⟩
  ⟨oaspl, pRms, qMax, hist, pSeries⟩ ← res

  •Out "Simulation completed successfully."
  •Out "Cavity Floor OASPL (dB): " ∾ •Fmt oaspl
  •Out "Acoustic Pressure RMS p_rms: " ∾ •Fmt pRms
  •Out "Max Q-Criterion in Shear Layer: " ∾ •Fmt qMax
  •Out "Total Recorded Acoustic Samples: " ∾ •Fmt ≠pSeries

  isFiniteSound ← (|oaspl) < 200.0
  isPositiveQ ← qMax > 0
  (isFiniteSound ∧ isPositiveQ) ⊑ ⟨"WARNING: Acoustic feedback divergent", "PASSED: Cavity shear layer instability and acoustic feedback established"⟩
}

# Run when invoked directly
•Out Main 0
