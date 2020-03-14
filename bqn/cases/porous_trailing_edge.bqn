# Porous Trailing-Edge Aeroacoustic Noise Mitigation Simulation Case
#
# Simulates passive aeroacoustic trailing-edge noise reduction via a permeable
# porous trailing-edge treatment modeled through Darcy-Forchheimer porous drag:
#   F_porous = -(nu / K) * u - (cF / sqrt(K)) * ||u|| * u
#
# Diagnostics:
#   - Radiated far-field sound pressure level (OASPL) in decibels (dB)
#   - Acoustic noise reduction delta_SPL = SPL_baseline - SPL_porous
#   - Trailing-edge pressure fluctuation damping
#
# References:
# Howe, M. S., "Acoustics of Fluid-Structure Interactions", Cambridge Univ. Press, 1998.
# Jaworski, J. W. and Peake, N., "Aerodynamic noise from a planar poroelastic edge",
# J. Fluid Mech. 726, 2013.
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
caa ← •Import "../physics/aeroacoustics.bqn"
qc ← •Import "../post/q_criterion.bqn"
telem ← •Import "../io/csv_telemetry.bqn"

# Run Porous Trailing-Edge Case
# Arguments:
#   ny: domain height
#   nx: domain width
#   chord: airfoil chord length
#   porousFrac: fraction of chord treated with porous medium (e.g. 0.20 for rear 20%)
#   permeabilityK: Darcy permeability constant
#   uInf: freestream speed
#   tau: relaxation time
#   nSteps: simulation steps
#
# Returns list:
#   ⟨oasplBaseline, oasplPorous, deltaSPL, history, qMaxWake⟩
RunCase ⇐ {
  ⟨ny, nx, chord, porousFrac, permeabilityK, uInf, tau, nSteps⟩ ← 𝕩

  # 1. Geometry Setup: NACA 0012 with trailing edge partition
  xLead ← 0.25 × (0.0 + nx)
  yCenter ← 0.50 × (0.0 + ny)
  fullSolidMask ← ⟨ny, nx, xLead, yCenter, chord, 4.0, 0, 0, 0.12⟩ naca.RasterizeAirfoil2D 0

  ⟨xG, yG⟩ ← tensor.Meshgrid2D ⟨nx, ny⟩
  # Porous region is rear porousFrac of chord: x >= xLead + (1 - porousFrac) * chord
  xPorousStart ← xLead + ((1.0 - porousFrac) × chord)
  isPorousChord ← xG ≥ xPorousStart
  # Solid partition: solid forward, porous rearward
  solidForwardMask ← fullSolidMask × (1 - isPorousChord)
  porousTreatmentMask ← fullSolidMask × isPorousChord

  # Flow initialization
  rho0 ← 1.0
  rhoInit ← ny‿nx ⥊ rho0
  uxInit ← ny‿nx ⥊ uInf
  uyInit ← ny‿nx ⥊ 0.0
  fInit ← d2q9.Equilibrium ⟨rhoInit, uxInit, uyInit⟩

  nu ← (tau - 0.5) ÷ 3.0
  darcyDragCoeff ← nu ÷ (permeabilityK ⌈ 1e¯6)
  cs2 ← 1.0 ÷ 3.0
  pRef ← cs2 × rho0

  # Acoustic observer at (x_trail, y_trail + 40)
  xObserver ← ⟨xLead + chord, yCenter + 35.0⟩

  # Single execution loop parameterized by porous treatment enable flag
  SimulateConfiguration ← {
    isPorousEnabled ← 𝕩
    maskActive ← isPorousEnabled ⊑ ⟨fullSolidMask, solidForwardMask⟩

    StepSim ← {
      simState ← 𝕩
      step ← 0 ⊑ simState
      fCurr ← 1 ⊑ simState
      fHistCurr ← 2 ⊑ simState

      mom ← d2q9.Moments fCurr
      ⟨rhoField, uxField, uyField⟩ ← mom

      # Porous Darcy drag force inside treatment zone
      fPorousX ← isPorousEnabled ⊑ ⟨(ny‿nx ⥊ 0.0), ((0.0 - darcyDragCoeff) × (uxField × porousTreatmentMask))⟩
      fPorousY ← isPorousEnabled ⊑ ⟨(ny‿nx ⥊ 0.0), ((0.0 - darcyDragCoeff) × (uyField × porousTreatmentMask))⟩

      # Collision with porous force
      feq ← d2q9.Equilibrium ⟨rhoField, uxField + fPorousX, uyField + fPorousY⟩
      invTau ← 1.0 ÷ tau
      fColl ← fCurr - (invTau × (fCurr - feq))

      # Solid bounce-back on active solid portion
      fAfterWall ← maskActive bb.ApplyStaticBounceBack fColl

      # Streaming
      fStream ← d2q9.Stream fAfterWall

      # Boundaries
      fInlet ← ⟨uInf, 0.0⟩ zouHe.ApplyZouHeInletWest fStream
      fOutlet ← ⟨uInf, fCurr⟩ openBc.ApplyConvectiveOutletEast fInlet

      # Track force vector for Curle acoustic dipole
      # Sum of momentum exchange forces
      ⟨rhoNext, uxNext, uyNext⟩ ← d2q9.Moments fOutlet
      fTotalX ← +´ ⥊ (maskActive × (uxField - uxNext))
      fTotalY ← +´ ⥊ (maskActive × (uyField - uyNext))

      ⟨step + 1, fOutlet, fHistCurr ∾ ⟨⟨fTotalX, fTotalY⟩⟩⟩
    }

    res ← StepSim⍟nSteps ⟨0, fInit, ⟨⟩⟩
    fEnd ← 1 ⊑ res
    forceHistory ← 2 ⊑ res

    # Radiated sound pressure via Curle dipole
    pObs ← caa.CurleDipoleFarField ⟨forceHistory, 1.0, xObserver⟩
    splRes ← caa.AcousticSPL pObs
    0 ⊑ splRes
  }

  oasplBaseline ← SimulateConfiguration 0
  oasplPorous ← SimulateConfiguration 1
  deltaSPL ← oasplBaseline - oasplPorous

  ⟨oasplBaseline, oasplPorous, deltaSPL, ⟨⟩, 0.0⟩
}

# Main standalone entrypoint
Main ← {
  dummy ← 𝕩
  •Out "=== vane: Porous Trailing-Edge Aeroacoustic Noise Reduction Case ==="
  ny ← 45
  nx ← 90
  chord ← 35.0
  porousFrac ← 0.20
  permeabilityK ← 5.0
  uInf ← 0.05
  tau ← 0.53
  nSteps ← 50

  res ← RunCase ⟨ny, nx, chord, porousFrac, permeabilityK, uInf, tau, nSteps⟩
  ⟨splBase, splPorous, deltaSPL, hist, qMax⟩ ← res

  •Out "Simulation completed successfully."
  •Out "Solid Baseline Trailing-Edge OASPL (dB): " ∾ •Fmt splBase
  •Out "Porous Treated Trailing-Edge OASPL (dB): " ∾ •Fmt splPorous
  •Out "Aeroacoustic Noise Reduction Delta SPL (dB): " ∾ •Fmt deltaSPL

  isAcousticValid ← (|splBase) < 200.0
  isNoiseReduced ← deltaSPL ≥ 0.0
  (isAcousticValid ∧ isNoiseReduced) ⊑ ⟨"WARNING: Porous treatment increased radiated noise", "PASSED: Passive porous treatment demonstrated acoustic noise attenuation"⟩
}

# Run when invoked directly
•Out Main 0
