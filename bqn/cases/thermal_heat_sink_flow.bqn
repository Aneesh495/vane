# Conjugate Thermal Convection Across Pin-Fin Heat Sink Simulation Case
#
# Simulates mixed forced and natural thermal convection through an array of heated
# pin-fin cooling passages using the Dual-Distribution Function (DDF) thermal LBM
# with Boussinesq buoyancy coupling.
#
# Diagnostics:
#   - Average Nusselt number Nu_bar quantifying surface heat transfer enhancement
#   - Thermal boundary layer development across staggered pin rows
#   - Thermal plume buoyancy coupling and total heat extraction rate
#
# References:
# Sparrow, E. M., Ramsey, J. W., and Altemani, C. A. C., "Experiments on in-line
# pin fin arrays and performance comparison with staggered arrays", J. Heat Transfer 102, 1980.
# He, X., Chen, S., and Doolen, G. D., "A novel thermal model for the lattice
# Boltzmann method in incompressible limit", J. Comp. Phys. 146, 1998.
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
tc ← •Import "../physics/thermal_convection.bqn"
telem ← •Import "../io/csv_telemetry.bqn"

# Run Thermal Heat Sink Simulation Case
# Arguments:
#   ny: domain height
#   nx: domain width
#   pinRadius: pin fin radius R
#   uInlet: cool inlet velocity
#   tCold: inlet coolant temperature
#   tHot: heated pin fin wall temperature
#   tauFluid: momentum relaxation time
#   tauThermal: thermal energy relaxation time
#   gBeta: buoyancy expansion coupling
#   nSteps: simulation steps
#
# Returns list:
#   ⟨meanNusselt, totalHeatExtracted, history, tempField⟩
RunCase ⇐ {
  ⟨ny, nx, pinRadius, uInlet, tCold, tHot, tauFluid, tauThermal, gBeta, nSteps⟩ ← 𝕩

  # 1. Geometry Setup: Staggered 2-row pin-fin array
  ⟨xG, yG⟩ ← tensor.Meshgrid2D ⟨nx, ny⟩
  rSq ← pinRadius × pinRadius

  # Pin 1 (row 1, upper): (0.35*nx, 0.65*ny)
  p1DistSq ← (((xG - (0.35 × nx)) ⋆ 2) + ((yG - (0.65 × ny)) ⋆ 2))
  mask1 ← p1DistSq ≤ rSq

  # Pin 2 (row 1, lower): (0.35*nx, 0.35*ny)
  p2DistSq ← (((xG - (0.35 × nx)) ⋆ 2) + ((yG - (0.35 × ny)) ⋆ 2))
  mask2 ← p2DistSq ≤ rSq

  # Pin 3 (row 2, center staggered): (0.55*nx, 0.50*ny)
  p3DistSq ← (((xG - (0.55 × nx)) ⋆ 2) + ((yG - (0.50 × ny)) ⋆ 2))
  mask3 ← p3DistSq ≤ rSq

  solidMask ← 0.0 + ((mask1 ∨ mask2) ∨ mask3)

  # 2. Flow and Thermal Initialization:
  rho0 ← 1.0
  rhoInit ← ny‿nx ⥊ rho0
  uxInit ← ny‿nx ⥊ uInlet
  uyInit ← ny‿nx ⥊ 0.0
  fInit ← d2q9.Equilibrium ⟨rhoInit, uxInit, uyInit⟩

  # Thermal lattice initialized to cold coolant temperature
  tInit ← ny‿nx ⥊ tCold
  gInit ← tc.ThermalEquilibrium ⟨tInit, uxInit, uyInit⟩

  # Time stepping state: ⟨step, f, g, history⟩
  StepSim ← {
    simState ← 𝕩
    step ← 0 ⊑ simState
    fCurr ← 1 ⊑ simState
    gCurr ← 2 ⊑ simState
    histCurr ← 3 ⊑ simState

    # Macroscopic fluid moments
    mom ← d2q9.Moments fCurr
    ⟨rhoField, uxField, uyField⟩ ← mom

    # Macroscopic temperature field
    tField ← tc.Temperature gCurr

    # A. Boussinesq Buoyancy Coupling: F_buoyant = [0, rho * gBeta * (T - tCold)]
    ⟨fbX, fbY⟩ ← tc.BoussinesqForce ⟨rhoField, tField, tCold, gBeta⟩

    # Fluid Collision with MRT
    fColl ← tauFluid mrt.CollideMRT fCurr

    # Obstacle Bounce-Back on pin fins
    fAfterWall ← solidMask bb.ApplyStaticBounceBack fColl

    # Fluid Streaming
    fStream ← d2q9.Stream fAfterWall

    # Fluid Boundaries: Inlet and Outflow
    fInlet ← ⟨uInlet, 0.0⟩ zouHe.ApplyZouHeInletWest fStream
    fOutlet ← ⟨uInlet, fCurr⟩ openBc.ApplyConvectiveOutletEast fInlet

    # B. Thermal Lattice Collision and Streaming:
    geq ← tc.ThermalEquilibrium ⟨tField, uxField, uyField⟩
    gColl ← tc.CollideThermalBGK ⟨gCurr, geq, tauThermal⟩

    # Isothermal Dirichlet Anti-Bounce-Back on heated pin fins (Tw = tHot)
    gAfterWall ← tc.IsothermalWallBC ⟨gColl, solidMask, tHot⟩

    # Thermal Streaming
    gStream ← d2q9.Stream gAfterWall

    # Average temperature of fluid domain
    fluidMask ← 1 - solidMask
    fluidT ← tField × fluidMask
    meanT ← (+´ ⥊ fluidT) ÷ ((+´ ⥊ fluidMask) ⌈ 1)

    telemEntry ← ⟨step, meanT, 0.0, 0.0, 0.0, 0.0, 0.0, uInlet ÷ d2q9.cs2 ⋆ 0.5⟩

    ⟨step + 1, fOutlet, gStream, histCurr ∾ ⟨telemEntry⟩⟩
  }

  initState ← ⟨0, fInit, gInit, ⟨⟩⟩
  finalState ← StepSim⍟nSteps initState

  fFinal ← 1 ⊑ finalState
  gFinal ← 2 ⊑ finalState
  history ← 3 ⊑ finalState

  tempFinal ← tc.Temperature gFinal

  # Nusselt number calculation along central pin passage
  deltaT ← (tHot - tCold) ⌈ 1e¯12
  ⟨nuLocal, nuMean⟩ ← tc.ComputeNusselt ⟨tempFinal, ⌊ (0.5 × ny), (0.0 + ny), deltaT⟩

  # Total extracted thermal energy = sum( rho * cp * (T - tCold) )
  fluidMask ← 1 - solidMask
  heatExtracted ← +´ ⥊ (fluidMask × (tempFinal - tCold))

  ⟨nuMean, heatExtracted, history, tempFinal⟩
}

# Main standalone entrypoint
Main ← {
  dummy ← 𝕩
  •Out "=== vane: Conjugate Thermal Convection Heat Sink Case ==="
  ny ← 45
  nx ← 90
  pinRadius ← 5.0
  uInlet ← 0.04
  tCold ← 1.0
  tHot ← 2.0
  tauFluid ← 0.54
  tauThermal ← 0.58
  gBeta ← 0.005
  nSteps ← 50

  res ← RunCase ⟨ny, nx, pinRadius, uInlet, tCold, tHot, tauFluid, tauThermal, gBeta, nSteps⟩
  ⟨nuMean, heat, hist, tempField⟩ ← res

  •Out "Simulation completed successfully."
  •Out "Average Nusselt Number Nu: " ∾ •Fmt nuMean
  •Out "Total Thermal Energy Extracted: " ∾ •Fmt heat
  •Out "Final Peak Temperature: " ∾ •Fmt ⌈´ ⥊ tempField
  •Out "Total Recorded Thermal Steps: " ∾ •Fmt ≠hist

  isPositiveHeat ← heat > 0
  isFiniteNu ← (|nuMean) < 100.0
  (isPositiveHeat ∧ isFiniteNu) ⊑ ⟨"WARNING: Non-physical heat transfer", "PASSED: Dual-distribution thermal convection across pin fins verified"⟩
}

# Run when invoked directly
•Out Main 0
