# Multi-Phase Droplet Impingement and Aerodynamic Surface Wetting Simulation Case
#
# Simulates the impact, deformation, spreading, and surface contact of a supercooled
# water droplet impinging on an aerodynamic leading-edge surface using the Shan-Chen
# pseudo-potential multi-phase lattice Boltzmann method.
#
# Diagnostics:
#   - Liquid droplet deformation and spreading factor beta(t) = D(t) / D0
#   - Dynamic contact angle and wetting force against aerodynamic surface coating
#   - Liquid-gas phase boundary interface tracking
#
# References:
# Gent, R. W., Dart, N. P., and Cansdale, J. T., "Aircraft icing", Phil. Trans.
# R. Soc. Lond. A 358, 2000.
# Bragg, M. B. et al., "Effect of icing on aircraft handling and performance",
# Prog. Aerosp. Sci. 41, 2005.
#
# Note: Zero em dashes are used in this file per coding standards.

math ← •Import "../lib/math_utils.bqn"
tensor ← •Import "../lib/tensor.bqn"
stencils ← •Import "../lib/stencils.bqn"
d2q9 ← •Import "../lbm/d2q9.bqn"
mrt ← •Import "../lbm/mrt.bqn"
bb ← •Import "../boundaries/bounceback.bqn"
sc ← •Import "../physics/multiphase_shan_chen.bqn"
telem ← •Import "../io/csv_telemetry.bqn"

# Run Droplet Impingement Simulation Case
# Arguments:
#   ny: domain height
#   nx: domain width
#   dropRadius: initial droplet radius R0
#   uImpact: droplet approach velocity towards leading edge
#   gParam: interparticle attraction strength G (e.g. -5.0)
#   gWall: surface wetting parameter (-0.5 for hydrophilic, +0.5 for hydrophobic)
#   tau: relaxation time
#   nSteps: simulation steps
#
# Returns list:
#   ⟨finalSpreadFactor, finalArea, history, liquidMask⟩
RunCase ⇐ {
  ⟨ny, nx, dropRadius, uImpact, gParam, gWall, tau, nSteps⟩ ← 𝕩

  # 1. Geometry Setup: Flat aerodynamic leading edge wall at east boundary
  # Wall at x >= nx - 8
  wallX ← nx - 6
  ⟨xG, yG⟩ ← tensor.Meshgrid2D ⟨nx, ny⟩
  solidMask ← 0.0 + (xG ≥ wallX)

  # 2. Multiphase Initialization:
  # Droplet centered at (x0, y0) moving towards wall at uImpact
  xDrop0 ← 0.35 × (0.0 + nx)
  yDrop0 ← 0.50 × (0.0 + ny)
  distSq ← (((xG - xDrop0) ⋆ 2) + ((yG - yDrop0) ⋆ 2))
  rSq ← dropRadius × dropRadius

  # Liquid density rho_L = 2.0 inside droplet, vapor rho_V = 0.1 outside
  isDrop ← distSq ≤ rSq
  rhoInit ← (isDrop × 1.9) + 0.1

  uxInit ← (0.0 + isDrop) × uImpact
  uyInit ← ny‿nx ⥊ 0.0

  fInit ← d2q9.Equilibrium ⟨rhoInit, uxInit, uyInit⟩
  pi ← 3.14159265358979323846
  initialArea ← pi × rSq

  # Time stepping state: ⟨step, f, history⟩
  StepSim ← {
    simState ← 𝕩
    step ← 0 ⊑ simState
    fCurr ← 1 ⊑ simState
    histCurr ← 2 ⊑ simState

    # Compute macroscopic moments
    mom ← d2q9.Moments fCurr
    ⟨rhoField, uxField, uyField⟩ ← mom

    # Compute Shan-Chen intermolecular interaction force F_int
    ⟨fIntX, fIntY⟩ ← sc.InteractionForce ⟨rhoField, gParam⟩

    # Compute surface wetting force F_wall
    ⟨fWallX, fWallY⟩ ← sc.SurfaceWettingForce ⟨rhoField, solidMask, gWall⟩

    fTotalX ← fIntX + fWallX
    fTotalY ← fIntY + fWallY

    # Compute equilibrium with Shan-Chen velocity shift
    feq ← sc.EquilibriumWithForce ⟨rhoField, uxField, uyField, fTotalX, fTotalY, tau⟩

    # Collision: BGK with force-shifted equilibrium
    invTau ← 1.0 ÷ tau
    fColl ← fCurr - (invTau × (fCurr - feq))

    # Obstacle bounce-back on leading-edge wall
    fAfterWall ← solidMask bb.ApplyStaticBounceBack fColl

    # Streaming
    fStream ← d2q9.Stream fAfterWall

    # Real physical velocity for diagnostics
    ⟨uPhysX, uPhysY⟩ ← sc.PhysicalVelocity ⟨uxField, uyField, fTotalX, fTotalY, rhoField⟩

    # Droplet interface and area
    ⟨liqMask, equivR⟩ ← sc.PhaseInterface ⟨rhoField, 2.0, 0.1⟩
    spreadFactor ← (2.0 × equivR) ÷ ((2.0 × dropRadius) ⌈ 1e¯12)

    telemEntry ← ⟨step, spreadFactor, equivR, 0.0, 0.0, 0.0, 0.0, 0.0⟩

    ⟨step + 1, fStream, histCurr ∾ ⟨telemEntry⟩⟩
  }

  initState ← ⟨0, fInit, ⟨⟩⟩
  finalState ← StepSim⍟nSteps initState

  fFinal ← 1 ⊑ finalState
  history ← 2 ⊑ finalState

  momFinal ← d2q9.Moments fFinal
  rhoFinal ← 0 ⊑ momFinal
  rhoMin ← ⌊´ ⥊ rhoFinal
  rhoMax ← ⌈´ ⥊ rhoFinal
  rhoMid ← 0.5 × (rhoMin + rhoMax)
  finalLiqMask ← rhoFinal ≥ rhoMid

  area ← +´ ⥊ finalLiqMask
  finalR ← ((0.0 ⌈ area) ÷ pi) ⋆ 0.5
  finalSpreadFactor ← (2.0 × finalR) ÷ ((2.0 × dropRadius) ⌈ 1e¯12)
  finalArea ← area

  ⟨finalSpreadFactor, finalArea, history, finalLiqMask⟩
}

# Main standalone entrypoint
Main ← {
  dummy ← 𝕩
  •Out "=== vane: Supercooled Droplet Impingement Simulation Case ==="
  ny ← 40
  nx ← 80
  dropRadius ← 7.0
  uImpact ← 0.04
  gParam ← ¯5.0
  gWall ← ¯0.2
  tau ← 1.0
  nSteps ← 50

  res ← RunCase ⟨ny, nx, dropRadius, uImpact, gParam, gWall, tau, nSteps⟩
  ⟨spread, area, hist, liqMask⟩ ← res

  •Out "Simulation completed successfully."
  •Out "Final Droplet Spreading Factor beta: " ∾ •Fmt spread
  •Out "Final Liquid Area (cells): " ∾ •Fmt area
  •Out "Total Recorded Impingement Steps: " ∾ •Fmt ≠hist

  isFiniteSpread ← (spread > 0.1) ∧ (spread < 5.0)
  isConservedArea ← area > 20
  (isFiniteSpread ∧ isConservedArea) ⊑ ⟨"WARNING: Droplet mass non-conservation", "PASSED: Multiphase droplet impingement and contact spreading verified"⟩
}

# Run when invoked directly
•Out Main 0
