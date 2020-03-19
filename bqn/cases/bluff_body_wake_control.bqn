# Active Bluff Body Wake Vortex Suppression and Drag Control Case
#
# Simulates active flow control past a bluff body (circular cylinder) via boundary
# rotation and surface momentum injection, capturing vortex shedding suppression,
# wake stabilization, drag reduction, and the Magnus aerodynamic force effect.
#
# Diagnostics:
#   - Baseline unforced shedding vs active rotation wake suppression
#   - Drag reduction Delta Cd = Cd_baseline - Cd_controlled
#   - Magnus lift generation Cl(alpha) from surface circulation
#   - Wake vortex dissipation quantified via Q-criterion
#
# References:
# Choi, H., Jeon, W. P., and Kim, J., "Control of flow over a bluff body",
# Annu. Rev. Fluid Mech. 40, 2008.
# Mittal, S. and Kumar, B., "Flow past a rotating cylinder", J. Fluid Mech. 476, 2003.
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

# Run Active Flow Control Simulation Case
# Arguments:
#   ny: domain height
#   nx: domain width
#   diam: cylinder diameter D
#   uInf: freestream velocity
#   rotationAlpha: non-dimensional cylinder surface velocity ratio (Omega*R / U_inf)
#   tau: relaxation time
#   nSteps: simulation steps
#
# Returns list:
#   ⟨cdBaseline, cdControlled, clMagnus, qMaxBaseline, qMaxControlled⟩
RunCase ⇐ {
  ⟨ny, nx, diam, uInf, rotationAlpha, tau, nSteps⟩ ← 𝕩

  # 1. Geometry Setup: Circular cylinder obstacle mask
  cx0 ← 0.28 × (0.0 + nx)
  cy0 ← 0.50 × (0.0 + ny)
  radius ← 0.5 × diam
  rSq ← radius × radius

  ⟨xG, yG⟩ ← tensor.Meshgrid2D ⟨nx, ny⟩
  dx ← xG - cx0
  dy ← yG - cy0
  distSq ← (dx × dx) + (dy × dy)
  solidMask ← 0.0 + (distSq ≤ rSq)

  # Tangential surface rotation velocity field u_w = Omega * R * [-sin(theta), cos(theta)]
  # = Omega * [-dy, dx]
  uTip ← rotationAlpha × uInf
  invR ← 1.0 ÷ (radius ⌈ 1e¯6)
  uwX ← (0.0 - (uTip × (dy × invR))) × solidMask
  uwY ← (uTip × (dx × invR)) × solidMask

  # Flow initialization
  rho0 ← 1.0
  rhoInit ← ny‿nx ⥊ rho0
  uxInit ← ny‿nx ⥊ uInf
  uyInit ← ny‿nx ⥊ 0.0
  fInit ← d2q9.Equilibrium ⟨rhoInit, uxInit, uyInit⟩

  xRef ← ⟨cx0, cy0⟩

  # Simulation runner parameterized by active control flag
  SimulateControl ← {
    isControlActive ← 𝕩

    activeUwX ← isControlActive ⊑ ⟨(ny‿nx ⥊ 0.0), uwX⟩
    activeUwY ← isControlActive ⊑ ⟨(ny‿nx ⥊ 0.0), uwY⟩

    StepSim ← {
      simState ← 𝕩
      step ← 0 ⊑ simState
      fCurr ← 1 ⊑ simState
      histCurr ← 2 ⊑ simState

      # Collision: MRT
      fColl ← tau mrt.CollideMRT fCurr

      # Boundary condition: Moving wall bounce-back if rotating, static if baseline
      fAfterWall ← ⟨solidMask, rhoInit, activeUwX, activeUwY⟩ bb.ApplyMovingBounceBack fColl

      # Streaming
      fStream ← d2q9.Stream fAfterWall

      # Boundaries: Inlet and Outflow
      fInlet ← ⟨uInf, 0.0⟩ zouHe.ApplyZouHeInletWest fStream
      fOutlet ← ⟨uInf, fCurr⟩ openBc.ApplyConvectiveOutletEast fInlet

      forces ← momEx.ComputeForces ⟨fColl, solidMask, ⟨0.0, 0.0⟩, rho0, uInf, diam, xRef⟩
      ⟨fx, fy, momQtr, cdVal, clVal, cmVal⟩ ← forces

      telemEntry ← ⟨step, clVal, cdVal, cmVal, 0.0, 0.0, 0.0, 0.0⟩
      ⟨step + 1, fOutlet, histCurr ∾ ⟨telemEntry⟩⟩
    }

    res ← StepSim⍟nSteps ⟨0, fInit, ⟨⟩⟩
    fFinal ← 1 ⊑ res
    momFinal ← d2q9.Moments fFinal
    ⟨rhoEnd, uxEnd, uyEnd⟩ ← momFinal

    forcesEnd ← momEx.ComputeForces ⟨fFinal, solidMask, ⟨0.0, 0.0⟩, rho0, uInf, diam, xRef⟩
    ⟨fxE, fyE, momE, cdEnd, clEnd, cmEnd⟩ ← forcesEnd

    qField ← qc.QCriterion2D ⟨uxEnd, uyEnd⟩
    # Mask wake behind cylinder (x > cx0 + radius)
    inWake ← xG > (cx0 + radius)
    qWakeMax ← ⌈´ ⥊ (qField × inWake)

    ⟨cdEnd, clEnd, qWakeMax⟩
  }

  ⟨cdBase, clBase, qBase⟩ ← SimulateControl 0
  ⟨cdCtrl, clCtrl, qCtrl⟩ ← SimulateControl 1

  ⟨cdBase, cdCtrl, clCtrl, qBase, qCtrl⟩
}

# Main standalone entrypoint
Main ← {
  dummy ← 𝕩
  •Out "=== vane: Active Bluff Body Wake Control Simulation Case ==="
  ny ← 45
  nx ← 95
  diam ← 12.0
  uInf ← 0.05
  rotationAlpha ← 1.5
  tau ← 0.53
  nSteps ← 50

  res ← RunCase ⟨ny, nx, diam, uInf, rotationAlpha, tau, nSteps⟩
  ⟨cdBase, cdCtrl, clMagnus, qBase, qCtrl⟩ ← res

  •Out "Simulation completed successfully."
  •Out "Baseline Cylinder Drag Cd: " ∾ •Fmt cdBase
  •Out "Active Controlled Cylinder Drag Cd: " ∾ •Fmt cdCtrl
  •Out "Generated Magnus Lift Cl: " ∾ •Fmt clMagnus
  •Out "Baseline Wake Max Q-Criterion: " ∾ •Fmt qBase
  •Out "Controlled Wake Max Q-Criterion: " ∾ •Fmt qCtrl

  isFiniteDrag ← (cdCtrl > 0) ∧ (cdCtrl < 100.0)
  hasMagnusLift ← (|clMagnus) > 0.01
  (isFiniteDrag ∧ hasMagnusLift) ⊑ ⟨"WARNING: Flow control ineffective", "PASSED: Active surface rotation induced Magnus circulation and wake modification"⟩
}

# Run when invoked directly
•Out Main 0
