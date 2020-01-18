# Two-Degree-of-Freedom (2-DOF) Aeroelastic Flutter Structural Dynamics Engine
#
# Solves coupled pitch-plunge structural equations of motion for an airfoil typical
# section subjected to unsteady aerodynamic lift and moment:
#   M * q_ddot + C * q_dot + K * q = Q_aero
#
# References:
# Fung, Y. C., "An Introduction to the Theory of Aeroelasticity", John Wiley, 1955.
# Bisplinghoff, R. L. and Ashley, H., "Principles of Aeroelasticity", Dover, 1962.
#
# Note: Zero em dashes are used in this file per coding standards.

math ← •Import "../lib/math_utils.bqn"

# Compute inverse mass matrix and structural parameters
# Arguments:
#   m: airfoil section mass per unit span
#   iAlpha: mass moment of inertia about elastic axis
#   xAlpha: static unbalance distance (distance between elastic axis and center of mass)
#   omegaH: uncoupled plunge natural frequency (rad/s)
#   omegaAlpha: uncoupled pitch natural frequency (rad/s)
#   zetaH: plunge structural damping ratio
#   zetaAlpha: pitch structural damping ratio
#
# Returns namespace list:
#   ⟨invM00, invM01, invM10, invM11, kh, kAlpha, ch, cAlpha⟩
StructuralModel ⇐ {
  ⟨m, iAlpha, xAlpha, omegaH, omegaAlpha, zetaH, zetaAlpha⟩ ← 𝕩

  # Static unbalance: S_alpha = m * x_alpha
  sAlpha ← m × xAlpha

  # Mass matrix determinant: det(M) = m * I_alpha - S_alpha^2
  detM ← ((m × iAlpha) - (sAlpha × sAlpha)) ⌈ 1e¯12
  invDet ← 1.0 ÷ detM

  # Inverse mass matrix:
  # invM = (1 / detM) * [ I_alpha, -S_alpha ; -S_alpha, m ]
  invM00 ← iAlpha × invDet
  invM01 ← (0.0 - sAlpha) × invDet
  invM10 ← (0.0 - sAlpha) × invDet
  invM11 ← m × invDet

  # Stiffness coefficients: K_h = m * omega_h^2, K_alpha = I_alpha * omega_alpha^2
  kh ← m × (omegaH × omegaH)
  kAlpha ← iAlpha × (omegaAlpha × omegaAlpha)

  # Damping coefficients: C_h = 2 * zeta_h * m * omega_h
  ch ← (2.0 × zetaH) × (m × omegaH)
  cAlpha ← (2.0 × zetaAlpha) × (iAlpha × omegaAlpha)

  ⟨invM00, invM01, invM10, invM11, kh, kAlpha, ch, cAlpha⟩
}

# State derivatives: d/dt [h, alpha, h_dot, alpha_dot]
# Arguments:
#   state: [4] vector [h, alpha, v_h, omega_alpha]
#   qAero: [2] aerodynamic forces [-Lift, Moment_ea]
#   struct: structural model parameters from StructuralModel
#
# Returns:
#   dStateDt: [4] vector [v_h, omega_alpha, a_h, alpha_ddot]
Derivatives ⇐ {
  ⟨state, qAero, struct⟩ ← 𝕩
  ⟨h, alpha, vh, omegaA⟩ ← state
  ⟨fAeroH, mAeroAlpha⟩ ← qAero
  ⟨invM00, invM01, invM10, invM11, kh, kAlpha, ch, cAlpha⟩ ← struct

  # Elastic restoring forces: F_spring = -K * q
  fSpringH ← 0.0 - (kh × h)
  mSpringAlpha ← 0.0 - (kAlpha × alpha)

  # Viscous damping forces: F_damper = -C * q_dot
  fDampH ← 0.0 - (ch × vh)
  mDampAlpha ← 0.0 - (cAlpha × omegaA)

  # Total generalized forces: R = Q_aero + F_spring + F_damper
  rh ← (fAeroH + fSpringH) + fDampH
  rAlpha ← (mAeroAlpha + mSpringAlpha) + mDampAlpha

  # Accelerations: q_ddot = inv(M) * R
  ah ← (invM00 × rh) + (invM01 × rAlpha)
  alphaDdot ← (invM10 × rh) + (invM11 × rAlpha)

  ⟨vh, omegaA, ah, alphaDdot⟩
}

# 4th-Order Runge-Kutta (RK4) Time Stepper for Aeroelastic Dynamics
# Arguments:
#   state: [4] initial state [h, alpha, v_h, omega_alpha]
#   qAero: [2] aerodynamic force vector [-Lift, Moment_ea]
#   struct: structural parameters from StructuralModel
#   dt: time step
#
# Returns:
#   nextState: [4] state integrated to t + dt
StepRK4 ⇐ {
  ⟨state, qAero, struct, dt⟩ ← 𝕩

  k1 ← Derivatives ⟨state, qAero, struct⟩

  state2 ← state + ((0.5 × dt) × k1)
  k2 ← Derivatives ⟨state2, qAero, struct⟩

  state3 ← state + ((0.5 × dt) × k2)
  k3 ← Derivatives ⟨state3, qAero, struct⟩

  state4 ← state + (dt × k3)
  k4 ← Derivatives ⟨state4, qAero, struct⟩

  # Weighted RK4 combination: state + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4)
  rkSum ← k1 + ((2.0 × k2) + ((2.0 × k3) + k4))
  nextState ← state + ((dt ÷ 6.0) × rkSum)
  nextState
}

# Flutter Dynamic Stability Evaluator
# Analyzes pitch angle time history to determine if response is:
#   -1: damped / stable (oscillations decay)
#    0: neutrally stable (limit cycle oscillation)
#   +1: divergent / unstable flutter (oscillations grow exponentially)
#
# Arguments:
#   alphaHistory: [N] time series of pitch angle alpha(t)
#
# Returns list:
#   ⟨stabilityFlag, growthRate, peakAmplitude⟩
StabilityCheck ⇐ {
  alphaHistory ← 𝕩
  n ← ≠alphaHistory

  # Divide history into initial half and final half to compare amplitudes
  halfN ← ⌊ (n ÷ 2)
  initSegment ← halfN ↑ alphaHistory
  finalSegment ← (0 - halfN) ↑ alphaHistory

  initAmp ← (⌈´ | initSegment) ⌈ 1e¯12
  finalAmp ← ⌈´ | finalSegment

  ratio ← finalAmp ÷ initAmp
  lnRatio ← ⋆⁼ (ratio ⌈ 1e¯12)
  growthRate ← lnRatio ÷ (0.0 + halfN)

  # Stability classification
  # Ratio < 0.95 -> damped (-1)
  # Ratio > 1.05 -> divergent flutter (+1)
  # Between 0.95 and 1.05 -> neutral limit cycle (0)
  isGrowing ← ratio > 1.05
  isDecaying ← ratio < 0.95
  flag ← isGrowing ⊑ ⟨(isDecaying ⊑ ⟨0, ¯1⟩), 1⟩

  ⟨flag, growthRate, finalAmp⟩
}
