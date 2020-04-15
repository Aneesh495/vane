# vane - Multi-Physics Engine Test Suite
# Tests: Aeroacoustics, Multiphase Shan-Chen, Thermal Convection, Aeroelastic Flutter.
# Zero em dashes used in comments or strings.

ac  ← •Import "../../bqn/physics/aeroacoustics.bqn"
sc  ← •Import "../../bqn/physics/multiphase_shan_chen.bqn"
th  ← •Import "../../bqn/physics/thermal_convection.bqn"
fl  ← •Import "../../bqn/physics/aeroelastic_flutter.bqn"
d2q9 ← •Import "../../bqn/lbm/d2q9.bqn"

# Assertion helpers
Assert ← { 𝕨 ! 𝕩 }
AssertNear ← {
  diff ← | 𝕩 - 𝕨
  "Assertion near-equality failed" ! diff < 0.05
}

TestAeroacoustics ← {
  dummy ← 𝕩
  # Sound speed cs = 1 / sqrt(3) ~ 0.57735
  (1.0 ÷ √3.0) AssertNear ac.c0

  # Acoustic SPL of a sinusoidal pressure signal:
  # p'(t) = 1e-4 * sin(2*pi*i/50)
  pi ← 3.141592653589793
  pSignal ← { 1e¯4 × •math.Sin (𝕩 × 2.0 × pi ÷ 50.0) }¨ ↕100
  ⟨oaspl, rmsPressure⟩ ← ac.AcousticSPL pSignal

  # RMS of sine with amplitude A = 1e-4 is A / sqrt(2) ~ 7.071e-5
  expectedRMS ← 1e¯4 ÷ √2.0
  expectedRMS AssertNear rmsPressure

  # OASPL = 20 * log10(pRMS / pRef) with pRef = 1e-5:
  # ratio ~ 7.071 -> 20 * log10(7.071) ~ 16.99 dB
  expectedSPL ← 20.0 × ((⋆⁼ (expectedRMS ÷ ac.pRef)) ÷ ⋆⁼ 10.0)
  expectedSPL AssertNear oaspl

  •Out "  PASS: Computational aeroacoustics engine validated."
  1
}

TestShanChenMultiphase ← {
  dummy ← 𝕩
  # Test pseudo-potential: psi(0) = 1 - exp(0) = 0
  0.0 AssertNear (sc.Psi 0.0)

  # At high density rho -> inf: psi -> 1.0
  1.0 AssertNear (sc.Psi 10.0)

  # Uniform density field has zero interaction force due to symmetry
  ny ← 10
  nx ← 10
  rhoUniform ← ny‿nx ⥊ 1.2
  ⟨fX, fY⟩ ← sc.InteractionForce ⟨rhoUniform, ¯5.0⟩
  maxFX ← ⌈´ ⥊ | fX
  maxFY ← ⌈´ ⥊ | fY
  0.0 AssertNear maxFX
  0.0 AssertNear maxFY

  •Out "  PASS: Shan-Chen multiphase pseudo-potential validated."
  1
}

TestThermalConvection ← {
  dummy ← 𝕩
  # Macroscopic temperature recovers sum of thermal distribution
  # For 9 equal components of 0.2, T = 1.8
  gUniform ← { 4‿4 ⥊ 0.2 }¨ ↕9
  temp ← th.Temperature gUniform
  1.8 AssertNear (0‿0 ⊑ temp)

  •Out "  PASS: Dual-distribution thermal convection validated."
  1
}

TestAeroelasticFlutter ← {
  dummy ← 𝕩
  # Set up symmetric uncoupled 2-DOF wing section:
  # mass m = 10.0, iAlpha = 1.0, xAlpha = 0.0 (elastic axis on CG)
  # omegaH = 10 rad/s, omegaAlpha = 20 rad/s, zetaH = 0.05, zetaAlpha = 0.05
  struct ← fl.StructuralModel ⟨10.0, 1.0, 0.0, 10.0, 20.0, 0.05, 0.05⟩

  # Step forward from initial pitch angle alpha = 0.05 rad, zero velocity
  initState ← ⟨0.0, 0.05, 0.0, 0.0⟩
  qAeroZero ← ⟨0.0, 0.0⟩
  dt ← 0.01

  # Run 150 time steps using BQN repeat operator ⍟
  Step ← {
    ⟨s, h⟩ ← 𝕩
    sNext ← fl.StepRK4 ⟨s, qAeroZero, struct, dt⟩
    ⟨sNext, h ∾ ⟨1 ⊑ sNext⟩⟩
  }
  simRes ← Step⍟150 ⟨initState, ⟨0.05⟩⟩
  alphaHistory ← 1 ⊑ simRes

  ⟨flag, growthRate, peakAmp⟩ ← fl.StabilityCheck alphaHistory
  # Damped free vibration should yield stabilityFlag = -1
  "Aeroelastic structural damping decay" Assert (flag = ¯1)

  •Out "  PASS: 2-DOF aeroelastic flutter RK4 structural dynamics validated."
  1
}

# Run all test cases in this file
•Out "Running Multi-Physics Engine Tests:"
TestAeroacoustics ⟨⟩
TestShanChenMultiphase ⟨⟩
TestThermalConvection ⟨⟩
TestAeroelasticFlutter ⟨⟩
•Out "All Physics tests passed successfully."
