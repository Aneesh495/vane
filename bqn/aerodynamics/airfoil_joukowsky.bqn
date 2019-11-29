# vane/bqn/aerodynamics/airfoil_joukowsky.bqn
# Conformal Joukowsky transformation and analytical potential flow airfoil benchmark.
# Note: Zero em dashes are used in this file per coding standards.

t ← •Import "../lib/tensor.bqn"
m ← •Import "../lib/math_utils.bqn"
Linspace ← t.Linspace
pi ← m.pi

# Generate Joukowsky Airfoil Coordinates
# aParam: conformal mapping scale parameter (e.g. 1.0)
# muX: center x-offset (< 0 controls thickness, e.g. -0.1)
# muY: center y-offset (> 0 controls camber, e.g. 0.08)
# aoaDeg: angle of attack in degrees
# nPoints: number of circumferential sample points
GenerateJoukowsky ⇐ {
  # 𝕩: ⟨aParam, muX, muY, aoaDeg, nPoints⟩
  ⟨aParam, muX, muY, aoaDeg, nPoints⟩ ← 𝕩
  
  # Radius of cylinder in zeta plane: R = sqrt((a - muX)^2 + muY^2)
  dxA ← aParam - muX
  rCyl ← √ ((dxA × dxA) + (muY × muY))
  
  # Angle array around circle: theta in [0, 2*pi]
  twoPi ← 2.0 × pi
  theta ← Linspace ⟨0.0, twoPi, nPoints⟩
  
  # Circle coordinates in zeta plane: xi + i*eta
  xi ← muX + (rCyl × •math.Cos theta)
  eta ← muY + (rCyl × •math.Sin theta)
  
  # Conformal transformation: z = zeta + a^2 / zeta
  # x = xi * (1 + a^2 / (xi^2 + eta^2))
  # y = eta * (1 - a^2 / (xi^2 + eta^2))
  rSq ← (xi × xi) + (eta × eta)
  safeRSq ← 1e¯12 ⌈ rSq
  factorA ← (aParam × aParam) ÷ safeRSq
  
  xPhys ← xi × (1.0 + factorA)
  yPhys ← eta × (1.0 - factorA)
  
  # Rotate by angle of attack alpha
  alphaRad ← aoaDeg × (pi ÷ 180.0)
  cosA ← •math.Cos alphaRad
  sinA ← •math.Sin alphaRad
  
  xRot ← (xPhys × cosA) + (yPhys × sinA)
  yRot ← (-(xPhys) × sinA) + (yPhys × cosA)
  
  ⟨xRot, yRot⟩
}

# Analytical Kutta-Joukowsky Theoretical Lift Coefficient
# C_L = 2 * pi * (1 + |muX|/a) * sin(alpha + beta)
TheoreticalLiftCoeff ⇐ {
  # 𝕩: ⟨aParam, muX, muY, aoaDeg⟩
  ⟨aParam, muX, muY, aoaDeg⟩ ← 𝕩
  alphaRad ← aoaDeg × (pi ÷ 180.0)
  
  # Camber angle beta = atan(muY / (a - muX))
  dxA ← aParam - muX
  betaRad ← •math.Atan (muY ÷ dxA)
  
  # Theoretical lift slope
  thickFactor ← 1.0 + ((|muX) ÷ aParam)
  clTheoretical ← 2.0 × pi × thickFactor × •math.Sin (alphaRad + betaRad)
  clTheoretical
}
