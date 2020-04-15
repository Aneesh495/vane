# vane/bqn/lbm/trt.bqn
# Two-Relaxation-Time (TRT) collision operator with magic parameter boundary tuning.
# Note: Zero em dashes are used in this file per coding standards.

d2q9 ← •Import "d2q9.bqn"
opp ← d2q9.opp

# Symmetric and Anti-Symmetric Population Decomposition
DecomposeSymAnti ⇐ {
  # 𝕩: populations f (list of 9 matrices)
  f ← 𝕩
  fSym ← {
    i ← 𝕩
    iOpp ← i ⊑ opp
    fi ← i ⊑ f
    fiOpp ← iOpp ⊑ f
    0.5 × (fi + fiOpp)
  }¨ ↕9
  
  fAnti ← {
    i ← 𝕩
    iOpp ← i ⊑ opp
    fi ← i ⊑ f
    fiOpp ← iOpp ⊑ f
    0.5 × (fi - fiOpp)
  }¨ ↕9
  
  ⟨fSym, fAnti⟩
}

# Calculate Anti-Symmetric Relaxation Time from Magic Parameter Lambda
# Lambda = (tauPlus - 0.5) * (tauMinus - 0.5)
# Standard magic parameters:
# Lambda = 1/4 (optimal for bounce-back location stability)
# Lambda = 3/16 (optimal for steady-state advection-diffusion)
TauMinusFromLambda ⇐ {
  # 𝕨: magic parameter Lambda (default 0.25)
  # 𝕩: symmetric relaxation time tauPlus
  lambdaVal ← 𝕨
  tauPlus ← 𝕩
  diff ← tauPlus - 0.5
  den ← ((|diff) < 1e¯12) ⊑ ⟨diff, 1e¯6⟩
  0.5 + (lambdaVal ÷ den)
}

# Two-Relaxation-Time Collision Operator
CollideTRT ⇐ {
  # 𝕨: ⟨tauPlus, tauMinus⟩ or tauPlus
  isScalar ← (≠𝕨) ≤ 1
  wIn ← 𝕨 ⋄ tauPlus ← isScalar ◶ ⟨{0 ⊑ wIn}, {wIn}⟩ 0
  tauMinus ← isScalar ◶ ⟨{1 ⊑ wIn}, {0.5 + 0.25 ÷ (0.0001 ⌈ (wIn - 0.5))}⟩ 0
  f ← 𝕩 ⋄ mom ← d2q9.Moments f
  feq ← d2q9.Equilibrium mom
  
  ⟨fSym, fAnti⟩ ← DecomposeSymAnti f
  ⟨feqSym, feqAnti⟩ ← DecomposeSymAnti feq
  
  invTauPlus ← 1.0 ÷ tauPlus
  invTauMinus ← 1.0 ÷ tauMinus
  
  # Collide symmetric and anti-symmetric modes independently:
  # f* = f - (1/tauPlus)*(fSym - feqSym) - (1/tauMinus)*(fAnti - feqAnti)
  fPost ← {
    i ← 𝕩
    fi ← i ⊑ f
    fsi ← i ⊑ fSym
    feqsi ← i ⊑ feqSym
    fai ← i ⊑ fAnti
    feqai ← i ⊑ feqAnti
    
    fi - (invTauPlus × (fsi - feqsi)) - (invTauMinus × (fai - feqai))
  }¨ ↕9
  
  fPost
}

# Full D2Q9 TRT Time Step: TRT Collision followed by Streaming
StepTRT ⇐ {
  # 𝕨: relaxation parameters (tauPlus or [tauPlus, tauMinus])
  # 𝕩: populations f
  rates ← 𝕨
  f ← 𝕩
  fPost ← rates CollideTRT f
  d2q9.Stream fPost
}
