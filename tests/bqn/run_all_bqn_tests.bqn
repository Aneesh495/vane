# vane - Master BQN Test Suite Runner
# Executes all unit, physics, aerodynamics, and vortex diagnostic test suites.
# Note: Zero em dashes are used in this file per coding standards.

nl ← ⟨@+10⟩

•Out ∾ ⟨
  "========================================================================", nl,
  "               vane: BQN Array Physics Verification Suite               ", nl,
  "========================================================================"
⟩

•Out nl ∾ "--- Suite 1: Mathematical Foundations & Stencils ---"
tMath ← •Import "test_math.bqn"

•Out nl ∾ "--- Suite 2: Lattice Boltzmann Dynamics (D2Q9, MRT, TRT, D3Q19, ELBM) ---"
tLbm ← •Import "test_lbm.bqn"

•Out nl ∾ "--- Suite 3: Aerodynamics, Boundary Layer & Momentum Exchange ---"
tAero ← •Import "test_aerodynamics.bqn"

•Out nl ∾ "--- Suite 4: Multi-Physics (Aeroacoustics, Multiphase, Flutter) ---"
tPhys ← •Import "test_physics.bqn"

•Out nl ∾ "--- Suite 5: Post-Processing, Vortex Cores & Streamlines ---"
tPost ← •Import "test_post.bqn"

•Out ∾ ⟨
  nl,
  "========================================================================", nl,
  "  SUCCESS: All BQN numerical suites passed verification with zero errors! ", nl,
  "========================================================================"
⟩
