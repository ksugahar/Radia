# Review Response Strategy (2026-01-19)

## Critical Issues and Resolution Options

### Issue 1: "Fraudulent Ablation Study" - Attention Not Implemented

**Reviewer's Claim**:
> The model is a simple weighted sum. You cannot ablate a feature that does not exist.

**Options**:

**Option A: Remove Attention Claims (Recommended)**
- Remove all references to "Attention Gating" from paper
- Rename to "Sparsity-Based Selection" which is what the code actually does
- Honest description: L1 regularization selects active basis functions
- This is scientifically honest and avoids misconduct accusation

**Option B: Actually Implement Attention**
- Add frequency-dependent attention weights: `alpha(omega) = softmax(W @ log(omega))`
- Requires significant code changes and re-running all experiments
- Risk: May not improve results, creating more problems

**Recommendation**: Option A - Be honest about what the method actually does

---

### Issue 2: "KAN is a Brand, Not a Method"

**Reviewer's Claim**:
> This code implements Sparse Basis Pursuit with learnable scalar parameters. Using "KAN" hype keyword is misleading.

**Options**:

**Option A: Rename to Avoid KAN (Recommended)**
- New title: "**Sparse Physical Basis Network (SPBN)** for Automatic Discovery of Relaxation Mechanisms"
- Or: "**Circuit-Compatible Basis Expansion (CCBE)** with Sparse Selection"
- Remove all KAN references from abstract and introduction
- Acknowledge in Related Work: "Unlike KAN which uses learnable splines, our approach uses fixed physical basis functions with learnable parameters"

**Option B: Justify "KAN-inspired"**
- Keep "KAN-inspired" but add extensive justification
- Claim: "KAN philosophy of interpretable building blocks" not "KAN architecture"
- Risk: Reviewer may still reject as misleading

**Recommendation**: Option A - Rename to avoid controversy

---

### Issue 3: n_restarts=5 is Insufficient

**Reviewer's Claim**:
> 5 random restarts using Adam is statistically guaranteed to fail finding global optimum.

**Response Options**:

**Option A: Increase n_restarts and Document**
- Change default to `n_restarts=20`
- Add convergence probability analysis in paper
- Show empirically that 95% convergence is achieved

**Option B: Add Better Initialization (Recommended)**
- Use logarithmic spacing for initial tau values (cover full frequency range)
- Add "warm start" from simplified VF poles
- This addresses the root cause, not just symptom

**Option C: Use Global Optimizer**
- Replace Adam with Differential Evolution for initial sweep
- Then fine-tune with Adam
- More expensive but more robust

**Recommendation**: Option B + increase n_restarts to 10

---

### Issue 4: Strawman VF Comparison

**Reviewer's Claim**:
> The benchmark results are scientifically invalid [using home-brewed VF].

**Already Addressed**:
- `benchmark_urn_vs_skrf_vf.py` uses scikit-rf VectorFitting
- `validate_real_world.py` includes VF comparison

**Additional Action**:
- Ensure all paper figures use scikit-rf VF results
- Add explicit statement: "VF comparison uses scikit-rf implementation of vectfit3"
- Show VF performs well on synthetic data (honest comparison)

---

## Proposed Paper Changes

### Title Change
**Before**: "KAN-inspired Universal Relaxation Network for Automatic Discovery..."
**After**: "Sparse Physical Basis Network for Automatic Discovery of Relaxation Mechanisms with SPICE-Compatible Circuit Synthesis"

### Abstract Changes
- Remove "KAN-inspired"
- Remove "attention-based"
- Emphasize: "sparse selection of circuit-compatible basis functions"
- Emphasize: "guaranteed passive synthesis"

### Section Changes

1. **Section 1 (Introduction)**
   - Remove KAN claims
   - Frame as "sparse basis expansion with physical constraints"
   - Clearly state: "Unlike neural networks, our basis functions have direct circuit equivalents"

2. **Section 3 (Methodology)**
   - Remove "Attention Gating" subsection entirely
   - Rename to "Sparse Basis Selection via L1 Regularization"
   - Honest description of what the algorithm actually does

3. **Section 6 (Ablation Study)**
   - Remove fake attention ablation
   - Keep only: (a) Learnable exponents, (b) Sparsity weight sensitivity, (c) n_restarts analysis
   - These are real, measurable effects

4. **Section 7 (Comparison with VF)**
   - Use only scikit-rf VF results
   - Show honest comparison where VF may win on some metrics
   - Emphasize URN advantages: passivity, interpretability, SPICE output

---

## Summary: Honest Paper Positioning

**What URN Actually Is**:
1. A parametric basis expansion with 29 circuit-compatible basis functions
2. L1 regularization for automatic sparse selection
3. Gradient-based optimization (Adam with multi-restart)
4. SPICE netlist generation from fitted parameters

**What URN Is NOT**:
1. Not a Kolmogorov-Arnold Network (no learnable activation functions)
2. Not attention-based (no softmax, no query-key-value)
3. Not guaranteed to find global optimum (local optimizer)

**URN's Actual Advantages**:
1. Passive by construction (positive RLC elements)
2. Interpretable (named physical mechanisms)
3. SPICE-compatible (direct circuit synthesis)
4. Handles broadband data (multi-decade frequency range)

---

## Next Steps

1. [ ] Revise title and abstract (remove KAN, remove attention)
2. [ ] Revise Section 3 methodology (honest description)
3. [ ] Revise Section 6 ablation (remove fake attention study)
4. [ ] Ensure all VF comparisons use scikit-rf
5. [ ] Increase n_restarts default to 10-20
6. [ ] Add convergence probability analysis
7. [ ] Recompile paper and resubmit
