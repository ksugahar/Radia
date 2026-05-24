# Sampled Pattern Matching (SPM) — Saotome 1995

Reference: H. Saotome, J.-L. Coulomb, Y. Saito, J.-C. Sabonnadiere,
"Magnetic Core Shape Design by the Sampled Pattern Matching Method",
IEEE Trans. Magnetics 31(3):1976-1979, May 1995.
DOI: 10.1109/20.376431.

Lab copy: W:/.../04_機械学習と最適化/03_最適化_モンテカルロ/
Magnetic Core Shape Design by the Sampled Pattern Matching Method.pdf

Affiliations: Chiba Univ + INPG CNRS Grenoble + Hosei Univ.

## The new objective function

Classical least squares for inverse problems (electromagnetic device
design from a target field):

    E_LS = ||x_t - x_e||^2                 (least-squares error)

where `x_t` is the target field vector and `x_e` is the evaluated
(simulated) field vector.

SPM replaces it with the **cosine similarity**:

    r = (x_t . x_e) / (||x_t|| ||x_e||)    (in [-1, 1], maximize)

## Why cosine similarity is better

The least-squares `E_LS` depends on:
- The shape/pattern of `x_e` vs `x_t`
- AND on the magnitude of `x_e`

The cosine `r`:
- Depends ONLY on the **shape/pattern** — magnitude is normalized
  out
- Equivalent to minimizing the NORMALIZED error:

      E_r = || x_t/||x_t|| - x_e/||x_e|| ||^2

In design problems where the field source can be re-scaled (DC coil
current, magnet remanence, etc.), the SHAPE is what matters; the
magnitude is set by the source. Cosine similarity directly targets
shape.

## Important properties

- **Source-magnitude invariance**: If `x_e = K * x_b` for any scalar
  K, then `r` is unchanged. The user does not need to know the
  required current density.
- **Bounded in [-1, 1]**: r=1 is perfect match; r=-1 is anti-aligned.
- **r=1 implies perfect pattern**, BUT magnitude must be set
  separately (multiply x_b by a scalar to match a calibration point).

## Algorithm: best-element selection on contour

For a 2D magnetic core shape design:

1. Mesh the design region with N small triangular elements.
2. Initialize: an "initial shape" of magnetic material (e.g. a few
   accumulated elements).
3. **For each candidate element on the current contour**:
   a. Add the element to the magnetic material (set mu_r = 500).
   b. Solve the FEM forward problem → get field x_e at target.
   c. Compute `r(candidate) = cos(x_t, x_e)`.
4. Select the candidate that maximizes `r`. Add to material; new
   contour.
5. Repeat until target tolerance, or all candidates rejected.

## Comparison: CLS vs SPM trajectory

For a typical target (uniform B_y on a target surface):
- CLS path: zigzags around the optimum because magnitude is mis-set
- SPM path: monotonically improves pattern, magnitude scales out

Saotome shows: CLS error reduces 18% → 8%; SPM error reduces 18% →
2% on the same example. ★ 4x better.

## Where SPM applies

| Domain | Reference |
|--------|-----------|
| Magnetic core shape design | Saotome 1995 (this paper) |
| Biomagnetic source identification | refs [5-7] in Saotome 1995 |
| Non-destructive testing (NDT) | ref [8] |
| Proton/heavy ion therapy dose | Mizukami 2010 (lab copy) |
| Radioactive resin shape estimation | Mizukami 2010 (lab copy) |
| EM tomography | adaptable |

## Why MCMC-shaped (and why it lives in radia_mcp.mcmc)

The "sample" in "sampled pattern matching" is the sampled vector at
mesh nodes / sensors → SPM is a forward sampling method (not MCMC in
the strict sense). It is here because:

1. **Same lab folder** as MCMC papers (`03_最適化_モンテカルロ/`)
   — the lab groups SPM with Monte Carlo methods.
2. **Same problem class**: inverse problem for EM design from
   sampled measurements.
3. **Complementary to MCMC inverse**: SPM gives a deterministic
   point estimate; MCMC gives a posterior. Use SPM as initialization
   for MCMC.

## Python prototype

```python
import numpy as np
from radia.panels.calc_inductance import run_fem_field

def spm_objective(x_t, x_e):
    # Cosine similarity (Saotome SPM rate r).
    return np.dot(x_t, x_e) / (np.linalg.norm(x_t) * np.linalg.norm(x_e))

def spm_design_iteration(mesh, current_material_mask, target_field,
                          target_locations):
    # One SPM iteration: pick best contour element to add.
    contour_elements = find_contour(mesh, current_material_mask)
    best_r, best_elem = -np.inf, None
    for elem in contour_elements:
        trial_mask = current_material_mask.copy()
        trial_mask[elem] = True
        x_e = run_fem_field(mesh, mu_r=500*trial_mask,
                             obs=target_locations)
        r = spm_objective(target_field, x_e)
        if r > best_r:
            best_r, best_elem = r, elem
    return best_elem, best_r

def spm_design(mesh, initial_mask, target_field, target_locations,
                max_iter=100, r_target=0.99):
    mask = initial_mask.copy()
    for it in range(max_iter):
        elem, r = spm_design_iteration(mesh, mask, target_field,
                                        target_locations)
        if elem is None or r < spm_objective(target_field,
                                              run_fem_field(...)):
            break
        mask[elem] = True
        if r >= r_target:
            break
    return mask, r
```

This is a greedy add-element search. Modern variants:
- Add-or-remove: at each step, consider also removing an element
- Multi-element batches: add the best K simultaneously
- Combined with topology gradient (`radia_mcp.topology_optimization`)
