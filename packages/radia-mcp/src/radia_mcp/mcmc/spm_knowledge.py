"""Sampled Pattern Matching (SPM) method for EM design — Saotome 1995."""

from ..common import load_prompt

SPM_OVERVIEW = load_prompt("mcmc", "spm_overview")

_SPM_OVERVIEW_INLINE_DEPRECATED = r"""
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
"""


SPM_VS_LEAST_SQUARES = r"""
# SPM vs Least Squares: when does cosine similarity win?

## Decision: prefer SPM when ...

1. **Source magnitude is a free design choice**
   - Coil current `I` can be adjusted after the design is fixed
   - Magnet remanence `Br` will be selected from material catalog
   - The DESIGN sets the field SHAPE; magnitude is downstream

2. **Target is a pattern, not absolute values**
   - "Uniform B over the target surface" — shape matters, scale doesn't
   - "Linear gradient field for MRI" — shape matters
   - "Quadrupole field for accelerator focusing" — shape matters

3. **Sensor calibration is unknown**
   - Biomagnetic measurements: sensor gain may be unknown
   - NDT: probe coupling varies; cosine similarity is gain-invariant

## Prefer least-squares when ...

1. **Absolute calibration is meaningful**
   - "B = 1.5 T at the gap" — magnitude IS the spec
   - "P_wp = 5 kW heating power" — absolute energy delivery

2. **No source-magnitude degree of freedom**
   - Permanent magnet with fixed Br (can't tune)
   - Geometry fully determined; no current to scale

3. **Want unbiased moment estimator** (Bayesian MAP / ML)
   - Likelihood naturally Gaussian → least squares
   - SPM has no direct probabilistic interpretation

## Hybrid approach (lab recommendation)

1. **Stage 1: SPM** to find the right shape (gets pattern right
   quickly, magnitude-invariant)
2. **Stage 2: Least squares / NUTS** to calibrate magnitude
   (scale + offset to match absolute targets)

This mirrors the Sato 2023 MCTS pattern: MCTS for global structure,
NSGA-II for shape refinement.

## Mathematical relationship

For unit-norm x_t and x_e, the two are equivalent:

    ||x_t - x_e||^2 = ||x_t||^2 + ||x_e||^2 - 2 (x_t.x_e)
                    = 2 - 2 cos(x_t, x_e)
                    = 2 (1 - r)

So minimizing LS on normalized vectors = maximizing r. The difference
appears when |x_e| varies during optimization.

## Bibliographic note

The cosine similarity has many names:
- Saotome 1995 calls it "SPM rate" (sampled pattern matching rate)
- Information retrieval: "cosine similarity"
- Signal processing: "normalized cross-correlation"
- Statistics: "Pearson correlation" (after subtracting means)

All are the same formula `r = (a · b) / (||a|| ||b||)`. Saotome's
contribution is the application to EM design with FEM forward model
and best-element-on-contour search.
"""


SPM_LAB_APPLICATIONS = r"""
# Lab applications of SPM

The lab corpus has 3 papers using SPM (all in
`W:/.../04_機械学習と最適化/03_最適化_モンテカルロ/`):

## 1. Magnetic Core Shape Design (Saotome 1995)

Original paper. 2D magnetic core (DC excited) → uniform B_y on a
target surface. SPM achieves 2% error vs CLS 8% in same iterations.

## 2. Proton/Heavy Ion Therapy Dose Optimization

Mizukami 2010, "Dose optimization of proton and heavy ion therapy
using generalized sampled pattern matching". Generalizes SPM to
**non-magnetic** dose distributions.

Application: tune beam intensity profile such that delivered dose
matches the target tumor shape, while sparing surrounding tissue.

## 3. Radioactive Spent Resin Shape Estimation

Mizukami 2010, "Shape Estimation of Radioactive Spent Resins in a
Storage Tank Using the Sampled Pattern Matching". Inverse problem:
gamma-ray sensor measurements → estimate spatial distribution of
radioactive material inside a tank.

## Why these matter for EM

Lab specializes in EM, but the same SPM math applies to:
- **Coil shape for target B field** (synthesis of arbitrary B via
  shaped coil)
- **Magnetic shield design** (minimize external B leakage with min
  shielding mass)
- **WPT coil pattern** (target mutual inductance distribution)
- **NMR shim coil** (target field uniformity over imaging volume)

## Lab-extension recipe

For a new EM design problem:

1. Define **target field vector** `x_t` (e.g. B_y at 100 grid points
   over a target surface).
2. Define **forward model**: given material/coil mask M, compute the
   simulated field `x_e(M)` at the same 100 grid points.
3. Initialize `M_0` (e.g. uniform thin layer, or empty).
4. Iterate:
   - Find candidate moves (add/remove material element, swap, ...)
   - For each candidate, recompute `x_e(M_candidate)` and `r`
   - Accept best `r`
5. Stop at convergence or `r >= 0.99`.

For 2D problems with N ~ 1000 candidate moves, each forward solve is
seconds → full optimization is hours-to-days.

## Cross-links

- `radia_mcp.topology_optimization` — gradient-based shape opt
  (faster than SPM if gradient is available)
- `radia_mcp.optuna` — outer BBO over hyperparameters of SPM
  (e.g. step size, initialization)
- `radia_mcp.mcmc('bayesian_inverse')` — Bayesian alternative giving
  posterior over the design

## Lab paper list (full names)

In `W:/.../04_機械学習と最適化/03_最適化_モンテカルロ/`:
- `Magnetic Core Shape Design by the Sampled Pattern Matching Method.pdf`
- `Dose optimization of proton and heavy ion therapy using generalized sampled pattern matching.pdf`
- `Shape Estimation of Radioactive Spent Resins in a Storage Tank Using the Sampled Pattern Matching.pdf`
- `A Wavelet Transform Approach to Inverse Problems of Vandermonde Type Systems.pdf`
  (related: another inverse problem method using wavelets instead of SPM)
- `MAXIMUM-ENTROPY-BASED TOMOGRAPHIC RECONSTRUCTION.pdf`
  (related: MaxEnt for inverse — see `radia_mcp.mcmc('maximum_entropy')`)
"""


MAXIMUM_ENTROPY = r"""
# Maximum Entropy Method (MaxEnt) — sister of MCMC

MaxEnt finds the probability distribution that:
- Matches given moments (mean, variance, correlation with data)
- Maximizes Shannon entropy subject to those constraints

For inverse problems / image reconstruction, MaxEnt picks the
"least committal" image consistent with the data.

## Why MaxEnt belongs in the MCMC subpackage

MaxEnt and MCMC share the **same Boltzmann / Gibbs distribution** as
their target:

    p(x) ∝ exp(-beta * E(x))

MCMC samples from this distribution; MaxEnt finds the mode
analytically when E(x) is a known constraint structure.

Combined workflow:
1. MaxEnt for prior → reduces to a smooth low-dim subspace
2. MCMC for posterior → samples within that subspace

## Lab paper: "MAXIMUM-ENTROPY-BASED TOMOGRAPHIC RECONSTRUCTION"

In `04_機械学習と最適化/03_最適化_モンテカルロ/`. Application:
reconstruct internal density distribution from projection
measurements (CT-like). MaxEnt regularization prevents over-fitting.

For EM applications:
- **MEG/EEG source reconstruction**: many sources can produce same
  external field; MaxEnt picks the smoothest one
- **Magnetic permeability tomography**: estimate mu_r(x, y, z) from
  external coil measurements
- **Image deblurring after MRI**: MaxEnt removes ringing

## MaxEnt algorithm sketch

Solve:
    minimize -H(p) = -∫ p(x) log p(x) dx
    subject to E_p[f_i(x)] = mu_i for each known constraint

By Lagrange:
    p(x) = exp(-sum_i lambda_i * f_i(x)) / Z

The Lagrange multipliers `lambda_i` are found by Newton iteration
to match the constraints. Z is the partition function.

## Python implementation pointers

- `scipy.optimize` for the dual Lagrange problem
- `scikit-learn` has `MaxEnt` in `LogisticRegression` (different
  context: classification)
- For tomography specifically: `scikit-image` `radon` + custom
  MaxEnt regularization

## Cross-references

- `radia_mcp.mcmc('bayesian_inverse')` — alternative inverse method
  with posterior (more informative but more expensive)
- `radia_mcp.mcmc('spm')` — alternative deterministic inverse via
  cosine similarity
- `radia_mcp.mcmc('algorithms')` — Boltzmann/Gibbs distribution is
  the MCMC target as well
"""


def get_spm_documentation(topic: str = "all") -> str:
    """Dispatch SPM topics.

    Topics:
        overview            - Saotome 1995 SPM + cosine similarity
        spm_vs_ls           - When SPM beats least-squares
        lab_applications    - 3 lab papers: magnetic core, proton therapy, resin
        maximum_entropy     - Sister method MaxEnt + lab paper
        all                 - Everything
    """
    topic = topic.lower().strip()
    if topic in ("overview", "intro", "saotome", "spm_overview"):
        return SPM_OVERVIEW
    if topic in ("spm_vs_ls", "spm_vs_least_squares", "vs_ls",
                 "cosine_similarity", "comparison"):
        return SPM_VS_LEAST_SQUARES
    if topic in ("lab_applications", "applications", "lab"):
        return SPM_LAB_APPLICATIONS
    if topic in ("maximum_entropy", "maxent", "entropy"):
        return MAXIMUM_ENTROPY
    if topic == "all":
        return "\n\n".join([SPM_OVERVIEW, SPM_VS_LEAST_SQUARES,
                            SPM_LAB_APPLICATIONS, MAXIMUM_ENTROPY])
    return (f"Unknown topic '{topic}'. Available: overview, spm_vs_ls, "
            "lab_applications, maximum_entropy, all.")
