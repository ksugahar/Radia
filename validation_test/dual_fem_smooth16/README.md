# Dual A/Omega FEM on smooth periodic problems: reference-first 1:1 and 2:1 blends (2026-11-08)

- Each job ran on one logical CPU (group 0, chosen with `--logical-cpu`), at BelowNormal priority, inside a job object with a 16 GiB cap.
- Runtime: NGSolve 6.2.2608.
- Result JSON files are kept outside the repository.

## Setup
- Periodic structured meshes on the unit torus (2D/3D). There are no outer boundaries, corners, material interfaces or geometry errors.
- Every case is an analytic manufactured solution.
  - H\* = −∇φ\*, B\* = curl A\* (rot ψ\* in 2D), S = B\* + μ∇φ\*.
  - μ is a smooth, positive periodic function.
- Discretization:
  - Omega: periodic H1(p).
  - A: a periodic H1(p) stream function in 2D; in 3D, periodic HCurl(p, nograds) with an H1(1) multiplier and three harmonic constants.
  - Both fields lie in P_{p−1}.
  - Assembly uses 8 extra integration orders; post-processing reads out at two integration orders, 2p+12 and 2p+18.
- Q_A = C − 2J_A ≤ Q\* ≤ Q_Ω = C + 2J_Ω holds. This is a corrected bracket for any conforming discrete fields.
- The cases, ladders and blend test were declared in `cases.py` before any run. A later correction (R1) changed only the 3D reference level. It left the coarse level, the compared weights and the decision rule unchanged.

## Procedure (set by the user)
1. Converge A and Omega separately. Freeze the bracket [Q_A, Q_Ω]: its midpoint becomes the reference Q_ref and its half-width the uncertainty.
   - 2D: n=4, p=6. 3D: n=6, p=4 (R1).
   - The reference counts as converged when the change in the last refinement step (drift) is at most 10 % of |midpoint − 2:1| at the coarse level.
2. At a coarse, unconverged level (2D: n=4, p=2; 3D: n=3, p=2), compare the 1:1 (midpoint) and 2:1 blends with Q_ref.
   - When the two absolute errors differ by no more than the bracket width (twice the half-width), the comparison is inconclusive.
3. The analytic value Q\* is used only to check Q_ref independently.

## Results (coarse-level errors relative to Q_ref)

| Case | Reference half-width | Q_ref−Q\* | 1:1 error | 2:1 error | r (coarse) | Closer |
|---|---:|---:|---:|---:|---:|---|
| F1-2D | 3e-11 | −4e-16 | −1.1e-15 | −2.6e-3 | 1.000 | midpoint (exact by rotational symmetry; not evidence) |
| F2-2D | 2e-6 | −2e-6 | −1.41e-1 | −1.89e-1 | 0.011 | midpoint |
| F3-2D | 3e-11 | +1e-13 | +3.8e-5 | −2.6e-3 | 1.010 | midpoint (nearly symmetric) |
| F4-2D | 1e-8 | −8e-9 | −1.59e-2 | −2.45e-2 | 0.232 | midpoint |
| F5-2D | 6e-8 | +5e-8 | +8.78e-2 | +4.43e-2 | 5.12 | 2:1 |
| F6-2D | 4e-8 | +4e-8 | +5.06e-2 | +3.13e-2 | 14.7 | 2:1 |
| F7-2D | 3e-5 | +3e-5 | +2.24 | +1.47 | 61.4 | 2:1 |
| F8-2D | 2e-9 | −2e-9 | −7.99e-3 | −1.12e-2 | 0.100 | midpoint |
| F1-3D | 2e-6 | +2e-6 | +8.64e-2 | +5.30e-2 | 13.4 | 2:1 |
| F2-3D | 4e-5 | −1e-5 | −4.45e-2 | −1.30e-1 | 0.705 | midpoint |
| F3-3D | 1e-6 | −7e-8 | −2.76e-3 | −2.09e-2 | 0.903 | midpoint |
| F4-3D | 3e-5 | −3e-5 | −1.44e-1 | −1.95e-1 | 0.028 | midpoint |
| F5-3D | 5e-5 | +5e-5 | +4.77e-1 | +2.95e-1 | 15.0 | 2:1 |
| F6-3D | 1e-6 | −2e-7 | −1.59e-2 | −3.15e-2 | 0.494 | midpoint |
| F7-3D | 7e-2 | +7e-2 | — | — | 85 | **inconclusive** (the Omega reference has not converged; the last step changed it by 100 %) |
| F8-3D | 3e-5 | −3e-5 | −1.52e-1 | −2.05e-1 | 0.024 | midpoint |

- 15 of the 16 cases reached a verdict:
  - midpoint closer: 10 (5 in 2D, 5 in 3D);
  - 2:1 closer: 5 (3 in 2D, 2 in 3D);
  - inconclusive: 1 (F7-3D).
- In all 16 cases, the analytic value lies inside the reference bracket.
- Convergence rates: under h-refinement, the field-norm order is about 2 for both A and Omega in every case.
  - The one exception is F7-3D, where Omega has not reached the asymptotic range (1.50).
  - Cases with different p-decay constants (8 mixed and 1 different in the summary) reflect the different wavenumbers of the approximated functions.

## Interpretation and limits
- The winner depends only on r = E_Ω/E_A.
  - The midpoint error is E_A(r−1)/2 and the 2:1 error is E_A(r−2)/3, so the midpoint is closer for r < 1.4 and 2:1 for r > 1.4.
  - Every result here followed this rule.
- In these manufactured solutions, Omega approximates φ\* and A approximates ψ\* (A\*), two independently chosen functions linked only through the artificial source S. The value of r (0.01 to 85) and the winner therefore reflect the chosen functions, not the formulations.
- The study shows only that, on smooth periodic problems, A and Omega share the same h-convergence order. It does not establish that either blend is better in general.
- Comparing physically meaningful weights needs a setting where A and Omega approximate the same physical state (no free S, B = μH). One candidate is a homogenization cell problem: μ varies periodically under a uniform applied field, and tilted laminations have an analytic solution. This should be declared in advance as a separate category.

## Reproduce
```
python -I run.py --out <results dir> [--logical-cpu N]
python -I amend_r1.py --out <results dir> [--logical-cpu N]
python -I aggregate.py <results dir>
```
There is no MCP route yet; the study runs as a standalone benchmark.
