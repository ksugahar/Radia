# Nonlinear ESIM-coupled hierarchical Cauer

Output data (.npz) from the ESIM Karl iteration pipeline for nonlinear
steel workpieces (Paper 2 application). Scripts live with the
production ESIM solver under `src/radia/` (`esim_cell_problem.py`,
`esim_coupled_solver.py`) — these examples store the verification /
demo numerical outputs.

| File | Origin |
|---|---|
| `esim_operating_points.npz` | ESIM cell-problem operating points (\|H_t\|, T) |
| `karl_iteration_data.npz` | Karl iteration convergence (Phase 3) |
| `full_nonlinear_lti_data.npz` | Full nonlinear LTI + Karl coupled (Phase 4) |
| `extension_3d_steel.npz` | 3D steel workpiece extension (Phase 5) |

The σ(T) perturbative CLN coupling (Phase 6) was abandoned 2026-05-24
because $\omega \tau_{\rm thermal} \gg 1$ at IH frequencies makes the
bilinear coupling physically meaningless; precomputed
$q_{\rm surf}(T, |H_t|)$ table lookup is the production path
(see memory `project_sigma_T_perturbative_cln_abandoned.md`).
