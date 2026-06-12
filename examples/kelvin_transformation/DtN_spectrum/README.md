# DtN-spectrum view of Kelvin-transformation coarse-mesh accuracy

These standalone scripts are the verified experiments behind the
`dtn_coarse_mesh` knowledge module
(`packages/radia-mcp/src/radia_mcp/radia_ngsolve/knowledge/dtn_coarse_mesh.py`).

**Idea.** Every open-boundary closure (Kelvin / BEM / PML / Robin) approximates the
one exterior Dirichlet-to-Neumann (Steklov–Poincaré) operator `Λ_ext`, whose sphere
eigenvalues are the mesh-independent ladder `λ_n = −(n+1)/R` (3D), `−n/R` (2D). The
discrete matrix `Λ_h` reproduces the LOW-degree eigenvalues accurately and almost
independently of mesh size `h`; the per-mode defect grows with degree `n`. Because a
compact source's field is dominated by low multipoles, a coarse Kelvin mesh already
resolves everything that matters — the empirical coarse-mesh accuracy, read off the
spectrum.

Every script PRINTS its result (no files written) and depends only on
`numpy` / `ngsolve` / `netgen` and the open `radia_mcp.radia_ngsolve` helpers
(`bem_integral`, `fem_bem_coupling`).

| Script | Question | Knowledge topic | Verified result |
|---|---|---|---|
| `sufficient_mesh.py` | How fine must the truncation mesh be? | `numerics` | sufficient-mesh criterion `N_surf ≳ n_src/√ε`; accurate band ≈ 0.12–0.18 of angular Nyquist; per-degree defect `~ n²(h/R)⁴` |
| `p_vs_h_study.py` | Is Kelvin an h-method or a p-method? | `p_method` | order `p ≥ n` kills the polynomial error (then a curved-geometry floor: ~5–6 digits in 3D = Kameari's result, deeper in 2D); raising `p` beats refining `h` by ~20–80× DOF |
| `poly_vs_sphere.py` | Does a faceted (polyhedron) truncation hurt? | `p_method` | dipole robust flat-or-curved (2–13×); quadrupole faceting error ~369× → `p` and geometry-order must rise together |
| `demo1_hp_lshape.py` | h vs p when the interior has a corner singularity | `p_method` | L-shape `r^{2/3}` corner: measured `α_h = 0.357`, `α_p = 0.661` (`α_p/α_h = 1.85`) → regularity decides; Kelvin region is the p-favourable analytic part |
| `demo2_dual_bracket.py` | Can we CERTIFY the open-BC error, not just estimate it? | `formulation` | complementary (Prager–Synge) bracket `E_primal ≤ E(f_h) ≤ E_comp`, gap `O(h²)`, equilibration residual ~1e-15 |
| `demo3_A_dtn_gradient.py` | Is the DtN gradient block formulation-dependent (Ω vs A)? | `formulation` | the gradient block `−(n+1)/R` is the SAME for the vector-potential A formulation (dipole `−2/R`; `B·n = curl_Γ A_t`, rel-L2 3.2e-5) |
| `demo_d_multipole_spectrum.py` | Source factor `c_n` for a magnetised square (Kameari Q-d) | `datasheet` | edge-charge spectrum has only `n ≡ 1 (mod 4)`; `n=3` FORBIDDEN by symmetry (single-resolution false-positive trap); leading correction `a₅/a₁ = (4/15)(a/R)⁴` |
| `demo_e_optimal_R.py` | Optimal Kelvin/truncation radius (Kameari Q-e) | `p_method` | disk/sphere = monotone → smallest R; square/cuboid → interior optimum `R/a ≈ 3` (DOF proxy `(R/a)²·p(R)²`) |
| `kelvin_exterior_mesh.py` | Does refining the exterior VOLUME mesh help? | `p_method` | `‖u_h−P_n‖≈1.5e-15` on every volume mesh (p≥n): the volume solve is Galerkin-exact → refining the exterior interior does nothing; only Γ's geometry (a surface effect) moves λ |
| `kelvin_exterior_mesh3.py` | Isolate Γ-surface vs interior volume | `p_method` | with the Γ surface fixed, λ is set by `∫_Ω|∇P_n|²/∮_Γ P_n²` — a fixed polynomial over a fixed domain; the exterior volume mesh enters nowhere. p=1<n control shows refinement only helps when the order is deficient |
| `floor_vs_curve.py` | Is the 5–6 digit floor really *geometry*? | `p_method` | fix FE order `p≥n` + mesh, raise only the isoparametric Curve order `k`: error drops ~1000× (`k=1` flat ~1% → `k=3` ~1e-5) → the floor is the curved-sphere **geometry**, not the multipole/method |

## Running

```bash
python demo_d_multipole_spectrum.py     # pure numpy, instant
python poly_vs_sphere.py                # needs ngsolve + netgen.occ
```

The `radia_mcp.radia_ngsolve` import path resolves when `radia-mcp` is installed
(`pip install -e packages/radia-mcp`).

## Manuscript figures

The figure-generation scripts for the companion write-up (`gen_fig*.py`) are kept
out of the repo as venue-specific plotting plumbing; the eight scripts above are the
reusable, self-checking experiments. Ask if the figure generators should be added.
