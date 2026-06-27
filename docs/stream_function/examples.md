# Demo ledger — what each example does

All under [`examples/stream_function/`](../../examples/stream_function/).

## Core demos (axisymmetric / fingerprint / planar)

| File | Topology | Best RMS | Notes |
|------|----------|----------|-------|
| `demo_coil_design_gz.py` | Cylinder, Gz gradient (1D) | 0.5 % | Original SA-25-020 demo |
| `demo_coil_design_gx.py` | Cylinder, Gx fingerprint (2D) | 0.8 % | Continuous SF solve, no chain |
| `demo_sf_to_peec_gz.py` | Cylinder, Gz + single-stroke + PEEC + CAD | — | Full pipeline incl. STEP export |
| `demo_sf_to_peec_gx.py` | Cylinder, Gx + chain methods + PEEC + CAD + sheet-metal distortion | 9.3 % (field_aware) → **1.4 %** (`--distort`, 1 current) | `--chain-method {field_aware, kuijpers, lobe, greedy}`; `--distort` = single-current radial bend; `--regularize {tsvd, tikhonov, h1}` |
| `demo_planar_uniform_coil.py` | Plane source, uniform Bz target | 0.58 % (+ Path-A) | Basis-loop pipeline |
| `demo_planar_uniform_fem_psi.py` | Plane, FE-direct H¹ ψ + Path-A + LS-OMP shim + sheet-metal distortion | **183 ppm** (10 feeds) / **605 ppm** (1 current, `--distort`) | `--order 3 --nlevels 30 --shim-tol-ppm 200`; or `--distort` = single-current 3D wire bend (no extra feeds); honest `--eval-n` grid |
| `demo_planar_uniform_fem_psi_aca.py` | Plane, FE matrix + HACApK ACA+TSVD | 0.67 % (+ Path-A) | Validates ACA+ on FE matrix |
| `demo_sphere_fe_direct.py` | **Sphere former** (ANY curved surface), FE-direct H¹ ψ + single-stroke + sheet-metal distort | uniform **0.24 %** / Z2 shim 4.3 % → **0.36 %** (`--distort`, 1 current) | `--target {uniform, z2}`; the case basis-loop CANNOT grid (no structured (φ,z) grid on a sphere); NMR shim on a curved former; manufacturability spacing gate |
| `demo_planar_uniform_fem_psi_advanced.py` | Plane, 6 regs + deformation + Optuna (RMS / constrained / Pareto) | 0.58 % NSGA-II Pareto accuracy end | `--regularize {...}`, `--deform`, `--minimize-reg --eps-rms ε`, `--pareto` |
| `demo_regularized_aca.py` | Plane, 5-mode sweep through cached `RegularizedTSVD` | 1.12 % (linf IRLS) | Single ACA+ factorisation reused across all 5 regularisations |
| `demo_reg_hyperparam_aca.py` | Plane, Optuna over σ(x,y) shape; ACA+ base reused across trials | **0.73 %** (vs 2.09 % uniform σ) | Fixed surface → A constant → ACA+TSVD computed ONCE, only fold rebuilt per trial |
| `demo_pareto_tikhonov_aca.py` | Plane, **(homogeneity, peak-J) Pareto front** via folded Tikhonov + ACA+TSVD | front, not a single RMS | `--front {energy, peak, both}`; ONE ACA factorisation, α-sweep ≈ **50 µs/Pareto-point**; L∞-IRLS seminorm pushes the peak front **−18 %** (median) vs the H¹ L-curve at matched homogeneity |
| `demo_pareto_geometry_nsga.py` | Plane, **push the front with the GEOMETRY lever + NSGA-II** | front | `--mode {geometry, nsga, both}`; former-size envelope (bigger former → **−34 %** exact-homog peak, diminishing returns); NSGA-II over (former, α) traces the 3-objective (homogeneity, peak, **former size**) surface — low peak ⟺ large former |
| `demo_pareto_cylinder.py` | **Cylinder Gx fingerprint** (MRI/shim), folded Tikhonov + ACA front | front | geometry lever = cylinder **length**; reuses the cylinder basis-loop kernel from `demo_sf_to_peec_gx.py`; length has an **OPTIMUM** (≈ 50 cm, **−37 %**) — longer stops helping once it covers the DSV (unlike the monotone planar former) |
| `demo_pareto_deform.py` | **Planar sheet-metal (板金) deformation** pushes the front | front | FORM the surface `z=f(x,y)` (folded-Tikhonov inner solve), CMA-ES optimise the shape per homogeneity level; `--zero-mean` (default) = genuine bending at FIXED standoff (**−17 %** exact-homog, whole front −5…−18 %); `--allow-standoff` ≈ −53 % but standoff-dominated |
| `demo_pareto_cylinder_deform.py` | **Cylinder in-surface sheet-metal (板金)** — the dominant cylinder lever | front | length-preserving axial reparametrisation (radius FIXED → 100 % genuine, no standoff); spacing-weighted seminorm + non-uniform-spacing peak; whole Gx-fingerprint front **−10…−25 %**. Radial forming is WEAK (~−3 %): the lever **flips** vs the plane |
| `demo_cmaes_magnet_design.py` | 16 magnetisation angles | — | SA-25-020 CMA-ES outer loop |
| `demo_magnet_array.py` | MMM/MSC magnet array | — | Validates kernel-agnostic for magnetic materials |
| `demo_coil_field_synthesis.py` | Simple coil field synthesis | — | First-principles demo |

## Visualisation

| File | Purpose |
|------|---------|
| `view_sf_coil_gx_gmsh.py` | Open Gx coil in GMSH (3 modes: `contours`, `chain`, `step`).  Has off-screen window prevention via explicit `General.GraphicsPositionX/Y` + `Width/Height`. |
| `bench_aca_vs_dense.py` | (ACA+)+TSVD timing benchmark vs full SVD |

## Benchmarks

Under [`examples/stream_function/benchmarks/`](../../examples/stream_function/benchmarks/):

| File | Status | Reference |
|------|--------|-----------|
| `bench_helmholtz_pair.py` | shipped | Analytical Maxwell pair (uniform Bz) |

Future benchmark targets are tracked in [benchmarks.md](benchmarks.md), not
as non-runnable `.py` stubs under `examples`.

See [benchmarks.md](benchmarks.md) for the validation strategy.

## Recommended order of demo runs (= learning path)

1. **`demo_coil_design_gz.py`** — see the SFM math in 1D.  Verify the
   wire pattern matches a saddle helix.
2. **`demo_sf_to_peec_gz.py --with-peec`** — full pipeline: SF design →
   single-stroke helix → CAD STEP → PEEC inductance → Bz verification.
3. **`demo_planar_uniform_coil.py --compensated-iter 30`** — planar
   case + Path-A iteration in the easy mode (basis-loop).  See the
   best-effort tracking find ~0.6 % RMS.
4. **`demo_planar_uniform_fem_psi.py --regularize h1 --compensated-iter 100 --compensated-step 0.05`**
   — FE-direct ψ + H¹ min-seminorm + Path-A.  Watch iter 40-47 drop
   monotonically 0.62 % → 0.49 %.  This is the proof Path-A *can*
   converge.
5. **`demo_planar_uniform_fem_psi_advanced.py --regularize h1 --order 3 --deform --deform-params bump --deform-trials 20`**
   — combine high-order FE + surface deformation outer loop.
6. **`demo_regularized_aca.py`** — see all 5 regularisations (L², H¹,
   σ-weighted H¹, inductance-diagonal, IRLS L∞) routed through the
   same cached `RegularizedTSVD` factorisation.  Per-mode `solve(B)`
   cost is sub-ms after the one-shot ~40 ms `S⁻¹V` + `W_inv` fold.
   See [regularization.md](regularization.md) for the closed-form
   derivation `ψ = S⁻¹V · W⁻¹ · Σ⁻¹ · UᵀB`.
7. **`demo_sf_to_peec_gx.py --chain-method kuijpers`** — hard tier:
   cylindrical fingerprint coil + 3 chain methods.  See why Path-A
   doesn't help here (tier-bounded).
8. **`view_sf_coil_gx_gmsh.py --mode contours`** vs `--mode chain` —
   visualise the SF design vs the manufactured single-stroke chain.
9. **`benchmarks/bench_helmholtz_pair.py`** — see how our planar SF
   compares to the analytical Maxwell pair baseline.

## Cross-reference

  - API per-function: [api.md](api.md)
  - Theory: [theory.md](theory.md)
  - Reproducing published benchmarks: [benchmarks.md](benchmarks.md)
