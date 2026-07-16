# NGSolve-style FEM basis functions — Mathematica symbolic reference

General-purpose, symbolic (Mathematica/wls) implementation of the
**NGSolve high-order hierarchical FEM shape functions** (H1 / H(curl) /
H(div)), ported from NGSolve's C++ element library.  Kept here (not in an
application directory such as `examples/CLN/`) because the basis is a
**reusable building block**: the `mathematica_*` MCP tools and any future
work that needs the *analytical* form of a shape function — e.g. the
analytical solid-angle / `1/r` integral for an HDiv-VIM integral-equation
kernel built on the FEEC basis — draws on it here.

## Why this matters (the loopless HDiv-VIM solver)

The discrete de Rham structure of these bases is exactly what removes the
"loop" problem from a surface-charge / integral-equation magnetostatic or
eddy-current solver:

- `curl(H(curl))` is **div-free to machine precision** (verified on
  distorted hexes, NGSolve, order 1–3).  So the div-free magnetization /
  current circulations ("loops", `ker N`) are **machine-zero field-null by
  construction** — no hand-built null vectors, no over-counting, no
  per-element ad-hoc element engineering (the thing that hand-crafted
  "distorted-element" formulations were doing to suppress the spurious
  loop-star component).
- **High order is one line** (`order = p`), so the low-order
  Raviart–Thomas / Nédélec accuracy loss on distorted hexes (Naff–Russell–
  Wilson) is removed by raising `p`, not by engineering the element.

To build the integral-equation **field operator** on this basis, the only
new ingredient is the field of each shape function: the solid-angle / `1/r`
integral of its (rational, Piola-mapped) density over the volume/face
element.  Two routes:

1. **Analytical** — integrate the *symbolic* shape function (this package)
   against the `1/r` kernel in Mathematica.  Where it closes, it is exact
   and fast.
2. **Gaussian quadrature** — needs no symbolic basis at all: evaluate the
   NGSolve shape function at quadrature points (NGSolve Python API) and sum
   against `1/r`.  Always available; the safe fallback.

(The CLN/Foster work used **arbitrary-precision adaptive integration**,
above machine double, when converting a transfer function to Foster form —
the same `NIntegrate[..., WorkingPrecision -> n, Method -> "GlobalAdaptive"]`
discipline applies when the analytical route is taken numerically.)

## Construction (Zaglmayr hierarchical)

```
l[i, x]        = LegendreP[i, x]                       (* Legendre        *)
L[i, x]        = Integrate[l[i-1, z], {z, -1, x}]      (* integrated Legendre *)
lS[n, s, t]    = t^n  l[n, s/t]                        (* scaled Legendre  *)
LS[n, s, t]    = t^n  L[n, s/t]                        (* scaled int. Leg. *)
```
then barycentric coordinates `lamda[i]` per element, and the
vertex / edge / face / interior shape functions assembled so the
H(curl) gradient block (`ker curl`) is explicit at the shape-function level
(`nograds` drops it).  Ref: Zaglmayr PhD (JKU 2006) §5–6; Schöberl &
Zaglmayr, COMPEL 24 (2005).  See `../notes_fem_hcurl.md`.

## Coverage (current `.wls`)

| space | high-order (any p) — Hex/Quad | high-order (any p) — Tet/Trig | rational |
|-------|-------------------------------|-------------------------------|----------|
| H1      | Segm, Trig, Tet, Quad, Hex, Prism (hierarchical) | (= H1 columns)   | **Pyramid** (p=1) |
| H(curl) | **Hex, Quad** (tensor de Rham) | **Tet, Trig** (classical Nédélec `N_k`) | Pyramid TODO |
| H(div)  | **Hex, Quad** (shared tensor de Rham; explicit RT/BDM aliases) | **Tet, Trig** (distinct `BDM_p` and Raviart–Thomas `RT_p`) | Pyramid TODO |

- **Hex / Quad** vector spaces: tensor de Rham construction (`Q`=Legendre/L2,
  `W`=integrated-Legendre/H1, `d: W→Q`) — spans NGSolve's order-`p`
  `H(curl)`/`H(div)` exactly; clean explicit polynomials for the VIM `1/r`
  integral.  NGSolve's `RT` selector does not change the Hex/Quad local space,
  so `HDivHexRT[p]`/`HDivHexBDM[p]` and
  `HDivQuadRT[p]`/`HDivQuadBDM[p]` are intentional aliases.
- **Tet / Trig** vector spaces (`simplex_ho.wls`): the classical
  `BDM_p = [P_p]^d`, `RT_p = [P_p]^d ⊕ x̃·P̃_p`, and
  `N_k = [P_k]^d ⊕ {p∈[P̃_{k+1}]^d : p·x̃=0}`.
  NGSolve uses a Zaglmayr integrated-Jacobi *hierarchical basis* of these same
  spaces; the complete-polynomial form is basis-equivalent and what the VIM
  field operator needs (it integrates whatever spans the space).
- **RT/BDM study invariant**: on simplices, `div(BDM_p)=P_(p-1)` and
  `div(RT_p)=P_p`, while their `ker(div)` dimensions are equal.  RT's extra
  equal-order DoFs are therefore charge-carrying modes, not additional loop
  modes.  `HDivTetFamilyLedger[p]` and `HDivTrigFamilyLedger[p]` return the
  symbolic DoF split used by the HDiv-MMM study.

The local dimension contract shared with NGSolve is:

| element | BDM p=1/2/3 | RT p=1/2/3 |
|---------|-------------|------------|
| Trig | 6 / 12 / 20 | 8 / 15 / 24 |
| Tet | 12 / 30 / 60 | 15 / 36 / 70 |
| Quad | 12 / 24 / 40 | 12 / 24 / 40 |
| Hex | 36 / 108 / 240 | 36 / 108 / 240 |

`tests/test_hdiv_family_dimensions.py` locks this family/order mapping against
the installed NGSolve API; the `.wls` self-tests independently prove the
polynomial ranks.
- **Pyramid `p≥2`** (rational edge/face bubbles) is the remaining gap.

## Source of truth

- **NGSolve C++**: `fem/recursive_pol.hpp` (Legendre/Jacobi recursions),
  `fem/h1hofe*`, `fem/hcurlhofe*`, `fem/hdivhofe*` (per-element high-order
  shape functions).  Re-derivable from these directly.
- **Original working notebooks (LAB-local, not shipped)**:
  `public-safe curated corpus functions\*.nb`
  (`H1 Shape functions.nb`, `Hcurl Shape Functions.nb`,
  `High Order Nédélec Elements 3D.nb`, ...; 2026-02 vintage).  These are
  large binary `.nb` and are **not** committed here (PyPI size + diff
  hygiene); the clean `.wls` in this directory are the version-controlled,
  shippable form.

## NGSolve C++ source (the porting reference)

The authoritative per-element shape-function source is NGSolve's
`fem/h1hofe_*.cpp`, `fem/hcurlhofe_*.cpp`, `fem/hdivhofe_impl.hpp`,
`fem/recursive_pol.hpp` — covering **all** element types incl.
`h1hofe_pyramid.cpp`, `h1hofe_prism.cpp`, `hcurlhofe_pyramid.cpp`,
`hcurlhofe_prism.cpp` (the hard rational pyramid/prism cases).  A local
copy is kept at **`C:/temp/ngsolve_fem_src/`** (recovered LAB-local; NOT
committed — it is NGSolve's own GPL source, kept only as the porting
reference).  Each clean `.wls` here is authored from these, then
self-tested.

## Status / plan

| file | element(s) / content | status |
|------|----------------------|--------|
| `RadiaBasis.m`      | Phase-1 Mathematica package for triangle/tetrahedron H1 Lagrange, RT0/RWG, L2; canonical source mirrored by `tests/basis/test_basis_functions.py` | done, promoted from the old examples tree so MCP knowledge owns the symbolic reference |
| `recursive_pol.wls` | Legendre / integrated-Legendre / scaled | done, self-test PASS |
| `h1.wls`            | H1 Trig, Tet, Quad, Hex, Prism          | done, self-test PASS (dim, PoU, edge-vanish, independence) |
| `h1.wls`            | H1 **Pyramid** (rational, p=1 vertices) | done, self-test PASS (PoU, rationality, independence); p>=2 edge/face bubbles TODO |
| `derham.wls`        | Whitney de Rham complex on tet (W0/W1/W2/W3) | done, self-test PASS (`d∘d=0`, `div(curl W1)=0` = loops div-free, Whitney dof) — the cohomology / loop foundation |
| `hcurl.wls`         | H(curl) Nedelec: **Hex/Quad tensor (high-order p)** + tet/trig Whitney W1 + grad-block (ker curl) + SZ helpers | done, self-test PASS (dims, `grad H1 c H(curl)` curl-free, `curl H(curl) c H(div)` div-free, Whitney tangential dof) |
| `hdiv.wls`          | H(div): **Hex/Quad shared tensor space** with explicit RT/BDM aliases + tet/trig Whitney RT0 + SZ covariant helpers | done, self-test PASS (dims, RT/BDM tensor alias equality, `div(H(div)) = L2`, RT0 face dof, **de Rham `curl H(curl) c H(div)` showcase**) |
| `simplex_ho.wls`    | high-order **Tet/Trig** H(curl)/H(div): distinct BDM `BDM_p`, Raviart–Thomas `RT_p`, and Nédélec `N_k` | done, self-test PASS (dimensions, independence, `BDM_p c RT_p`, normal-trace rank, `div BDM_p=P_(p-1)`, `div RT_p=P_p`, equal `ker(div)`, `curl N_p c BDM_p c RT_p`) |
| `cohomology.wls`    | discrete de Rham complex (d0/d1 incidence), Betti b0/b1/b2, **harmonic H^1 generator** (tree-cotree, `b1=E-V+C`), **metric-independence of Betti** (harmonic-1-form dim = `b1` under ANY SPD metric M1 — the Hodge/material can't change topology) | done, self-test PASS (filled-disk/cycle/annulus/figure-eight; generator is a cycle but not a gradient; harmonic dim = b1 under random SPD metrics) — the GLOBAL loops, and *topology is not in the Hodge* |
| `evrs_derham_bridge.wls` | **EVRS / EID T-method bridge**: `T in HCurl`, `J = curl T`; arbitrary retained EVRS combinations map to HDiv-compatible solenoidal currents; the discrete chain `phi --G--> T --C--> J --D--> rho` satisfies `C G = 0`, `D C = 0`; flat surface-Omega has `K = n x grad_Gamma Omega`, `div_Gamma K = 0` | done, self-test PASS (`curl grad = 0`, `div curl(T)=0`, `curl(W1) c W2`, `C G = 0`, `D C = 0`, `R_T=C^T M_R C` gauge-null/symmetric, EVRS Galerkin identity, port RHS gauge-invariant, surface-Omega tangential/solenoidal) |
| `eddy_topology_reduction.wls` | **topology-aware HCurl-VIM eddy-bubble reduction**: face role is decided by the two neighboring material labels; conductor-air/exterior faces become SIBC/surface-Omega trace class; conductor-insulator faces become non-SIBC boundary traces; conductor-conductor faces become bridge class; the bridge class is reduced by conductive-graph cycle coordinates `j_bridge = Z gamma`, with `B Z = 0` | done, self-test PASS (face role classification, SIBC-only exterior faces, conductor-insulator is not SIBC, air-air ignored, `components=V-rank(B)`, `cycle_rank=E-rank(B)`, bridge-cycle currents have zero cell divergence, surface-Omega tangential/solenoidal, parent order `p` is an admissibility condition rather than an empirical constant) |
| `eddy_parent_order_admissibility.wls` | **parent HCurl order ledger**: `p` is selected by spatial requirements, not by the CLN/Krylov rank; `p >= max(r_bulk, r_bridge, r_surface)`, surface-Omega scalar degree is `r_surface+1`, and enriched bridge traces have `cycle_rank * dim P_r(face)` modes | done, self-test PASS (`dim P_r(face)=Binomial[r+2,2]`, constant bridge trace gives one mode per cycle, quadratic trace multiplies by 6, p=6 can be admissible/conservative without being optimal, high-order bridge trace can force p=6 symbolically, current coarse count `24+7+3=34`) |
| `vim_field.wls`     | **VIM field operator**: charge extraction (sigma=M.n, rho=-div M), 1/r field, van Oosterom-Strackee solid angle (solid-angle face-charge kernel) | done, self-test PASS (div-free M -> rho=0 = loop field-null; VIM -> dipole far-field O(1/R^2); VOS solid angle = quadrature) |
| `vim_loopfree.wls`  | **THE loop-mode question**: is the FEEC VIM formulated so no spurious loop modes arise? | done, self-test PASS -- loops = curl(interior H(curl)) (+) cohomology are CHARGE-FREE (div=0 AND M.n=0) => field-null **by construction on any element** (Piola preserves div + normal-trace); constant-M misses them (avg=0) -> its tree-cotree loops only approx ker(N) on distorted hexes = the defect |
| `infinite_element_derham.wls` | **de Rham (exact-sequence) INFINITE ELEMENT** on the spherical EXTERIOR (open boundary) -- Mathematica twin of the maintained IE notes in `docs/open_boundary/INFINITE_ELEMENT_SOTA.md` | done, self-test PASS -- via Mathematica's built-in ORTHONORMAL spherical `Grad`/`Curl`/`Div`: the radial decay families shift **+1 per form degree** (S0={n+1..n+P}, S1=S0+1, S2=S0+2, S3=S0+3) so grad/curl/div COMMUTE (grad(V0) subset V1, curl(V1) subset V2, div(V2) subset V3 with explicit structure constants; curl.grad=0, div.curl=0; toroidal/div closure via the Legendre eig -n(n+1), shown m-independent). Demkowicz-Pal (CMAME 164, 1998), STATIC/low-freq; the 0-form tower is the scalar Bettess IE. NOTE the shipped de Rham open boundary is instead the coordinate-mapping family (Kelvin / coordinate-scaling IE), de Rham inherited free; high-freq needs an oscillatory exp(ikr) basis (Astley/Demkowicz-Pal radiating) |

**The de Rham complex `H1 →grad→ H(curl) →curl→ H(div) →div→ L2` is now verified
symbolically** (both maps exact: `curl∘grad=0`, `div∘curl=0`, each image lands in
the next space), and the **cohomology `H^1`** (global loops on multiply-connected
bodies) is the tree-cotree cycle count. Loop-free magnetic-material solving is HDiv-VIM's domain.

**What is settled at the `.wls` level** (the right gate BEFORE any C++): the FEEC VIM is
**loop-mode-free by construction** (`vim_loopfree.wls`) — the field-null loops are exactly
the charge-free space `curl(interior H(curl)) ⊕ cohomology H^1`, and the Piola map keeps
them charge-free (div=0 AND M·n=0) on *any* distorted element, so no spurious loop modes
arise and no per-geometry numerical null-vector patch is needed.  This is the formulation
question that had to be answered symbolically first.  The full de Rham complex + cohomology
+ VIM field operator are covered for the FEEC element families Radia needs
(tet/hex/wedge plus quad/trig/prism) — 9 files, 100+ self-test assertions, all PASS —
now including the **EXTERIOR open-boundary** de Rham element (`infinite_element_derham.wls`),
the FEEC counterpart on the unbounded side (the same `H1→H(curl)→H(div)→L2` exactness, with
decay families that commute under grad/curl/div).

Still to settle at the `.wls`/formulation level before C++:
- the loop/star (Hodge) SPLIT of the assembled high-order system — confirm the star
  (charge-carrying) block is well-conditioned once the loops are removed;
- the VIM right-hand side / collocation correspondence (the evaluation points must match
  the high-order basis, per the lab's "Yano element" experience);
- a small end-to-end VIM solve on a distorted multi-element patch, loops removed, vs a
  trusted reference — entirely in `.wls`/Python before committing to a C++ kernel.

Deferred (NOT on the solver's critical path):
- **pyramid** p>=2 (rational edge/face bubbles): Radia has **no pyramid element**
  (legacy face-charge hex/wedge and tetra moment paths), so this is a completeness item only.  The space is the
  rational collapsed-coordinate bubbles (`xt=x/(1-z)`, scaled integrated-Legendre +
  triangle bubbles); NGSolve's `EdgeOrthoPol`/`TrigOrthoPol` are dual-shape
  optimizations, not needed for the primal space.
- (optional) the NGSolve Zaglmayr integrated-Jacobi *hierarchical* simplex basis,
  if a bit-exact match to NGSolve's dof ordering is ever needed (the spaces are
  already covered by `simplex_ho.wls`).
