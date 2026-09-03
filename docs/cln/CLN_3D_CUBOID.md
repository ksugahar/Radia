# 3D Cauer Ladder Network for the Cu Cuboid Benchmark

This document describes the path from the (working) 2D axisymmetric Cu disk
Cauer cross-validation to the 3D Cu cuboid analysis that is the actual
target of the Nagamine 2026 paper [^Nagamine2026].

The 2D disk benchmark (see [`docs/axifem/AXIFEM.md`](../axifem/AXIFEM.md))
showed that two completely different formulations -- BEM-Foster integral
equation (Mathematica + 50-digit mpmath Nagamine pipeline) and Henrotte
axisymmetric FE (`radia.axifem` in radia-core, `axihenrotte` order=1/2) +
Hiruma 3-term --
agree on the leading Cauer rung `τ_pair[0] = L_1/R_0` to **0.28 %**. The
3D cuboid is the next step.

---

## 1. Why HCurl, not H1, for 3D

In the 2D axisymmetric Cu disk problem, the unknown `A_φ(r, z)` is a
**scalar** (only the φ-component of the vector potential is non-zero by
symmetry), so an H1 (or Henrotte H1) FE space suffices.

In 3D, the magnetic vector potential `A(x, y, z)` is genuinely **3-component
vector**. The natural FE space is **HCurl** (curl-conforming) because:

* `B = ∇×A` lives in HDiv (divergence-conforming) by construction;
* the bilinear form `∫ ν |∇×A|² dV` integrates over HCurl basis functions
  cleanly;
* tangential continuity at element interfaces is preserved (the right
  conformity for `A` in a magnetostatic / eddy-current weak form);
* gauge fixing by tree-cotree edge masking removes the kernel of `curl`
  without polluting the spectrum.

A scalar reduction with three separate H1 components for `(A_x, A_y, A_z)`
does **not** enforce `∇·A = 0` (Coulomb gauge) or `A·n = 0` on conductor
surfaces, and so admits spurious gradient modes. The
`cuboid_521_kameari_kelvin_v23_hiruma_3term.py` script is one such scalar
reduction (H1 Order=3, single-component, exploiting the cuboid's
1-direction symmetry of one mode at a time); it is **not** the right
generic 3D framework. **HCurl** is.

---

## 2. Existing 3D HCurl assets

The imported research assets are tracked under
[`validation_test/maglev/research_cln/ngsolve_validation/`](../../validation_test/maglev/research_cln/ngsolve_validation/):

| File | Purpose |
|---|---|
| [`cuboid_521_3dir_full_foster.py`](../../validation_test/maglev/research_cln/ngsolve_validation/cuboid_521_3dir_full_foster.py) | HCurl Order=2 Foster eigsh decomposition for Cu cuboid 5×2×1 mm in 3 directions (B_x, B_y, B_z), tree-cotree gauge, 25-30 modes per direction |
| `bem_foster_cauer.py` | 3D Cauer extraction from BEM Foster (cuboid version of `disk_bem_cauer.py`) |
| `bem_foster_cauer_highprec.wls` | Mathematica BEM Foster eigvals + amplitudes for cuboid (input to `bem_foster_cauer.py`) |
| `foster_cln_python_pure.py` | Pure-Python BEM K-matrix + Cauer pipeline for cuboid (mp.quad Tanh-Sinh integration; ~17 digits) |

The HCurl mesh + matrices in `cuboid_521_3dir_full_foster.py` are the right
starting point:

```python
from ngsolve import HCurl, BilinearForm, curl, dx
fes = HCurl(mesh, order=ORDER, dirichlet="conductor_surface", nograds=True)

a_form = BilinearForm(fes)
a_form += (1.0/mu0) * curl(u) * curl(v) * dx       # K (stiffness)
m_form = BilinearForm(fes)
m_form += sigma_Cu * u * v * dx                     # M (sigma mass)
```

Tree-cotree gauge (`build_spanning_tree_interior` in the same file) masks
out the gradient kernel of curl so the symmetric eigenproblem is
non-singular.

---

## 3. Roadmap: 3-way validation for the Cu cuboid

The same logical structure as the disk validation:

```
              .                           .
              .  3-way cross-validation   .
              .                           .

  (A) BEM Foster ──┐
       Mathematica │
       elliptic    │      Nagamine 2026 paper:
       integrals   │      summation -> Taylor moments alpha_n -> CFE
                   ├──>   -> Cauer ladder R_{2k}, L_{2k+1}
                   │
                   │      mpmath 50-digit (no interval arithmetic;
                   │       structurally Nagamine, not "verified")
                   │
  (B) HCurl FE ────┤
       NGSolve     │
       Order=2 +   │      Hiruma 3-term Lanczos:
       tree gauge  │      lambda_{2k+1} = 1/R_{2k}
                   │      lambda_{2k+2} = L_{2k+1}
                   │      tau_pair[k]   = lambda_{2k+1} * lambda_{2k+2}
                   │
                   ▼
                tau_pair[k] comparison (normalisation-invariant)
                per direction: B_x, B_y, B_z
```

Three directions because for a 5×2×1 mm cuboid the eddy-current modes
split by which component of an applied **uniform B** drives them:

* `B_x` drives currents in the y-z cross-section (effective dimensions 2×1 mm);
* `B_y` drives currents in the x-z cross-section (5×1 mm);
* `B_z` drives currents in the x-y cross-section (5×2 mm).

Each direction has its own Cauer ladder. The `cuboid_521_3dir_full_foster.py`
script already enumerates analytical TE-mode targets per direction (eq.
`tau_TE = μ_0 σ / (π² (m²/L₁² + n²/L₂²))`) and computes Foster R, L per
mode; we now need the **Cauer extraction** (Nagamine pipeline) on the BEM
side and the **Hiruma 3-term** on the HCurl side.

---

## 3.5 Pre-flight: 3D HCurl validation against axifem p=2 on the Cu disk

Before tackling the cuboid (which has no axisymmetric reference solution),
the 3D HCurl + Hiruma 3-term pipeline is **first calibrated on the Cu disk**
(R = 10 mm, t = 2 mm) for which `radia.axifem` in radia-core
(`axihenrotte p=2`) already
gives the canonical answer:

```
axihenrotte p=2 fine, B_z drive:
  tau_pair[0] = L_1/R_0 = 218.71 us         (BEM Cauer: 219.32 us, gap 0.28%)
  tau_pair[1]           =  78.12 us
  tau_pair[2]           =  39.54 us
  ...
```

We now mesh the **same disk in 3D**, solve with HCurl, and verify that
`tau_pair[0]` reproduces 218.7 µs.

### Key experimental constraint: matched curve order and FE order

For the comparison to be a clean test of the FE formulation rather than
a measurement of the geometric-discretisation error, the **mesh curve
order must equal the FE polynomial order**:

```python
mesh = Mesh(OCCGeometry(disk).GenerateMesh(maxh=h))
mesh.Curve(ORDER)                                # geometry curve order
fes = HCurl(mesh, order=ORDER, dirichlet="...", nograds=True)
```

For HCurl `order=2`, set `mesh.Curve(2)` so that the cylindrical conductor
boundary is represented by quadratic surface patches rather than flat
facets. With unmatched orders the boundary geometry approximation
(`mesh.Curve(1)`, flat tets) introduces an `O(h²)` error on the curved
side surface that swamps the FE convergence.

### Direction: only B_z is needed for the disk

By rotational symmetry, the three directions B_x, B_y, B_z that the
cuboid script enumerates collapse on the disk: B_x and B_y give modes
that vanish identically (the disk has no preferred radial direction in
the x-y plane), and B_z is the only physically meaningful drive. This
matches the axisymmetric reduction in `radia.axifem` (which only treats
B_z).

### Expected outcome and what failure means

* If the 3D HCurl + Hiruma 3-term result matches axifem p=2 to ≲ 1 % at
  `tau_pair[0]`: the 3D HCurl pipeline is correctly implemented, and we
  can confidently move to the cuboid (where no axisymmetric reference
  exists).
* If they disagree: the 3D HCurl side has a bug (most likely tree-cotree
  gauge masking, A_ext projection, or BC application) — the
  axisymmetric Q2 result is the gold standard since it has been
  cross-validated against BEM-Foster (Nagamine pipeline) at 0.28 %.

This pre-flight is a **necessary** step before the cuboid claim is
credible: without the disk validation, a 5-10 % cuboid disagreement
between BEM and HCurl could be FE error or BEM error, but with the disk
validation pinned, any cuboid disagreement is unambiguously attributable
to the cuboid's specific challenges (anisotropy, mesh quality at sharp
corners).

---

## 4. Implementation tasks

### 4.1 HCurl + Hiruma 3-term (FE side)

A new script `cuboid_521_3dir_hiruma_3term.py` next to
`cuboid_521_3dir_full_foster.py`. Reuse:

* the HCurl FESpace with `nograds=True` and `dirichlet="conductor_surface"`,
* the tree-cotree gauge masking,
* `K` and `M` matrices,
* the `A_ext` dictionary keyed by direction.

For each direction `k_dir ∈ {x, y, z}`:

1. project `A_ext[k_dir]` onto HCurl via `set` (or L² projection) to get
   the starting GridFunction `A_ext_h`;
2. compute `b = M @ A_ext_h_vec` (the Hiruma RHS — the σ-mass-weighted
   external A);
3. apply the Hiruma 3-term recursion (port from
   [`validation_test/axifem/research/verification/test_hiruma_disk_q1.py:hiruma_3term`](../../validation_test/axifem/research/verification/test_hiruma_disk_q1.py))
   on `(K_red, M_red, b_red)`, where `_red` = restriction to FreeDofs;
4. read off Nagamine `R_{2k}, L_{2k+1}, tau_pair[k]` per stage as in the
   disk script.

Output: `cuboid_521_3dir_hiruma_results.json` with the same shape as the
disk-Q2 results (one `stages` list per direction).

### 4.2 BEM Cauer (Nagamine pipeline)

The 3D BEM scripts (`bem_foster_cauer*.py/wls`) already produce Foster
eigenvalues and amplitudes per direction. To get Cauer `R_{2k}, L_{2k+1}`,
extend them to:

1. compute Taylor moments `α_n = Σ_k g_k² τ_k^(n+1) · (-1)^n` per direction;
2. apply `cauer_extract` (50-digit mpmath, classical Cauer extraction) →
   `p_1, p_2, p_3, ...`;
3. map to Nagamine: `R_{2k} = p_{2k+1}`, `L_{2k+1} = 1/p_{2k+2}`,
   `tau_pair[k] = L_{2k+1}/R_{2k}`.

This is a direct copy of the disk pipeline
([`disk_bem_cauer.py`](../../../W%3A/30_CauerLadderNetwork/2026_04_01_長方形CLN/ngsolve_validation/disk_bem_cauer.py)).

### 4.3 3-way validation test

Mirror [`tests/test_3way_cauer_cross_validation.py`](../../validation_test/axifem/research/verification/test_3way_cauer_cross_validation.py)
in `cuboid_521_3way_cauer_cross_validation.py`: load BEM Cauer + HCurl
Hiruma JSON files, print the per-direction `τ_pair[k]` table, and assert
the leading rung agrees to ≲ 1 %.

---

## 5. Open questions and notes

* **Anisotropy**: the cuboid 5×2×1 mm is anisotropic (aspect ratio 5:2:1).
  The three directions give three different R, L sequences. Does Nagamine
  treat them separately, or as a 3×3 admittance tensor `Y_{ij}(s)`? The
  paper analyses each direction separately (eq 11 has scalar Y for one
  applied direction at a time); we follow the same convention.
* **Verified bounds**: Nagamine's paper uses MPFR + interval arithmetic
  for *rigorous* error bounds. Our pipeline uses 50-digit `mpmath`
  floats — high precision but not interval-rigorous. Adding `mpfi` or
  similar is a future enhancement.
* **Q-element generalisation**: the `axihenrotte p=k` Q-element trick
  (polynomial basis in `s = r²`) is **specific to the axisymmetric `1/r`
  weight**. In 3D Cartesian there is no such weight, so standard NGSolve
  HCurl-Order=k is already the right tool — no further basis engineering
  needed for the 3D case.
* **v23 status**: the existing `cuboid_521_kameari_kelvin_v23_hiruma_3term.py`
  uses scalar H1 + Kameari reduction. It is *not* the canonical 3D solver
  -- it works on a single-direction-aligned reduction and was used as a
  cross-check in 2D. The canonical path going forward is HCurl as
  described above.

---

## 6. Connection to radia-core axifem

The `radia-core` `radia.axifem` module
([`docs/axifem/AXIFEM.md`](../axifem/AXIFEM.md)) is **specific to
axisymmetric problems** -- its Henrotte basis polynomial in `s = r²` only
makes sense when the integrand has the axisymmetric `1/r` weight. The 3D
cuboid problem does not need a special basis; it uses standard NGSolve
HCurl. The two share:

* the **Hiruma 3-term Lanczos recurrence** (same algorithm, applied to
  whichever K, M pair the geometry produces);
* the **Nagamine Cauer-ladder convention** R_{2k}, L_{2k+1};
* the **3-way cross-validation methodology** (BEM-Foster + Nagamine
  pipeline vs FE Hiruma 3-term).

Result: the 3D extension is a direct application of the 2D framework
with the FE side swapped from `axihenrotte` to NGSolve's stock HCurl
on a tetrahedral/hex mesh of the cuboid. No new theory is needed for
the 3D case beyond what Nagamine 2026 already covers.

[^Nagamine2026]: H. Nagamine, T. Yamaguchi, K. Sugahara, S. Hiruma, T.
    Mifune, T. Matsuo, "Verified Numerical Computations of the Cauer
    Network Representation of a Square Prism Conductor", manuscript
    2026-05-04 (Japan Journal of Industrial and Applied Mathematics
    submission).
