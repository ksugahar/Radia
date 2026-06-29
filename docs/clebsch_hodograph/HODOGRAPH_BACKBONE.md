# The coordinate-transform backbone — one diagram, every variant a cell

This is the **spine** for the `docs/clebsch_hodograph/demos/` line. The two
companion notes are specializations of it:

- [`DESIGN_METHODOLOGY.md`](DESIGN_METHODOLOGY.md) — the A–φ design loop (iron
  pole = equipotential), the φ-side application of the spine.
- [`HDIV_VIM_CLEBSCH_BRIDGE.md`](HDIV_VIM_CLEBSCH_BRIDGE.md) — the HDiv-VIM /
  Clebsch / flux-line face, the de Rham application of the spine.
- [`DIFFERENTIAL_GEOMETRY_WLS.md`](DIFFERENTIAL_GEOMETRY_WLS.md) — the
  executable Mathematica guardrail for the spine (`weakform_hodge`,
  `hodograph`, `canonical`, `surface_derham`, `dtn_geometry`).

The public result-saved demo entry point is
[`public_demo.ipynb`](public_demo.ipynb). The full source/result catalog is
[`examples_catalog.ipynb`](examples_catalog.ipynb). Together they keep the
runnable example sources, saved figures/JSON, and protected validation
references visible beside this theory map.

The point of this note (the user's *"微分幾何で理解するべき"* — understand it
through differential geometry): the Kelvin transform, the hodograph (both
kinds), the von Mises partial transform, the bidirectional pullback, and the
Chaplygin saturation linearization are **not five tricks**. They are five
**cells of one diagram**, indexed by three independent axes. Reading the diagram
tells you (i) why each cell behaves as it does, (ii) which pieces are classical,
and (iii) which cells are built and which are open. It is a working map of the
hodograph machinery, not an originality claim.

---

## 1. The one diagram

Magnetostatics is a **closed 2-form with a (possibly nonlinear) Hodge-star
constitutive law**:

```
  B ∈ Ω²   (flux 2-form)          dB = 0           (de Rham: B = dA, "no monopoles")
  H ∈ Ω¹   (field 1-form)         dH = J           (Ampère, source)
  H = ⋆_ν B                       (constitutive: the Hodge star carries the material ν)
```

Every transform in this line is a **diffeomorphism** `φ : 𝓜 → 𝓝` (a change of
coordinates, or a change of which variables are the coordinates). Under `φ`, the
two halves behave completely differently — and that split is the whole story:

| half | what it is | how it transforms | consequence |
|---|---|---|---|
| `dB = 0` | topology / de Rham | **natural** — pulls back with no metric data | preserved by *any* `φ` |
| `H = ⋆_ν B` | geometry / metric | the **Hodge star** pulls back with a **weight** `W` | this is the only thing that moves |

For the scalar/Dirichlet energy `∫ du ∧ ⋆dv = ∫ ∇u·∇v dx`, the pulled-back
Hodge star is the explicit weight

```
  ∫_physical ∇u·∇v dx  =  ∫_computational (∇u)ᵀ W (∇v) dξ ,
  W = |det J| (Jᵀ J)⁻¹ ,        J = ∂(physical)/∂(computational).
```

`W` **is** the discrete Hodge star written in computational coordinates. This
identity is verified to machine precision (deformed-mesh assembly == reference +
`W`, `~1e-16`, conformal and non-conformal `φ`) in
[`bidirectional_coordinate_transform_2d.py`](demos/bidirectional_coordinate_transform_2d.py)
(golden `test_bidirectional_coordinate_transform_2d`). Everything below is what
`W` and `⋆_ν` do in each cell.

> **Symbolic (Mathematica) twin of this section.** The Hodge / material half is
> reproduced symbolically in
> [`packages/radia-mcp/.../mathematica/differential_geometry/`](../../packages/radia-mcp/src/radia_mcp/mathematica/differential_geometry/):
> `weakform_hodge.wls` derives the material modulation `ν' = ν |det P|(PᵀP)⁻¹`
> (the **same** formula as `W = |det J|(JᵀJ)⁻¹`) and lands the pullback-Kelvin
> `ν'` — spherical isotropic `(r'/R)²`, cylindrical anisotropic
> `diag(1,1,(ρ'/R)⁴)` — plus the SPD nonlinear `⋆_ν` tangent (no fold);
> `hodograph.wls` does the 3-axis backbone (cochain-map split, conformal
> weight-freeness, Clebsch + helicity, `A_z`-as-Clebsch-potential, Chaplygin
> linearisation); `canonical.wls` records the Hamiltonian / Legendre reading of
> the hodograph; `surface_derham.wls` splits HOIBC into surface topology plus
> analytic DtN / Steklov geometry; and `dtn_geometry.wls` records the DtN /
> Steklov operator itself as condensed exterior Hodge geometry. All self-tests
> `ALL PASS`. The human
> index is
> [`DIFFERENTIAL_GEOMETRY_WLS.md`](DIFFERENTIAL_GEOMETRY_WLS.md); the MCP mirror
> is `differential_forms_mathematica_recipes(topic="differential_geometry" |
> "weakform_hodge" | "hodograph" | "canonical" | "surface_derham" |
> "dtn_geometry")`.

---

## 2. The three axes

Each cell of the diagram is a choice on three independent axes.

### Axis 1 — which map `φ` (and *how much* of the field it moves)

| map | new coordinates | PDE in new coords | ∞ | single-valued? | weight `W` |
|---|---|---|---|---|---|
| **geometric Kelvin** `z ↦ R²/z̄` | base space, inverted | unchanged (same `⋆`) | → **point** | yes (geometry) | conformal in 2-D ⇒ `I`; `(R/ρ′)²` in 3-D/axisym |
| **potential plane** `(Φ, A_z)` | two potentials (0-forms) | **trivial** (CR / `W=I`) | **unbounded** (`|W_pot|→∞`) | yes | `I` (2-D conformal) |
| **field plane** `(θ, q=|B|)` | two field comps | non-trivial (Chaplygin) | → **point** | can **fold** | variable |
| **partial / von Mises** `(x, A_z)` | keep one base coord + **one** Clebsch potential | reduced | → **bounded** (`A_z` range finite) | **yes** (no fold) | reduced |

The last row is the **partial (von Mises) transform** this line works with:
*keep one coordinate, transform one Clebsch potential.* Among the four it is the
one that is both single-valued / fold-free **and** compactifies the open boundary
(the `A_z` range is bounded because flux is conserved), at the price of only a
*partial* PDE simplification (not the full Chaplygin linearization). The full
field-plane `(θ,q)` buys the full linearization but can fold; the potential plane
buys a trivial PDE but leaves ∞ unbounded. No single 2-variable hodograph gives
all three — so "Kelvin in the hodograph" is a *composition* (a partial / Clebsch
representation **plus** a geometric Kelvin inversion), not a single map.

### Axis 2 — dimension (what survives in `W`)

The 2-D miracle is purely metric: **the Hodge star on 1-forms is conformally
invariant in 2-D**, so a conformal `φ` leaves `⋆` (hence `W`) untouched.

| dim | conformal `φ` ⇒ `W` | weight that survives | reason |
|---|---|---|---|
| **2-D** | `J = sR ⇒ W = s²(s²I)⁻¹ = I` | none (**weight-free**) | `⋆` conformally invariant on `Ω¹(ℝ²)` |
| **axisym** | — | `2πr` | the revolution Jacobian is a genuine metric factor |
| **3-D** | `J = sI ⇒ W = sI` | scale `s` (and more) | `⋆` on `Ω¹(ℝ³)` is **not** conformally invariant |

This is exactly why
[`hodograph_kelvin_2d.py`](demos/hodograph_kelvin_2d.py)
is weight-free (`μ′ = μ₀`) while
[`hodograph_kelvin_axisym.py`](demos/hodograph_kelvin_axisym.py)
and
[`clebsch_kelvin_3d.py`](demos/clebsch_kelvin_3d.py)
carry the `(R/ρ′)²` / `2πr` Kelvin weight (the `symbolic_pullback_check()` in the
bidirectional example derives all three from the same `W`).

Dimension also decides whether the **field-side** transform even exists globally:
in 3-D the Clebsch representation `B = ∇α×∇β` needs **two** potentials for
**three** coordinates, and a global pair exists iff the **helicity**
`h = ∫A·B` vanishes (Moffatt 1969). Verified in
[`clebsch_3d_closing_condition.py`](demos/clebsch_3d_closing_condition.py)
(Clebsch `h = −7e-17`, ABC Beltrami `h = 3(2π)³ ≠ 0`, chaotic). **This is why 3-D
does not auto-linearize** — the full field-plane interchange has a topological
obstruction; only the *partial* (one-potential) transform is always available.

### Axis 3 — material (linear `⋆` vs nonlinear `⋆_ν`)

Saturation makes the Hodge star **field-dependent**: `H = ν(|B|) B`, i.e.
`⋆_ν` is a metric that depends on the form it acts on. Its tangent (the
linearization the Newton/Picard step sees) is

```
  ∂H/∂B  =  ν(|B|) I  +  ν′(|B|) (B ⊗ B)/|B| .
```

For ordinary magnetic saturation `ν` and `ν′` are **≥ 0** (reluctivity rises
monotonically), so the tangent metric stays **symmetric positive-definite**. The
problem therefore stays **elliptic** — `⋆_ν` never changes type. **This is the
differential-geometric reason magnetic saturation does not fold the hodograph**:
unlike transonic gas dynamics (where the analogous operator changes type at the
sonic line and a *limiting line* appears), the magnetic `⋆_ν` has no limiting
line. Saturation can only deform the metric, never flip its signature; *only the
geometry (axis 1, the field plane) can fold.*

The same SPD fact picks the well-conditioned formulation: the convex coenergy
`∫ W(|B|)` (`W′ = ν` monotone ⇒ `W` convex) is the **B-input A-formulation**,
the de Rham dual where the solve is stable at the knee — see
[`saturation_loop_2d.py`](demos/saturation_loop_2d.py)
and the conditioning discussion in
[`HDIV_VIM_CLEBSCH_BRIDGE.md`](HDIV_VIM_CLEBSCH_BRIDGE.md) §2.

---

## 3. The cell map (built vs open — honest)

Reading the axes against the realized examples. `✓` = verified + golden-locked
(`tests/feec/test_clebsch_hodograph_research.py`); `→` = open / planned.

| map (axis 1) | dim | material | cell artifact | status |
|---|---|---|---|---|
| geometric Kelvin | 2-D | linear | `hodograph_kelvin_2d.py` (err `~2e-8` vs air-box `3e-3`) | ✓ |
| geometric Kelvin | axisym | linear | `hodograph_kelvin_axisym.py` | ✓ |
| geometric Kelvin | 3-D | linear | `clebsch_kelvin_3d.py` | ✓ |
| geometric Kelvin | 3-D | **saturable** | `clebsch_kelvin_nonlinear_3d.py` (merged geometry+material Picard, rung 3) | ✓ (stable regime) |
| potential plane `(Φ, A_z)` | 2-D | linear | `a_method_clebsch_2d.py`, `hdiv_vim_clebsch_2d_az.py` (`A_z` *is* the Clebsch potential, recovered `4e-7`) | ✓ |
| field plane `(θ, q)` | 2-D | **saturable** | `chaplygin_hodograph_2d.py` (1-shot quadrature == full nonlinear loop), `chaplygin_turning_guide_2d.py` | ✓ |
| field plane `(θ, q)`, free boundary | 2-D | linear / saturable | `chaplygin_free_boundary_2d.py`, `chaplygin_inverse_vonmises_2d.py` (`(Φ,A)` rectangle), `chaplygin_inverse_nonlinear_2d.py` | ✓ |
| **partial `(x, A_z)`** | 2-D | **saturable** + open boundary | the partial-transform saturable cell | **→ task #50, not yet built** |
| partial / Clebsch | 3-D | saturable | `clebsch_kelvin_nonlinear_3d.py` is the *merged* surrogate; a genuine field-coordinate 3-D reduction is the prize | **→ task #51 / open (helicity)** |

Two honest reads fall straight out of the table:

1. **The main unbuilt hodograph piece is the partial `(x, A_z)` saturable open
   boundary.** The partial transform is realized in the *linear* `A_z`-is-Clebsch
   sense (`hdiv_vim_clebsch_2d_az.py`), and saturation is in hand via the *full*
   field-plane Chaplygin (`chaplygin_hodograph_2d.py`), but the two have not been
   combined: the partial `(x, A_z)` *saturable* open-boundary solve is not yet
   assembled. The "saturation = one linear solve" result in hand is the full
   `(θ,q)` double-swap; the partial-map version (task #50) is the hodograph work to
   do — the open boundary compactifies there because the `A_z` range is bounded, a
   different mechanism from the field-plane `q → H₀` point.

2. **3-D saturable open boundary is genuinely open**, gated by axis 2's helicity
   obstruction (no global Clebsch pair where `h ≠ 0`). The current 3-D capstone
   (`clebsch_kelvin_nonlinear_3d.py`) sidesteps it with a single merged
   geometry+material Picard rather than a field-coordinate linearization. Task #51
   (air-box-vs-Kelvin on a 3-D saturable body) measures the *open-boundary* win in
   that regime; the *linearization* prize remains open.

---

## 4. Prior art — the pieces are classical

Every constituent transform is classical; this note is a working map of the
machinery, not an originality claim. Stated in the diagram's own language:

| piece | differential-geometry statement | source |
|---|---|---|
| complex-potential pole design `(Φ, A_z)` | conformal `φ` ⇒ trivial pullback of `⋆` | 1900s |
| 2-D conformal weight-freeness | `⋆` conformally invariant on `Ω¹(ℝ²)` | classical |
| geometric Kelvin FE open boundary | conformal compactification of the *base* | Remacle/Lowther/Imhoff/Brunotte/Nicolet (~30 yr) |
| Chaplygin / von Mises hodograph | linearize `⋆_ν` in field coordinates | gas dynamics, classical |
| Clebsch + helicity obstruction | global `B = dα∧dβ` ⇔ `h = 0` | Robert 1991, Moffatt 1969 |

The work in this line is the **hodograph machinery itself** and its
finite-element realization — composing a partial / Clebsch representation with a
geometric Kelvin inversion, keeping the nonlinear `⋆_ν`, and pushing the open
cells in §3.

**The modern differential-forms framing of these transforms is itself recent and
already occupied — it is NOT an open gap.** Two CEFC 2026 papers do precisely the
exterior-calculus / pullback re-expression of the machinery in §1:

- **Dervisha, Marjamäki, Rasilo, Tarhasaari (Tampere), CEFC 2026** — *Bidirectional
  Coordinate Transformation and Its Application to 2-D Magnetic Field Problems*: the
  exterior-calculus derivation of the bidirectional `(x,y) ↔ (A, φ)` map (both
  potentials as the new coordinates — the full potential-plane chart), with
  `⋆dA = µ dφ`, `⋆dφ = −(1/µ)dA`, and `µ = µ(A,φ)` permitted; it even redraws the
  material geometry in `(A,φ)`. Example is a *bounded* core-limb + airgap — **no open
  boundary, no partial `(x, A_z)` chart.**
- **Nagamine, Yamaguchi, Sugahara, CEFC 2026** — *A Pullback-Based Formulation of
  Kelvin Transformation in Electromagnetic Field Analysis*: the pullback `k*`
  formulation of the Kelvin open boundary, with the material law `ν' = (r'/R)² ν`
  (spherical, conformal — isotropic) and `ν' = diag(1,1,(ρ'/R)⁴) ν` (cylindrical,
  non-conformal — anisotropic) derived by **equating the bilinear energy functionals**
  — i.e. exactly the `W`-as-pulled-back-`⋆` of §1, and the conformal/non-conformal
  split of axis 2. **Linear material; no saturation.** (`symbolic_pullback_check()`
  in `bidirectional_coordinate_transform_2d.py` re-derives this same `ν'`.)

So the geometric machinery (§1) is established prior art. What **neither** paper does —
the *combination* of a partial / field chart **with** the Kelvin open boundary **and** a
nonlinear `⋆_ν` (the §3 open cells, tasks #50/#51) — is simply not yet built. Recorded
here as "what is and isn't done," not as an originality claim.

---

## 5. How to use this spine

- **Designing a new variant?** Pick a cell `(map, dim, material)`. The axes tell
  you the weight (`W` from axis 2), whether it can fold (axis 1), and whether it
  stays elliptic (axis 3) *before* you write the solve.
- **Debugging a transform?** Check the two halves separately: `dB = 0` is
  `φ`-independent (if it breaks, the representation leaked solenoidal content —
  trace a flux line, `HDIV_VIM_CLEBSCH_BRIDGE.md` §dynamical face); the metric is
  in `W` / `⋆_ν` only.
- **Tempted to stake an originality claim?** Don't, by default — the pieces are
  classical (§4). Put the effort into the hodograph cell instead.

---

## References

Classical foundations are cited in the companion notes; the spine-specific ones:

- **P. Robert**, "Clebsch Potentials and the Visualization of Three-Dimensional
  Solenoidal Vector Fields," *IEEE Trans. Magn.* **27**(5), 1991 — field `H ∈ Ω¹`,
  flux `B ∈ Ω²`, global Clebsch as a de Rham question.
- **H. K. Moffatt**, "The degree of knottedness of tangled vortex lines,"
  *J. Fluid Mech.* **35**, 1969 — helicity as the obstruction (axis 2).
- **von Mises / Chaplygin / Molenbroek** hodograph (gas dynamics) — the field-plane
  linearization and the limiting-line (fold) that magnetic `⋆_ν` does *not* have
  (axis 3), and the source of the partial (von Mises) transform itself.
- Kelvin-transformation open-boundary FE: Remacle 1995, Lowther 1989,
  Imhoff/Brunotte/Nicolet — the geometric (base-space) compactification (axis 1);
  lineage Wong & Ciric (*COMPEL* 4(3), 1985), Freeman & Lowther (*IEEE TMag* 24(6),
  1988), Sugahara (*IEEE TMag* 58(9), 2022).
- **Differential-forms framing of the transforms themselves (recent prior art):**
  A. Dervisha, A. Marjamäki, P. Rasilo, T. Tarhasaari (Tampere), "Bidirectional
  Coordinate Transformation and Its Application to 2-D Magnetic Field Problems,"
  *CEFC 2026* — exterior-calculus `(x,y) ↔ (A, φ)`, `µ(A,φ)` allowed, bounded domain;
  H. Nagamine, T. Yamaguchi, K. Sugahara, "A Pullback-Based Formulation of Kelvin
  Transformation in Electromagnetic Field Analysis," *CEFC 2026* — pullback Kelvin,
  material law `ν'` via bilinear energy equivalence, linear material. The §3 open
  cells are the un-built *combination* (partial chart + Kelvin + nonlinear `⋆_ν`).
