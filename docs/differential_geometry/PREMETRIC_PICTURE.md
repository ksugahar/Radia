# The premetric picture — one split (`d` vs `⋆`) behind weak forms, Kelvin, IE, PML, HOIBC, hodograph

A single differential-geometry split organises the lab's whole open-boundary /
material / coil-design machinery. It is **premetric electromagnetism**
(Hehl–Obukhov; Kottler–Cartan–van Dantzig): Maxwell splits into a **topological**
half that is metric-free, and **one** constitutive operator — the **Hodge star** —
that carries *all* the metric, material, coordinate map, and open-boundary content.

Every claim below is locked by a self-testing Mathematica `.wls` in
[`packages/radia-mcp/src/radia_mcp/mathematica/differential_geometry/`](../../packages/radia-mcp/src/radia_mcp/mathematica/differential_geometry/)
(+ `../basis_functions/` for the FEEC bases and cohomology). This note is the map;
the `.wls` are the proof.

---

## 1. The split

```
  dF = 0 ,  dG = J        TOPOLOGY  —  metric-free, invariant under ANY chart   (the "d" half)
  G  = ⋆ F                GEOMETRY  —  the Hodge star carries metric+material   (the "⋆" half)
```

- **`d` (the exterior derivative)** is the *basis / topology* half. Its discrete
  form is the **de Rham complex** `H1 →grad→ H(curl) →curl→ H(div) →div→ L2`
  (Whitney / Nédélec / Raviart–Thomas elements). `d∘d = 0` is built in; its
  **cohomology** (Betti numbers `b_k`) counts the holes / global loops. **`d` knows
  nothing about the metric.**
- **`⋆` (the Hodge star)** is the *operator / geometry* half. It carries the
  **material** `ν, μ, ε`, the **metric** of a coordinate map, the **open-boundary
  DtN**, and — when complex — **absorption (PML)**. Everything that *moves* when you
  change coordinates, material, or truncation is in `⋆`.

The weak (energy) form `a(w,A) = ∫⟨dw, ν dA⟩` is a **Hodge pairing**: `d` supplies
the topology, `⋆_ν` supplies the metric + material. (`weakform_hodge.wls`.)

---

## 2. The master table — domain × (`d`, `⋆`)

| domain | basis = `d` (de Rham / FEEC) | operator = `⋆` (Hodge / DtN) | exact / approx |
|---|---|---|---|
| **volume** (interior FEM) | Whitney / Nédélec / RT, any order | scalar or tensor **material** `ν(x)` | exact (the PDE) |
| **exterior** (open boundary) | the **same** de Rham basis on a finite layer | **tensor Hodge**: **Kelvin** (real metric) / **PML** (complex metric) = the **infinite element** | Kelvin/IE **exact** DtN; PML approximate |
| **surface** (BEM / impedance) | **surface** de Rham (RWG / Nédélec-surface) | **surface Hodge = DtN/Steklov**: exact (Kelvin/IE) or **HOIBC** (local Padé) | HOIBC **approximate** |

**Reading the table:** the basis (`d`) is the *same* de Rham family everywhere; only
the operator (`⋆`) changes. The infinite element is **not a new element** — it is the
standard de Rham basis with a **tensor Hodge** (the Kelvin/PML metric) on a finite
layer. HOIBC is **not a new basis** — it is the **surface** de Rham basis with a
**local approximation of the surface DtN**.

---

## 3. Transformation optics = "push the geometry onto the Hodge"

Because only `⋆` carries the metric, a coordinate transformation can be traded
**entirely** for a material — this is **transformation optics**, and in premetric
language it is the slogan *push everything onto the Hodge*. The one map is the
**metric → material** relation

```
  χ(g) = √(det g) · g⁻¹       (the Hodge star written as a material tensor)
```

- our **Kelvin** material rule is `ν' = χ(P^T P)` (pullback metric), reproducing the
  **Nagamine–Yamaguchi–Sugahara** pullback law (spherical isotropic `(r'/R)²`,
  cylindrical anisotropic `diag(1,1,(ρ'/R)⁴)`);
- **Pendry's** vacuum TO material `ε' = Λ Λ^T / det Λ` is the **same** `χ` on the
  pushforward metric `(Λ Λ^T)⁻¹`. One map, two coordinate expressions.
(`weakform_hodge.wls §8`.)

**But only the GEOMETRY moves onto the Hodge — the TOPOLOGY cannot.** Two precise
limits:

1. **Curvature stays zero.** A coordinate map / Kelvin / material changes the metric
   but a diffeomorphism cannot curve flat space: the pullback metric `J^T J` is flat
   (`Riemann = 0`, conformal *and* non-conformal), while a genuine sphere is curved
   (`R = 2`). The metric is the *geometry*, not the *curvature*. (`weakform_hodge.wls §7,§9`.)
2. **Betti is metric-independent.** The harmonic-form count = `b_k` for **any** SPD
   metric (Hodge theorem): the Hodge / material **cannot change a hole**. *Topology is
   not in the Hodge.* (`cohomology.wls`.)

So: **push the geometry onto `⋆`; the topology stays in `d`.** This is exactly why
Kelvin / PML / transformation optics can never change the cohomology (the coil
windows, the yoke loops) — they only move the metric.

---

## 4. What the one split explains

| phenomenon | in the picture | evidence |
|---|---|---|
| infinite element (de Rham IE) | de Rham basis + exterior **tensor Hodge**; `= ` standard element ∘ Kelvin map | `infinite_element_derham.wls`, `weakform_hodge.wls` |
| Kelvin open boundary | **real** tensor Hodge; reproduces the exact DtN | `weakform_hodge.wls §3,§8` |
| PML | **complex** tensor Hodge (geometry-flexible, approximate) | (same `χ`, complex `g`) |
| HOIBC | **surface** de Rham basis + **local Padé** of the surface DtN | `surface_derham.wls [A],[B]` |
| two vector impedance ladders | surface Hodge–Helmholtz → H(curl) tangential `n/R`, H(div) normal `(n+1)/R` | `surface_derham.wls [C]` (cf. act7_30) |
| **DtN / Steklov operator itself** | the **exterior `⋆` condensed** to Γ: Schur complement / Riccati fixed point; `= ` shifted `√(-Δ_Γ)` (nonlocal, order +1); self-adjoint positive (the boundary `H^{1/2}` metric) | `dtn_geometry.wls` |
| nonlinear material (saturation) | `⋆_ν` field-dependent; tangent SPD ⇒ elliptic, no fold | `weakform_hodge.wls §6`, `hodograph.wls [5]` |
| Chaplygin linearisation | a canonical transformation that makes `μ(q)` a coefficient | `hodograph.wls [5]` |
| weak form = Hamiltonian's shadow | `δ(action)=0`; Legendre `B↔H`, flux line = Hamilton flow (`A_z`) | `canonical.wls` |
| global loops / cuts | the cohomology `H¹` (`d` side), `b1` metric-independent | `cohomology.wls` |

---

## 5. Honest limits

- **Metric ≠ curvature.** Our spaces are **flat** `ℝ³` (MQS / Laplace kernel); the
  Hodge looks non-trivial only because of the chart/material. Genuine curvature needs
  curved spacetime (GR) — same formalism, different regime.
- **HOIBC is approximate; Kelvin/IE are exact.** The surface DtN `Λ(n)=(n+1)/R` (and
  the vector `n/R`) is **nonlocal**; a finite-order HOIBC (a local Laplace–Beltrami
  polynomial) matches finitely many modes and then deviates at the *knee*. Higher
  order = more matched modes, never exact.
- **Discrete needs resolution.** At the continuous level `⋆` carries everything, but
  the de Rham basis must be **rich enough** (order `p`, isoparametric curving `k`) to
  resolve what a sharp/anisotropic `⋆` demands — the `2·min(p,k)` floor law and the
  curved-truncation requirement.
- **Real vs complex Hodge.** Static / MQS = **real** tensor Hodge (radia's regime);
  radiating / full-wave = **complex** Hodge (PML, Astley/Demkowicz–Pal wave-envelope).

---

## 6. References

- **Premetric EM:** F. W. Hehl, Yu. N. Obukhov, *Foundations of Classical
  Electrodynamics* (Birkhäuser, 2003); Kottler / Cartan / van Dantzig.
- **Differential forms in EM / FEEC:** A. Bossavit, *Computational
  Electromagnetism* (1998); D. Arnold, R. Falk, R. Winther, *Finite element
  exterior calculus* (Acta Numerica 2006).
- **Transformation optics:** J. B. Pendry, D. Schurig, D. R. Smith, *Science* 312
  (2006); U. Leonhardt, T. Philbin, *Geometry and Light* (2010).
- **Kelvin open boundary (pullback):** H. Nagamine, T. Yamaguchi, K. Sugahara,
  *A Pullback-Based Formulation of Kelvin Transformation*, CEFC 2026 (lineage
  Wong–Ciric 1985, Freeman–Lowther 1988, Sugahara 2022).
- **Infinite element (de Rham):** L. Demkowicz, J. T. Oden / Demkowicz–Pal, CMAME
  164 (1998); P. Bettess, *Infinite Elements* (1992).
- **Hodograph / Clebsch / helicity:** see
  [`../clebsch_hodograph/HODOGRAPH_BACKBONE.md`](../clebsch_hodograph/HODOGRAPH_BACKBONE.md);
  Robert (IEEE TMag 1991), Moffatt (JFM 1969).
- **DtN spectrum / IE-vs-Kelvin:** `docs/open_boundary/`,
  `examples/kelvin_transformation/DtN_spectrum/` (act7_* incl. act7_30 two ladders).
