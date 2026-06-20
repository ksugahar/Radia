# Infinite elements for open-boundary problems — state of the art

Where the infinite-element (IE) method stands today, and how it sits against the closures Radia
already ships (**Kelvin** inversion-FE) or can reach (**PML**, **BEM**). Companion to
[`OPEN_BOUNDARY_MAP.md`](OPEN_BOUNDARY_MAP.md) and the measured
[`DTN_SPECTRUM_COMPARISON.md`](DTN_SPECTRUM_COMPARISON.md). Citations are Crossref-verified (DOIs
below); foundational works without a DOI here are named.

## 0. TL;DR positioning

| closure | implementation ease | geometry | regime | maturity / status |
|---|---|---|---|---|
| **Kelvin** (inversion / shell transform) | **easiest by far** — one coordinate map + standard FE, parameter-free | **sphere-locked** (Liouville: 3-D conformal maps are Möbius) | static / quasi-static (Laplace) | classical, ubiquitous in low-freq magnetostatics |
| **PML / CFS-PML** | easy — a complex-stretched layer; NGSolve ships it | flexible (`Cartesian`/`Brick`/radial) | waves + (CFS) low-freq | **the de-facto default** for general open-boundary EM |
| **infinite element** | **hard** — special radial basis, ∞-integration, conditioning-aware basis; NGSolve has none | flexible (wraps any convex / spheroidal surface) | static → waves (regime-specific basis) | **mature in acoustics**, rigorous-modern via Hardy space; niche in EM |
| **BEM** | medium (dense; `ngsolve.bem` FMM) | any | all | exact, dense `N²` |

The honest one-liner: **on the sphere the IE and Kelvin are the *same method*** (measured,
[`act7_28`](../../examples/kelvin_transformation/DtN_spectrum/act7_28_ie_vs_kelvin_fair_dtn.py));
the IE's reason to exist over Kelvin is **geometry** (no sphere-lock), and its cost is
**implementation + a conditioning-aware basis**.

## 1. The families

1. **Mapped / decay-function IE** (the origin). Zienkiewicz & Bettess mapped infinite elements;
   Bettess, *Infinite Elements* (Penshaw, 1992). In EM the low-frequency "ballooning" open
   boundary (Silvester–Lowther lineage) is the same idea. Radial shape functions decay as `(a/r)^k`
   — exactly the basis of [`act7_25`](../../examples/kelvin_transformation/DtN_spectrum/act7_25_infinite_element_dtn.py).

2. **Acoustic conjugated IE — the production workhorse.** For exterior Helmholtz the radial basis
   must carry the outgoing phase `e^{ikr}`; the **Astley–Leis "mapped wave-envelope" / conjugated**
   formulation uses a complex-conjugate test weight so the ∞-radial integral converges (the
   dominant production IE for exterior acoustics). **Burnett's** prolate/oblate-spheroidal multipole
   IE (DOI [10.1121/1.411286](https://doi.org/10.1121/1.411286), 1994; 3-D follow-up
   [10.1121/1.419883](https://doi.org/10.1121/1.419883)) is the separable spheroidal variant — the
   natural IE for an **elongated** body. Review: Gerdes, *J. Comput. Acoust.* 8 (2000).

3. **Conditioning is the practical crux** (and it is solved). The naive monomial/`(a/r)^k` basis is
   Hilbert/Cauchy-ill-conditioned — **measured in-repo**: `cond ≈ 10 → 4e9` for `P=2..8`
   ([`act7_27`](../../examples/kelvin_transformation/DtN_spectrum/act7_27_ie_vs_kelvin_vs_pml_gate1.py)).
   The fix is a **Jacobi-polynomial radial basis**: Dreyer & von Estorff, "Improved conditioning of
   infinite elements for exterior acoustics" (DOI [10.1002/nme.804](https://doi.org/10.1002/nme.804),
   2003; robustness [10.1016/j.cma.2005.01.019](https://doi.org/10.1016/j.cma.2005.01.019), 2006;
   "Efficient Infinite Elements based on Jacobi Polynomials"
   [10.1007/978-3-540-77448-8_9](https://doi.org/10.1007/978-3-540-77448-8_9), 2008). Our
   [`act7_28`](../../examples/kelvin_transformation/DtN_spectrum/act7_28_ie_vs_kelvin_fair_dtn.py)
   reproduces the lesson: an orthogonal basis for the *same* space drops `cond` to `≈ 2 → 339`.

4. **Maxwell / vector (de Rham) IE.** Demkowicz & Pal, "An infinite element for Maxwell's
   equations" (*CMAME* 164, 1998) — the exact-sequence vector IE that
   [`act7_26`](../../examples/kelvin_transformation/DtN_spectrum/act7_26_derham_infinite_element.py)
   verifies in the static limit; Cecot–Demkowicz–Rachowicz built the 3-D hp version.

## 2. The modern state of the art — Hardy-space infinite elements

The rigorous current SotA is the **Hardy-space IE** (pole-condition radiation): Hohage & Nannen,
"Hardy Space Infinite Elements for Scattering and Resonance Problems"
(*SINUM*, DOI [10.1137/070708044](https://doi.org/10.1137/070708044), 2009). Properties that make
it the reference today:

- **Proven convergence** — super-algebraic for separable problems; the error decays **exponentially**
  in the number of Hardy modes.
- **Eigenstructure-preserving** — it does not pollute the spectrum, so it is the method of choice for
  **resonance / scattering-pole** problems (nano-optics, open resonators).
- **Handles non-compact inhomogeneities / waveguides** (arXiv:1004.1025) and **dispersive media with
  opposite phase/group-velocity signs** (Halla et al., *Numer. Math.*,
  DOI [10.1007/s00211-015-0739-0](https://doi.org/10.1007/s00211-015-0739-0), 2015).
- **de Rham / Maxwell, high order:** Nannen et al., **"Exact Sequences of High-Order Hardy Space
  Infinite Elements for Exterior Maxwell Problems"** (*SISC*,
  DOI [10.1137/110860148](https://doi.org/10.1137/110860148), 2013) — the modern, convergent,
  high-order, exact-sequence **vector** IE. This is the rigorous descendant of Demkowicz–Pal and the
  production-grade target if a de Rham IE is ever built; `act7_26` is its static, low-order kernel.

## 3. Electromagnetics specifically

Gómez-Revuelto, García-Castillo & Demkowicz, "A comparison between PML, infinite elements and an
iterative BEM as mesh truncation methods for hp self-adaptive procedures in electromagnetics"
(*PIER*, DOI [10.2528/pier12020201](https://doi.org/10.2528/pier12020201), 2012): all three are
exact at the continuous level and hp-compatible; **no universal winner — choice by application**.
In practice **PML is the common EM default** (geometry-flexible, wave-native), with FEM-BEM where
exactness is required; the IE is the niche choice. (For low-frequency magnetostatics the inversion /
Kelvin shell is the classical workhorse.)

## 4. What our in-repo measurements contribute

- **IE ≡ Kelvin on the sphere** (`act7_28`): the Kelvin map sends `r^{−(n+1)}` → polynomial
  `ξ^{n+1}`, and the IE matrix `A_kl = a(kl+n(n+1))/((k+l)−1)` *is* the Kelvin-mapped monomial energy
  → identical DtN at matched DOF. They differ only in basis conditioning (item 3 above).
- **Geometry is the IE's edge** (`act7_27a`): exterior DOF scales `~AR²` for Kelvin (sphere-lock)
  vs `~AR` for IE/box-PML — so the IE earns its keep on **elongated / planar** bodies.
- **The de Rham static kernel is verified twice** — sympy (`act7_26`) and a Mathematica twin
  (`radia_mcp/.../basis_functions/infinite_element_derham.wls`) — the low-order static analog of
  Nannen 2013.
- **The static de Rham HIGH-ORDER IE is completed** (`act7_29`): the exact sequence holds at
  arbitrary radial order, the orthogonal basis bounds conditioning in *every* form end (scalar-H1
  `10→4e9`→`2→339`, vector-L2 `19→1.5e10`→`3→15` for `P=2..8`), and the DtN is spectral
  (`n≤P−1`) — a usable, arbitrary-order, exact-sequence exterior element (the static kernel of
  Nannen 2013; C++/3-D port is downstream).

## 5. Build decision for Radia (Gate 2 / Gate 3 spec)

**If an IE is built, it is a re-implementation of a mature published method, not novelty.** The spec:

- **Basis (static):** Jacobi / integrated-Legendre radial basis (Dreyer–von Estorff) — never the
  naive monomials (`act7_27` conditioning).
- **Basis (waves):** Astley–Leis conjugated `e^{ikr}` envelope, or the Hardy-space modes
  (Hohage–Nannen) for proven convergence + resonance use.
- **Vector / de Rham:** the exact-sequence Hardy-space Maxwell IE (Nannen 2013); `act7_26` is its
  kernel.
- **The honest case to build it:** a **homogeneous-exterior, static/low-freq, elongated or planar**
  magnetostatics problem where Kelvin wastes air (sphere-lock) *and* the exterior carries no iron
  (Kelvin's material capability is then moot). On a compact / spherical body, or with an iron
  exterior, **Kelvin is easier and at least as good** — build nothing.

**Gates** (verify-first, "settle in Python before C++"):
- **Gate 1 — done** (`act7_27`): geometry (Kelvin out `~AR²`) + the naive-basis conditioning trap.
- **Gate 1-fair — done** (`act7_28`): on the DtN yardstick IE≡Kelvin; deficit is basis-only, fixable.
- **Gate 2 — next:** a real 3-D IE (Jacobi basis) vs **box-PML** vs Kelvin on an elongated body —
  does the geometry edge survive in a full solve, and beat the already-shipped box-PML?
- **Gate 3:** a pybind11 C++ IE **only if** Gate 2 confirms a decisive, recurring need.

## 6. References (Crossref-verified DOIs)

- Burnett 1994, *JASA*, [10.1121/1.411286](https://doi.org/10.1121/1.411286) (+ 1997
  [10.1121/1.419883](https://doi.org/10.1121/1.419883)) — prolate-spheroidal multipole IE.
- Dreyer & von Estorff 2003, *IJNME*, [10.1002/nme.804](https://doi.org/10.1002/nme.804);
  2006, *CMAME*, [10.1016/j.cma.2005.01.019](https://doi.org/10.1016/j.cma.2005.01.019);
  von Estorff 2008, [10.1007/978-3-540-77448-8_9](https://doi.org/10.1007/978-3-540-77448-8_9) —
  Jacobi-polynomial conditioning.
- Hohage & Nannen 2009, *SINUM*, [10.1137/070708044](https://doi.org/10.1137/070708044) — Hardy-space IE.
- Nannen et al. 2013, *SISC*, [10.1137/110860148](https://doi.org/10.1137/110860148) — exact-sequence
  high-order Hardy-space IE for **exterior Maxwell** (de Rham).
- Halla et al. 2015, *Numer. Math.*, [10.1007/s00211-015-0739-0](https://doi.org/10.1007/s00211-015-0739-0)
  — Hardy-space IE, dispersive (phase/group sign).
- Gómez-Revuelto, García-Castillo & Demkowicz 2012, *PIER*,
  [10.2528/pier12020201](https://doi.org/10.2528/pier12020201) — PML vs IE vs BEM for hp EM.
- Demkowicz & Pal 1998, *CMAME* 164 — infinite element for Maxwell's equations (DOI not verified here).
- Bettess 1992, *Infinite Elements* (book); Astley–Leis mapped wave-envelope / conjugated IE; Gerdes
  2000, *J. Comput. Acoust.* 8 (IE review) — foundational, named (DOIs not verified here).
