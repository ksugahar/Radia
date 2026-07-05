"""MoM foundations: Harrington 1968, RWG, basis functions, collocation vs Galerkin."""

INTRODUCTION = r"""
# Method of Moments (MoM) — foundations

Reference: R.F. Harrington, "Field Computation by Moment Methods",
Macmillan 1968 (reissued IEEE Press 1993).  THE foundational text.

## The Method in 3 lines

1. Discretize unknown as basis-function sum: f(x) = Σ_n a_n φ_n(x)
2. Substitute into integral equation: L f = g  →  Σ_n a_n L φ_n = g
3. Test with test functions ψ_m: ⟨ψ_m, Σ_n a_n L φ_n⟩ = ⟨ψ_m, g⟩
   → matrix system Z·a = b

## Basis × test function choice = "MoM variant"

| Variant | Basis φ_n | Test ψ_m | Notes |
|---------|-----------|----------|-------|
| **Collocation** | local (e.g. pulse) | Dirac δ | Cheap, simple to code; less accurate |
| **Galerkin** | local | same as basis | Symmetric Z; better convergence |
| **Petrov-Galerkin** | local | different (e.g. flat pulses test) | For non-symmetric IE |
| **Subdomain pulse** | pulses on segments | match | Original wire-grid NEC style |
| **Entire-domain** | sin/cos modes | same | For small N, periodic geom |

Modern surface BEM = Galerkin RWG; modern wire MoM = Pocklington/Hallen
with sinusoidal basis (NEC).

## Choice of formulation (kernel)

For PEC scattering:
- **EFIE**: integrate over surface current J on PEC → unique except interior res.
- **MFIE**: integrate over surface current J via dual → unique except interior res.
- **CFIE = α·EFIE + (1-α)·η·MFIE**: avoids interior resonance, well-conditioned
  at high freq, breaks down at low freq.

For dielectric / penetrable scattering:
- **PMCHWT** (Poggio-Miller-Chang-Harrington-Wu-Tsai): coupled E/H on both
  sides of interface.  Well-posed but ill-conditioned at low freq.

For eddy current:
- **SIBC** (Leontovich): assume H_tangential is small inside conductor;
  enforce Z·H_t = E_t on surface.
- **Bíró-Preis A-V**: full coupled (A, V) inside conductor.
- **HOIBC** (high-order IBC): Mitzner curvature correction + Dong-Di Rienzo.

See `surface_ie` topic for details on each.
"""


RWG_BASIS = r"""
# RWG Basis Functions — Rao-Wilton-Glisson 1982

Reference: S.M. Rao, D.R. Wilton, A.W. Glisson, "Electromagnetic
Scattering by Surfaces of Arbitrary Shape", IEEE Trans. Antennas Propag.
30(3):409-418, 1982.  DOI: 10.1109/TAP.1982.1142818.

THE foundation of all modern surface BEM.  Every commercial EM tool
(HFSS, FEKO, CST surface) and every open-source library (Bembel,
ngsolve.bem, scuff-EM) uses RWG or RWG-derived bases.

## The RWG basis function

For each interior edge e = T_+ ∩ T_- between two triangles T_+ and T_-:

```
f_e(r) = (l_e / 2A_+) (r - r_+)    for r in T_+
       = (l_e / 2A_-) (r_- - r)    for r in T_-
       = 0                          otherwise
```

where:
- l_e = length of edge e
- A_± = area of triangle T_±
- r_± = vertex opposite to edge e in T_±

## Key properties

| Property | Why it matters |
|----------|----------------|
| Tangential component continuous across e | RWG ∈ H(div,Γ) → physical surface currents |
| Normal component zero on free edges | No spurious charge accumulation |
| div_Γ f_e = ±l_e/A_± (piecewise constant) | Charge density is constant per triangle |
| Star/loop decomposition possible | Loop-star basis for low-freq stabilization |

## In NGSolve

```python
from ngsolve import HDivSurface, Mesh
mesh = Mesh("model.vol")
fes = HDivSurface(mesh, order=1)   # RWG = HDivSurface order 1
u, v = fes.TnT()
# ... assemble EFIE / MFIE
```

In `ngsolve.bem`:
```python
from ngsolve.bem import SingleLayerPotentialOperator
S = SingleLayerPotentialOperator(fes, ...)
```

## Higher-order RWG

| Order | Name | DOFs per triangle |
|-------|------|-------------------|
| 0.5 (constant tangential) | RT0 / Raviart-Thomas | 3 |
| 1.5 | BDM1 | 6 |
| 2 | RT1 | 8 |
| ... | hp-version | (p+1)(p+2)/2 + (p+1) |

NGSolve `HDivSurface(mesh, order=p)` supports arbitrary p.
Modern Isogeometric BEM (Bembel) uses splines instead of polynomials —
even higher accuracy per DOF.
"""


WIRE_GRID = r"""
# Wire-grid modeling — legacy MoM (NEC heritage)

Before RWG, surface MoM used WIRE-GRID approximations: cover the
conductor surface with a grid of thin wires and solve the thin-wire
EFIE on each segment.

## NEC (Numerical Electromagnetics Code)

- Burke-Poggio 1981 LLNL, Fortran code
- Solves Pocklington / Hallén integral equation on wire segments
- Sinusoidal basis with point matching
- Still widely used for antenna design (HF, VHF, UHF)

## Issues with wire-grid

1. **Wire radius vs grid spacing** rule: l_wire ≥ 2a, h_grid ≈ 0.1λ
2. **Aliasing**: triangular surface currents get mapped to wire currents
   with errors O(h²)
3. **Near-field accuracy poor** within ~1 grid spacing
4. **Cannot model dielectric** scatterers (PEC only without extensions)

## Why RWG superseded wire-grid

RWG triangular basis:
- Models any topology (curved surfaces) cleanly
- Natural for FEM meshing tools
- Surface charge / current divergence well-defined
- 10-100x fewer DOFs than equivalent wire-grid

## When wire-grid is still used

- **Legacy NEC compatibility** — antenna design community
- **Quick rough estimates** for HF antennas where surface mesh overkill
- **Education** — easier to derive Pocklington than EFIE-RWG

The lab does NOT use wire-grid.  All BEM workflows go through RWG
surface meshes (or HDiv-VIM volume elements).

References (this folder):
- `01_mom_basics/NEC/nec2prt1.pdf` — NEC-2 manual
- `01_mom_basics/WireGrid/On_WireGrid_Representation_Solid_Metallic.pdf`
- `01_mom_basics/モーメント法_入門.pdf` — JP intro covering wire-grid
"""


def get_mom_foundations_knowledge(topic: str = "introduction") -> str:
    """Dispatch MoM foundations topics.

    Topics:
        introduction      - The Method in 3 lines (DEFAULT)
        rwg_basis         - Rao-Wilton-Glisson basis (1982)
        wire_grid         - Legacy NEC-style wire-grid MoM
        all               - Everything
    """
    topic = topic.lower().strip()
    if topic in ("introduction", "intro", "overview"):
        return INTRODUCTION
    if topic in ("rwg_basis", "rwg", "basis"):
        return RWG_BASIS
    if topic in ("wire_grid", "wiregrid", "nec"):
        return WIRE_GRID
    if topic == "all":
        return "\n\n".join([INTRODUCTION, RWG_BASIS, WIRE_GRID])
    return (f"Unknown topic '{topic}'. Available: introduction, rwg_basis, "
            "wire_grid, all.")
