"""Additional EM force topics (2026-05-22):
  - permanent_magnet_force (M-fixed body, bound surface current)
  - energy_method_classical (constant flux vs constant current)
  - sensitivity_vwp (proper VWP via shape derivative / sensitivity equations)
  - lorentz_canonical (Lorentz as a first-class canonical method)
  - meissner_force (Type-I superconductor, mu_r -> 0)

These complement the 7-method catalog in forces_knowledge.py and the
practical recipe in em_force_ngsolve_recipe_knowledge.py.
"""

PERMANENT_MAGNET_FORCE = r"""
# Permanent magnet (PM) force computation

## Distinction from soft iron

| Soft iron | Permanent magnet |
|-----------|------------------|
| M induced from H via nonlinear M(H) | M fixed by material (residual mag.) |
| Needs Solve() to find M | No Solve needed for M |
| Force ∝ grad(B²)/μ (Kelvin) | Force ∝ (M · ∇)B (Kelvin on dipole) |
| F vanishes if H = 0 | F nonzero in any external field |

## Three equivalent PM force formulations

### (A) Bound-current model
A PM with magnetization M is equivalent to a body with:
  - Volume bound current density J_M = curl(M)
  - Surface bound current density K = M × n
For UNIFORM M: curl(M) = 0 inside, K = M × n on boundary only.

Then Lorentz force: F_PM = ∮_PM_surf (K × B) dA

### (B) Magnetic-dipole / Kelvin form
For a body with magnetization M(x):
  F_PM = ∫_PM (M · ∇) B_ext dV
where B_ext is the field from EVERYTHING ELSE (other PMs + coils + iron).

For UNIFORM M and slowly-varying B_ext:
  F_PM ≈ V_PM · (M · ∇) B_ext  (at PM center)

### (C) Equivalent surface charge (scalar potential)
For magnetostatic with no currents in PM region:
  rho_M = -div(M) (volume magnetic charge density)
  sigma_M = M · n   (surface magnetic charge)
Then F_PM = ∫ rho_M H_ext dV + ∮ sigma_M H_ext dA

## Sign convention

For a bar magnet aligned in +x in an external field gradient:
  - If grad(B_x) > 0 in +x direction: F_PM is +x (toward stronger field)
  - This is the "magnetic dipole attracted to gradient" rule
  - Same direction as soft iron in linear regime (both align with gradient)

## NGSolve implementation in calc_em_force.py

```bash
python calc_em_force.py --vol pm_assembly.vol \
    --pm-magnetization "magnet1:954930,0;magnet2:0,954930" \
    --source-materials "coil:1e6" \
    --methods maxwell,lorentz,vwp \
    --fes-order 3
```

- `--pm-magnetization "mat:Mx,My"`: M in A/m (e.g., 954930 = 1.2 T / mu_0)
- Multiple PMs separated by `;`
- Output: `F_pm_per_material` dict + `F_pm_target_N` for target body

The implementation uses formulation (A) bound surface current K = M × n
integrated on PM body boundary.  Uniform-M assumption inside body.

## When to use which method

- **Method (A) bound current**: cleanest for FE; surface integral on
  PM boundary; works for any geometry
- **Method (B) Kelvin form**: simplest for analytical estimates; requires
  grad(B_ext), so 2nd derivative of A in FE
- **Method (C) scalar potential**: useful when phi-formulation is used
  (Omega-reduced, FEMM-style)

## Reference

- Coey, "Magnetism and Magnetic Materials", Cambridge 2010 (Ch 2)
- du Trémolet de Lacheisserie et al., "Magnetism Vol 1", Springer 2005
"""


ENERGY_METHOD_CLASSICAL = r"""
# Energy method: constant flux vs constant current

The CLASSICAL textbook treatment (Bewley 1948, Hammond 1986) starts here.

## The two thermodynamic potentials

For a magnetostatic system with current sources i and flux linkages λ:

| Quantity | Definition | Independent var. |
|----------|------------|------------------|
| Magnetic energy W | ∫_0^λ i(λ', x) dλ' | flux λ |
| Co-energy W'     | ∫_0^i λ(i', x) di' | current i |
| Linear case      | W = W' = (1/2) Σ_jk L_jk i_j i_k | |

The Legendre transform: W + W' = Σ_k i_k λ_k.

## Force formulas (the SAME physical force, different bookkeeping)

For a body moving by dx_d in direction d:

| Constraint | Formula | Why |
|-----------|---------|-----|
| Flux constant (λ fixed) | F_d = -∂W/∂x_d \|_λ | external sources cannot exchange energy |
| Current constant (i fixed) | F_d = +∂W'/∂x_d \|_i | sources DO exchange energy at constant i |

Sign FLIPS between the two!  Both give the same F_d.

For linear systems (most FE solvers):
  F_d = (1/2) i_j i_k ∂L_jk/∂x_d   (mutual-inductance form)

## When to use which (practical guide)

**Constant current (use W', co-energy)** — DEFAULT in FE:
- Stator winding driven by a current source
- Field computed for given i
- ∂W'/∂x_d via finite difference: shift body, solve, finite-diff W'
- Easy to combine with circuit coupling (i known)

**Constant flux (use W, energy)** — RARE in FE:
- Superconducting coil (flux frozen by Meissner effect)
- Open-circuit measurement (flux preserved)
- Implementation: needs Lagrange multiplier on ∫ B·n dA = const
- Sometimes preferred for short-circuit transient simulations

## FE implementation of co-energy method

For 2D A_z formulation with current source J:
  W' = (1/2) ∫ A · J dV     (linear case)
  F_d = +∂W'/∂x_d via:
    1. Solve A(x_d=0)        → W'_0
    2. Solve A(x_d=+delta)   → W'_+
    3. Solve A(x_d=-delta)   → W'_-
    4. F_d ≈ (W'_+ - W'_-) / (2 delta)

In `calc_em_force.py`:
- `--methods vwp` reports the co-energy diagnostic W' (no shifting)
- `--methods vwp_true` implements (1)-(4) via IfPos virtual translation
  (EXPERIMENTAL — FE artifacts dominate)
- `--methods sensitivity_vwp` implements (4) via the sensitivity
  equation K (dA/dx) = -(dK/dx) A (EXPERIMENTAL)

## Common confusion: which sign?

Students often miss that:
  F_d > 0 (force in +d direction) means body MOVES toward INCREASING W' at constant i.
  F_d > 0 means body MOVES toward DECREASING W at constant λ.

For an iron piece attracted toward a coil:
  - Co-energy W' INCREASES as iron approaches coil  → F_attract = +∂W'/∂x > 0 ✓
  - Energy W DECREASES as iron approaches (sources do positive work)  → F_attract = -∂W/∂x > 0 ✓

Both formulas agree.  The sign confusion comes from the constraint.

## Reference

- Bewley, "Two-Dimensional Fields in Electrical Engineering", 1948 (Ch 7)
- Hammond, "Energy Methods in Electromagnetism", Oxford 1986
- White-Woodson, "Electromechanical Energy Conversion", 1959 (Ch 3)
"""


SENSITIVITY_VWP = r"""
# Sensitivity-method VWP — the PROPER replacement for the buggy IfPos translation

## The problem with naive virtual-translation

Naive Pile-2018 VWP:
  F_d = +(W_co(+delta) - W_co(-delta)) / (2 delta) | at constant i

Implementation in calc_em_force.py: shift body's material region via
IfPos hard-step indicator, re-solve A 5 times (0, ±delta_x, ±delta_y),
finite-diff W_co.

**THE BUG**: IfPos crosses mesh element boundaries.  The FE integration
error depends on HOW the iron-air interface aligns with the mesh.
This error (O(0.1 J) typically) is much LARGER than the actual
force-induced W_co change (O(1e-4 J)).  → numerical garbage.

## The fix: sensitivity equation, not energy FD

Instead of finite-differencing W_co (a SCALAR), differentiate A:
  K(x) A(x) = f(J)
  ⇒ K (dA/dx_d) = -(dK/dx_d) A

Then:
  F_d = +(1/2) ∫ J · (dA/dx_d) dV

This requires:
  1. Build dK/dx_d analytically (shape derivative of stiffness matrix)
  2. Solve sensitivity equation K (dA/dx_d) = rhs_sens
  3. Integrate J · dA over coil

Advantage: dK/dx_d is a SHAPE DERIVATIVE that can be computed via
boundary integral exactly, no FE noise:

  dK/dx_d  applied to (u, v) = -(nu_iron - nu_air) ∫_iron_surf
                                  n_d (grad u · grad v) dA

So:
  rhs_sens = +(nu_iron - nu_air) ∫_iron_surf n_d (grad u_test · grad A) dA

This is a CLEAN BOUNDARY INTEGRAL, no IfPos hard-step needed.

## Equivalent direct formula (Pironneau-Allaire)

The full shape derivative can be done analytically:
  F_d = -(1/2) A^T (dK/dx_d) A
      = +(1/2)(nu_iron - nu_air) ∫_iron_surf n_d |B|^2 dA

This is THE PRODUCTION-GRADE VWP for linear isotropic μ.
Compare with Maxwell stress on iron surface:
  F_d_maxwell = (1/(2 mu_0)) ∫_iron_surf [2 (B·n) B_d - |B|^2 n_d] dA

The shape-derivative form has ONLY the (1/2) |B|² n_d term, NOT the
(B·n) B_d term.  Mathematical reason: shape derivative computes the
energy gradient, not the stress integral; the (B·n) B term integrates
to 0 for a closed surface (∫ (B·n) B dA = 0 by divergence theorem on
divergence-free B).

## NGSolve implementation status (2026-05-22)

`--methods sensitivity_vwp` in calc_em_force.py:
- Uses FD on bbox shift to approximate dK/dx_d
- Returns 0 if shift < element_size (FE noise floor)
- Production version SHOULD use the analytical shape derivative
  boundary form (NOT YET IMPLEMENTED)

TODO: replace FD with analytical shape derivative.

## Reference

- Pironneau, "Optimal Shape Design for Elliptic Systems", Springer 1984
- Allaire-Jouve-Toader 2004 (level-set shape opt for elliptic problems)
- Aitchison-Mokhtarzadeh-Wexler 2009 (shape derivative for E&M)
"""


LORENTZ_CANONICAL = r"""
# Lorentz force as a first-class canonical method

## Why Lorentz deserves its own category

The Lorentz force on currents:
  F_Lorentz = ∫_conductor J × B dV

is often LUMPED with VWP/co-energy ("they're equivalent at constant i").
But operationally, Lorentz has unique properties:

1. **No body mu_r needed**: works for ANY current carrier (coil in air,
   PCB trace, plasma loop).
2. **Volume integral, no boundary trace**: immune to BND projection
   issues that affect Maxwell stress on iron surface.
3. **Independent cross-check**: by Newton's 3rd law, sum over all
   bodies = 0; checking Lorentz on coil + Maxwell on iron should give
   zero sum.
4. **Decomposable per-conductor**: F_each_coil is well-defined and
   doesn't require partitioning ambiguity.

## When Lorentz is the RIGHT method

- **Force on a coil** (single or array): always Lorentz, no ambiguity
- **Force on a busbar / power-electronics conductor**: Lorentz
- **Plasma stability**: Lorentz on plasma current
- **Magnetorheological fluid / ionic liquid**: Lorentz on charged carriers
- **Electric field response of conductors**: Lorentz + electric force

## When Lorentz is the WRONG method

- **Force on a passive iron body** (no current): Lorentz gives 0
  → use Maxwell volume / Coulomb / Henrotte instead
- **Force on a permanent magnet** (no free current): Lorentz on free
  current gives 0 → use bound-current Lorentz (= surface integral
  of (M × n) × B) or Kelvin (M · ∇) B form

## Equivalence with co-energy / VWP

For a single coil in a fixed external field:
  W_co_coil = (1/2) L i² + i ψ_ext   (linear case)
  F = +∂W_co/∂x | i const = +i ∂ψ_ext/∂x
       = i · ∂(∫ A·dl)/∂x          (Stokes: ψ = ∮ A·dl)
       = ∫ (i × B) · dl             (after some algebra)
       = ∫ J × B dV                  (volume form for distributed J)

So Lorentz = co-energy gradient for current-driven systems.
The reason BOTH should be reported in cross-check: they expose
DIFFERENT FE error modes.

## NGSolve implementation in calc_em_force.py

`--methods lorentz`: direct F = ∫ J × B dV on the target body.
Per-material breakdown reported as `F_lorentz_per_material`.
Returns (0, 0) if target_body has no current.

`--methods vwp`: same Lorentz but with Newton-3 inversion when
target_body is passive (sums -F on other bodies).

Use BOTH for cross-validation: Lorentz on a current-carrying body +
Maxwell on the passive partner should give equal-and-opposite forces.

## Reference

- Griffiths "Introduction to Electrodynamics" 4th ed, §8.2.2
- Jackson "Classical Electrodynamics" 3rd ed, §6.7 (Maxwell stress)
"""


MEISSNER_FORCE = r"""
# Meissner force (Type-I superconductor, mu_r → 0)

## The physics

A Type-I superconductor expels magnetic field from its interior
(Meissner-Ochsenfeld 1933).  This is equivalent to:
  - Perfect diamagnet: mu_r → 0
  - Perfect conductor in flux-frozen state
  - Surface persistent current that exactly cancels external B inside

For force calculation, the superconductor experiences a REPULSIVE
force from any external B source — opposite to soft iron's attraction.

## Image-current method (analytical)

For a Type-I SC half-plane at x=0 and a current I at distance d:
  Image current: -I at distance d on the SC side (antiparallel)
  Force per unit length: F = +mu_0 I² / (4 pi · 2d)   (REPULSIVE)

Sign: the SC pushes the current AWAY from itself (sign opposite to
soft-iron half-plane which ATTRACTS).

## NGSolve implementation in calc_em_force.py

```bash
python calc_em_force.py --vol sc_levitation.vol \
    --target-body iron --meissner --methods maxwell,nodal,vwp
```

The `--meissner` flag overrides `--mu-r-iron` to 1e-6, making the
body act as a perfect diamagnet (B ≈ 0 inside).  All standard methods
(maxwell, nodal, vwp) work; the result has REVERSED sign vs
ferromagnetic case (repulsion not attraction).

## Validation: levitation force

For a small magnet levitating over a SC plane:
  F_levitation = magnet weight  (force-balance equilibrium height)

Test setup:
  - PM material with M=954930 A/m (1.2 T equivalent)
  - SC body with mu_r ≈ 0 (--meissner flag)
  - Measure F_PM as function of height → equilibrium where F = mg

## Limitations

The mu_r → 0 model captures ONLY the Meissner state.  Real SCs have:
  - Type-II vortex penetration above H_c1 (need flux-pinning model)
  - Critical current J_c (Bean model for force in flux-pinned state)
  - Temperature-dependent λ_London penetration depth

For HTS levitation (most modern applications), use the Bean model:
  - Constant J_c inside vortex-penetrated region
  - F = ∫ J_c × B dV (Lorentz on critical current)

This is NOT yet implemented; use mu_r → 0 only for Meissner-state
Type-I or rough HTS first-cut.

## Reference

- Hellman et al., "Magnetic levitation by high-Tc superconductors",
  Rev. Mod. Phys. 1995
- Brandt, "Levitation in physics", Science 1989
- Tinkham, "Introduction to Superconductivity" 2nd ed, McGraw-Hill 1996
"""


def get_em_force_extras(topic: str = "all") -> str:
    """Dispatch by topic.

    Topics:
        permanent_magnet_force  - PM force formulations (bound current,
                                    Kelvin, scalar charge)
        energy_method_classical - Constant flux vs constant current,
                                    Legendre duality, FE recipes
        sensitivity_vwp         - Proper shape-derivative VWP (replaces
                                    IfPos buggy translation)
        lorentz_canonical       - Lorentz as a first-class method,
                                    decomposable per-conductor
        meissner_force          - Type-I superconductor mu_r → 0,
                                    image method, levitation
        all                     - All five
    """
    topic = topic.lower().strip()
    if topic in ("permanent_magnet_force", "pm", "permanent_magnet",
                  "magnet"):
        return PERMANENT_MAGNET_FORCE
    if topic in ("energy_method_classical", "energy_method", "energy",
                  "co_energy", "flux_vs_current"):
        return ENERGY_METHOD_CLASSICAL
    if topic in ("sensitivity_vwp", "sensitivity", "shape_derivative",
                  "pironneau"):
        return SENSITIVITY_VWP
    if topic in ("lorentz_canonical", "lorentz", "jxb"):
        return LORENTZ_CANONICAL
    if topic in ("meissner_force", "meissner", "superconductor",
                  "type_i", "diamagnet"):
        return MEISSNER_FORCE
    if topic == "all":
        return "\n\n".join([
            PERMANENT_MAGNET_FORCE,
            ENERGY_METHOD_CLASSICAL,
            SENSITIVITY_VWP,
            LORENTZ_CANONICAL,
            MEISSNER_FORCE,
        ])
    return (
        f"Unknown topic '{topic}'. Available: permanent_magnet_force, "
        "energy_method_classical, sensitivity_vwp, lorentz_canonical, "
        "meissner_force, all."
    )
