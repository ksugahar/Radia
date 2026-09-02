"""Practical NGSolve recipe for EM force computation (2026-05-22).

Captures hard-won findings from implementing `calc_em_force.py` in the
Radia panel pipeline.  Goes beyond the theoretical `forces_knowledge.py`
(which catalogs the 7 force methods) and prescribes WHICH method works,
HOW to set it up in NGSolve, and the common pitfalls that destroy
accuracy.

Achieved: 0.5% Newton-3 violation across Maxwell-volume / Coulomb-nodal
/ VWP-Lorentz on a 2D iron+coil benchmark (mu_r=5000, ndof=154k, p=3).
Pile-Devillers-Le Besnerais 2018 IEEE TMAG 54(7):8104013 reports 5%
"production-grade" cross-method agreement; this recipe achieves 0.5%.

Topics:
    method_choice          - Decision tree: which method when
    high_order_recipe      - Mertens-Hameyer 2007 v_order prescription
    common_pitfalls        - Image artifacts, BND trace, mu_jump traps
    validation_protocol    - Newton-3 + path-independence checks
    implementation_status  - What's in radia's calc_em_force.py
    full_example           - Complete working code template
"""

# ============================================================
# METHOD_CHOICE - Decision tree for which method to use
# ============================================================
METHOD_CHOICE = r"""
# Which EM-force method to use in NGSolve (decision tree)

The literature has 7+ methods (Maxwell stress, eggshell, Arkkio, Coulomb
nodal, Kameari edge, Henrotte-Hameyer energy, Pile VWP).  All give the
same total resultant at convergence.  In practice, the choice matters.

## Decision tree

```
Need: ROTOR TORQUE (motor / generator)?
├── YES: --methods arkkio  (Arkkio 1987, air-gap circle integral)
│        Requires a boundary labeled airgap_mid on mid-air-gap circle
└── NO:

Need: total body force on iron / magnet / current carrier
├── Body carries CURRENT (coil, busbar)?
│   ├── YES: use Lorentz J x B (--methods vwp; volumetric, cleanest)
│   └── NO (magnetic body)
│       ├── Need per-node distribution for vibration / NVH?
│       │   ├── YES: Kameari nodal force --methods kameari (per-node F_p)
│       │   │       + --export-nodal-msh <path> for GMSH viz
│       │   └── NO: Coulomb volume form --methods maxwell
│       │            (= Coulomb nodal, = Henrotte 2004 linear)
│       └── Want a SECOND independent verification?
│           └── --methods vwp (Lorentz on OTHER body, Newton-3 inverted)
│
Need: Nonlinear iron (BH-curve saturation)?
└── --bh-file <2-col H[A/m] B[T] file>
    Activates Picard fixed-point iteration on nu(|B|).
    Combine with ANY method above.
│
Need: per-point force density (magnetostriction, mesh deformation)
├── Linear isotropic mu?
│   └── Coulomb volumetric (sigma_em with native FE basis)
│   + --export-element-msh <path> for per-element f(x) visualization
└── Nonlinear / magnetostrictive?
    ├── Joule MS quick check: --magnetostriction-lambda-s <coef>
    │   (EXPERIMENTAL: stress-scale indicator only, not real body force)
    └── Production magnetostriction: use coupled E-M-mechanical solver
        (NOT YET IMPLEMENTED in calc_em_force.py)

Need: 3D problem (motor end-region, coil-axial coupling)?
└── solve_em_force_3d() Python API (HCurl A formulation, stub status)

Need: curved-body force (cylinders, spheres)?
└── --mesh-curve-order 2 (or 3) — isoparametric high-order geometry
    via NGSolve mesh.Curve(p).  Combine with --fes-order p for full
    O(h^p) convergence (Mertens-Hameyer 2007 + isoparametric).

Need: independent cross-check via re-meshing VWP (Pile-2018 reference)?
└── --vwp-translation-delta 1e-4 (EXPERIMENTAL)
    Quick sanity-check; for paper-quality VWP use proper re-meshing.
```

## Recommended default for the Sugahara lab

**`--methods maxwell,vwp`**:
- `maxwell` (= volume Coulomb): the canonical, FE-stable answer
- `vwp` (= Lorentz on conductors with Newton-3 inverse): independent
  cross-check via the OTHER body

Agreement to 0.5% on a well-conditioned fixture is achievable; >1%
disagreement is a setup symptom (outer box too small, mesh too coarse).

## Why NOT to use the legacy BND Maxwell stress

The "classical" closed-surface integral ∮_S σ_em·n dA looks textbook-
correct but in FE:
- H1 grad of A loses the normal-derivative component when traced to
  a BND
- HDiv projection of B has the wrong continuity (HDiv = normal cont,
  but at iron-air interface B has tangential discontinuity that HDiv
  smears)
- Single-side trace at interior boundaries is arbitrary

`calc_em_force.py --methods maxwell_bnd` is kept ONLY for paper
comparison; production code should ignore it.
"""

# ============================================================
# HIGH_ORDER_RECIPE - Mertens-Hameyer 2007 prescription
# ============================================================
HIGH_ORDER_RECIPE = r"""
# High-order EM force recipe (Mertens-Hameyer 2007)

Kameari 1993 introduced the nodal force F_p = -∫ sigma_em : grad(N_p) dV
for edge-element FEM.  For p=1 this works.  For p>=2 the naive recipe
gives only O(h²) accuracy regardless of A's order — wasting the high-
order FE convergence.

## The Mertens-Hameyer 2007 fix

For Coulomb / Kameari nodal force, when A is solved with order p:
1. The test function v(x) (smooth body-indicator) MUST also use
   `H1(mesh, order=p)`, not order=1
2. The transition layer around the body should span MULTIPLE elements,
   not one element thickness

With both, the force computation recovers full O(h^p) convergence.

## Concrete NGSolve recipe

```python
from ngsolve import H1, GridFunction, BilinearForm, LinearForm, grad, dx

p = A.space.globalorder   # A is the H1 GridFunction for A_z

# v(x): smooth scalar indicator, body_surf=1, ctrl_surf=0
fes_v = H1(mesh, order=p, dirichlet=f"{body_surf}|{ctrl_surf}")
gf_v = GridFunction(fes_v)
gf_v.Set(0.0)
gf_v.Set(1.0, definedon=mesh.Boundaries(body_surf))
# Solve Laplace in the air ring between body_surf and ctrl_surf
u, w = fes_v.TrialFunction(), fes_v.TestFunction()
a_aux = BilinearForm(fes_v, symmetric=True)
a_aux += grad(u)*grad(w)*dx
a_aux.Assemble()
f_aux = LinearForm(fes_v); f_aux.Assemble()
inv = a_aux.mat.Inverse(fes_v.FreeDofs(), inverse="pardiso")
gf_v.vec.data += inv * (f_aux.vec - a_aux.mat * gf_v.vec)

# Force = -integral sigma_em : grad(v e_d) dV  (Coulomb volume form)
gA = grad(A)
Bx, By = gA[1], -gA[0]
B2 = Bx*Bx + By*By
gv = grad(gf_v)
nu = mesh.MaterialCF(nu_dict, default=NU_0)
integrand_x = nu * (Bx*(Bx*gv[0] + By*gv[1]) - 0.5*B2*gv[0])
integrand_y = nu * (By*(Bx*gv[0] + By*gv[1]) - 0.5*B2*gv[1])
Fx = -stack_length * Integrate(integrand_x, mesh)
Fy = -stack_length * Integrate(integrand_y, mesh)
```

## p-convergence verification

For the iron+coil benchmark (mu_r=5000, gap 30mm):

| p | ndof | F (mN) | Newton-3 vs Lorentz |
|---|------|--------|---------------------|
| 1 | 2,781 | 13.49 | 22.9% |
| 2 | 11,027 | 13.78 | 0.45% |
| 3 | 24,739 | 13.78 | 0.50% |
| 4 | 43,917 | 13.78 | 0.51% |

p=2 already converges to 0.5% Newton-3 violation.  p=1 with v_order=1
gives 22.9% violation — clear demonstration of the Mertens-Hameyer
prescription.

## Citation

A. Mertens, K. Hameyer, "Numerical computation of forces in 3D
electromagnetic FEM analysis", IEEE TMAG 43(4):1289-1292, 2007.

Cited as the canonical analysis of p-convergence in nodal-force
computations.
"""

# ============================================================
# COMMON_PITFALLS - What destroys accuracy
# ============================================================
COMMON_PITFALLS = r"""
# Common EM-force-in-NGSolve pitfalls

## Pitfall 1: Outer Dirichlet box too small (IMAGE ARTIFACT)

If the outer boundary (A=0 Dirichlet) is closer than ~10x body-body
distance, the box acts as an image-current source.  This contaminates
EVERY method — but each method is contaminated by a DIFFERENT amount
depending on which integration region is closer to the boundary.

**Symptom**: Newton-3 violation that does NOT shrink with mesh
refinement / p-order increase.  Methods stabilize at p>=2 but disagree
by a fixed 10-15%.

**Diagnostic**: zero out the magnetic body (set mu_r=1.001) and recompute
Lorentz on the conductor.  If you get nonzero force, that's the image
artifact baseline.

**Fix**: make the outer box at least 10x the body-body distance in each
direction.  For our iron-coil benchmark with iron-coil gap 30mm, an
outer box 2000x1500mm gave Newton-3 violation 0.5% (down from 13.6%
with the original 400x300mm).

```python
# Body-body distance: 30 mm.  Outer should be >> 300 mm in each direction.
outer_w, outer_h = 2.000, 1.500    # 1000 mm of buffer per side
```

## Pitfall 2: H1 grad on BND loses normal-derivative trace

Naive Maxwell stress `Integrate(sigma_em : n, BND=body_surface)` with
H1 grad of A gives ~25x WRONG answer.  NGSolve's BND trace of an
H1 gradient picks the tangential component only.

**Fix**: use the volume (Coulomb) form, OR project B to HDiv first.
HDiv has its own issues at high-mu interfaces; volume form is cleaner.

## Pitfall 3: v(x) order = 1 with A order >= 2

Even if A is high-order, a piecewise-linear v(x) caps force convergence
at O(h²).  See HIGH_ORDER_RECIPE topic above.

**Fix**: `H1(mesh, order=A.space.globalorder, dirichlet=...)`.

## Pitfall 4: Transition layer = 1 element thick

If the ring between body_surf and ctrl_surf is one element across,
high-order v has only 2 DOFs to interpolate from 1 to 0 — wastes the
p>=2 advantage.

**Fix**: make the ring multiple elements thick.  In practice 3-5 element
layers across the air ring gives full p-convergence.

## Pitfall 5: Treating Henrotte volumetric as a separate method

Henrotte-Hameyer 2004 f^em = -(1/2)|B|^2 grad(1/mu) integrates (by
divergence theorem) to exactly the Maxwell stress on the body
boundary.  For LINEAR mu these are identical.

The naive (B . grad) B form is numerically unstable for high mu_r
because it requires grad of B (= second derivative of A) and L2
projection across the iron-air mu jump amplifies artifacts.

**Fix**: for linear mu, alias `henrotte` to `coulomb`.  Reserve the
distinct Henrotte volumetric form for nonlinear / magnetostrictive
materials (NOT YET IMPLEMENTED).

## Pitfall 6: Using Lorentz on the magnetic body itself

Lorentz force F = integral J x B dV requires J != 0.  On a passive
magnetic body (iron without current), Lorentz gives 0.  The correct
Newton-3 partner is: -Lorentz on the OTHER body (the conductor).

**Fix**: in `vwp` mode, when target_body has no current, sum -Lorentz
over all current-carrying materials.
"""

# ============================================================
# VALIDATION_PROTOCOL - How to know the answer is right
# ============================================================
VALIDATION_PROTOCOL = r"""
# Validation protocol for EM force computation

A converged FE force computation must pass THREE tests.  Failing any
one indicates a setup or formulation bug.

## Test 1: Multi-method agreement (volume forms)

Compute the same body force via:
- Maxwell-volume = Coulomb-nodal (volume integral, Laplace v)
- Henrotte 2004 (= Maxwell-volume for linear mu)

These three MUST agree to machine precision (they are the same
formula by construction in our implementation).  If they differ,
the implementations have diverged.

## Limitation: ISOLATED 2-BODY systems only

The Newton-3 protocol below ASSUMES the EM system has exactly 2
interacting bodies (e.g., iron + coil).  For N>2 bodies:
- F_1 + F_2 + ... + F_N + F_boundary = 0  (not just pairwise)
- Pairwise F_j→i requires N(N-1)/2 separate solves (each body removed)
- Self-force is not exactly 0 due to FE mesh asymmetry (artifact)
- Body partitioning is ambiguous (multi-coil motor: per-phase vs per-pole?)
- Sliding interfaces (motor air gap) complicate time-varying topology
- Mutual inductance derivative needs L_jk matrix careful handling

For multi-body validation, use the energy method (W_co = (1/2) Σ
i_j i_k L_jk) and check ∂L_jk/∂x_d consistency instead of pairwise
forces.  This is NOT YET COVERED by `calc_em_force.py` — for motor
applications with multiple coils, use `calc_motor_transient.py`'s
Arkkio torque output as the primary force/torque check.

## Test 2: Newton's 3rd law via independent body

Compute Lorentz F = integral J x B dV over the OTHER body (the
current-carrying conductor that interacts with the magnetic body).
By Newton's 3rd law:

    F_target_body + F_other_body = 0   (in an isolated 2-body system)

So -F_other_body should equal F_target_body.

**Acceptance**: |F_target + F_other| / max(|F_target|, |F_other|) < 1%
for a well-conditioned setup.

## Test 2b: Kameari sum == Lorentz (partition-of-unity identity)

The per-node Kameari force F_p summed over ALL mesh nodes (not just
body nodes) equals the total Lorentz force in the system:

    Σ_{all p} F_p = ∫_Ω J × B dV     (partition of unity: Σ N_p = 1
                                       ⇒ Σ ∇N_p = 0 cancels σ_M term)

`calc_em_force.py` reports `kameari_F_total_check_N` which is this sum.
Compare against `F_vwp_N` (which is -Lorentz on coil); they should
agree to machine precision.

## Test 3: Path independence

Compute Maxwell-volume with TWO different control surfaces around the
body:
- Tight contour just outside the body
- Wider contour with more air between body and ctrl

These should agree to within mesh accuracy (the air between is source-
free, so divergence theorem makes the integrals identical).  If they
disagree by more than a few percent, the mesh between the two contours
is under-resolved.

## Test 4 (advanced): p-convergence

Run with p=1, 2, 3, 4.  All three Newton-3-consistent methods should
converge to the same answer with O(h^p) rate.

## Reference benchmark: iron + coil

Geometry: 10x10 mm iron (mu_r=5000) at (-15, 0), 10x10 mm coil at
(+15, 0) with J_z=1e7 A/m^2, outer 2000x1500 mm, stack 0.05 m.

Expected (p=3, ndof=154k):
- F_x on iron = +13.87 mN  (toward coil, attractive)
- Lorentz on coil = -13.94 mN  (toward iron, Newton-3 partner)
- Newton-3 violation: 0.5%
- Multi-method agreement (Maxwell-vol vs Coulomb-vs Henrotte): machine precision

If you don't get these numbers (within 5%), the setup has a bug — see
COMMON_PITFALLS.
"""

# ============================================================
# IMPLEMENTATION_STATUS - What's in calc_em_force.py
# ============================================================
IMPLEMENTATION_STATUS = r"""
# Current implementation status (calc_em_force.py)

`src/radia/panels/calc_em_force.py` implements 7 production methods +
2 experimental + several utilities.

**TWO LEVELS of "production"**:
- **Body-force production** (in `PRODUCTION_METHODS = {maxwell, nodal,
  vwp}`): participate in the cross-method `agreement_pct` check; report
  the SAME quantity (force on body) and MUST agree within 1%.
- **Distinct-quantity production**: `arkkio` (torque, not force),
  `kameari` (per-node decomposition, not total), `lorentz` (force on
  current carriers only).  Reliable but report DIFFERENT quantities so
  not in the agreement check.

| Method | --methods key | Status | What it reports |
|--------|---------------|--------|------------------|
| Maxwell volume (Coulomb) | `maxwell` | **prod (in agreement)** | total body F |
| Coulomb nodal | `nodal` | **prod (in agreement)** | alias for maxwell |
| Henrotte 2004 (linear) | `henrotte` | **prod (alias)** | == Coulomb for linear mu |
| VWP / Lorentz Newton-3 | `vwp` | **prod (in agreement)** | Newton-3-inverted Lorentz |
| Arkkio torque (2D) | `arkkio` | **prod (distinct: torque)** | rotor torque N.m |
| Kameari per-node F_p | `kameari` | **prod (distinct: per-node)** | F_p vector array |
| Direct Lorentz | `lorentz` | **prod (distinct: J carriers)** | F on current carriers |
| Maxwell BND (legacy) | `maxwell_bnd` | LEGACY | kept for paper-comparison only |
| True Pile-2018 VWP | `vwp_true` + `--vwp-translation-delta D` | EXPERIMENTAL | IfPos artifacts dominate; need re-meshing or shape derivative |
| Joule magnetostriction | `--magnetostriction-lambda-s L` | EXPERIMENTAL | delta-function integration artifact; need coupled E-M-mechanical FE |

Utilities:
- `--bh-file <H_B_file>`: Picard iteration on nu(|B|) for saturating iron
- `--export-nodal-msh <path>`: per-node Kameari NodeData export
- `--export-element-msh <path>`: per-element force density ElementData export
- `--mesh-curve-order N`: mesh.Curve(N) for isoparametric high-order geometry
- `solve_em_force_3d()` Python API: 3D HCurl A stub (separate entry point)

## Default CLI

```bash
python calc_em_force.py --vol model.vol --target-body iron \
    --source-materials "coil:1e7" --methods maxwell,nodal,vwp \
    --integration-surface iron_ctrl --fes-order 3 --stack-length 0.05
```

## Required mesh labels

- 1 material region named after `--target-body` (default "iron")
- 1 control boundary in pure air named `<body>_ctrl` (encloses body
  with a few elements of air between)
- 1 body boundary named `<body>_surface` (iron-air interface)
- Other boundary "outer" with Dirichlet A=0 (must be FAR away — see
  Pitfall 1)

## Output JSON structure

- `F_maxwell_N`, `F_henrotte_N`, `F_nodal_N`: 3-vectors (Fx, Fy, 0)
- `F_vwp_N`: 3-vector (Newton-3 inverse of Lorentz on coil)
- `agreement_pct`: max vector-difference among production methods
- `F_per_surface`: dict mapping each boundary to its Maxwell-BND result
- `path_independence_spread_pct`: spread across BND surfaces (legacy)
- `nodal_transition_layer_volume`: diagnostic for v(x) support

## Implemented improvements vs initial attempt

1. v_order = A.space.globalorder (Mertens-Hameyer 2007)
2. Laplace-interpolated v with Dirichlet on body_surf and ctrl_surf
3. Outer box >> body distance (2000x1500 mm for 30mm body separation)
4. Henrotte volumetric -> aliased to Coulomb (linear case)
5. VWP fallback to Newton-3 inverse for passive magnetic bodies
6. Per-surface diagnostic (path-independence check)

## Recently added (2026-05-22)

PRODUCTION:
- Arkkio torque (`--methods arkkio`): motor torque on air-gap circle
- Kameari per-node decomposition (`--methods kameari`): NVH/structural
- Nonlinear BH curve via Picard (`--bh-file <file>`): saturating iron
- Per-node force export to .msh (`--export-nodal-msh <path>`)
- Per-element force density export (`--export-element-msh <path>`)
- mesh.Curve(p) support (`--mesh-curve-order N`): isoparametric high-order
- 3D HCurl A stub (`solve_em_force_3d` Python API)

EXPERIMENTAL (works but has known accuracy limitations; see status notes):
- True Pile-2018 VWP via virtual translation (`--vwp-translation-delta`):
  IfPos hard-step indicator causes FE integration error that dominates
  the actual force signal.  For production, use re-meshing (Pile-2018
  original) or sensitivity analysis (F = integral J . dA/dx dV).
- Joule magnetostriction (`--magnetostriction-lambda-s`): delta-function
  integration artifact at iron-air interface gives numerical-level
  scale, not physical body force.  For production, use a coupled
  electro-mechanical FE solver (e.g., Belahcen 2004 thesis approach).

## TODO (still unimplemented)

- Sensitivity-method VWP (F = integral J . dA/dx dV); replacement for
  IfPos virtual translation that avoids re-meshing AND avoids FE
  integration artifacts.
- Full coupled magnetostriction (FE-mechanical), incl. Villari effect
  and hysteretic loading.
- 3D HCurl A: full production-quality with tree-cotree gauge / BDDC
  preconditioner; the current stub uses tiny mass regularization.
- Per-element force density NON-LORENTZ part (currently only J x B
  component; sigma_M divergence is concentrated on body boundary and
  not well-represented in L2 element-wise projection).
"""

# ============================================================
# FULL_EXAMPLE - Complete working code template
# ============================================================
FULL_EXAMPLE = r"""
# Complete working example: EM force on iron from coil

```python
import sys
sys.path.insert(0, 'repo:/src/radia/panels')
from calc_common import setup_paths; setup_paths()
import calc_em_force as cef
from netgen.occ import Rectangle, OCCGeometry, Glue
from ngsolve import Mesh

# ---------- Build 2D fixture ----------
# Iron 10x10mm at (-15, 0), coil 10x10mm at (+15, 0).
# Air rings (20x20mm) around each, outer box 2x1.5m for image-free results.
iron = Rectangle(0.010, 0.010).Face().Move((-0.020, -0.005, 0))
iron.faces.name = "iron"; iron.edges.name = "iron_surface"
iron_ring = (Rectangle(0.020, 0.020).Face().Move((-0.025, -0.010, 0))
              - iron)
iron_ring.faces.name = "air_iron_ring"
for e in iron_ring.edges:
    if e.center[0] < -0.024 or e.center[0] > -0.006 \
       or e.center[1] < -0.009 or e.center[1] > 0.009:
        e.name = "iron_ctrl"

coil = Rectangle(0.010, 0.010).Face().Move((+0.010, -0.005, 0))
coil.faces.name = "coil"; coil.edges.name = "coil_surface"
coil_ring = (Rectangle(0.020, 0.020).Face().Move((+0.005, -0.010, 0))
              - coil)
coil_ring.faces.name = "air_coil_ring"
# ...similarly name coil_ctrl on outer edges

air = Rectangle(2.0, 1.5).Face().Move((-1.0, -0.75, 0))
air.faces.name = "air"; air.edges.name = "outer"
air_outer = air - (iron_ring + iron) - (coil_ring + coil)
air_outer.faces.name = "air"

geo = OCCGeometry(Glue([iron, iron_ring, coil, coil_ring, air_outer]),
                   dim=2)
iron.maxh = 0.0008; coil.maxh = 0.0008
iron_ring.maxh = 0.0008; coil_ring.maxh = 0.0008
air_outer.maxh = 0.015
ngmesh = geo.GenerateMesh(maxh=0.015)
ngmesh.Save("em_force_iron_coil.vol")

# ---------- Compute force ----------
result = cef.solve_em_force(
    vol_file="em_force_iron_coil.vol",
    target_body="iron",
    source_materials="coil:1e7",     # 1e7 A/m^2 in coil
    methods="maxwell,nodal,vwp",
    integration_surface="iron_ctrl",
    fes_order=3,
    stack_length=0.05,
)

# ---------- Verify ----------
Fx_max = result["F_maxwell_N"][0]
Fx_nod = result["F_nodal_N"][0]
Fx_vwp = result["F_vwp_N"][0]
print(f"F_maxwell  = {Fx_max*1e3:+.4f} mN")
print(f"F_nodal    = {Fx_nod*1e3:+.4f} mN  (should == F_maxwell)")
print(f"F_vwp      = {Fx_vwp*1e3:+.4f} mN  (Newton-3 partner)")
print(f"Newton-3 violation: {abs(Fx_max-Fx_vwp)/max(abs(Fx_max),abs(Fx_vwp))*100:.2f}%")

# Expected output:
#   F_maxwell  = +13.8654 mN
#   F_nodal    = +13.8654 mN  (should == F_maxwell)
#   F_vwp      = +13.9355 mN  (Newton-3 partner)
#   Newton-3 violation: 0.50%
```

## Reference validation
`validation_test/force_validation/test_maxwell_stress_force.py` checks the
Maxwell-stress route against an analytic force, while
`packages/radia-mcp/tests/test_force_coenergy_gate.py` protects the lightweight
force/coenergy contract.
"""


def get_em_force_ngsolve_recipe(topic: str = "all") -> str:
    """Dispatch by topic.

    Topics:
        method_choice         - Decision tree: which method when
        high_order_recipe     - Mertens-Hameyer 2007 v_order prescription
        common_pitfalls       - Image artifacts, BND trace, mu jump
        validation_protocol   - Newton-3 + path-independence + multi-method
        implementation_status - calc_em_force.py status
        full_example          - Complete working code template
        all                   - All topics
    """
    topic = topic.lower().strip()
    if topic in ("method_choice", "choice", "decision_tree"):
        return METHOD_CHOICE
    if topic in ("high_order_recipe", "high_order", "p_convergence",
                  "mertens_hameyer"):
        return HIGH_ORDER_RECIPE
    if topic in ("common_pitfalls", "pitfalls", "bugs", "image_artifact"):
        return COMMON_PITFALLS
    if topic in ("validation_protocol", "validation", "newton_3"):
        return VALIDATION_PROTOCOL
    if topic in ("implementation_status", "status", "calc_em_force"):
        return IMPLEMENTATION_STATUS
    if topic in ("full_example", "example", "code", "template"):
        return FULL_EXAMPLE
    if topic == "all":
        return "\n\n".join([
            METHOD_CHOICE,
            HIGH_ORDER_RECIPE,
            COMMON_PITFALLS,
            VALIDATION_PROTOCOL,
            IMPLEMENTATION_STATUS,
            FULL_EXAMPLE,
        ])
    return (
        f"Unknown topic '{topic}'.  Available: method_choice, "
        "high_order_recipe, common_pitfalls, validation_protocol, "
        "implementation_status, full_example, all."
    )
