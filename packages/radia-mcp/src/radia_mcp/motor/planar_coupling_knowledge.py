"""radia_mcp.motor knowledge: the 2D PLANAR machine-modelling stack in radia --
MMMM / HDiv-VIM soft-iron demag + the SHARED planar postprocessing + the staggered
eddy-current coupling (PM motor / induction machine / eddy-current brake).

Public-safe, analytic-led (Bessel / dipole / monolithic-FEM gated).  This is the queryable
companion to the goldens under validation_test/feec/test_moment2d_*.py,
test_planar_eddy_coupling.py, test_hdiv_vim_2d_magnets.py.
"""
from __future__ import annotations

SECTIONS: dict[str, str] = {
    "overview": """\
# 2D planar machine modelling in radia (per-unit-length motor cross-section)

Two soft-iron demag methods, ONE shared postprocessing + coupling layer:

- **collocation MMMM** -- `radia.mmmm2d` (C++ `rad_moment2d`): per-element edge line-charge DOF,
  the constitutive law imposed on the field MOMENTS about each centroid (1 monopole + 2 dipole +
  (nEdge-3) quad rows).  2D Laplace kernel G = -ln(r)/(2 pi).  Dense LU (few-element method).
- **HDiv-VIM** -- `radia.vim._vim2d` (RT1 charge Gram): loop-free H(div) flux formulation.

Both produce a per-element magnetisation M, so the postprocessing + coupling is METHOD-AGNOSTIC:

- **`radia.planar_charges`** -- shared exterior field / A_z / complex (eddy phasor) field / Maxwell
  torque + force / magnet_field (the M.n edge-charge cloud, log kernel; C++ hot loops).
- **`radia.planar_eddy`** -- shared staggered eddy-current coupling (maglev / IM / ECB).

Use MMMM for per-element accuracy on a few-element cross-section; HDiv-VIM for loop-free convergence.
See `planar_2d` on the mmm_core tool for the MMMM element formulation itself.
""",
    "eddy_coupling": """\
# Staggered eddy-current coupling -- radia.planar_eddy (maglev / IM / eddy-current brake)

Weak (sequential) coupling of the ANALYTIC soft-iron demag (MMMM or HDiv-VIM) with an NGSolve
REDUCED-POTENTIAL complex A_z eddy FEM in a SEPARATE conductor.  The lab maglev method (Chadebec 2006
IEM open boundary + Biro 2000 reduced potential): the iron field is analytic (NO air mesh -- a magnet
may move without re-meshing); NGSolve solves ONLY the eddy reaction in the conductor.

- METHOD-AGNOSTIC: `iron_solve` is a callback `H_ext_complex (nEl,2) -> M (nEl,2)` -- MMMM and the
  HDiv-VIM both fit.
- The iron field is injected as the shared M.n log-charge cloud rendered as an NGSolve **atan2 CF**
  (exact, no interpolation).  Branch cut: place the iron so the conductor sees it from one side.
- The eddy system matrix (grad.grad + j w mu0 sigma mass) is iteration-INDEPENDENT -> assembled +
  factored ONCE (pardiso), RHS-only per staggered step.  Typically ~4 staggered iters.
- COMPLEX iron (chi real): M = solve(Re H) + j solve(Im H).

Gates: standalone eddy == analytic conducting-cylinder Bessel  <Bx>/B0 = 2 I1(z)/(z I0(z)),
z = (1+j) a/delta (NOT the center value 1/I0(z) -- the interior field is non-uniform at finite freq,
match the AVERAGED quantity); the full staggered solve reproduces a MONOLITHIC FEM to 6e-4 (iron
M_avg) / 1.6e-4 (conductor <Bx>).
""",
    "pm_motor": """\
# Permanent magnets: separate body (design A) vs embedded region (design B), unified rotor

A hard PM has a RIGID magnetisation -> a ONE-WAY source (it does not demagnetise; no iteration): its
field (shared `exterior_field` / `magnet_field`) is added to the applied field the soft iron reacts
to.  Two placements:

- **design A (separate body)** -- `solve_planar_demag(iron, mu_r=, magnets=[(pm_mesh, M_fixed)])`
  (MMMM) or `hdiv_demag_solve(iron, mu_r, H_ext, magnets=[...])` (HDiv-VIM).  The PM is a SEPARATE mesh.
- **design B (embedded region)** -- `solve_planar_demag(mesh, mu_r={"iron": ..}, pm={"pm": [Mx,My]})`
  (MMMM): a PM SEGMENT is a REGION of the SAME mesh as the iron (a real PM-motor rotor: magnets
  embedded in the iron).  The mesh is partitioned soft/hard; only the soft subsystem is solved =
  design A on one partitioned mesh (fields superpose; no new C++).  (design B for HDiv-VIM is a
  follow-up -- it needs a soft/hard PlanarDemagBody partition.)

**Unified rotor (PM + iron + eddy)** -- `planar_eddy.couple_mmmm(rotor, fem, sigma, freq, mu_r=,
pm={..})`: a rotor = soft iron + embedded PM coupled to a conductor eddy (PM-motor / eddy-current
brake).  The PM is treated as an in-phase phasor source at omega (moving-rotor / ECB convention): its
field magnetises the iron AND drives the conductor eddy.  Gated vs a monolithic AC+PM FEM (2.9e-3).
""",
    "nonlinear": """\
# Nonlinear soft iron + eddy (effective-chi AC)

`couple_mmmm(..., bh_table=[[H,B],..])`: nonlinear soft iron in the staggered eddy coupling via an
amplitude-based EFFECTIVE-chi Picard -- the 1st-harmonic AC approximation:
    chi_eff = M(|H|) / |H|,   |H| = sqrt(|Hx|^2 + |Hy|^2)  (the PHASOR magnitude).
|H| reduces to |H| for a real / DC field, so the sigma->0 limit recovers the DC nonlinear demag.
It captures amplitude-dependent saturation, NOT harmonic generation (that needs time stepping).  The
chi state warm-starts across staggered calls.

Gates: sigma->0 recovers the standalone DC nonlinear MMMM demag (2e-3); a low drive recovers the
linear chi0 result (1e-3).  DC nonlinear demag itself: MMMM/HDiv scalar-chi Picard + safeguarded
Anderson(1), M = M_of(H0 - D M) fixed point (disk D=1/2).
""",
    "api": """\
# API quick reference (all inside `with ngsolve.TaskManager():`)

    import radia.mmmm2d as m2, radia.planar_charges as pc, radia.planar_eddy as pe

    # soft-iron demag (linear or nonlinear; scalar or {region: value} dict)
    r = m2.solve_planar_demag(mesh, mu_r=1000.0, H_ext=(H0, 0.0))          # or bh_table=[[H,B],..]
    r = m2.solve_planar_demag(mesh, mu_r={"stator": 4000, "rotor": 2000}, H_ext=cf)

    # permanent magnets: separate body (A) or embedded region (B)
    r = m2.solve_planar_demag(iron, mu_r=200, H_ext=(0,0), magnets=[(pm_mesh, M_fixed)])
    r = m2.solve_planar_demag(mesh, mu_r={"iron": 200}, H_ext=(0,0), pm={"pm": [Mx, My]})

    # shared postprocessing (either method's M feeds these)
    H  = pc.exterior_field(mesh, r["M"], P)         # H at points P (n,2), in air
    T  = pc.maxwell_torque(mesh, r["M"], Rc, H_ext=(H0,0))     # reluctance torque, unit length
    F  = pc.force_between([(mesh_a, Ma), (mesh_b, Mb)], Rc, center)     # inter-body (maglev)
    sw = m2.torque_angle_sweep(mesh, H0, angles, Rc, mu_r=...)  # LINEAR: matrix factored ONCE

    # staggered eddy coupling (maglev / IM / PM-motor / eddy brake)
    res = pe.couple_mmmm(iron_or_rotor, fem_mesh, sigma=, freq=, mu_r=|bh_table=, pm=, B0=)
    res["M"]      # complex per-element M      res["gfu"]   # eddy reaction potential A_r
""",
    "validation": """\
# Validation ladder (all gated, self-contained; validation_test/feec/)

- MMMM element: disk D=1/2, ellipse a:b -> Dx=b/(a+b) (0.005%), chi-sweep, quad==triangulated.
- MMMM vs HDiv-VIM vs analytic scattered dipole (linear cylinder): M 0.000%, scattered H 0.061%.
- Maxwell torque circle == mu0 A (M x H0); force: uniform field -> 0, Newton F(A)=-F(B).
- Eddy FEM == analytic conducting-cylinder Bessel 2 I1(z)/(z I0(z)) ~1e-4.
- Staggered coupling == monolithic FEM: 6e-4 (iron) / 1.6e-4 (conductor).
- Unified PM+iron+eddy rotor == monolithic AC+PM FEM: 2.9e-3.
- PM: rigid disk exterior == 2D dipole a^2 M/(2 r^2); design B == design A (1e-9); PM+iron ==
  monolithic magnetostatic FEM (2e-2); HDiv-VIM magnets= == MMMM magnets= (3e-2).
- Nonlinear+eddy: sigma->0 -> DC nonlinear demag (2e-3); low-drive -> linear chi0 (1e-3).
""",
}


def get_planar_coupling(topic: str = "overview") -> str:
    """2D planar MMMM / HDiv-VIM demag + shared eddy/PM coupling knowledge."""
    t = (topic or "overview").strip().lower()
    if t in ("all", "*"):
        return "\n\n---\n\n".join(SECTIONS[k] for k in SECTIONS)
    if t not in SECTIONS:
        valid = ", ".join(sorted(SECTIONS))
        return "Unknown topic %r. Valid topics: %s (or 'all')." % (topic, valid)
    return SECTIONS[t]
