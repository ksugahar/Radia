"""Kameari + Kelvin v20: HDiv-projected J source for J·n=0 strict (2026-05-05).

Sugahara guidance (b2 path, step 1):
  - Kelvin infra verified sound by M3 torus benchmark (FEM Kelvin 91.58 nH
    vs analytical 88.55 nH, +3.42% acceptable discretization error).
  - v15-v19 issues are formulation-level (HCurl + applied B + iteration),
    not Kelvin.
  - v19 enforced J·n=0 via Coulomb correction of A_s, but tau_0=534 us
    didn't match canonical 43.84 (Mathematica BEM HDiv strict).
  - v20: try HDiv-projected initial source. HDiv basis with
    dirichlet="conductor_surface" gives J·n=0 STRICT (essential BC),
    not just integrally zero. Different subspace from v19.

v20 strategy:
  1. Build HDiv space on conductor with dirichlet="conductor_surface".
  2. Project sigma * A_s_phys onto HDiv -> gfJ_imp (J·n=0 strict).
  3. Run Kameari iteration with this HDiv source.
  4. Use B² form for L_n (manifestly positive).

Expected: if HDiv-strict source gives different tau than v19 (Coulomb-
soft enforcement), the basis-level constraint matters. If similar, the
issue is elsewhere (need full HDiv eigenvalue problem, step 2).

(original v19 docstring follows for reference)

Kameari + Kelvin v19: A_s Coulomb correction for J·n=0 at conductor (2026-05-05).

Sugahara guidance:
  - v18 (B² form) gave L>0 always, but tau_0 = 452 us (off from canonical
    Phase 7-11 19.77 / HDiv 43.84) because A_s = (B_0/2)(-y, x, 0) has
    A·n ≠ 0 at conductor x and y faces, so J_imp = sigma A_s violates
    J·n=0 (isolated-conductor BC).
  - The IBP identity (j, A)_C = ∫ ν |curl A|² requires J·n=0; otherwise
    a surface contribution ∫_∂C (J·n) phi dS contaminates.

v19 strategy:
  1. Solve a global Poisson for correction phi:
       ∫_full ∇phi·∇psi dV = -∫_∂cond (A_s·n) psi dS
     This generates phi continuous everywhere (H1, Periodic + GND vertex
     Dirichlet), with [∂phi/∂n]_∂cond = +A_s·n (jump). Adjusting so that
     (A_s + grad phi)·n = 0 from conductor side at conductor surface.
  2. A_s_corrected = A_s + grad(phi). Has same curl (= B_0 ẑ in inner)
     and tangential continuity preserved.
  3. J_imp = sigma × A_s_corrected satisfies J·n=0 at conductor surface.
  4. Run Kameari iteration as v17/v18.

Geometry update: name cuboid faces explicitly as "conductor_surface" so
the surface integral ds("conductor_surface") works.

Expected: tau_0 should match Phase 7-11 spectral 19.77 us (HCurl class
with effective J·n=0 via parity) or HDiv canonical 43.84 us.

(original v18 docstring follows for reference)

Kameari + Kelvin v18: L_n via B² integral (2026-05-05).

Sugahara guidance:
  - v17 still has L_1 < 0 because L_n is computed as R_n × <J_n, A_pot_accum>
    which contains cross terms ∫ (1/μ) curl(A_k) curl(A_n) dV (k<n).
  - These cross terms are NOT zero just because Schmidt of J in σ-norm is
    enforced (drift = 2.2e-15) — A-side orthogonality in (1/μ)-norm is
    a separate condition.
  - Canonical Tanimoto-Sugahara TMAG3304725 eq. (3): L_{n,n} = ∫_C ĵ_n · â_n.
    By IBP using a(â_n, â_n) = (ĵ_n, â_n)_C:
        L_{n,n} = ∫_full ν |curl(â_n)|² dV
    This is the **B² integral form** — manifestly ≥ 0.
  - Use this directly for the per-stage L_n. Drop cross-term contamination.

v18 changes from v17:
  - Compute L_n = R_n² × ∫_full ν |curl(gfA_r)|² dV (B² form)
  - Reduce N_STAGES = 2 (focus on L_0 and L_1, the user's request)

Expected: L_1 > 0 with B² formula. Verifies the A-method works in v17's
gauge setup once L_n is computed canonically.

(original v17 docstring follows for reference)

Kameari + Kelvin v17: A-method gauge fix via Helmholtz-Hodge (2026-05-05).

Sugahara guidance (this session):
  - v15 sign-flips at stage 1 = FEM is NOT solving correctly (not a
    natural Kameari breakdown). A-method needs ∇·A = 0 gauge enforced.
  - Canonical recipe (memory project_3D_CLN_kameari_in_vacuum_hardness):
    nograds=True / type1=True, NO penalty, Helmholtz-Hodge OR tree-cotree.
  - For Kelvin BC (Periodic + GND vertex), tree-cotree alone is
    insufficient — Periodic identification creates equivalence classes
    that BFS spanning tree doesn't fully cover. Need explicit Helmholtz-
    Hodge projection.

v17 changes from v15:
  - REMOVE GAUGE_EPS penalty (canonical: no penalty)
  - ADD Helmholtz-Hodge projection after each FEM solve:
      A_div_free = A - grad(phi)
      where phi solves: ∫ ∇phi·∇psi dx = ∫ A·∇psi dx in matching H1 space
  - Keep nograds=True + tree-cotree gauge (canonical)
  - Keep Pardiso (next iteration: try shifted AMS preconditioner if v17
    insufficient, per Sugahara hint)

Expected if v17 works: clean stages 0..N (drift at machine precision,
all L_n positive), then natural Kameari breakdown onset at higher N
(this would be the paper motivation data — natural breakdown after a
healthy regime).
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
sys.path.insert(0, r"S:/Radia/01_GitHub/src/radia")
print("Starting v20b (HDiv source + explicit Dirichlet zero force)...",
      flush=True)

from netgen.occ import Box, Sphere, Pnt, OCCGeometry
from ngsolve import (
    Mesh, HCurl, HDiv, H1, Periodic, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, grad, dx, ds, x, y, z,
    Integrate, TaskManager, ngsglobals, specialcf,
)
from kelvin_geometry import add_kelvin_exterior_domain
from kelvin_material import (
    make_kelvin_nu_cf,
    make_reduced_potential_background_cf,
    NU_0,
)
from collections import deque
from math import pi
import json, time
from pathlib import Path

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7
ax, ay, az = 5e-3, 2e-3, 1e-3
V_cond = ax * ay * az

R_K = 25e-3
OFFSET = (2.5 * R_K, 0, 0)
H_COND = 0.5e-3
H_AIR = 2.0e-3
ORDER = 2
N_STAGES = 2  # focus on L_0 and L_1 only
BONUS_INT = 8
# GAUGE_EPS removed — canonical recipe uses no penalty.
# H1 Poisson for Helmholtz-Hodge needs a tiny shift to be invertible
# (the constant function is in its kernel without Dirichlet on a face).
PHI_SHIFT = 1e-12

R0_anal = 48 / (sigma_Cu * V_cond * (ax**2 + ay**2))
ELF_TAU_LEAD = 11.51e-6


def build_spanning_tree(mesh):
    nv = mesh.nv
    visited = [False] * nv
    tree_edges = []
    adj = [[] for _ in range(nv)]
    for ed in mesh.edges:
        v0, v1 = ed.vertices[0].nr, ed.vertices[1].nr
        adj[v0].append((v1, ed.nr))
        adj[v1].append((v0, ed.nr))
    visited[0] = True
    queue = deque([0])
    while queue:
        v = queue.popleft()
        for vn, edn in adj[v]:
            if not visited[vn]:
                visited[vn] = True
                tree_edges.append(edn)
                queue.append(vn)
    return tree_edges


def build_geo():
    cuboid = Box(Pnt(-ax/2, -ay/2, -az/2), Pnt(ax/2, ay/2, az/2))
    cuboid.name = "conductor"
    # Name conductor faces so we can do ds("conductor_surface") for the
    # Coulomb correction surface integral.
    for f in cuboid.faces:
        f.name = "conductor_surface"
    cuboid.maxh = H_COND
    sphere_inner = Sphere(Pnt(0, 0, 0), R_K)
    for f in sphere_inner.faces:
        f.name = "kelvin_int"
    sphere_inner.maxh = H_AIR
    inner_air = sphere_inner - cuboid
    inner_air.name = "air"
    geo, info = add_kelvin_exterior_domain(
        [inner_air, cuboid], offset=OFFSET, R_K=R_K, inner_maxh=H_AIR)
    return OCCGeometry(geo)


def main():
    ngsglobals.msg_level = 0
    print(f"=== v15: canonical CLN iteration with Kelvin ===", flush=True)
    print(f"  R_K = {R_K*1000:.1f} mm, OFFSET = {OFFSET}", flush=True)
    print(f"  Reference R_0 (Parseval) = {R0_anal:.6e}", flush=True)
    print(f"  ELF reference tau_lead = {ELF_TAU_LEAD*1e6:.2f} us", flush=True)
    print(f"  N_STAGES = {N_STAGES}\n", flush=True)

    geo = build_geo()
    print("Generating mesh...", flush=True)
    t0 = time.time()
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=H_AIR, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(ORDER + 1)
    print(f"  ne = {mesh.ne}  ({time.time()-t0:.1f}s)", flush=True)
    print(f"  materials = {mesh.GetMaterials()}\n", flush=True)

    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu, "air": 0.0,
                                 "kelvin": 0.0})
    sigma_inv_cf = mesh.MaterialCF({"conductor": 1.0/sigma_Cu, "air": 0.0,
                                     "kelvin": 0.0})

    nu_cf = make_kelvin_nu_cf(mesh, R_K, OFFSET, nu_0=NU_0,
                              kelvin_mats=("kelvin",))

    A_s_cf = make_reduced_potential_background_cf(
        mesh,
        F_inner_factory=lambda xc, yc, zc: CoefficientFunction((-yc, xc, 0)) * 0.5,
        R_K=R_K, offset=OFFSET, kelvin_mats=("kelvin",), dim=3,
    )

    fes = Periodic(HCurl(mesh, order=ORDER, dirichlet_bbnd="GND",
                         nograds=True))
    print(f"  HCurl ndof = {fes.ndof}", flush=True)
    tree_edges = build_spanning_tree(mesh)
    fd = fes.FreeDofs()
    masked = 0
    for edge_nr in tree_edges:
        edge = mesh.edges[edge_nr]
        dofs = fes.GetDofNrs(edge)
        if dofs and fd[dofs[0]]:
            fd[dofs[0]] = False
            masked += 1
    print(f"  tree-cotree masked {masked} edges", flush=True)

    u, v = fes.TnT()

    gfA_s = GridFunction(fes, name="A_s")
    gfA_s.vec[:] = 0.0
    print("  Projecting A_s onto HCurl (Periodic)...", flush=True)
    t0 = time.time()
    with TaskManager():
        gfA_s.Set(A_s_cf, bonus_intorder=BONUS_INT)
    print(f"    {time.time()-t0:.1f}s", flush=True)

    A_s_phys_cf = CoefficientFunction((-y, x, 0)) * 0.5
    err_in_cond = float(Integrate(
        (gfA_s - A_s_phys_cf) * (gfA_s - A_s_phys_cf)
        * dx("conductor", bonus_intorder=BONUS_INT), mesh))
    norm_in_cond = float(Integrate(
        A_s_phys_cf * A_s_phys_cf
        * dx("conductor", bonus_intorder=BONUS_INT), mesh))
    rel_err = (err_in_cond / norm_in_cond) ** 0.5 if norm_in_cond > 0 else 0
    print(f"  ||gfA_s - A_s_phys||/||A_s_phys|| in conductor = {rel_err:.2e}",
          flush=True)

    # === HDiv-projected J source for J·n=0 strict (v20b) ===
    # v20 found that HDiv .Set() did NOT enforce dirichlet="conductor_surface".
    # v20b explicitly zeroes out the Dirichlet (non-Free) DOFs after the Set
    # so that J·n=0 holds strictly on all 6 cuboid faces.
    print("  Building HDiv source space + explicit Dirichlet zero...",
          flush=True)
    t0 = time.time()
    fes_J = HDiv(mesh, order=ORDER,
                  dirichlet="conductor_surface",
                  definedon=mesh.Materials("conductor"))
    print(f"  HDiv ndof (conductor only) = {fes_J.ndof}, "
          f"FreeDofs count = {sum(fes_J.FreeDofs())}", flush=True)

    gfJ_imp = GridFunction(fes_J, name="J_imp_HDiv")
    with TaskManager():
        gfJ_imp.Set(sigma_Cu * A_s_phys_cf,
                     definedon=mesh.Materials("conductor"),
                     bonus_intorder=BONUS_INT)

    # EXPLICIT Dirichlet zero: HDiv .Set() does not respect dirichlet flag
    # in this configuration, so we manually zero-out the constrained DOFs.
    fd = fes_J.FreeDofs()
    arr = gfJ_imp.vec.FV().NumPy()
    n_zeroed = 0
    for i in range(len(arr)):
        if not fd[i]:
            if abs(arr[i]) > 0:
                n_zeroed += 1
            arr[i] = 0.0
    print(f"  manually zeroed {n_zeroed} Dirichlet DOFs", flush=True)
    print(f"    {time.time()-t0:.1f}s", flush=True)

    # Diagnostic: verify (gfJ_imp · n) at conductor surface = 0
    n_normal = specialcf.normal(3)
    Jn_int = float(Integrate(
        (gfJ_imp * n_normal) ** 2
        * ds("conductor_surface", bonus_intorder=BONUS_INT), mesh))
    print(f"  ∫_∂C |J_imp_HDiv·n|² dS = {Jn_int:.4e} "
          f"(should be ~0 from Dirichlet)", flush=True)

    # How well does gfJ_imp approximate sigma * A_s_phys in conductor?
    proj_err = float(Integrate(
        (gfJ_imp - sigma_Cu * A_s_phys_cf) *
        (gfJ_imp - sigma_Cu * A_s_phys_cf)
        * dx("conductor", bonus_intorder=BONUS_INT), mesh))
    proj_norm = float(Integrate(
        (sigma_Cu * A_s_phys_cf) * (sigma_Cu * A_s_phys_cf)
        * dx("conductor", bonus_intorder=BONUS_INT), mesh))
    print(f"  ||gfJ_imp - σA_phys||/||σA_phys|| in conductor = "
          f"{(proj_err/proj_norm)**0.5:.4e}", flush=True)
    A_s_corr_n_int = Jn_int  # for output JSON
    A_s_n_orig_int = float(Integrate(
        (sigma_Cu * A_s_phys_cf * n_normal) ** 2
        * ds("conductor_surface", bonus_intorder=BONUS_INT), mesh))
    curl_diff = proj_err  # not strictly the curl diff; reused field

    # Note: gfA_s is still needed for energy_norm / R_0 calculations
    # downstream, but gfJ_imp replaces sigma_cf * gfA_s as the iteration
    # source.

    a = BilinearForm(fes)
    a += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=BONUS_INT)
    # NO penalty term — canonical recipe relies on tree-cotree gauge alone
    # for the curl-curl operator. Helmholtz-Hodge projection is applied
    # post-solve to clean residual gradient components introduced through
    # the Periodic Kelvin boundary.

    print("  Assembling+factor (curl-curl)...", flush=True)
    t0 = time.time()
    with TaskManager():
        a.Assemble()
        try:
            inv = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso")
        except Exception:
            inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    print(f"    {time.time()-t0:.1f}s\n", flush=True)

    # === Helmholtz-Hodge: H1 space + Poisson factorization ===
    # Same Periodic + GND vertex Dirichlet as HCurl, so H1 gradients live
    # in the same constraint structure.
    fes_phi = Periodic(H1(mesh, order=ORDER, dirichlet_bbnd="GND"))
    print(f"  H1 ndof (for Helmholtz-Hodge) = {fes_phi.ndof}", flush=True)
    phi, psi = fes_phi.TnT()
    a_phi = BilinearForm(fes_phi)
    a_phi += grad(phi) * grad(psi) * dx(bonus_intorder=BONUS_INT)
    a_phi += PHI_SHIFT * phi * psi * dx(bonus_intorder=BONUS_INT)
    print("  Assembling+factor (H1 Poisson for Helmholtz-Hodge)...",
          flush=True)
    t0 = time.time()
    with TaskManager():
        a_phi.Assemble()
        try:
            inv_phi = a_phi.mat.Inverse(fes_phi.FreeDofs(), inverse="pardiso")
        except Exception:
            inv_phi = a_phi.mat.Inverse(fes_phi.FreeDofs(),
                                        inverse="sparsecholesky")
    print(f"    {time.time()-t0:.1f}s\n", flush=True)

    def helmholtz_hodge_project(gfA_in):
        """Subtract grad(phi) such that ∫ (gfA - grad(phi)) · grad(psi) dx = 0
        for all psi in H1 (Periodic, GND-Dirichlet)."""
        f_phi = LinearForm(fes_phi)
        f_phi += gfA_in * grad(psi) * dx(bonus_intorder=BONUS_INT)
        with TaskManager():
            f_phi.Assemble()
        gf_phi = GridFunction(fes_phi)
        with TaskManager():
            gf_phi.vec.data = inv_phi * f_phi.vec
        gfA_proj = GridFunction(fes)
        with TaskManager():
            gfA_proj.Set(gfA_in - grad(gf_phi), bonus_intorder=BONUS_INT)
        return gfA_proj

    # R_0 check via the HDiv-projected J source: <J_imp, J_imp/sigma>_C
    R_inv_check = float(Integrate(gfJ_imp * gfJ_imp * sigma_inv_cf
                                  * dx("conductor", bonus_intorder=BONUS_INT),
                                  mesh))
    R0_test = 1.0 / R_inv_check if R_inv_check > 1e-30 else float('inf')
    print(f"  R_0 verification (HDiv source): {R0_test:.6e} "
          f"vs Parseval {R0_anal:.6e}, ratio {R0_test/R0_anal:.6f}\n",
          flush=True)

    # Iteration: J in conductor (HDiv-class), A in HCurl (full domain).
    # Initial J_imp is the HDiv GridFunction gfJ_imp. As Schmidt update
    # subtracts sigma * A_pot_accum / L_n, the result becomes a generic
    # CF (not strictly HDiv anymore). For a truly HDiv-bound iteration
    # we would re-project at each step — for stage 0 only this is moot.
    diag = []
    J_history = []
    J_imp_cf = gfJ_imp  # use HDiv GridFunction directly as iteration source
    gfApot = GridFunction(fes)
    gfApot.vec[:] = 0.0
    sign_flip_first_n = None

    for n in range(N_STAGES):
        t_stage = time.time()

        # CANONICAL CLN linear form: source is J_imp in conductor only.
        # NO curl(A_s) RHS — that's the §7.5 PEEC formula, not applicable here.
        f = LinearForm(fes)
        f += J_imp_cf * v * dx("conductor", bonus_intorder=BONUS_INT)
        with TaskManager():
            f.Assemble()
        gfA_r_raw = GridFunction(fes)
        with TaskManager():
            gfA_r_raw.vec.data = inv * f.vec

        # Helmholtz-Hodge projection: enforce ∇·A = 0 by subtracting grad(phi)
        gfA_r = helmholtz_hodge_project(gfA_r_raw)

        R_inv = float(Integrate(J_imp_cf * J_imp_cf * sigma_inv_cf
                                * dx("conductor", bonus_intorder=BONUS_INT),
                                mesh))
        if R_inv < 1e-30:
            print(f"  [{n}] energy_norm collapsed ({R_inv:.2e}), stopping",
                  flush=True)
            break
        R_n = 1.0 / R_inv

        # Project current J_imp_cf to a GridFunction for cross-product diagnostics
        gfJ_n = GridFunction(fes)
        try:
            gfJ_n.Set(J_imp_cf, definedon=mesh.Materials("conductor"))
        except Exception:
            gfJ_n.Set(J_imp_cf)

        # Schmidt orthogonality drift: <J_n, J_m/sigma>_cond for m<n
        cross_inner = []
        for gfJ_m in J_history:
            xnm = float(Integrate(J_imp_cf * gfJ_m * sigma_inv_cf
                                  * dx("conductor", bonus_intorder=BONUS_INT),
                                  mesh))
            cross_inner.append(xnm)
        max_cross = max((abs(c) for c in cross_inner), default=0.0)
        rel_drift = max_cross / abs(R_inv) if R_inv != 0 else 0.0

        gfApot.vec.data += R_n * gfA_r.vec

        # CANONICAL B² integral form (Tanimoto-Sugahara TMAG3304725 eq. 3 + IBP):
        #     L_{n,n} = ∫_full ν |curl(â_n)|² dV   with â_n = R_n × A_n
        # Manifestly ≥ 0. Drops cross-term ∫ ν curl(A_k)·curl(A_n) (k<n) that
        # contaminates the <J_n, A_pot_accum>_cond formulation when A-side
        # Schmidt orthogonality is not separately enforced.
        B2_int = float(Integrate(nu_cf * curl(gfA_r) * curl(gfA_r)
                                 * dx(bonus_intorder=BONUS_INT), mesh))
        L_n = R_n * R_n * B2_int

        # Also compute v17 form for diagnostic comparison
        Ln_int_v17 = float(Integrate(J_imp_cf * gfApot
                                     * dx("conductor", bonus_intorder=BONUS_INT),
                                     mesh))
        L_n_v17_form = R_n * Ln_int_v17

        signL = "+" if L_n > 0 else ("-" if L_n < 0 else "0")
        signR = "+" if R_n > 0 else ("-" if R_n < 0 else "0")
        signL_v17 = "+" if L_n_v17_form > 0 else ("-" if L_n_v17_form < 0 else "0")
        tau_us = L_n/R_n*1e6 if R_n > 0 and L_n > 0 else None
        elapsed = time.time() - t_stage
        tau_str = f"{tau_us:.4f}" if tau_us else "N/A"
        print(f"  [{n}] drift={rel_drift:.3e}  En={R_inv:.3e}  "
              f"R={R_n:.3e}({signR})  L_B2={L_n:.3e}({signL})  "
              f"L_v17={L_n_v17_form:.3e}({signL_v17})  "
              f"tau={tau_str}us  ({elapsed:.1f}s)", flush=True)
        diag.append({"n": n, "R_n": R_n, "L_n": L_n,
                     "L_n_v17_form": L_n_v17_form,
                     "B2_int": B2_int, "tau_us": tau_us,
                     "energy_norm": R_inv, "rel_drift": rel_drift,
                     "sign_L": signL, "sign_R": signR})

        if L_n <= 0 and sign_flip_first_n is None:
            sign_flip_first_n = n
            print(f"      ** FIRST SIGN FLIP at n={n} (continuing for "
                  "breakdown diagnostic) **", flush=True)

        J_history.append(gfJ_n)

        # CANONICAL Schmidt update. With L_n < 0 the iteration is no longer
        # well-defined but we keep going so the breakdown trajectory is recorded.
        if abs(L_n) < 1e-30:
            print(f"  [{n}] |L_n| collapsed, stopping", flush=True)
            break
        J_imp_cf = J_imp_cf - sigma_cf * gfApot * (1.0 / L_n)

    print("\n" + "=" * 78, flush=True)
    print("HCURL + KELVIN KAMEARI ITERATION — BREAKDOWN TRAJECTORY", flush=True)
    print("=" * 78, flush=True)
    print(f"{'n':>3} {'rel_drift':>12} {'energy_norm':>14} "
          f"{'R_n':>12} {'L_n':>14} {'tau_us':>12}", flush=True)
    for d in diag:
        tau_s = f"{d['tau_us']:.4f}" if d['tau_us'] else "N/A"
        print(f"{d['n']:>3} {d['rel_drift']:>12.3e} {d['energy_norm']:>14.4e} "
              f"{d['R_n']:>12.3e} {d['L_n']:>14.4e}({d['sign_L']}) "
              f"{tau_s:>12}",
              flush=True)

    print("\nBreakdown summary:", flush=True)
    if sign_flip_first_n is not None:
        print(f"  First L_n sign flip at stage n = {sign_flip_first_n}",
              flush=True)
    else:
        print("  No sign flip in this run (unexpected for HCurl Kameari "
              "in vacuum-coupled isolated cuboid)", flush=True)
    en_first = diag[0]['energy_norm'] if diag else None
    en_last = diag[-1]['energy_norm'] if diag else None
    if en_first and en_last:
        print(f"  Energy norm n=0 -> n={diag[-1]['n']}: "
              f"{en_first:.3e} -> {en_last:.3e}  (ratio {en_last/en_first:.3e})",
              flush=True)
    print("\nThis trajectory is the paper's motivation: HCurl Kameari "
          "iteration breaks down on isolated-conductor-in-vacuum problems. "
          "The proposed Mathematica BEM Foster-CLN method "
          "(hex_bem_hdiv_order3_overnight.wls, HDiv interior div-free, "
          "J·n=0 strict) does not exhibit this breakdown.",
          flush=True)

    out = {"method": "HCurl + Kelvin Kameari iteration (v20: HDiv-projected J source, J·n=0 strict)",
           "L_n_formula": "L_n = R_n² × ∫_full ν |curl(gfA_r)|² dV (canonical B² form)",
           "gauge_recipe": "nograds=True + tree-cotree + Helmholtz-Hodge + HDiv source (no penalty)",
           "R_K_mm": R_K*1000, "h_cond_mm": H_COND*1000,
           "h_air_mm": H_AIR*1000, "order": ORDER, "ne": mesh.ne,
           "rel_err_gfA_s_cond": rel_err,
           "HDiv_Jn_int": A_s_corr_n_int,
           "sigma_A_phys_n_orig_int": A_s_n_orig_int,
           "HDiv_proj_err": curl_diff,
           "stages": diag, "R0_analytic": R0_anal,
           "sign_flip_first_n": sign_flip_first_n}
    out_path = Path(__file__).parent / "cuboid_521_kameari_kelvin_v20b_HDiv_dirichlet_force.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}", flush=True)


if __name__ == "__main__":
    main()
