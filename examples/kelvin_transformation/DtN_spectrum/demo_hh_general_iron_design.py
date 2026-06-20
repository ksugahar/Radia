# -*- coding: utf-8 -*-
# DEMO (hh): GENERAL coil inverse design with ARBITRARY (non-concentric) iron, via the
# material-aware Kelvin-FEM transfer matrix M -- the Track-B next step after demo_ff.
#
# demo_ff was a CONCENTRIC, MODAL toy: M[k,n] = R_n(r_t) P_n(cos th_k) with R_n the analytic
# layered-sphere radial transfer (3 spherical-harmonic mode amplitudes).  This promotes it to:
#   * a REAL winding-surface stream function psi -- the piecewise-(order p) H1 trace on the coil
#     surface (hundreds of nodal DoFs), NOT 3 modal amplitudes;
#   * ARBITRARY iron -- a non-concentric blob where NO closed-form layered/Sommerfeld Green
#     function exists, so the material-aware kernel can ONLY come from the Kelvin-FEM;
#   * the transfer matrix M[target_k, psi_dof_j] built DIRECTLY from ONE sparse factorisation of
#     the Kelvin-FEM (one back-substitution per coil DoF -- the demo_v Schur idea, specialised to
#     the Dirichlet(coil)->field(target) map the design inverts);
#   * an INVERSE design psi = M^+ B_target (folded TSVD, the radia.streamfunction regularisation);
#   * a FRESH full Kelvin-FEM forward solve of the designed psi (gf.Set(psi)+solve -- it does NOT
#     touch M) confirming the iron-aware design HITS the target while the free-space-designed coil
#     MISSES it in the real iron system.
#
# Coil = the inner Dirichlet boundary (r=a); the design variable psi = the magnetic-scalar-potential
# trace there.  Targets sit in the PHYSICAL vacuum region (a<r<R_out, outside the iron) and are read
# DIRECTLY as gf(target) -- so the verification does NOT depend on the inverse-Kelvin map convention
# (only the open BC at infinity uses the Kelvin ball).  Physics is anchored on the CONCENTRIC sub-case
# against the analytic layered transfer (Scenario A); the non-concentric design is the new result.
#
# Kelvin-FEM gotchas honoured (each is a ~1e7x blow-up if missed): single-face periodic glue (one
# Identify); mu carried as an FE coefficient incl. the kelvin-ball (R/rho)^2 weight; point ground +
# coil-Dirichlet pin the constant; n>=1 only (no spurious monopole -- magnetostatics has none anyway).
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import time
import numpy as np
import ngsolve as ng
import netgen.occ as occ
from netgen.occ import Sphere, Pnt, Vec, IdentificationType, OCCGeometry
from radia_mcp.radia_ngsolve.fem_bem_coupling import _solid_harmonic

a, R_out, offset = 0.5, 1.0, 3.0                 # coil r=a, Kelvin truncation/ball R_out, ball offset
r_t = 0.95                                       # target radius in the physical vacuum region (a<r_t<R_out)
INTORDER = 8


# ----------------------------------------------------------------------------------------------------
# analytic layered-sphere transfer (CONCENTRIC anchor only): external value R_n(r) for a unit Y_n
# coefficient on r=a, magnetic shell mu_s in [b,c], vacuum elsewhere, decay at infinity.
# ----------------------------------------------------------------------------------------------------
def analytic_R(n, mu_s, r, b, c):
    M = np.zeros((5, 5)); rhs = np.zeros(5)
    M[0, 0] = a**n; M[0, 1] = a**-(n + 1); rhs[0] = 1.0
    M[1, 0] = b**n; M[1, 1] = b**-(n + 1); M[1, 2] = -b**n; M[1, 3] = -b**-(n + 1)
    M[2, 0] = n * b**(n - 1); M[2, 1] = -(n + 1) * b**-(n + 2)
    M[2, 2] = -mu_s * n * b**(n - 1); M[2, 3] = -mu_s * (-(n + 1)) * b**-(n + 2)
    M[3, 2] = c**n; M[3, 3] = c**-(n + 1); M[3, 4] = -c**-(n + 1)
    M[4, 2] = mu_s * n * c**(n - 1); M[4, 3] = mu_s * (-(n + 1)) * c**-(n + 2)
    M[4, 4] = -(-(n + 1)) * c**-(n + 2)
    return np.linalg.solve(M, rhs)[4] * r**-(n + 1)


# c_n: the probe trace _solid_harmonic(n)/a^n equals c_n * P_n on r=a.  The current zonal solid
# harmonic is _solid_harmonic(n) = r^n P_n(cos th) (recurrence build, fem_bem_coupling), so /a^n
# traces to P_n EXACTLY => c_n = 1 for all n.  (Was {1,2,2} for an earlier _solid_harmonic
# normalisation; that became stale when the recurrence build landed -- the energy-quotient results
# are scale-invariant and unchanged, but this ABSOLUTE transfer comparison is not, so c_n is now
# MEASURED at runtime in scenario_A and asserted against this dict to fail loud on any future drift.)
_CN = {1: 1.0, 2: 1.0, 3: 1.0}


# ----------------------------------------------------------------------------------------------------
# build the Kelvin-FEM (coil = inner Dirichlet boundary; iron = "shell" material; open via Kelvin ball)
#   iron_spec: None (no iron) | ("concentric", b, c) | ("blob", (cx,cy,cz), r_i)
# ----------------------------------------------------------------------------------------------------
def build_fem(iron_spec, order, maxh):
    """Build the Kelvin-FEM mesh + space ONCE; return (mesh, fes, make_A) where make_A(mu_shell)
    assembles A for a given shell permeability.  The SAME mesh/DoFs serve both the material-aware
    operator (mu_shell=mu_r) and the free-space operator (mu_shell=1): the iron region is always
    MESHED (a distinct 'shell' material), so the two operators share the coil discretisation and a
    free-space-designed psi is directly realisable in the iron forward solve (no cross-mesh map)."""
    inner = Sphere(Pnt(0, 0, 0), a); outer = Sphere(Pnt(0, 0, 0), R_out)
    for f in inner.faces: f.name = "inner"
    for f in outer.faces: f.name = "kelvin_int"
    if iron_spec is None or iron_spec[0] == "none":
        vac = (outer - inner); vac.mat("vac"); solids = [vac]
    elif iron_spec[0] == "concentric":
        _, b, c = iron_spec
        s_b = Sphere(Pnt(0, 0, 0), b); s_c = Sphere(Pnt(0, 0, 0), c)
        sh1 = (s_b - inner); sh1.mat("vac"); shm = (s_c - s_b); shm.mat("shell")
        sh3 = (outer - s_c); sh3.mat("vac"); solids = [sh1, shm, sh3]
    elif iron_spec[0] == "blob":
        _, ctr, r_i = iron_spec
        iron = Sphere(Pnt(*ctr), r_i); iron.mat("shell")
        vac = (outer - inner) - iron; vac.mat("vac"); solids = [vac, iron]
    else:
        raise ValueError("bad iron_spec %r" % (iron_spec,))
    kball = Sphere(Pnt(offset, 0, 0), R_out)
    for f in kball.faces: f.name = "kelvin_ext"
    kball.mat("kelvin"); gnd = occ.Vertex(Pnt(offset, 0, 0)); gnd.name = "GND"
    fi = [f for s in solids for f in s.faces if f.name == "kelvin_int"][0]
    fe = [f for f in kball.faces if f.name == "kelvin_ext"][0]
    fi.Identify(fe, "kelvin", IdentificationType.PERIODIC, occ.gp_Trsf.Translation(Vec(offset, 0, 0)))
    mesh = ng.Mesh(OCCGeometry(occ.Glue(solids + [kball, gnd])).GenerateMesh(maxh=maxh)).Curve(min(order + 1, 4))
    x, y, z = ng.x, ng.y, ng.z; rp2 = (x - offset)**2 + y * y + z * z + 1e-20
    has_shell = "shell" in mesh.GetMaterials()
    fes = ng.Periodic(ng.H1(mesh, order=order, dirichlet="inner|GND"))
    u, v = fes.TnT()

    def make_A(mu_shell):
        matvals = {"vac": 1.0, "kelvin": R_out**2 / rp2}
        if has_shell:
            matvals["shell"] = mu_shell
        mu = mesh.MaterialCF(matvals, default=1.0)
        A = ng.BilinearForm(mu * ng.grad(u) * ng.grad(v) * ng.dx(bonus_intorder=INTORDER)); A.Assemble()
        return A
    return mesh, fes, make_A


def coil_dofs(mesh, fes):
    g = fes.GetDofs(mesh.Boundaries("inner"))
    return [i for i in range(fes.ndof) if g[i]]


def targets_on_sphere(mesh, thetas, phis):
    pts, dirs = [], []
    for th in thetas:
        for ph in phis:
            d = (np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th))
            dirs.append(d); pts.append(mesh(r_t * d[0], r_t * d[1], r_t * d[2]))
    return pts, np.array(dirs)


def solve_trace(fes, A, Ainv, gf):
    """In-place: complete the Dirichlet lift (gf already carries the coil trace) by the interior solve."""
    rr = gf.vec.CreateVector(); rr.data = -(A.mat * gf.vec); gf.vec.data += Ainv * rr


def build_transfer(mesh, fes, A, Ainv, cdofs, targets):
    """M[k,j] = scalar potential at target k from a unit Dirichlet bump on coil DoF j.
    ONE factorisation (Ainv) reused; one back-substitution + target eval per coil DoF."""
    M = np.zeros((len(targets), len(cdofs)))
    for col, idof in enumerate(cdofs):
        g = ng.GridFunction(fes); g.vec[:] = 0.0; g.vec[idof] = 1.0
        solve_trace(fes, A, Ainv, g)
        M[:, col] = [g(t) for t in targets]
    return M


def tsvd_pinv(M, rcond=1e-10):
    U, s, Vt = np.linalg.svd(M, full_matrices=False)
    keep = s > rcond * s[0]
    return (Vt[keep].T * (1.0 / s[keep])) @ U[:, keep].T


def fresh_forward(mesh, fes, A, Ainv, cdofs, psi, targets):
    """A FRESH full Kelvin-FEM solve: set the designed nodal psi on the coil, solve, read targets.
    Does NOT use M -- so it independently checks the M-assembly + the design."""
    gf = ng.GridFunction(fes); gf.vec[:] = 0.0
    for j, idof in enumerate(cdofs):
        gf.vec[idof] = psi[j]
    solve_trace(fes, A, Ainv, gf)
    return np.array([gf(t) for t in targets])


# ====================================================================================================
# SCENARIO A -- CONCENTRIC anchor: the FEM transfer matrix M matches the analytic layered transfer.
# ====================================================================================================
def scenario_A(order=2, maxh=0.25):
    print("=" * 96)
    print("SCENARIO A -- CONCENTRIC anchor: FEM transfer M vs analytic layered-sphere transfer")
    print("=" * 96)
    b, c, mu_r = 0.7, 0.9, 50.0
    mesh, fes, make_A = build_fem(("concentric", b, c), order, maxh)
    A = make_A(mu_r)
    Ainv = A.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    cdofs = coil_dofs(mesh, fes)
    thetas = np.deg2rad([10, 35, 60, 90, 120, 150]); phis = np.deg2rad([0.0])
    targets, dirs = targets_on_sphere(mesh, thetas, phis)
    t0 = time.time(); M = build_transfer(mesh, fes, A, Ainv, cdofs, targets); tM = time.time() - t0
    print("  iron shell mu_r=%.0f in [%.2f,%.2f]; coil DoFs=%d; %d targets at r_t=%.2f; M built %.2fs"
          % (mu_r, b, c, len(cdofs), len(targets), r_t, tM))
    from scipy.special import eval_legendre
    cth = np.cos(thetas)
    print("\n   n   FEM R_n(r_t)   analytic*c_n   rel.err   (c_n)   M@psi==fresh-solve")
    okA = True
    for n in (1, 2, 3):
        gf = ng.GridFunction(fes)
        gf.Set(_solid_harmonic(n) / a**n, ng.BND, definedon=mesh.Boundaries("inner"))
        psi_n = np.array([gf.vec[i] for i in cdofs])
        v_M = M @ psi_n                                            # matrix prediction
        v_fresh = fresh_forward(mesh, fes, A, Ainv, cdofs, psi_n, targets)   # fresh solve (no M)
        book = np.linalg.norm(v_M - v_fresh) / (np.linalg.norm(v_fresh) + 1e-300)
        P = eval_legendre(n, cth)
        fem_R = float(np.dot(v_fresh, P) / np.dot(P, P))           # best-fit R_n from the angular pattern
        # measure the probe trace's P_n coefficient on r=a (robust to the _solid_harmonic normalisation;
        # asserts against _CN so a future normalisation change fails LOUD instead of silently rescaling).
        c_meas = float(np.dot([gf(mesh(a * np.sin(t), 0.0, a * np.cos(t))) for t in thetas], P) / np.dot(P, P))
        assert abs(c_meas - _CN[n]) < 0.05, (
            "probe normalisation drifted: measured c_%d=%.3f != _CN=%.1f -- update _CN to match "
            "_solid_harmonic (fem_bem_coupling)" % (n, c_meas, _CN[n]))
        an_R = analytic_R(n, mu_r, r_t, b, c) * c_meas
        rel = abs(fem_R - an_R) / abs(an_R)
        okA = okA and (rel < 3e-2) and (book < 1e-9)
        print("   %d   %12.6e   %12.6e   %.2e   (%.2f)   %.1e"
              % (n, fem_R, an_R, rel, c_meas, book))
    print("\n  => the directly-assembled material-aware transfer M reproduces the analytic layered")
    print("     transfer (physics) AND M@psi == a fresh FEM solve (assembly).  [%s]"
          % ("PASS" if okA else "CHECK"))
    return okA


# ====================================================================================================
# SCENARIO B -- the new result: GENERAL design with NON-CONCENTRIC iron; design HITS, free-space MISSES.
# ====================================================================================================
def _design_once(blob, mu_r, order, maxh):
    """One full design + fresh-forward at a given resolution.  Returns (err_iron, err_free, n_coil)."""
    # ONE mesh/space; the iron-aware (mu_shell=mu_r) and free-space (mu_shell=1) operators share the
    # coil discretisation, so a free-space-designed psi transplants directly into the iron forward solve.
    mesh, fes, make_A = build_fem(blob, order, maxh)
    cdofs = coil_dofs(mesh, fes)
    A_iron = make_A(mu_r); Ai = A_iron.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    A_free = make_A(1.0);  Af = A_free.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    thetas = np.deg2rad(np.linspace(12, 168, 8)); phis = np.deg2rad([0, 90, 180, 270])
    targets, _ = targets_on_sphere(mesh, thetas, phis)
    M_iron = build_transfer(mesh, fes, A_iron, Ai, cdofs, targets)
    M_free = build_transfer(mesh, fes, A_free, Af, cdofs, targets)
    # PRODUCIBLE target: a tilted dipole trace on the coil (engages the off-axis iron), realised in
    # the IRON system -> target potential samples the design must reproduce.
    gfw = ng.GridFunction(fes)
    gfw.Set((ng.z + 0.6 * ng.x) / a, ng.BND, definedon=mesh.Boundaries("inner"))
    target = M_iron @ np.array([gfw.vec[i] for i in cdofs])
    psi_iron = tsvd_pinv(M_iron) @ target
    psi_free = tsvd_pinv(M_free) @ target
    # FRESH full Kelvin-FEM forward solves in the REAL IRON system (do NOT use M) ----------------
    ach_iron = fresh_forward(mesh, fes, A_iron, Ai, cdofs, psi_iron, targets)
    ach_free = fresh_forward(mesh, fes, A_iron, Ai, cdofs, psi_free, targets)
    nt = np.linalg.norm(target)
    return (np.linalg.norm(ach_iron - target) / nt, np.linalg.norm(ach_free - target) / nt,
            len(cdofs), len(targets))


def scenario_B(order=2):
    print("\n" + "=" * 96)
    print("SCENARIO B -- NON-CONCENTRIC iron: real-surface psi design HITS; free-space design MISSES")
    print("=" * 96)
    blob = ("blob", (0.0, 0.0, 0.68), 0.16); mu_r = 50.0
    print("  iron blob r=%.2f @ %s (off-centre on +z, between coil r=%.1f and R_out=%.1f); mu_r=%.0f"
          % (blob[2], blob[1], a, R_out, mu_r))
    t0 = time.time()
    err_iron, err_free, n_coil, n_tgt = _design_once(blob, mu_r, order, 0.22)
    print("  coil DoFs=%d; %d targets at r_t=%.2f; design built in %.2fs" % (n_coil, n_tgt, r_t, time.time() - t0))
    print("\n  DESIGN realised in a FRESH full Kelvin-FEM solve with the iron (relative target error):")
    print("    design WITH iron  (invert material-aware M): %.2e   <- HITS (design inverts the exact M)"
          % err_iron)
    print("    design IGNORING iron (invert free-space M) : %.2e   <- MISSES by ~%.0f%% (iron not in kernel)"
          % (err_free, 100 * err_free))

    # mesh-robustness: the free-space MISS is a PHYSICAL shield effect, not a coarse-mesh artifact ---
    print("\n  mesh-refinement (is the free-space miss physical, not a discretisation artifact?):")
    print("    maxh   coil DoFs   err_iron (hit)   err_free (miss)")
    free_misses = []
    for mh in (0.30, 0.24, 0.19):
        ei, ef, nc, _ = _design_once(blob, mu_r, order, mh)
        free_misses.append(ef)
        print("    %.2f      %4d        %.2e        %.3f" % (mh, nc, ei, ef))
    # the robust statement is "the free-space design ALWAYS misses by a lot" (the exact % is a
    # physical quantity that shifts with the discretised producible target -- here 31-39%).
    okB = (err_iron < 1e-2) and (min(free_misses) > 0.25)
    print("    -> free-space miss is consistently large (%.0f-%.0f%%): a real shield effect, not a"
          % (100 * min(free_misses), 100 * max(free_misses)))
    print("       coarse-mesh artifact (the iron-aware design stays exact to machine precision)")
    print("\n  => with NON-CONCENTRIC iron (no closed-form Green function) the Kelvin-FEM transfer M is")
    print("     the correct design kernel: the iron-aware design reproduces the target in an independent")
    print("     forward solve; the free-space-Biot-Savart design is off by ~%.0f%%.   [%s]"
          % (100 * err_free, "PASS" if okB else "CHECK"))
    return okB, dict(err_iron=float(err_iron), err_free=float(err_free),
                     n_coil=n_coil, n_target=n_tgt,
                     free_misses=[float(x) for x in free_misses])


if __name__ == "__main__":
    okA = scenario_A()
    okB, statB = scenario_B()
    print("\n" + "=" * 96)
    print("RESULT: Scenario A (concentric anchor) %s ; Scenario B (non-concentric design) %s"
          % ("PASS" if okA else "CHECK", "PASS" if okB else "CHECK"))
    print("  Track B next step done: real-surface stream-function coil design with arbitrary iron via the")
    print("  material-aware Kelvin-FEM transfer matrix M -- design inverts M, free-space kernel misses by")
    print("  %.0f%%, iron-aware design hits in an independent forward solve." % (100 * statB["err_free"]))
