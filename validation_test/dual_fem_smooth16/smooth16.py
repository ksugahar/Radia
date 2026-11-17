"""Smooth periodic A/Omega dual FEM rate study (manufactured truth, no boundary).

Omega: periodic H1(p) scalar potential, solving  int mu grad(phi).grad(v) = int S.grad(v).
A:     2D periodic H1(p) stream function psi with B = rot psi;
       3D periodic HCurl(p, nograds=True) with B = curl A, an H1(1) Lagrange multiplier
       for the Whitney gradients, and three scalar multipliers for the harmonic constants.
       It solves  int nu (B - S).curl(W) = 0.
Both formulations give fields in P_{p-1}.  The truth derivatives come from
NGSolve CoefficientFunction.Diff (symbolic forward differentiation).  NGSolve
owns the spaces, periodicity, quadrature and assembly.
"""

from math import log, sqrt
from time import perf_counter

import numpy as np


def _grad(ng, f, dim):
    return ng.CF(tuple(f.Diff(v) for v in (ng.x, ng.y, ng.z)[:dim]))


def truth(ng, dim, fields):
    mu = fields["mu"]
    phi = fields["phi"]
    grad_phi = _grad(ng, phi, dim)
    if dim == 2:
        psi = fields["psi"]
        b = ng.CF((psi.Diff(ng.y), -psi.Diff(ng.x)))
    else:
        a = fields["A"]
        b = ng.CF((a[2].Diff(ng.y) - a[1].Diff(ng.z),
                   a[0].Diff(ng.z) - a[2].Diff(ng.x),
                   a[1].Diff(ng.x) - a[0].Diff(ng.y)))
    h = -grad_phi
    s = b + mu * grad_phi
    return dict(mu=mu, nu=1 / mu, phi=phi, grad_phi=grad_phi, h=h, b=b, s=s)


def torus(ng, dim, n):
    from ngsolve.meshes import MakeStructured2DMesh, MakeStructured3DMesh
    if dim == 2:
        return MakeStructured2DMesh(quads=False, nx=n, ny=n, periodic_x=True, periodic_y=True)
    return MakeStructured3DMesh(hexes=False, nx=n, ny=n, nz=n,
                                periodic_x=True, periodic_y=True, periodic_z=True)


def verify_truth(ng, mesh, dim, t, order=24):
    """Check the analytic equations on the torus by sampling and integration."""
    rng = np.random.default_rng(20261108)
    pts = rng.random((40, dim))
    div_b = sum(t["b"][i].Diff(v) for i, v in enumerate((ng.x, ng.y, ng.z)[:dim]))
    resid_div = max(abs(div_b(mesh(*p))) for p in pts)
    # periodicity of phi, B, mu, S: compare x with x + e_i (kept inside [0,1) by shifting back)
    per = 0.0
    for p in pts:
        for i in range(dim):
            q = p.copy()
            q[i] = (q[i] + 0.5) % 1.0
            r = q.copy()
            r[i] = (r[i] + 0.5) % 1.0  # back to p modulo 1 via two half steps
            for key in ("phi", "b", "mu", "s"):
                a = np.atleast_1d(t[key](mesh(*p)))
                c = np.atleast_1d(t[key](mesh(*r)))
                per = max(per, float(np.max(np.abs(a - c))))
    means = {key: np.atleast_1d(ng.Integrate(t[key], mesh, order=order)).tolist()
             for key in ("phi", "h", "b")}
    scale_b = sqrt(ng.Integrate(t["b"] * t["b"], mesh, order=order))
    scale_h = sqrt(ng.Integrate(t["h"] * t["h"], mesh, order=order))
    return dict(div_b_max=resid_div, periodic_wrap_max=per, means=means,
                b_l2=scale_b, h_l2=scale_h)


def _pin_first(fd):
    for i in range(len(fd)):
        if fd[i]:
            fd.Clear(i)
            return i
    raise RuntimeError("no free dof to pin")


def solve_omega(ng, mesh, t, p, bonus):
    space = ng.Periodic(ng.H1(mesh, order=p))
    u, v = space.TnT()
    a = ng.BilinearForm(space, symmetric=True)
    a += t["mu"] * ng.grad(u) * ng.grad(v) * ng.dx(bonus_intorder=bonus)
    f = ng.LinearForm(space)
    f += t["s"] * ng.grad(v) * ng.dx(bonus_intorder=bonus)
    started = perf_counter()
    a.Assemble(); f.Assemble()
    fd = ng.BitArray(space.FreeDofs()); _pin_first(fd)
    gf = ng.GridFunction(space)
    gf.vec.data = a.mat.Inverse(fd, inverse="sparsecholesky") * f.vec
    seconds = perf_counter() - started
    res = f.vec.CreateVector(); res.data = f.vec - a.mat * gf.vec
    free = np.array([fd[i] for i in range(len(fd))])
    rel = float(np.linalg.norm(np.array(res)[free]) / max(np.linalg.norm(np.array(f.vec)[free]), 1e-300))
    grad_phi = ng.grad(gf)
    h = -grad_phi
    return dict(h=h, b=t["mu"] * h + t["s"], grad_phi=grad_phi, ndof=space.ndof,
                seconds=seconds, linear_residual=rel)


def solve_a(ng, mesh, dim, t, p, bonus):
    started = perf_counter()
    if dim == 2:
        space = ng.Periodic(ng.H1(mesh, order=p))
        u, w = space.TnT()
        rot = lambda g: ng.CF((ng.grad(g)[1], -ng.grad(g)[0]))
        a = ng.BilinearForm(space, symmetric=True)
        a += t["nu"] * rot(u) * rot(w) * ng.dx(bonus_intorder=bonus)
        f = ng.LinearForm(space)
        f += t["nu"] * t["s"] * rot(w) * ng.dx(bonus_intorder=bonus)
        a.Assemble(); f.Assemble()
        fd = ng.BitArray(space.FreeDofs()); _pin_first(fd)
        gf = ng.GridFunction(space)
        gf.vec.data = a.mat.Inverse(fd, inverse="sparsecholesky") * f.vec
        b = rot(gf)
        ndof = space.ndof
    else:
        vspace = ng.Periodic(ng.HCurl(mesh, order=p, nograds=True))
        qspace = ng.Periodic(ng.H1(mesh, order=1))
        nspace = ng.NumberSpace(mesh)
        space = vspace * qspace * nspace * nspace * nspace
        (u, q, c0, c1, c2), (w, r, d0, d1, d2) = space.TnT()
        dx = ng.dx(bonus_intorder=bonus)
        a = ng.BilinearForm(space)
        a += t["nu"] * ng.curl(u) * ng.curl(w) * dx
        a += (ng.grad(q) * w + ng.grad(r) * u) * dx
        a += (c0 * w[0] + c1 * w[1] + c2 * w[2] + d0 * u[0] + d1 * u[1] + d2 * u[2]) * dx
        f = ng.LinearForm(space)
        f += t["nu"] * t["s"] * ng.curl(w) * dx
        a.Assemble(); f.Assemble()
        fd = ng.BitArray(space.FreeDofs())
        offset = vspace.ndof
        for i in range(qspace.ndof):  # pin the multiplier constant
            if fd[offset + i]:
                fd.Clear(offset + i)
                break
        gf = ng.GridFunction(space)
        gf.vec.data = a.mat.Inverse(fd, inverse="umfpack") * f.vec
        b = ng.curl(gf.components[0])
        ndof = vspace.ndof
    seconds = perf_counter() - started
    res = f.vec.CreateVector(); res.data = f.vec - a.mat * gf.vec
    free = np.array([fd[i] for i in range(len(fd))])
    rel = float(np.linalg.norm(np.array(res)[free]) / max(np.linalg.norm(np.array(f.vec)[free]), 1e-300))
    return dict(b=b, h=t["nu"] * (b - t["s"]), ndof=ndof, seconds=seconds, linear_residual=rel)


def metrics(ng, mesh, dim, t, om, av, order):
    I = lambda cf: float(ng.Integrate(cf, mesh, order=order))
    eo = om["grad_phi"] - t["grad_phi"]
    eb = av["b"] - t["b"]
    err_o = I(t["mu"] * eo * eo)
    err_a = I(t["nu"] * eb * eb)
    c = I(t["nu"] * t["s"] * t["s"])
    j_o = 0.5 * I(t["mu"] * om["grad_phi"] * om["grad_phi"]) - I(t["s"] * om["grad_phi"])
    j_a = 0.5 * I(t["nu"] * (av["b"] - t["s"]) * (av["b"] - t["s"]))
    q_o, q_a = c + 2 * j_o, c - 2 * j_a
    energy = I(t["mu"] * t["grad_phi"] * t["grad_phi"])
    q_star = c - energy
    # kappa: load minus self-energy of each discrete solution (zero for exact Galerkin)
    kappa_o = I(t["s"] * om["grad_phi"]) - I(t["mu"] * om["grad_phi"] * om["grad_phi"])
    kappa_a = I(t["nu"] * t["s"] * av["b"]) - I(t["nu"] * av["b"] * av["b"])
    cross = I(eo * eb)  # Prager-Synge: grad error against curl error, 0 on the torus
    dh = om["h"] - av["h"]
    gap_norm = I(t["mu"] * dh * dh)
    mean_h_a = np.atleast_1d(ng.Integrate(av["h"], mesh, order=order)).tolist()
    mean_b_o = np.atleast_1d(ng.Integrate(om["b"], mesh, order=order)).tolist()
    mid, two = 0.5 * (q_a + q_o), (2 * q_a + q_o) / 3
    return dict(
        readout_order=order, energy_star=energy, q_star=q_star, q_a=q_a, q_omega=q_o,
        err2_omega=err_o, err2_a=err_a,
        rel_field_err_omega=sqrt(max(err_o, 0) / energy), rel_field_err_a=sqrt(max(err_a, 0) / energy),
        r=err_o / err_a if err_a > 0 else None,
        identity_omega=(q_o - q_star) - err_o, identity_a=(q_star - q_a) - err_a,
        kappa_omega=kappa_o, kappa_a=kappa_a, cross_orthogonality=cross,
        gap=q_o - q_a, gap_minus_norm=(q_o - q_a) - gap_norm,
        mean_h_a=mean_h_a, mean_b_omega=mean_b_o,
        midpoint_rel_err=(mid - q_star) / q_star, two_to_one_rel_err=(two - q_star) / q_star,
        a_rel_err=(q_a - q_star) / q_star, omega_rel_err=(q_o - q_star) / q_star,
    )


def run_level(ng, case, n, p, bonus, extras):
    case_id, dim, label, fields = case
    mesh = torus(ng, dim, n)
    t = truth(ng, dim, fields)
    with ng.TaskManager():  # one study level owns the CPU parallel region
        om = solve_omega(ng, mesh, t, p, bonus)
        av = solve_a(ng, mesh, dim, t, p, bonus)
        reads = [metrics(ng, mesh, dim, t, om, av, 2 * p + e) for e in extras]
    fine, coarse = reads[-1], reads[0]
    gap = abs(fine["gap"]) if fine["gap"] else 0.0
    kappa_shift = max(abs(fine["kappa_a"] - coarse["kappa_a"]), abs(fine["kappa_omega"] - coarse["kappa_omega"]),
                      abs(fine["q_a"] - coarse["q_a"]), abs(fine["q_omega"] - coarse["q_omega"]))
    return dict(case=case_id, dim=dim, n=n, p=p, elements=mesh.ne,
                ndof_omega=om["ndof"], ndof_a=av["ndof"],
                seconds_omega=om["seconds"], seconds_a=av["seconds"],
                linear_residual_omega=om["linear_residual"], linear_residual_a=av["linear_residual"],
                readouts=reads,
                kappa_readout_shift=kappa_shift,
                kappa_stable=bool(gap > 0 and kappa_shift <= 1e-3 * gap),
                floored=bool(min(fine["rel_field_err_omega"], fine["rel_field_err_a"]) < 1e-11))


def _slope(xs, ys):
    xs, ys = np.asarray(xs, float), np.asarray(ys, float)
    return float(np.polyfit(xs, ys, 1)[0])


def summarize(case_id, p_rows, h_rows):
    def usable(rows):
        return [r for r in rows if not r["floored"] and r["kappa_stable"]]
    out = dict(case=case_id)
    # p ladder: ln(err) vs p, p >= 2
    pr = [r for r in usable(p_rows) if r["p"] >= 2]
    if len(pr) >= 2:
        ps = [r["p"] for r in pr]
        so = -_slope(ps, [log(r["readouts"][-1]["rel_field_err_omega"]) for r in pr])
        sa = -_slope(ps, [log(r["readouts"][-1]["rel_field_err_a"]) for r in pr])
        out.update(p_decay_omega=so, p_decay_a=sa,
                   p_equal=bool(abs(so - sa) <= 0.10 * max(abs(so), abs(sa))))
    else:
        out.update(p_decay_omega=None, p_decay_a=None, p_equal=None)
    hr = usable(h_rows)
    seq_o, seq_a = [], []
    for c, f in zip(hr, hr[1:]):
        ratio = f["n"] / c["n"]
        seq_o.append(log(c["readouts"][-1]["rel_field_err_omega"] / f["readouts"][-1]["rel_field_err_omega"]) / log(ratio))
        seq_a.append(log(c["readouts"][-1]["rel_field_err_a"] / f["readouts"][-1]["rel_field_err_a"]) / log(ratio))
    out.update(h_orders_omega=seq_o, h_orders_a=seq_a,
               h_equal=bool(abs(seq_o[-1] - seq_a[-1]) <= 0.15) if seq_o else None)
    flags = [out["p_equal"], out["h_equal"]]
    if None in flags:
        out["verdict"] = "undetermined"
    else:
        out["verdict"] = {2: "equal", 1: "mixed", 0: "different"}[sum(flags)]
    out["excluded_levels"] = [dict(n=r["n"], p=r["p"], floored=r["floored"], kappa_stable=r["kappa_stable"])
                              for r in p_rows + h_rows if r["floored"] or not r["kappa_stable"]]
    return out


def reference_blend_test(p_rows, h_rows, coarse, reference):
    """User procedure: freeze Q_ref from the converged bracket, then judge blends
    on the predeclared coarse level.  The analytic Q* only checks Q_ref."""
    by = {(r["n"], r["p"]): r for r in p_rows + h_rows}
    ref_row = by[(reference["n"], reference["p"])]
    ref = ref_row["readouts"][-1]
    q_ref = 0.5 * (ref["q_a"] + ref["q_omega"])
    hw = 0.5 * (ref["q_omega"] - ref["q_a"])
    prev = by.get((reference["n"], reference["p"] - 1))
    drift = None
    if prev:
        pv = prev["readouts"][-1]
        drift = dict(a=ref["q_a"] - pv["q_a"], omega=ref["q_omega"] - pv["q_omega"])
    h_top = max(h_rows, key=lambda r: r["n"])["readouts"][-1]
    frozen = dict(q_ref=q_ref, half_width=hw, endpoints=[ref["q_a"], ref["q_omega"]],
                  last_step_drift=drift, readout_stable=ref_row["kappa_stable"],
                  h_top_bracket=[h_top["q_a"], h_top["q_omega"]],
                  h_top_contains_q_ref=bool(h_top["q_a"] - hw <= q_ref <= h_top["q_omega"] + hw),
                  analytic_check=dict(q_star=ref["q_star"], q_ref_minus_q_star=q_ref - ref["q_star"],
                                      q_star_inside=bool(ref["q_a"] <= ref["q_star"] <= ref["q_omega"])))
    c = by[(coarse["n"], coarse["p"])]["readouts"][-1]
    est = dict(A=c["q_a"], Omega=c["q_omega"], midpoint=0.5 * (c["q_a"] + c["q_omega"]),
               two_to_one=(2 * c["q_a"] + c["q_omega"]) / 3)
    errors = {k: dict(signed=v - q_ref, relative=(v - q_ref) / q_ref, abs=abs(v - q_ref),
                      signed_vs_analytic=(v - ref["q_star"]) / ref["q_star"]) for k, v in est.items()}
    diff = errors["two_to_one"]["abs"] - errors["midpoint"]["abs"]
    if abs(diff) <= 2 * hw:
        closer = "inconclusive"
    else:
        closer = "midpoint" if diff > 0 else "2:1"
    converged = bool(drift is not None and ref_row["kappa_stable"] and max(abs(drift["a"]), abs(drift["omega"]))
                     <= 0.1 * abs(est["midpoint"] - est["two_to_one"]))
    if not converged and closer != "inconclusive":
        closer_label = closer + " (reference convergence not established)"
    else:
        closer_label = closer
    den = q_ref - c["q_a"]
    return dict(reference=frozen, reference_converged=converged, coarse=coarse,
                coarse_errors=errors, closer=closer_label,
                abs_error_difference=diff,
                r_coarse=(c["q_omega"] - q_ref) / den if den > 0 else None)
