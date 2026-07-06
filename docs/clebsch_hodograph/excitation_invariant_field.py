r"""Self-contained compute helpers for excitation_invariant_field.ipynb.

Excitation-INVARIANT flux lines: keep the SAME field-line pattern as the drive current
rises (NOT a cyclotron, where the field is meant to change).  This module is imported by
the docs notebook and by the golden test tests/feec/test_excitation_invariant_field.py;
it carries its own small (s,z) end-pack geometry so the docs artifact stands alone.

Physics: below the iron knee the magnet is a LINEAR magnetostatic system, so scaling the
excitation by alpha scales B everywhere by alpha -- the flux-LINE pattern (streamlines
b_hat = B/|B|) is IDENTICAL, only the amplitude grows.  So "same flux lines as the current
rises" is AUTOMATIC in the linear regime; the ONLY thing that can move the pattern is the
nonlinearity, i.e. iron SATURATION (mu(|B|) dropping first at the pole-tip corner).  The
measurable quantity is the flux-line DIRECTION drift over an air region,

    D_dir(I) = rms_x || b_hat(x; I) - b_hat(x; I_lin) ||,

which is 0 while linear and grows once the iron saturates.  The end chamfer that relieves
the corner keeps the flux lines invariant several-fold DEEPER into saturation.

Relation to the hodograph design method (docs/clebsch_hodograph/DESIGN_METHODOLOGY.md):
the hodograph reads the pole as an equipotential (level set) of a prescribed LINEAR
(harmonic) potential -- so in the linear regime it produces exactly the designed field,
whose flux lines are excitation-invariant BY LINEARITY.  This module is the forward
characterization of how deep into saturation that designed pattern stays invariant, and
the corner-relief lever that extends it (the same lever the 2D Chaplygin hodograph
linearises analytically).
"""
import math
import os

import numpy as np

MU0 = 4.0e-7 * math.pi
MUR0 = 2000.0            # Froehlich unsaturated mu_r
BK = 1.2                 # Froehlich knee (T)

# ---- (s, z) end-pack geometry (meters); s = beam, z = gap (upper half) ----
GAP = 0.024              # SMALL gap -> the iron reluctance matters (saturation bites)
G2 = GAP / 2.0
S_POLE = 0.120           # iron half-length along the beam (terminates at +-s_pole)
S_BODY = 0.085           # flat body half-length (chamfer acts on s_body..s_pole)
POLE_T = 0.035           # iron thickness above the gap face
S_AIR = 0.230            # air half-extent along the beam (past the ends)
Z_AIR = 0.130            # air height (z, upper half)


def build_endpack(depth=0.0, exponent=1.0, maxh=None):
    """Upper-half (s, z) end pack: an air gap under a Froehlich-iron pole that TERMINATES
    at |s| = s_pole, with a parametrized end chamfer
        z_face(s) = g/2 + depth * ((|s| - s_body)/(s_pole - s_body))^exponent   (|s|>s_body)
    exponent < 1 = concave (early lift), > 1 = convex (late lift).  Boundaries: 'median'
    (z=0, the dipole antisymmetry plane, phi=0), 'irontop' (the driven pole back, phi=mmf).
    """
    import ngsolve as ng
    from netgen.occ import WorkPlane, OCCGeometry, Glue

    s_body = S_BODY
    z_back = G2 + POLE_T
    n_face = 61

    def face(s):
        a = abs(s)
        if a <= s_body or depth <= 0:
            return G2
        u = (a - s_body) / (S_POLE - s_body)
        return G2 + depth * (min(max(u, 0.0), 1.0)) ** exponent

    ss = np.linspace(-S_POLE, S_POLE, n_face)
    zf = [face(s) for s in ss]

    box = WorkPlane().MoveTo(-S_AIR, 0.0).Rectangle(2 * S_AIR, Z_AIR).Face()
    wp = WorkPlane().MoveTo(float(ss[0]), float(zf[0]))
    for s, z in zip(ss[1:], zf[1:]):
        wp.LineTo(float(s), float(z))
    wp.LineTo(float(ss[-1]), z_back)
    wp.LineTo(float(ss[0]), z_back)
    wp.Close()
    iron = wp.Face()
    air = box - iron
    air.faces.name = "air"
    iron.faces.name = "iron"
    shape = Glue([air, iron])
    for e in shape.edges:
        c = e.center
        if abs(c.y) < 1e-9:
            e.name = "median"
        elif abs(c.y - z_back) < 1e-9 and abs(c.x) < S_POLE + 1e-9:
            e.name = "irontop"
        elif abs(abs(c.x) - S_AIR) < 1e-9:
            e.name = "s_wall"
        elif abs(c.y - Z_AIR) < 1e-9:
            e.name = "top"
        else:
            e.name = "poleface"
    if maxh is None:
        maxh = GAP / 2.2
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=maxh))
    return mesh


def air_grid():
    """Fixed (s,z) points GUARANTEED in air for any (depth, exponent) in the search range:
    under the gap face (z <= 0.78*G2 stays below the >=G2 pole face) plus a fringe patch
    just past the pole end (where the flux lines bend most, so where saturation rotates
    them most)."""
    pts = []
    for s in np.linspace(0.0, 1.32 * S_POLE, 34):
        for z in np.linspace(0.30 * G2, 0.78 * G2, 4):
            pts.append((s, z))
    for s in np.linspace(1.02 * S_POLE, 1.34 * S_POLE, 8):
        for z in np.linspace(0.30 * G2, 1.5 * G2, 5):
            pts.append((s, z))
    return np.array(pts)


def corner_kappa(mesh, B_drive, order=3, relax=0.5, tol=1e-4, maxit=60):
    """Pole-tip corner concentration kappa = peak iron |B| / body |B| at the saturated
    drive -- a smooth L^10 norm over the END iron (|s| > 0.75*s_body), a differentiable
    proxy for the corner hot spot.  Context for the flux-line drift (the corner is what
    rotates the fringe field lines)."""
    from ngsolve import (H1, L2, BilinearForm, GridFunction, grad, dx, CF, Norm,
                         Integrate, IfPos, TaskManager, x as s_coord)
    mmf = B_drive * GAP / (2.0 * MU0)
    z_probe = 0.12 * G2
    ss = np.linspace(-1.4 * S_POLE, 1.4 * S_POLE, 241)
    with TaskManager():
        fes = H1(mesh, order=order, dirichlet="median|irontop")
        u, v = fes.TnT()
        gfu = GridFunction(fes)
        bccf = mesh.BoundaryCF({"irontop": mmf, "median": 0.0}, default=0.0)
        fes_mu = L2(mesh, order=0)
        mu_gf = GridFunction(fes_mu)
        mu_gf.Set(mesh.MaterialCF({"iron": MUR0}, default=1.0))
        iron_ind = mesh.MaterialCF({"iron": 1.0}, default=0.0)
        for _ in range(maxit):
            a = BilinearForm(fes)
            a += mu_gf * grad(u) * grad(v) * dx
            a.Assemble()
            gfu.Set(bccf, definedon=mesh.Boundaries("median|irontop"))
            r = gfu.vec.CreateVector()
            r.data = -a.mat * gfu.vec
            gfu.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * r
            B = MU0 * mu_gf * Norm(grad(gfu))
            froh = 1.0 + (MUR0 - 1.0) / (1.0 + (B / BK) ** 2)
            mu_t = (1.0 - iron_ind) * CF(1.0) + iron_ind * froh
            mu_n = GridFunction(fes_mu)
            mu_n.Set(mu_t)
            d = mu_n.vec.CreateVector()
            d.data = mu_n.vec - mu_gf.vec
            resid = d.Norm() / (mu_gf.vec.Norm() or 1.0)
            mu_gf.vec.data += relax * d
            if resid < tol:
                break
        g = grad(gfu)
        Bz = np.array([-MU0 * g(mesh(float(s), z_probe))[1] for s in ss])
        Bmag = MU0 * mu_gf * Norm(grad(gfu))
        end_ind = iron_ind * IfPos(s_coord * s_coord - (0.75 * S_BODY) ** 2, 1.0, 0.0)
        ev = float(Integrate(end_ind, mesh))
        n_pk = 10.0
        peak = (float(Integrate(end_ind * Bmag ** n_pk, mesh)) / ev) ** (1.0 / n_pk) \
            if ev > 0 else 0.0
    body = np.abs(ss) < 0.5 * S_BODY
    Bz_body = float(np.mean(Bz[body]))
    return float(peak / abs(Bz_body)) if abs(Bz_body) > 1e-30 else 0.0


def solve_bhat(mesh, B_drive, pts, order=3, saturate=True, relax=0.5, tol=1e-4, maxit=60):
    """Nonlinear scalar-potential Picard solve (Froehlich mu(|B|)).  Return the unit
    field-line direction b_hat = B/|B| = -grad(psi)/|grad(psi)| (mu_r=1 in air) at the air
    sampling points, plus the iron <mu_r>.  saturate=False forces mu = MUR0 (the linear
    reference, whose b_hat is DRIVE-INDEPENDENT by linearity)."""
    from ngsolve import (H1, L2, BilinearForm, GridFunction, grad, dx, CF, Norm,
                         Integrate, TaskManager)
    mmf = B_drive * GAP / (2.0 * MU0)
    with TaskManager():
        fes = H1(mesh, order=order, dirichlet="median|irontop")
        u, v = fes.TnT()
        gfu = GridFunction(fes)
        bccf = mesh.BoundaryCF({"irontop": mmf, "median": 0.0}, default=0.0)
        fes_mu = L2(mesh, order=0)
        mu_gf = GridFunction(fes_mu)
        mu_gf.Set(mesh.MaterialCF({"iron": MUR0}, default=1.0))
        iron_ind = mesh.MaterialCF({"iron": 1.0}, default=0.0)
        for _ in range(maxit):
            a = BilinearForm(fes)
            a += mu_gf * grad(u) * grad(v) * dx
            a.Assemble()
            gfu.Set(bccf, definedon=mesh.Boundaries("median|irontop"))
            r = gfu.vec.CreateVector()
            r.data = -a.mat * gfu.vec
            gfu.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * r
            if not saturate:
                break
            B = MU0 * mu_gf * Norm(grad(gfu))
            froh = 1.0 + (MUR0 - 1.0) / (1.0 + (B / BK) ** 2)
            mu_t = (1.0 - iron_ind) * CF(1.0) + iron_ind * froh
            mu_n = GridFunction(fes_mu)
            mu_n.Set(mu_t)
            d = mu_n.vec.CreateVector()
            d.data = mu_n.vec - mu_gf.vec
            resid = d.Norm() / (mu_gf.vec.Norm() or 1.0)
            mu_gf.vec.data += relax * d
            if resid < tol:
                break
        iv = float(Integrate(iron_ind, mesh))
        mur_mean = float(Integrate(iron_ind * mu_gf, mesh)) / iv if iv > 0 else 0.0
        g = grad(gfu)
        bh = np.zeros((len(pts), 2))
        for i, (s, z) in enumerate(pts):
            gg = g(mesh(float(s), float(z)))
            b = np.array([-gg[0], -gg[1]])                 # B ~ -grad(psi) in air
            n = np.hypot(b[0], b[1])
            bh[i] = b / n if n > 1e-30 else 0.0
    return bh, float(mur_mean)


def _drift(bh, ref):
    """rms flux-line direction difference (radians, small-angle ~ chord length)."""
    return float(np.sqrt(np.mean(np.sum((bh - ref) ** 2, axis=1))))


def flux_line_drift(depth, exponent, B_lin=0.15, B_sat=1.70, order=3, maxh=None, pts=None):
    """The OBJECTIVE: the flux-line DIRECTION drift between a low (linear-regime) drive and
    a high (saturated) drive.  2 solves.  Small = the flux lines keep their shape as the
    current is turned up.  Also returns the saturated iron <mu_r> and the pole-tip corner
    kappa (context: the corner is what drives the drift)."""
    mesh = build_endpack(depth=depth, exponent=exponent, maxh=maxh)
    if pts is None:
        pts = air_grid()
    bh_ref, mur_ref = solve_bhat(mesh, B_lin, pts, order=order, saturate=True)
    bh_sat, mur_sat = solve_bhat(mesh, B_sat, pts, order=order, saturate=True)
    kap = corner_kappa(mesh, B_sat, order=order)
    return {
        "depth_m": float(depth), "exponent": float(exponent),
        "D_dir": _drift(bh_sat, bh_ref),
        "mur_ref": mur_ref, "mur_sat": mur_sat, "corner_kappa_sat": kap,
    }


def invariance_curve(depth, exponent, drives=(0.15, 0.45, 0.75, 1.05, 1.35, 1.70),
                     order=3, maxh=None):
    """The full D_dir(B_drive) SATURATED sweep (reporting / figure), plus the LINEAR
    control (mu forced constant -> D_dir ~ 0 at every drive: linearity => invariant flux
    lines; saturation is the sole breaker).  Reference = the lowest drive."""
    mesh = build_endpack(depth=depth, exponent=exponent, maxh=maxh)
    pts = air_grid()
    drives = list(drives)
    bh_ref_sat, _ = solve_bhat(mesh, drives[0], pts, order=order, saturate=True)
    bh_ref_lin, _ = solve_bhat(mesh, drives[0], pts, order=order, saturate=False)
    D_sat, D_lin, murs = [], [], []
    for B in drives:
        bh_s, mur_s = solve_bhat(mesh, B, pts, order=order, saturate=True)
        bh_l, _ = solve_bhat(mesh, B, pts, order=order, saturate=False)
        D_sat.append(_drift(bh_s, bh_ref_sat))
        D_lin.append(_drift(bh_l, bh_ref_lin))
        murs.append(mur_s)
    return {
        "depth_m": float(depth), "exponent": float(exponent),
        "drives_T": drives, "D_dir_saturated": D_sat, "D_dir_linear_control": D_lin,
        "mur_sat_sweep": murs, "D_dir_max": float(max(D_sat)),
        "linear_control_max": float(max(D_lin)),
    }


def optimize(trials=24, B_lin=0.15, B_sat=1.70, order=3, seed=0):
    """Minimize the saturated flux-line DIRECTION drift D_dir(B_lin -> B_sat) over the end
    chamfer (depth, exponent) -- the shape that keeps the flux lines invariant deepest into
    saturation -- vs the flat-cut baseline.  Optuna (TPE) if available, else a coarse grid
    + local refine.  D_dir is monotone in the drive, so the top-drive drift is the
    worst-case invariance error."""
    depth_lo, depth_hi = 0.0, 0.95 * G2
    exp_lo, exp_hi = 0.4, 3.0
    pts = air_grid()
    flat = flux_line_drift(0.0, 1.0, B_lin=B_lin, B_sat=B_sat, order=order, pts=pts)
    history = []

    def obj(depth, exponent):
        r = flux_line_drift(depth, exponent, B_lin=B_lin, B_sat=B_sat, order=order, pts=pts)
        history.append({"depth_m": r["depth_m"], "exponent": r["exponent"],
                        "D_dir": r["D_dir"], "corner_kappa_sat": r["corner_kappa_sat"]})
        return r

    def cost(r):
        return r["D_dir"]

    best, used = None, "grid"
    try:
        import optuna
        optuna.logging.set_verbosity(optuna.logging.WARNING)
        used = "optuna_tpe"

        def _obj(trial):
            d = trial.suggest_float("depth", depth_lo, depth_hi)
            e = trial.suggest_float("exponent", exp_lo, exp_hi)
            return cost(obj(d, e))

        study = optuna.create_study(direction="minimize",
                                    sampler=optuna.samplers.TPESampler(seed=seed))
        study.optimize(_obj, n_trials=trials, show_progress_bar=False)
        bp = study.best_params
        best = flux_line_drift(bp["depth"], bp["exponent"], B_lin=B_lin, B_sat=B_sat,
                               order=order, pts=pts)
    except Exception:
        used = "grid_refine"
        depths = np.linspace(depth_lo, depth_hi, 5)
        exps = np.linspace(exp_lo, exp_hi, 5)
        grid = [obj(float(d), float(e)) for d in depths for e in exps]
        best = min(grid, key=cost)
        d0, e0 = best["depth_m"], best["exponent"]
        for dd in (d0 - 0.15 * G2, d0, d0 + 0.15 * G2):
            for ee in (e0 - 0.4, e0, e0 + 0.4):
                if depth_lo <= dd <= depth_hi and exp_lo <= ee <= exp_hi:
                    r = obj(float(dd), float(ee))
                    if cost(r) < cost(best):
                        best = r

    inv_factor = (flat["D_dir"] / best["D_dir"]) if best["D_dir"] > 0 else float("inf")
    return {
        "optimizer": used, "n_evals": len(history),
        "B_lin_T": float(B_lin), "B_sat_T": float(B_sat), "B_K_iron_T": float(BK),
        "flat_cut": flat, "optimized": best,
        "invariance_factor": float(inv_factor),
        "flux_lines_more_invariant": bool(best["D_dir"] < flat["D_dir"]),
        "flat_flux_lines_already_invariant": bool(flat["D_dir"] < 1e-2),
        "history": history,
    }


def make_figure(res, flat_curve, best_curve, png_path):
    """3-panel figure: (L) D_dir vs excitation flat/optimized + linear control; (M) the
    top-drive invariance bar; (R) the optimization scatter."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(15.6, 4.3), dpi=140)
    fc, oc = res["flat_cut"], res["optimized"]

    dr = np.array(flat_curve["drives_T"])
    ax[0].plot(dr, np.array(flat_curve["D_dir_saturated"]) * 1e3, "C3-o", lw=1.6, ms=4,
               label="flat cut (saturated)")
    ax[0].plot(dr, np.array(best_curve["D_dir_saturated"]) * 1e3, "C0-o", lw=1.6, ms=4,
               label="optimized (saturated)")
    ax[0].plot(dr, np.array(flat_curve["D_dir_linear_control"]) * 1e3, "k:", lw=1.2,
               label="linear control (mu=const)")
    ax[0].axvline(BK, color="0.6", lw=0.8, ls="--")
    ax[0].text(BK, ax[0].get_ylim()[1] * 0.92, " iron knee", fontsize=7.5, color="0.4")
    ax[0].set_xlabel("drive  $B_{body}$ [T]  (excitation)")
    ax[0].set_ylabel("flux-line direction drift  $D_{dir}$  [mrad]")
    ax[0].set_title("Flux lines stay invariant while linear,\n"
                    "drift once the iron saturates (control ~ 0)")
    ax[0].legend(fontsize=7.5)

    x = np.arange(2)
    dd = [fc["D_dir"] * 1e3, oc["D_dir"] * 1e3]
    bars = ax[1].bar(x, dd, color=["C3", "C2"])
    for b_, d_ in zip(bars, dd):
        ax[1].text(b_.get_x() + b_.get_width() / 2, d_, f"{d_:.2f}\nmrad",
                   ha="center", va="bottom", fontsize=9)
    ax[1].set_xticks(x); ax[1].set_xticklabels(["FLAT cut", "OPTIMIZED"])
    ax[1].set_ylabel(f"$D_{{dir}}$ at $B$={res['B_sat_T']:.1f} T  [mrad]")
    ax[1].set_ylim(0, max(dd) * 1.28)
    ax[1].set_title(f"End chamfer keeps flux lines invariant deeper\n"
                    f"({res['invariance_factor']:.1f}x smaller drift; "
                    f"d={oc['depth_m']*1e3:.1f} mm, p={oc['exponent']:.2f})")

    h = res["history"]
    d = np.array([x["depth_m"] for x in h]) * 1e3
    e = np.array([x["exponent"] for x in h])
    yy = np.array([x["D_dir"] for x in h]) * 1e3
    sc = ax[2].scatter(d, yy, c=e, cmap="viridis", s=28)
    ax[2].axhline(fc["D_dir"] * 1e3, color="C3", lw=1, ls="--",
                  label=f"flat cut {fc['D_dir']*1e3:.2f} mrad")
    ax[2].scatter([oc["depth_m"] * 1e3], [oc["D_dir"] * 1e3], marker="*", s=220,
                  color="C1", edgecolor="k", zorder=5, label="optimum")
    ax[2].set_xlabel("chamfer depth [mm]")
    ax[2].set_ylabel("$D_{dir}$ at top drive [mrad] (lower = more invariant)")
    ax[2].set_title(f"Optimization ({res['optimizer']}, {res['n_evals']} evals):\n"
                    "the flux-line-invariance-maximizing end chamfer")
    ax[2].legend(fontsize=8)
    fig.colorbar(sc, ax=ax[2], label="chamfer exponent")

    fig.tight_layout()
    fig.savefig(png_path, bbox_inches="tight")
    plt.close(fig)
    return os.path.abspath(png_path)
