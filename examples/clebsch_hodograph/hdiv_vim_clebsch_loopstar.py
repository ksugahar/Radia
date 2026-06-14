r"""HDiv-VIM <-> Clebsch: the loop modes of the demag operator ARE Clebsch magnetizations.

The de Rham CAPSTONE bridging the two threads of this directory's research line:
the HDiv-VIM demag SOLVER (operator side, `radia.hdiv_vim`) and the Clebsch
hodograph DESIGN line (potential side).  The Hodge / Helmholtz split of a
magnetization is

    M  =  grad(phi)     (+)     grad(alpha) x grad(beta)
          \________/             \___________________/
          gradient / "star"       Clebsch / "loop"
          carries the charge      charge-free  ->  FIELD-NULL in N = B^T G B
          pole-forming            flux-guiding (yoke return)

The HDiv-VIM is built so that its demag operator `N = B^T G B` (B = the magnetic
charge map M -> (rho = -div M, sigma = M.n), G = the charge Gram) annihilates the
divergence-free RT modes -- "loop modes are field-null by construction" (de Rham:
a curl carries no charge).  That property IS the statement that a **Clebsch**
magnetization `grad(alpha) x grad(beta) = d(alpha d beta)` (an exact, hence
closed, 2-form) makes **no demagnetizing field**.  So the Clebsch potentials
(alpha, beta) are the coordinates of the HDiv-VIM kernel; the loop-star split is
the Clebsch-gradient (Hodge) split is the pole / flux-guide split.

Verified here (unit sphere, HDiv order 1, exact analytic charge Gram):

  (1) FIELD-NULL   D(Clebsch) ~ 0   vs   D(gradient) ~ 1,   D(uniform) = 1/3.
                   The azimuthal Clebsch field M = (y,-x,0) = grad((x^2+y^2)/2) x grad(z)
                   has machine-zero discrete divergence; its only residual demag is
                   the faceting surface charge (-> 0 under mesh.Curve).

  (2) NO STRAY FIELD   the EXTERNAL field of the Clebsch M is ~0 (no charge -> no
                   field), tiny vs the gradient M's external field.

  (3) GAUGE / SUPERPOSITION   the external field of `grad(phi) + t * Clebsch` is
                   independent of t -- you can add any amount of flux-circulation
                   to a magnetization without changing the field it produces.

This is the linear (kinematic) level of the connection.  The nonlinear payoff
(the saturable HDiv-VIM solve reformulated in Clebsch/hodograph coordinates,
where the Chaplygin hodograph linearises saturation) is the next frontier; see
`saturation_loop_2d.py` / `chaplygin_hodograph_2d.py` for the hodograph side.

The de Rham view of this bridge is classical: P. Robert, "Clebsch Potentials and
the Visualization of Three-Dimensional Solenoidal Vector Fields", IEEE Trans.
Magn. 27(5), 1991 -- B is a 2-form, the loop modes are the closed 2-forms, and a
GLOBAL Clebsch pair B = grad(alpha) x grad(beta) exists iff the helicity vanishes.
Full references in docs/clebsch_hodograph/HDIV_VIM_CLEBSCH_BRIDGE.md.

run:  python hdiv_vim_clebsch_loopstar.py
"""
import os

import numpy as np
import ngsolve as ng
from netgen.csg import CSGeometry, Sphere, Pnt

from radia.hdiv_vim import DemagOperator, reconstruct_field_polynomial

x, y, z = ng.x, ng.y, ng.z

# the three test magnetizations (constructed, not solved):
CF_CLEBSCH = ng.CoefficientFunction((y, -x, 0))     # grad((x^2+y^2)/2) x grad(z) -- azimuthal, solenoidal
CF_GRAD = ng.CoefficientFunction((x, y, z))          # grad(|r|^2/2) -- radial, pure gradient (charged)
CF_UNIFORM = ng.CoefficientFunction((0, 0, 1))       # uniform -- surface charge only (sphere D = 1/3)


def _charges(fes, cf):
    """Discrete charge content of a magnetization: ||div M|| (volume) and ||M.n|| (surface)."""
    gf = ng.GridFunction(fes); gf.Set(cf)
    dv = float(ng.sqrt(ng.Integrate(ng.div(gf) ** 2, gf.space.mesh)))
    nrm = ng.specialcf.normal(3)
    sv = float(ng.sqrt(ng.Integrate(ng.InnerProduct(gf.Trace(), nrm) ** 2, gf.space.mesh, ng.BND)))
    return dv, sv


def analyze(maxh=0.5, order=1, ext_pts=((0, 0, 2.0), (2.0, 0, 0), (0, 2.0, 1.0))):
    """Bridge analysis on the unit sphere.  Returns the demag factors, charges, external-field ratio,
    and the gauge (add-Clebsch) field invariance."""
    g = CSGeometry(); g.Add(Sphere(Pnt(0, 0, 0), 1.0))
    obs = np.asarray(ext_pts, float)
    with ng.TaskManager():
        mesh = ng.Mesh(g.GenerateMesh(maxh=maxh))
        fes = ng.HDiv(mesh, order=order)
        op = DemagOperator(fes)

        # (1) field-null: the demag factor (Rayleigh quotient) of each magnetization
        D_cl = op.DemagFactor(CF_CLEBSCH)
        D_gr = op.DemagFactor(CF_GRAD)
        D_un = op.DemagFactor(CF_UNIFORM)
        div_cl, mn_cl = _charges(fes, CF_CLEBSCH)
        div_gr, mn_gr = _charges(fes, CF_GRAD)

        # (2) no stray field: external H from the charges.  The "star" reference is the UNIFORM
        # magnetization M = grad(z) -- a genuine charged gradient field with a DIPOLE external field
        # (the radial M=(x,y,z) is spherically symmetric, so by the shell theorem it ALSO has ~no
        # external field -- it would be a misleading stray-field reference).  The Clebsch M is
        # charge-free -> ~0 external field.
        gf_cl = ng.GridFunction(fes); gf_cl.Set(CF_CLEBSCH)
        gf_un = ng.GridFunction(fes); gf_un.Set(CF_UNIFORM)
        H_cl = reconstruct_field_polynomial(mesh, gf_cl, obs, quad=4)
        H_un = reconstruct_field_polynomial(mesh, gf_un, obs, quad=4)
        ext_cl = float(np.linalg.norm(H_cl)); ext_un = float(np.linalg.norm(H_un))

        # (3) gauge / superposition: external field of (uniform star) + t*Clebsch vs t=0 -- unchanged.
        gauge = []
        for t in (0.5, 1.0, 2.0, 5.0):
            gf = ng.GridFunction(fes); gf.Set(CF_UNIFORM + t * CF_CLEBSCH)
            H = reconstruct_field_polynomial(mesh, gf, obs, quad=4)
            gauge.append((t, float(np.linalg.norm(H - H_un) / ext_un)))
    return {
        "ne": int(mesh.GetNE(ng.VOL)), "order": order, "maxh": maxh,
        "D_clebsch": D_cl, "D_gradient": D_gr, "D_uniform": D_un,
        "ratio_D": abs(D_cl) / abs(D_gr),
        "divM_clebsch": div_cl, "Mn_clebsch": mn_cl,
        "divM_gradient": div_gr, "Mn_gradient": mn_gr,
        "ext_clebsch": ext_cl, "ext_uniform": ext_un, "ratio_ext": ext_cl / ext_un,
        "gauge": gauge,
    }


def _plot(r):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.4, 4.0), dpi=150)
    labels = ["Clebsch\n$\\nabla a\\times\\nabla b$", "gradient\n$\\nabla\\phi$", "uniform"]
    Ds = [abs(r["D_clebsch"]), abs(r["D_gradient"]), abs(r["D_uniform"])]
    bars = ax.bar(labels, Ds, color=["C0", "C3", "0.6"])
    ax.set_yscale("log"); ax.set_ylabel("demag factor  $D = m^T N m / m^T M_{mass} m$")
    ax.set_title("HDiv-VIM: Clebsch (loop) is field-null")
    ax.axhline(1.0 / 3.0, color="0.4", ls="--", lw=0.8)
    for b, d in zip(bars, Ds):
        ax.text(b.get_x() + b.get_width() / 2, d, f"{d:.1e}", ha="center", va="bottom", fontsize=8)
    ts = [0.0] + [t for t, _ in r["gauge"]]
    dev = [0.0] + [d for _, d in r["gauge"]]
    ax2.plot(ts, dev, "o-", color="C0")
    ax2.set_xlabel("$t$  (amount of Clebsch added to $\\nabla\\phi$)")
    ax2.set_ylabel("rel. change in external field")
    ax2.set_title("gauge: adding flux-circulation\nchanges nothing outside")
    ax2.set_ylim(bottom=0)
    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout(); fig.savefig(png, bbox_inches="tight"); plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    print("HDiv-VIM <-> Clebsch: the loop modes ARE Clebsch magnetizations\n")
    r = analyze()
    print(f"  unit sphere  ne={r['ne']}  HDiv order {r['order']}")
    print(f"  (1) FIELD-NULL  demag factor D:")
    print(f"        Clebsch  M=(y,-x,0) : D={r['D_clebsch']: .3e}  ||divM||={r['divM_clebsch']:.1e} "
          f" ||M.n||={r['Mn_clebsch']:.2e}")
    print(f"        gradient M=(x,y,z)  : D={r['D_gradient']: .3e}  ||divM||={r['divM_gradient']:.1e} "
          f" ||M.n||={r['Mn_gradient']:.2e}")
    print(f"        uniform  M=(0,0,1)  : D={r['D_uniform']: .3e}  (sphere analytic 1/3)")
    print(f"        -> D_Clebsch / D_gradient = {r['ratio_D']:.2e}  (the loop mode carries ~no demag)")
    print(f"  (2) NO STRAY FIELD  external |H|: Clebsch={r['ext_clebsch']:.2e}  uniform(dipole)="
          f"{r['ext_uniform']:.2e}  ratio={r['ratio_ext']:.2e}")
    print(f"  (3) GAUGE  external field change adding t*Clebsch to the uniform star:")
    for t, dev in r["gauge"]:
        print(f"        t={t:>4}: rel field change = {dev:.2e}")
    print("\n  => loop-star = Clebsch-gradient = Hodge = pole / flux-guide.  The Clebsch potentials")
    print("     coordinatize the HDiv-VIM kernel: the operator side and the potential side of one")
    print("     de Rham complex.")
    _plot(r)


if __name__ == "__main__":
    main()
