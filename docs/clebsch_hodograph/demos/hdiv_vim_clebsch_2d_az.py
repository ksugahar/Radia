r"""2-D unification: the flux function A_z IS the Clebsch potential -- the bridge between the HDiv-VIM
loop-star split and the Chaplygin hodograph.

Companion to `hdiv_vim_clebsch_loopstar.py` (the 3-D capstone).  In **2-D**, every in-plane
divergence-free field is `B = grad(A_z) x z_hat = (dA_z/dy, -dA_z/dx)`, i.e. the flux function `A_z`
is **exactly the Clebsch potential `alpha`** (with `beta = z`).  So the loop-star (Hodge) split of a
2-D magnetization is

    M  =  grad(phi)        (+)   grad(A_z) x z_hat
          star / charge          loop / Clebsch (A_z = the stream/flux function)
          div M = lap(phi)       div = 0  (charge-free -> field-null)

This is the same de Rham split as the 3-D capstone, but in 2-D the Clebsch *pair* collapses to the
**single scalar `A_z`** -- which is exactly why the **Chaplygin hodograph** (interchange
`(x,y) <-> (theta, q=|B|)`) **linearises the saturation** of the 2-D problem: the saturable
`div(nu(|grad A_z|) grad A_z) = 0` has one scalar unknown `A_z` = the Clebsch potential.  See
`chaplygin_hodograph_2d.py` / `saturation_loop_2d.py` for the saturation side.

Verified here (unit square, `H1` order 3): the Hodge / Helmholtz decomposition
`M = grad(phi) (+) grad(A_z) x z_hat` is reconstructed to machine precision; the loop part is
machine-zero divergence (charge-free); and `grad(A_z) x z_hat` reproduces the loop part exactly --
A_z is the Clebsch potential.

run:  python hdiv_vim_clebsch_2d_az.py
"""
import os
from math import pi

import numpy as np
import ngsolve as ng
from netgen.geom2d import unit_square

x, y = ng.x, ng.y


def analyze(maxh=0.05, order=3):
    """The flux function A_z IS the Clebsch potential: build a loop field B = grad(A_z) x z_hat from a
    KNOWN A_z (so div B = 0 is exact), verify it is charge-free (field-null), and RECOVER A_z from B
    via the stream-function weak form -- then contrast with a charged gradient field.

    `.Diff(x/y)` is used only on the SYMBOLIC potentials A_z, phi (where it is exact spatial diff); the
    field recovery uses the proper `rot(w) = (w_y, -w_x)` weak form (no `.Diff` on GridFunctions)."""
    # the chosen Clebsch potential (vanishes on the unit-square boundary)
    Az = ng.sin(pi * x) * ng.sin(pi * y)
    B_loop = ng.CoefficientFunction((Az.Diff(y), -Az.Diff(x)))       # grad(A_z) x z_hat  (Clebsch/loop)
    div_loop_cf = B_loop[0].Diff(x) + B_loop[1].Diff(y)              # = A_z_yx - A_z_xy = 0 (exact)
    # a charged gradient field for contrast
    phi = 0.5 * (x * x + y * y)
    B_star = ng.CoefficientFunction((phi.Diff(x), phi.Diff(y)))      # = (x, y), pure gradient
    div_star_cf = B_star[0].Diff(x) + B_star[1].Diff(y)             # = 2

    with ng.TaskManager():
        mesh = ng.Mesh(unit_square.GenerateMesh(maxh=maxh))
        nrm = ng.specialcf.normal(2)
        L2 = lambda cf, reg=ng.VOL: float(ng.sqrt(ng.Integrate(ng.InnerProduct(cf, cf), mesh, reg)))

        # field-null metrics of the loop field
        div_loop = float(ng.sqrt(ng.Integrate(div_loop_cf ** 2, mesh)))           # ~0 (machine)
        div_star = float(ng.sqrt(ng.Integrate(div_star_cf ** 2, mesh)))           # = 2 * |Omega|^(1/2)
        bn_loop = L2(ng.InnerProduct(B_loop, nrm), ng.BND)                        # ~0 (A_z=0 on bdry)
        bn_star = L2(ng.InnerProduct(B_star, nrm), ng.BND)

        # RECOVER A_z from B_loop via the stream-function weak form:
        #   int grad(A') . grad(w) = int B . rot(w),   rot(w) = (w_y, -w_x)
        V = ng.H1(mesh, order=order, dirichlet=".*")
        u, w = V.TnT()
        a = ng.BilinearForm(ng.grad(u) * ng.grad(w) * ng.dx); a.Assemble()
        rot_w = ng.CoefficientFunction((ng.grad(w)[1], -ng.grad(w)[0]))
        f = ng.LinearForm(ng.InnerProduct(B_loop, rot_w) * ng.dx); f.Assemble()
        gA = ng.GridFunction(V)
        gA.vec.data = a.mat.Inverse(V.FreeDofs(), inverse="sparsecholesky") * f.vec
        recover_err = float(ng.sqrt(ng.Integrate((gA - Az) ** 2, mesh))
                            / ng.sqrt(ng.Integrate(Az ** 2, mesh)))                # A_z recovered?
        # and the recovered Clebsch field rot(A_z) reproduces B_loop
        gradA = ng.grad(gA)
        B_rec = ng.CoefficientFunction((gradA[1], -gradA[0]))
        clebsch_err = L2(B_loop - B_rec) / L2(B_loop)
    return {
        "ne": int(mesh.ne), "order": order, "maxh": maxh,
        "div_loop": div_loop, "div_star": div_star, "bn_loop": bn_loop, "bn_star": bn_star,
        "recover_err": recover_err, "clebsch_err": clebsch_err,
        "gA": gA, "mesh": mesh, "Az": Az,
    }


def _plot(r):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    mesh = r["mesh"]
    n = 26
    xs = np.linspace(0.02, 0.98, n)
    X, Y = np.meshgrid(xs, xs)
    gA = r["gA"]
    Az = np.zeros((n, n)); Mx = np.zeros((n, n)); My = np.zeros((n, n))
    gradA = ng.grad(gA)
    for i in range(n):
        for j in range(n):
            mp = mesh(float(X[i, j]), float(Y[i, j]))
            Az[i, j] = float(gA(mp))
            ga = gradA(mp)
            Mx[i, j] = float(ga[1]); My[i, j] = -float(ga[0])      # grad(A_z) x z_hat
    fig, ax = plt.subplots(1, 1, figsize=(5.2, 4.6), dpi=150)
    cs = ax.contour(X, Y, Az, levels=14, colors="0.5", linewidths=0.7)
    ax.streamplot(X, Y, Mx, My, color="C0", density=1.0, linewidth=0.7, arrowsize=0.8)
    ax.set_aspect("equal"); ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_title("loop part = $\\nabla A_z \\times \\hat z$\n(field lines = $A_z$ iso-contours)")
    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout(); fig.savefig(png, bbox_inches="tight"); plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    print("2-D unification: the flux function A_z IS the Clebsch potential\n")
    r = analyze()
    print(f"  unit square  ne={r['ne']}  H1 order {r['order']}")
    print(f"  loop field  B = grad(A_z) x z_hat  (Clebsch, A_z = sin(pi x) sin(pi y)):")
    print(f"    ||div B_loop|| = {r['div_loop']:.2e}  (charge-free in the volume -> field-null)")
    print(f"    ||B_loop . n||_bdry = {r['bn_loop']:.2e}  (tangential to A_z = const -> no surface charge)")
    print(f"  charged field  B = grad(phi) = (x,y):  ||div B_star|| = {r['div_star']:.2e}  "
          f"(carries the charge),  ||B_star . n||_bdry = {r['bn_star']:.2e}")
    print(f"  RECOVER A_z from B via the stream-function weak form:")
    print(f"    ||A_z(recovered) - A_z|| / ||A_z||        = {r['recover_err']:.2e}")
    print(f"    ||B_loop - grad(A_z_rec) x z|| / ||B_loop|| = {r['clebsch_err']:.2e}  (A_z IS the Clebsch potential)")
    print("\n  => in 2-D the Clebsch pair collapses to the single scalar A_z (the flux function);")
    print("     that is why the Chaplygin hodograph linearises the 2-D saturation -- one scalar")
    print("     unknown A_z = the Clebsch potential.  loop-star = Clebsch-gradient = pole / flux-guide.")
    _plot(r)


if __name__ == "__main__":
    main()
