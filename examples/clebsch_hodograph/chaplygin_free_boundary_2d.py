"""Frontier 2 (the turning-guide free boundary): the hodograph IMAGE of a
turning flux guide -- a rectangle for constant width (1-D self-linearising),
a theta-DEPENDENT region (a FREE BOUNDARY) once the guide tapers.

RESEARCH example (track A, the open frontier of chaplygin_turning_guide_2d.py).
The turning-guide forward construction solved the linear Chaplygin PDE on a
FIXED hodograph rectangle and back-mapped to a physical patch.  The INVERSE
direction -- prescribe the physical guide, find its hodograph image -- is the
free-boundary problem.  This file does the achievable, honest piece: it
COMPUTES the hodograph image of a prescribed turning flux guide and shows
concretely

  - constant-width circular bend  -> image is a RECTANGLE in (q, theta):
    q in [q(r_out), q(r_in)] at EVERY position angle, theta spanning the bend.
    The q-extent does NOT depend on theta -> the image bounds are fixed ->
    this is the 1-D self-linearising case (|H| ~ 1/r is forced by geometry).
  - tapering bend (outer wall spirals inward) -> the gap width varies along the
    bend, so the q-extent VARIES with theta -> the image bounds are theta-
    DEPENDENT -> a genuine FREE BOUNDARY (the hodograph image is unknown a
    priori; recovering it from the physical boundary is the open inverse
    problem).

The hodograph image's boundary shape is a GEOMETRIC property of the field, so a
LINEAR flux solve (nu = 1) already exhibits it; saturation shifts the q-values
but not the rectangle-vs-free-boundary distinction.  A turning flux guide
(iron walls = flux lines) is bounded by const-A (the walls) and const-Phi
(inlet/outlet) curves; only when those happen to be const-q / const-theta lines
(constant width) is the image a rectangle.

run:  python chaplygin_free_boundary_2d.py
"""
import math
import os

from numpy import pi
import numpy as np
from ngsolve import (Mesh, H1, GridFunction, grad, InnerProduct, dx, CF, x, y,
                     sqrt, BilinearForm, TaskManager, Integrate)
from netgen.occ import WorkPlane, OCCGeometry


def _sector(r_in, r_out0, taper, phimax, n=80):
    """First-quadrant annular sector, inner radius r_in (constant), outer wall
    r_out(phi) = r_out0 (1 - taper * phi/phimax) (taper=0 -> circular bend).
    Boundary: inner arc + inlet edge + outer (possibly spiral) wall + outlet."""
    wp = WorkPlane()
    phis = np.linspace(0.0, phimax, n)
    # start at inner-arc, phi=0
    wp.MoveTo(r_in, 0.0)
    # inlet edge phi=0: r_in -> r_out0 (taper at phi=0 is 0)
    wp.LineTo(r_out0, 0.0)
    # outer wall phi: 0 -> phimax
    for ph in phis[1:]:
        ro = r_out0 * (1.0 - taper * ph / phimax)
        wp.LineTo(ro * math.cos(ph), ro * math.sin(ph))
    # outlet edge phi=phimax: r_out(phimax) -> r_in
    rmax = r_out0 * (1.0 - taper)
    wp.LineTo(r_in * math.cos(phimax), r_in * math.sin(phimax))
    # inner arc phimax -> 0
    for ph in phis[::-1][1:]:
        wp.LineTo(r_in * math.cos(ph), r_in * math.sin(ph))
    wp.Close()
    face = wp.Face()
    face.faces.name = "guide"
    rmid = 0.5 * (r_in + r_out0)
    for e in face.edges:
        c = e.center
        rc = math.hypot(c[0], c[1])
        ang = math.atan2(c[1], c[0])
        if ang < 0.05 * phimax:
            e.name = "inlet"
        elif ang > 0.95 * phimax:
            e.name = "outlet"
        elif rc < rmid:
            e.name = "inner"
        else:
            e.name = "outer"
    return OCCGeometry(face, dim=2)


def solve_image(r_in=0.5, r_out0=1.0, taper=0.0, phimax=0.5 * pi, Psi=1.0,
                order=3, maxh=0.03):
    """Linear flux solve A on the sector (inner wall A=0, outer A=Psi, inlet/
    outlet the radial profile), then sample the hodograph image (q=|B|,
    theta=arg B) over the interior.  Returns the image points + a measure of
    how theta-DEPENDENT the q-extent is (0 = rectangle, >0 = free boundary)."""
    with TaskManager():
        mesh = Mesh(_sector(r_in, r_out0, taper, phimax).GenerateMesh(maxh=maxh))
        mesh.Curve(order)
        fes = H1(mesh, order=order, dirichlet="inner|outer|inlet|outlet")
        u, v = fes.TnT()
        rr = sqrt(x * x + y * y)
        ramp = (rr - r_in) / (r_out0 - r_in)             # radial flux profile
        a = BilinearForm(InnerProduct(grad(u), grad(v)) * dx)
        a.Assemble()
        gf = GridFunction(fes)
        bcf = mesh.BoundaryCF({"inner": CF(0.0), "outer": CF(Psi),
                               "inlet": Psi * ramp, "outlet": Psi * ramp},
                              default=0.0)
        gf.Set(bcf, definedon=mesh.Boundaries("inner|outer|inlet|outlet"))
        r = gf.vec.CreateVector()
        r.data = -a.mat * gf.vec
        gf.vec.data += a.mat.Inverse(fes.FreeDofs(),
                                     inverse="sparsecholesky") * r
        B = grad(gf)                                     # B=(A_y,-A_x): use grad rot
        # B_x = dA/dy, B_y = -dA/dx
        Bx = grad(gf)[1]
        By = -grad(gf)[0]

        # sample the interior on a (position-angle, radial) grid -> hodograph
        qs, ths, posang = [], [], []
        for ph in np.linspace(0.18 * phimax, 0.82 * phimax, 40):   # avoid port fringing
            ro = r_out0 * (1.0 - taper * ph / phimax)
            for rr_s in np.linspace(r_in + 0.10 * (ro - r_in),
                                    ro - 0.10 * (ro - r_in), 16):
                px, py = rr_s * math.cos(ph), rr_s * math.sin(ph)
                bx = float(Bx(mesh(px, py)))
                by = float(By(mesh(px, py)))
                q = math.hypot(bx, by)
                qs.append(q)
                ths.append(math.atan2(by, bx))
                posang.append(ph)
    qs = np.array(qs); ths = np.array(ths); posang = np.array(posang)

    # theta-dependence of the q-extent: bin by field-angle theta, look at how the
    # q-range per bin varies across bins (rectangle -> all bins same q-range).
    tb = np.linspace(ths.min(), ths.max(), 9)
    qspan_lo, qspan_hi = [], []
    for i in range(len(tb) - 1):
        m = (ths >= tb[i]) & (ths < tb[i + 1])
        if m.sum() > 3:
            qspan_lo.append(np.percentile(qs[m], 5))
            qspan_hi.append(np.percentile(qs[m], 95))
    qspan_lo = np.array(qspan_lo); qspan_hi = np.array(qspan_hi)
    # how much the low/high q-edges drift across theta (relative) = free-boundary measure
    qmid = 0.5 * (qs.min() + qs.max())
    free_measure = float((np.ptp(qspan_lo) + np.ptp(qspan_hi)) / (2 * qmid))
    return {
        "taper": taper, "ne": int(mesh.ne),
        "q": qs, "theta": ths, "posang": posang,
        "q_range": (float(qs.min()), float(qs.max())),
        "theta_range_deg": (float(np.degrees(ths.min())),
                            float(np.degrees(ths.max()))),
        "free_measure": free_measure,   # ~0 rectangle, >0 theta-dependent (free bdry)
    }


def _plot(rc, rt):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.2), dpi=150)
    for ax, r, ttl in ((axes[0], rc, "constant width: RECTANGLE image"),
                       (axes[1], rt, "tapered: FREE BOUNDARY (theta-dependent)")):
        sc = ax.scatter(np.degrees(r["theta"]), r["q"], c=np.degrees(r["posang"]),
                        s=8, cmap="twilight")
        ax.set_xlabel(r"field direction  $\theta_B$  [deg]")
        ax.set_ylabel("hodograph coordinate  $q=|B|$")
        ax.set_title(f"{ttl}\nfree-boundary measure = {r['free_measure']:.2f}")
        fig.colorbar(sc, ax=ax, label="position angle [deg]")
    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    print("Frontier 2: the hodograph image of a turning flux guide "
          "(rectangle vs free boundary)\n")
    rc = solve_image(taper=0.0)
    rt = solve_image(taper=0.5)
    print(f"  constant-width bend: ne={rc['ne']}, q in [{rc['q_range'][0]:.3f}, "
          f"{rc['q_range'][1]:.3f}], theta_B in "
          f"[{rc['theta_range_deg'][0]:.0f}, {rc['theta_range_deg'][1]:.0f}] deg")
    print(f"    -> free-boundary measure = {rc['free_measure']:.3f}  (~0 => the "
          f"image is a RECTANGLE: q-extent theta-independent = self-linearising)")
    print(f"  tapered bend (outer wall spirals in 50%): ne={rt['ne']}, "
          f"q in [{rt['q_range'][0]:.3f}, {rt['q_range'][1]:.3f}]")
    print(f"    -> free-boundary measure = {rt['free_measure']:.3f}  (>0 => the "
          f"q-extent VARIES with theta = a genuine FREE BOUNDARY)")
    print("\n  => constant width = rectangle hodograph image (the 1-D self-")
    print("     linearising case); a tapering turn = theta-dependent image = a")
    print("     free boundary.  Recovering that image from the prescribed")
    print("     physical guide (the inverse hodograph solve) is the open frontier.")
    _plot(rc, rt)


if __name__ == "__main__":
    main()
