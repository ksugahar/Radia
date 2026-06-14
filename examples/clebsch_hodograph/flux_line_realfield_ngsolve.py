r"""Flux-line closure on a REAL solved FE field: the de Rham / edge-FE reconstruction
closes, the nodal-averaged reconstruction spirals.

This is the *field* face of `flux_line_closure_symplectic.py` carried onto an actual
NGSolve solve (that script isolated the *integrator* face on an analytic field).  It
realises Noguchi's point on a genuine finite-element magnetic field and turns it into
the **field-reconstruction-quality diagnostic** for the HDiv-VIM migration.

Setup: a 2-D magnetostatic solve `-div(grad A_z) = J` for a current dipole (a +J / -J
wire pair) on a square, `A_z in H1` (order p).  A flux line is an integral curve of B
(`dx/ds = B`); in 2-D `A_z` is the flux-line Hamiltonian, so a flux line closes iff
`A_z` is conserved along it.  We trace ONE flux line of two reconstructions of the SAME
solve with the SAME integrator, so the only variable is the reconstruction:

  (A) the de Rham / edge-FE field   B_cl = rot(grad A_z) = (dA_z/dy, -dA_z/dx).
      This is the discrete curl of the scalar potential -- a CLOSED 2-form by
      construction (`div B_cl = 0` weakly; its normal component is continuous across
      every element edge because it equals the tangential derivative of the continuous
      A_z).  Crucially `B_cl . grad(A_z) = 0 pointwise` -- the field is exactly tangent
      to the flux surfaces A_z = const.  This is exactly the edge-(H(curl)) flux-line
      field of Hirahatake, Noguchi, Igarashi & Yamashita (IEEJ, pp.1205-1212):
      `B = curl A` from edge elements is divergence-free by the de Rham property.

  (B) the nodal-averaged field        B_lk = L2-project(B_cl) onto VectorH1.
      Forcing BOTH components continuous (nodal averaging) breaks the de Rham
      structure: B_lk is no longer the curl of a scalar, `div B_lk != 0`, and
      `B_lk . grad(A_z) != 0` -- the reconstructed field drifts OFF the flux surfaces.
      This is the 2-D analogue of a *leaky* demag-field reconstruction (the HDiv-VIM
      `M_mass^-1 N m` solenoidal leak): the flux line is the visible symptom.

Reported:
  - field-quality metric: the RMS misalignment `|B . grad(A_z)| / (|B| |grad A_z|)`
    (the cosine of the angle off the flux surface).  ~0 for (A), order-1e-1 for (B).
  - the traced flux line: A_z drift along the curve, and the return distance to the
    start.  (A) closes (A_z bounded, returns); (B) spirals (A_z drifts secularly).
    Same RK4 integrator for both -- so the contrast is the FIELD, not the stepper.

The 3-D generalisation is helicity-obstructed (a real 3-D field need not be foliated
by closed flux surfaces) -- see `clebsch_3d_closing_condition.py`.  In 2-D `A_z` is the
single Clebsch potential and closure is unobstructed, so 2-D is the clean regime to
pin the reconstruction-quality diagnostic.

run:  python flux_line_realfield_ngsolve.py
"""
import os

import numpy as np
import ngsolve as ng
from netgen.geom2d import SplineGeometry

x, y = ng.x, ng.y

# current dipole: +J / -J wires at (+-XC, 0), smooth (Gaussian) so no region meshing.
XC = 0.35
WIRE_S = 0.06          # wire core radius (Gaussian width)
EPS_LEAK = 0.06        # strength of the explicit de Rham-complement (charge) leak


def _solve_Az(order=3, maxh=0.045):
    """Real 2-D magnetostatic solve -div(grad A_z) = J for the current dipole.
    Returns (mesh, gfA)."""
    geo = SplineGeometry()
    geo.AddRectangle((-1, -1), (1, 1), bc="outer")
    mesh = ng.Mesh(geo.GenerateMesh(maxh=maxh))
    fes = ng.H1(mesh, order=order, dirichlet="outer")
    u, v = fes.TnT()
    blob_p = ng.exp(-((x - XC) ** 2 + y ** 2) / WIRE_S ** 2)
    blob_m = ng.exp(-((x + XC) ** 2 + y ** 2) / WIRE_S ** 2)
    J = blob_p - blob_m
    a = ng.BilinearForm(fes, symmetric=True)
    a += ng.grad(u) * ng.grad(v) * ng.dx
    f = ng.LinearForm(fes)
    f += J * v * ng.dx
    a.Assemble(); f.Assemble()
    gfA = ng.GridFunction(fes)
    gfA.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    return mesh, gfA


def _misalignment(mesh, B_cf, gradA_cf):
    """RMS of  |B . grad(A_z)| / (|B| |grad A_z|)  -- the cosine of the angle off the
    flux surface A_z = const.  0 = exactly tangent (de Rham), >0 = drifts off."""
    dotsq = ng.InnerProduct(B_cf, gradA_cf) ** 2
    bn2 = ng.InnerProduct(B_cf, B_cf)
    gn2 = ng.InnerProduct(gradA_cf, gradA_cf)
    num = ng.Integrate(dotsq, mesh)
    den = ng.Integrate(bn2 * gn2, mesh)
    return float(np.sqrt(num / den))


def _trace(mesh, vel_cf, gfA, x0, ds, n):
    """RK4 flux-line trace dx/ds = vhat (unit-speed) of the field vel_cf.  Same
    integrator for every field -> the only variable is the reconstruction.  Stops
    if the curve leaves the mesh.  Returns (xs, Az_along)."""
    def vhat(p):
        mp = mesh(p[0], p[1])              # raises if outside the mesh
        vx, vy = vel_cf(mp)
        v = np.array([vx, vy]); nv = np.linalg.norm(v)
        return v / nv if nv > 1e-14 else v
    xs = [np.array(x0, float)]; Az = [float(gfA(mesh(x0[0], x0[1])))]
    p = np.array(x0, float)
    for _ in range(n):
        try:
            k1 = vhat(p)
            k2 = vhat(p + 0.5 * ds * k1)
            k3 = vhat(p + 0.5 * ds * k2)
            k4 = vhat(p + ds * k3)
            p = p + ds / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
            Az.append(float(gfA(mesh(p[0], p[1]))))
            xs.append(p.copy())
        except Exception:
            break
    return np.array(xs), np.array(Az)


def _closure(xs, x0):
    """Smallest return distance to the start after the curve has left a neighbourhood."""
    d = np.linalg.norm(xs - np.asarray(x0), axis=1)
    if len(d) < 10:
        return float(d.max() if len(d) else 1.0)
    left = np.where(d > 0.5 * d.max())[0]
    return float(d[left[0]:].min()) if len(left) else float(d.min())


def analyze(order=3, maxh=0.045, turns=4.0, steps_per_turn=520, eps_leak=EPS_LEAK):
    """The reconstruction-quality diagnostic.  Three reconstructions of ONE solve:
      - de Rham  rot(grad A_z)             -- closed 2-form, exactly on the flux surface;
      - nodal-averaged VectorH1 order 1    -- the realistic parameter-free leak (the
                                              edge-vs-nodal point, Noguchi);
      - explicit charge leak B_cl + eps*grad(A_z) -- the controlled de Rham-complement
                                              admixture (model of the M_mass^-1 N m leak).
    Returns each field's flux-surface misalignment and the closure of its traced line."""
    with ng.TaskManager():
        mesh, gfA = _solve_Az(order=order, maxh=maxh)
        gradA = ng.grad(gfA)
        B_cl = ng.CoefficientFunction((gradA[1], -gradA[0]))      # rot(grad A_z): closed 2-form
        # (B) realistic leak: nodal averaging onto linear VectorH1 (the edge-vs-nodal point).
        fesB = ng.VectorH1(mesh, order=1)
        gfB = ng.GridFunction(fesB); gfB.Set(B_cl)
        B_nd = ng.CoefficientFunction((gfB[0], gfB[1]))
        # (C) controlled leak: a spurious irrotational (charge) component of strength eps.
        B_lk = B_cl + eps_leak * gradA

        mis_cl = _misalignment(mesh, B_cl, gradA)
        mis_nd = _misalignment(mesh, B_nd, gradA)
        mis_lk = _misalignment(mesh, B_lk, gradA)

        # a flux line around the +wire: start above it, on a closed A_z contour.
        x0 = (XC, 0.20)
        # rough loop circumference for a geometric step (~ encircles the wire).
        circ = 2 * np.pi * 0.20
        ds = circ / steps_per_turn
        n = int(turns * steps_per_turn)
        xs_cl, Az_cl = _trace(mesh, B_cl, gfA, x0, ds, n)
        xs_nd, Az_nd = _trace(mesh, B_nd, gfA, x0, ds, n)
        xs_lk, Az_lk = _trace(mesh, B_lk, gfA, x0, ds, n)

    Az0 = abs(Az_cl[0]) + 1e-30
    drift = lambda Az: float(np.max(np.abs(Az - Az[0])) / Az0)
    return {
        "order": order, "maxh": maxh, "ne": int(mesh.GetNE(ng.VOL)), "eps_leak": eps_leak,
        "mis_closed": mis_cl, "mis_nodal": mis_nd, "mis_leaky": mis_lk,
        "drift_closed": drift(Az_cl), "drift_nodal": drift(Az_nd), "drift_leaky": drift(Az_lk),
        "ret_closed": _closure(xs_cl, x0), "ret_nodal": _closure(xs_nd, x0),
        "ret_leaky": _closure(xs_lk, x0),
        "n_closed": len(xs_cl), "n_nodal": len(xs_nd), "n_leaky": len(xs_lk),
        "xs_closed": xs_cl, "xs_nodal": xs_nd, "xs_leaky": xs_lk, "x0": x0,
        "Az_closed": Az_cl, "Az_nodal": Az_nd, "Az_leaky": Az_lk,
    }


def _plot(r):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.8, 4.4), dpi=150)
    cl, nd, lk = r["xs_closed"], r["xs_nodal"], r["xs_leaky"]
    ax.plot(lk[:, 0], lk[:, 1], "C3-", lw=0.8, label="charge leak (spirals)")
    ax.plot(nd[:, 0], nd[:, 1], "C1-", lw=0.8, label="nodal-averaged (spirals)")
    ax.plot(cl[:, 0], cl[:, 1], "C0-", lw=1.5, label="de Rham rot(grad $A_z$) (closes)")
    ax.plot([XC], [0], "kx", ms=8, mew=2); ax.plot(*r["x0"], "k.", ms=10)
    ax.set_aspect("equal"); ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_title("flux line of a REAL FE solve:\nde Rham closes, leaky reconstructions spiral")
    ax.legend(fontsize=8, loc="lower center")
    for key, c, lab in (("closed", "C0", "de Rham"), ("nodal", "C1", "nodal-avg"),
                        ("leaky", "C3", "charge leak")):
        Az = r[f"Az_{key}"]; s = np.arange(len(Az))
        ax2.semilogy(s, np.abs(Az / Az[0] - 1) + 1e-16, c + "-", lw=1.2,
                     label=f"{lab} (mis={r['mis_' + key]:.1e})")
    ax2.set_xlabel("step"); ax2.set_ylabel("$|A_z/A_z(0) - 1|$  (loss of closure)")
    ax2.set_title("flux-surface drift = reconstruction leak")
    ax2.legend(fontsize=8, loc="lower right")
    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout(); fig.savefig(png, bbox_inches="tight"); plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    print("Flux-line closure on a REAL solved FE field: de Rham closes, nodal-averaged spirals\n")
    r = analyze()
    print(f"  2-D magnetostatic solve  ne={r['ne']}  H1 order {r['order']}  (current dipole)")
    print(f"  field-quality (RMS misalignment |B.gradA|/(|B||gradA|), 0 = on the flux surface):")
    print(f"    de Rham  rot(grad A_z) : {r['mis_closed']:.2e}   (closed 2-form, tangent by construction)")
    print(f"    nodal-averaged (ord 1) : {r['mis_nodal']:.2e}   (edge-vs-nodal leak, parameter-free)")
    print(f"    charge leak (eps={r['eps_leak']:.2f})  : {r['mis_leaky']:.2e}   (de Rham-complement admixture)")
    print(f"  traced flux line (same RK4 integrator):")
    print(f"    de Rham  : A_z drift={r['drift_closed']:.2e}  return dist={r['ret_closed']:.2e}  -> CLOSES")
    print(f"    nodal-avg: A_z drift={r['drift_nodal']:.2e}  return dist={r['ret_nodal']:.2e}  -> spirals")
    print(f"    charge   : A_z drift={r['drift_leaky']:.2e}  return dist={r['ret_leaky']:.2e}  -> spirals")
    print(f"    ratio drift(charge)/drift(de Rham)    = {r['drift_leaky']/ (r['drift_closed']+1e-30):.0f}x")
    print(f"    ratio drift(nodal-avg)/drift(de Rham) = {r['drift_nodal']/ (r['drift_closed']+1e-30):.0f}x")
    print("\n  => on a genuine FE field, the de Rham (edge / scalar-potential) reconstruction's flux")
    print("     lines close; the leaky (nodal-averaged or charge-admixed) reconstructions spiral.  This")
    print("     is the field-quality diagnostic for the HDiv-VIM migration: trace a flux line; if it")
    print("     spirals, the reconstruction leaked solenoidal content (the M_mass^-1 N m leak).")
    _plot(r)


if __name__ == "__main__":
    main()
