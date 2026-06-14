r"""Does raising the element ORDER make flux lines close?  No -- the de Rham
REPRESENTATION does.  An order sweep that answers two old questions.

Two questions motivated this script:

  Q1 (extending Noguchi).  Noguchi's edge-FE flux-line method works because
     `B = curl A` from edge (`H(curl)`) elements is exactly divergence-free (a
     closed 2-form).  Can that be *extended* -- "if the field is de Rham, a
     symplectic (structure-preserving) tracker can be built on it"?

  Q2 (Kameari's remark).  "second-order elements -- the flux lines don't close."
     Is that because they are *not de Rham*?  Would a de Rham 2nd-order element
     close?

The honest, measured answer (this script reproduces the numbers):

  * `B = curl A` is divergence-free for ANY conforming `A` -- edge (`H(curl)`) OR
    nodal Lagrange (`[H1]^3`) -- at EVERY order (machine zero).  So "edge vs nodal
    A" is *not* the flux-line-closure discriminator (it is the spurious-eigenmode /
    interface-continuity discriminator -- the classical edge-element motivation,
    Bossavit / Nedelec / Kameari).  A first guess that "nodal curl leaks through
    normal jumps" is WRONG: `(curl A).n` depends only on the *tangential* trace of
    `A`, which is continuous for nodal `A` too.

  * What breaks closure is leaving the de Rham representation -- **nodally SMOOTHING
    B** (the legacy "evaluate the field at the nodes and interpolate" post-processing,
    an L2 projection onto a continuous nodal space).  The smoothed field is no longer
    an exact curl: it acquires a spurious divergence.  That leak DECREASES with order
    but is NEVER zero -- so the flux lines do *not* close *even at 2nd order* (Q2:
    yes, that matches Kameari; raising the order does not fix it).  The de Rham
    `B = curl A`, kept in its native representation (unsmoothed), is EXACTLY
    divergence-free at every order -- a de Rham 2nd-order field *closes* (Q2: yes).

  * de Rham is therefore the **enabling precondition** for a structure-preserving
    tracker (Q1): only an exactly closed 2-form has the conserved structure a
    symplectic (2-D, Hamiltonian `A_z`) / volume-preserving (3-D, `div B = 0`)
    integrator preserves.  A smoothed (leaky) field has no closed 2-form for the
    integrator to keep, so it spirals regardless of the stepper.  Noguchi supplies
    the field (de Rham, closed 2-form); accelerator/Feng-Kang structure-preserving
    integration supplies the tracker -- their union is the natural extension.

What is reported:

  Part A (3-D, the de Rham exactness): the weak interior divergence of `B` for
    de Rham `curl(HCurl_p)`, nodal-A `curl([H1]^3_p)`, and nodally-smoothed `B`,
    swept over order p.  The two curls are machine-zero at every order; the
    smoothed field leaks (decreasing with p, never zero).

  Part B (2-D, the closure tie): a real magnetostatic solve `-div(grad A_z) = J`;
    a flux line traced (same RK4) for the de Rham `rot(grad A_z)` vs the
    nodally-smoothed reconstruction, swept over order p.  de Rham is exactly tangent
    to the flux surfaces (misalignment 0.0) and closes at EVERY order; the smoothed
    reconstruction's misalignment / A_z drift fall with p but stay far above de Rham
    -- it does not close even at 2nd order.

3-D caveat: de Rham gives `div B = 0`, which in 2-D forces closure but in 3-D only
gives a *volume-preserving* flow; genuine 3-D closure additionally needs the
helicity to vanish (a global Clebsch pair) -- see `clebsch_3d_closing_condition.py`.

Refs: P. Robert, IEEE Trans. Magn. 27(5), 1991 (Clebsch / de Rham, B as a 2-form);
Hirahatake-Noguchi-Igarashi-Yamashita, IEEJ pp.1205-1212 (edge-FE flux lines);
A. Bossavit / J.C. Nedelec (edge elements, the discrete de Rham sequence);
Moffatt, J. Fluid Mech. 35, 1969 (helicity).  Full list in
docs/clebsch_hodograph/HDIV_VIM_CLEBSCH_BRIDGE.md.

run:  python derham_closure_order_sweep.py
"""
import os

import numpy as np
import ngsolve as ng
from netgen.csg import unit_cube

import flux_line_realfield_ngsolve as fl   # reuse the 2-D solve + trace helpers

x, y, z = ng.x, ng.y, ng.z

# a smooth, genuinely 3-D vector potential (so curl A is a non-trivial field).
A3D = ng.CoefficientFunction((ng.sin(y) * z, ng.cos(z) * x, ng.sin(x) * y))


def _weak_interior_div(B, mesh, order):
    """RMS interior weak divergence  ||q -> INT B.grad(q)||  over the H1_0 test space
    (q = 0 on the outer boundary, so only the INTERIOR spurious charge is measured),
    normalized by ||B||.  0 = exact closed 2-form (no spurious magnetic charge)."""
    Q = ng.H1(mesh, order=order, dirichlet=".*")
    q = Q.TestFunction(); u = Q.TrialFunction()
    Lf = ng.LinearForm(Q); Lf += ng.InnerProduct(B, ng.grad(q)) * ng.dx; Lf.Assemble()
    M = ng.BilinearForm(Q); M += (ng.grad(u) * ng.grad(q) + u * q) * ng.dx; M.Assemble()
    c = ng.GridFunction(Q)
    c.vec.data = M.mat.Inverse(Q.FreeDofs(), inverse="sparsecholesky") * Lf.vec
    riesz = float(np.sqrt(abs(ng.InnerProduct(c.vec, Lf.vec))))
    bnorm = float(np.sqrt(ng.Integrate(ng.InnerProduct(B, B), mesh)))
    return riesz / (bnorm + 1e-30)


def divergence_sweep(orders=(1, 2, 3), maxh=0.25):
    """Part A: 3-D weak interior divergence of three reconstructions, swept over order.
    de Rham curl(HCurl_p) and nodal-A curl([H1]^3_p) are machine-zero at every order;
    the nodally-smoothed B leaks (decreasing with p, never zero)."""
    rows = []
    with ng.TaskManager():
        mesh = ng.Mesh(unit_cube.GenerateMesh(maxh=maxh))
        for p in orders:
            Hc = ng.HCurl(mesh, order=p); a = ng.GridFunction(Hc); a.Set(A3D)
            B_dr = ng.curl(a)                                          # de Rham exact curl
            Hn = ng.VectorH1(mesh, order=p); an = ng.GridFunction(Hn); an.Set(A3D)
            G = ng.grad(an)
            B_nd = ng.CoefficientFunction((G[2, 1] - G[1, 2], G[0, 2] - G[2, 0],
                                           G[1, 0] - G[0, 1]))         # nodal A's curl
            gfB = ng.GridFunction(Hn); gfB.Set(B_dr)                   # SMOOTH B onto nodes
            B_av = ng.CoefficientFunction((gfB[0], gfB[1], gfB[2]))
            rows.append((p,
                         _weak_interior_div(B_dr, mesh, p),
                         _weak_interior_div(B_nd, mesh, p),
                         _weak_interior_div(B_av, mesh, p)))
    return {"ne": int(mesh.GetNE(ng.VOL)), "rows": rows}


def closure_sweep(orders=(1, 2, 3), maxh=0.06, turns=3.0, steps_per_turn=360):
    """Part B: 2-D flux-line closure of the de Rham rot(grad A_z) vs the nodally-smoothed
    reconstruction, swept over order.  de Rham is exactly tangent (misalignment 0.0) and
    closes at every order; the smoothed reconstruction leaks (mis / drift fall with p but
    stay far above de Rham) -- it does not close even at 2nd order."""
    x0 = (fl.XC, 0.20); circ = 2 * np.pi * 0.20
    rows = []
    for p in orders:
        with ng.TaskManager():
            mesh, gfA = fl._solve_Az(order=p, maxh=maxh)
            gradA = ng.grad(gfA)
            B_cl = ng.CoefficientFunction((gradA[1], -gradA[0]))      # de Rham rot(grad A_z)
            fesB = ng.VectorH1(mesh, order=p); gfB = ng.GridFunction(fesB); gfB.Set(B_cl)
            B_av = ng.CoefficientFunction((gfB[0], gfB[1]))           # nodally-smoothed
            mis_cl = fl._misalignment(mesh, B_cl, gradA)
            mis_av = fl._misalignment(mesh, B_av, gradA)
            ds = circ / steps_per_turn; n = int(turns * steps_per_turn)
            xs_cl, Az_cl = fl._trace(mesh, B_cl, gfA, x0, ds, n)
            xs_av, Az_av = fl._trace(mesh, B_av, gfA, x0, ds, n)
        drift = lambda Az: float(np.max(np.abs(Az - Az[0])) / (abs(Az[0]) + 1e-30))
        rows.append({
            "order": p, "mis_derham": mis_cl, "mis_smoothed": mis_av,
            "drift_derham": drift(Az_cl), "drift_smoothed": drift(Az_av),
            "ret_derham": fl._closure(xs_cl, x0), "ret_smoothed": fl._closure(xs_av, x0),
            "xs_derham": xs_cl, "xs_smoothed": xs_av,
            "Az_derham": Az_cl, "Az_smoothed": Az_av,
        })
    return {"x0": x0, "rows": rows}


def analyze(orders=(1, 2, 3)):
    return {"divergence": divergence_sweep(orders=orders),
            "closure": closure_sweep(orders=orders), "orders": list(orders)}


def _plot(r):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(13.5, 4.2), dpi=150)
    # Part A: 3-D weak divergence vs order
    dv = r["divergence"]["rows"]
    ps = [row[0] for row in dv]
    axA.semilogy(ps, [max(row[1], 1e-17) for row in dv], "C0o-", label="de Rham curl($H$(curl))")
    axA.semilogy(ps, [max(row[2], 1e-17) for row in dv], "C2s--", label="nodal-A curl($[H^1]^3$)")
    axA.semilogy(ps, [row[3] for row in dv], "C3^-", label="nodally-smoothed $B$")
    axA.set_xlabel("element order $p$"); axA.set_ylabel("weak interior divergence / $\\Vert B\\Vert$")
    axA.set_title("3-D: $B=$curl$\\,A$ is div-free at every order;\nsmoothing is the leak")
    axA.set_xticks(ps); axA.legend(fontsize=8, loc="center right")
    # Part B: 2-D flux-surface misalignment vs order
    cl = r["closure"]["rows"]
    ps2 = [row["order"] for row in cl]
    axB.semilogy(ps2, [max(row["mis_derham"], 1e-17) for row in cl], "C0o-", label="de Rham (closes)")
    axB.semilogy(ps2, [row["mis_smoothed"] for row in cl], "C3^-", label="smoothed (spirals)")
    axB.set_xlabel("element order $p$")
    axB.set_ylabel("flux-surface misalignment $|B\\cdot\\nabla A_z|/(|B||\\nabla A_z|)$")
    axB.set_title("2-D: raising the order does NOT close it;\nthe de Rham representation does")
    axB.set_xticks(ps2); axB.legend(fontsize=8, loc="center right")
    # Part C: the tracking-closure signal -- A_z (which MUST be conserved for the flux
    # line to close) along the traced line.  de Rham keeps it flat at the integrator
    # floor at every order; the smoothed reconstruction lets it drift (falling with
    # order but never to the floor) = the line does not close, even at 2nd order.
    row2 = next(rw for rw in cl if rw["order"] == 2)
    s = np.arange(len(row2["Az_derham"]))
    axC.semilogy(s, np.abs(row2["Az_derham"] / row2["Az_derham"][0] - 1) + 1e-16, "C0-",
                 lw=1.6, label="de Rham $p=2$ (closes)")
    for row, c in zip(cl, ("C1", "C3", "C4")):
        sm = row["Az_smoothed"]; ss = np.arange(len(sm))
        axC.semilogy(ss, np.abs(sm / sm[0] - 1) + 1e-16, c + "-", lw=1.0,
                     label=f"smoothed $p={row['order']}$")
    axC.set_xlabel("step along the flux line")
    axC.set_ylabel("$|A_z/A_z(0) - 1|$  (loss of closure)")
    axC.set_title("the tracking signal: smoothed loses $A_z$\n(does not close) even at $p=2$")
    axC.legend(fontsize=8, loc="lower right")
    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout(); fig.savefig(png, bbox_inches="tight"); plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    print("Does raising the element ORDER make flux lines close?  No -- the de Rham "
          "REPRESENTATION does.\n")
    r = analyze()
    print(f"  Part A (3-D, ne={r['divergence']['ne']}): weak interior divergence of B / ||B||")
    print(f"    order |  de Rham curl(HCurl) |  nodal-A curl([H1]^3) |  nodally-SMOOTHED B")
    for p, drv, ndv, avv in r["divergence"]["rows"]:
        print(f"      {p}   |     {drv:.2e}        |      {ndv:.2e}        |     {avv:.2e}")
    print("    -> both curls are machine-zero at EVERY order (div curl = 0 is a de Rham")
    print("       property, not an order property); the SMOOTHING leaks, decreasing with")
    print("       order but never zero -> 'even 2nd order does not close'.")
    print(f"\n  Part B (2-D flux-line closure): de Rham rot(grad A_z) vs nodally-smoothed")
    print(f"    order | de Rham mis / A_z drift / return | smoothed mis / A_z drift / return")
    for row in r["closure"]["rows"]:
        print(f"      {row['order']}   |   {row['mis_derham']:.1e} / {row['drift_derham']:.2e} / "
              f"{row['ret_derham']:.2e}  |  {row['mis_smoothed']:.1e} / {row['drift_smoothed']:.2e} / "
              f"{row['ret_smoothed']:.2e}")
    print("    -> de Rham misalignment is EXACTLY 0 at every order (closed 2-form by")
    print("       construction) -> closes; the smoothed misalignment FALLS with order but")
    print("       stays far above de Rham -> spirals even at 2nd order.")
    print("\n  => closure is governed by the REPRESENTATION (de Rham vs nodal-smoothed), NOT")
    print("     the order.  de Rham is the closed-2-form precondition that makes a symplectic")
    print("     / volume-preserving flux-line tracker meaningful (the Noguchi extension).")
    _plot(r)


if __name__ == "__main__":
    main()
