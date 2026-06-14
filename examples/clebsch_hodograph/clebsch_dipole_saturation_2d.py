r"""Saturation in the dipole design: the iron flux RETURN path is a saturable flux
guide (a Chaplygin guide), so the gap-field operating curve B_gap(NI) is a magnetic-
circuit 1-shot -- computed at LINEAR cost, validated against the nonlinear FEM.

The linear dipole workflow (`clebsch_dipole_design_workflow.py`) designed the pole
SURFACE as a scalar-potential level set assuming the iron is an ideal equipotential
(unsaturated).  At high excitation the iron SATURATES, and the gap field stops rising
linearly with the drive NI -- the engineering operating limit.  This file integrates
that into the design via the magnetic CIRCUIT:

    NI  =  H_gap * gap   +   SUM_segments  H_iron(seg) * L_seg
        =  (B_gap/mu0) gap  +  SUM_seg  nu(B_seg) B_seg L_seg ,

with the flux conserved around the ring (`B_seg = Psi/A_seg`, `Psi = B_gap A_gap`).
The iron sum IS the Chaplygin reluctance integral `INT nu(|B|)|B| dl` of
`chaplygin_design_sweep_2d.py` -- the iron return path is a saturable flux guide, and
a deliberately NECKED segment (a THROAT, cross-section `A_throat < A_gap`) is the
Chaplygin throat embedded in the magnet: it carries `B_throat = B_gap (A_gap/A_throat)`
and SATURATES FIRST, setting the knee of `B_gap(NI)`.

Given NI the circuit root-finds `B_gap` in microseconds (no mesh, no Picard); the
nonlinear FEM (a real 2-D window-frame electromagnet, A-formulation `nu(|B|)` Picard)
is the truth it is validated against.  Here the lumped circuit matches the FEM to
~6-11% (the lumped-magnetic-circuit error: gap fringing + corner flux crowding), so a
whole design space -- the operating curve B_gap(NI) over a family of throat widths --
is obtained at LINEAR cost, where the nonlinear FEM would need a Picard loop per point.

Design content: the throat cross-section sets the saturation knee
`B_gap,knee ~ B_k (A_throat/A_gap)`; a thinner throat saturates the magnet EARLIER (the
gap field softens at a lower drive).  This is the saturation lever the (linear) level-set
pole design cannot see -- the iron flux-path saturation, made a one-line circuit by the
Chaplygin hodograph view.

Honest scope: the lumped series circuit (gap + iron segments) is the standard
engineering approximation (~10% here vs the FEM); the iron-path reluctance is the
Chaplygin 1-shot whose slender-guide form is validated in `chaplygin_hodograph_2d.py`.
The 2-D FEM here is the end-to-end check that the lumped circuit captures the real
B_gap(NI).

run:  python clebsch_dipole_saturation_2d.py            # circuit operating curves (fast)
      python clebsch_dipole_saturation_2d.py --fem      # + nonlinear FEM validation (slow)
"""
import math
import os
import sys
import time

import numpy as np

MU0 = 4 * math.pi * 1e-7

# ---- 2-D window-frame electromagnet (meters); flux circulates, crosses a top-leg gap ----
AX, AY = 0.10, 0.08        # outer half-width / half-height
T = 0.025                  # leg thickness (= gap cross-section A_gap)
GW = 0.006                 # gap width (flux-direction extent in the top leg)
CW, CH = 0.05, 0.05        # coil bundle (threads the window)
BOX = 0.30                 # air-box half-size
MUR0, BK = 2000.0, 1.5     # Froehlich iron: mu_r0 at B=0, knee B_k


def _mu_r(B, mur0=MUR0, Bk=BK):
    """Froehlich saturation in |B|: mu_r0 at B=0 -> 1 as B -> infinity."""
    return 1.0 + (mur0 - 1.0) / (1.0 + (B / Bk) ** 2)


# --------------------------------------------------------------------------- #
# the magnetic-circuit 1-shot: B_gap(NI) at LINEAR cost (root-find, no mesh)
# --------------------------------------------------------------------------- #
def _NI_of_Bgap(B_gap, T_throat):
    """The drive NI required for a given gap field, by the series magnetic circuit
    (gap + iron segments incl. the saturable THROAT).  Flux Psi = B_gap*A_gap, each
    segment B_seg = Psi/A_seg, iron MMF = nu(B_seg) B_seg L_seg (the Chaplygin term)."""
    A_gap = T
    Psi = B_gap * A_gap
    # iron segments: (length, cross-section).  Throat = the thin bottom bar.
    segs = [(2 * AX, T),              # top leg (gapped)
            (2 * AX, T_throat),       # bottom leg = the THROAT (thin -> saturates first)
            (4 * AY, T)]              # two side legs
    mmf_iron = 0.0
    for L, A in segs:
        B = Psi / A
        mmf_iron += (B / (MU0 * _mu_r(B))) * L          # H_seg * L_seg
    mmf_gap = (B_gap / MU0) * GW
    return mmf_gap + mmf_iron


def circuit_oneshot(NI, T_throat, B_hi=20.0):
    """Invert the circuit: given the drive NI, the gap field B_gap (bisection -- the
    NI(B_gap) map is monotone).  Microseconds, no mesh, no Picard.  B_hi is just the
    bracket upper bound (well above any achievable B_gap -- not a physical cap)."""
    lo, hi = 1e-5, B_hi
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if _NI_of_Bgap(mid, T_throat) < NI:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def design_map(T_throats=(0.006, 0.010, 0.016), NI_list=None):
    """The operating curves B_gap(NI) over a family of throat widths -- the saturation
    DESIGN MAP, all circuit 1-shots (linear cost).  Returns the curves + the knee per
    throat (where the gap-field permeance dB_gap/dNI has fallen to half its small-NI
    value -- the saturation onset)."""
    if NI_list is None:
        NI_list = np.geomspace(200.0, 1.2e5, 40)
    NI_list = np.asarray(NI_list, float)
    t0 = time.perf_counter()
    curves = []
    for Tt in T_throats:
        Bg = np.array([circuit_oneshot(float(NI), Tt) for NI in NI_list])
        perm = np.gradient(Bg, NI_list)
        p0 = float(perm[0])
        knee = int(np.argmin(np.abs(perm - 0.5 * p0)))
        curves.append({"T_throat": float(Tt), "NI": NI_list.tolist(), "B_gap": Bg.tolist(),
                       "knee_NI": float(NI_list[knee]), "knee_Bgap": float(Bg[knee]),
                       "Bgap_knee_pred": float(BK * Tt / T)})
    return {"curves": curves, "n_points": int(len(T_throats) * len(NI_list)),
            "map_seconds": float(time.perf_counter() - t0), "NI": NI_list.tolist()}


# --------------------------------------------------------------------------- #
# the nonlinear 2-D FEM (the truth the circuit is validated against)
# --------------------------------------------------------------------------- #
def _build_mesh(T_throat, maxh=0.04):
    import ngsolve as ng
    from netgen.occ import WorkPlane, OCCGeometry, Glue

    def rect(x0, x1, y0, y1):
        return WorkPlane().MoveTo(x0, y0).Rectangle(x1 - x0, y1 - y0).Face()

    top = rect(-AX, AX, AY - T, AY)
    bot = rect(-AX, AX, -AY, -AY + T_throat)            # thin throat
    left = rect(-AX, -AX + T, -AY, AY)
    right = rect(AX - T, AX, -AY, AY)
    gap = rect(-GW / 2, GW / 2, AY - T, AY)
    iron = (top + bot + left + right) - gap
    iron.faces.name = "iron"; iron.faces.maxh = 0.008
    coil = rect(-CW / 2, CW / 2, -CH / 2, CH / 2)
    coil.faces.name = "coil"; coil.faces.maxh = 0.01
    air = rect(-BOX, BOX, -BOX, BOX) - iron - coil
    air.faces.name = "air"
    for e in air.edges:
        c = e.center
        if max(abs(c[0]), abs(c[1])) > 0.9 * BOX:
            e.name = "outer"
    with ng.TaskManager():
        return ng.Mesh(OCCGeometry(Glue([iron, coil, air]), dim=2).GenerateMesh(maxh=maxh))


def nonlinear_fem(NI, T_throat, mesh=None, order=2, niter=250, tol=1e-6, relax=0.4):
    """A real 2-D nonlinear window-frame electromagnet: -div(nu(|B|) grad A) = J, the
    iron Froehlich-saturable, A-formulation under-relaxed Picard.  Returns the gap field
    B_gap and the Picard iteration count (the per-point cost the circuit avoids)."""
    import ngsolve as ng
    from ngsolve import H1, GridFunction, grad, dx, CF, sqrt, BilinearForm, LinearForm

    if mesh is None:
        mesh = _build_mesh(T_throat)
    Jcf = mesh.MaterialCF({"coil": NI / (CW * CH)}, default=0.0)
    fes = H1(mesh, order=order, dirichlet="outer")
    u, v = fes.TnT()

    def lin(nucf):
        with ng.TaskManager():
            a = BilinearForm(nucf * grad(u) * grad(v) * dx); a.Assemble()
            f = LinearForm(Jcf * v * dx); f.Assemble()
            g = GridFunction(fes)
            g.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
        return g

    gf = lin(CF(1.0 / MU0)); prev = np.array(gf.vec); nit = 0
    for it in range(niter):
        B = grad(gf); Bmag = sqrt(B[0] * B[0] + B[1] * B[1] + 1e-30)
        nucf = mesh.MaterialCF({"iron": 1.0 / (MU0 * _mu_r(Bmag))}, default=1.0 / MU0)
        gnew = lin(nucf)
        gf.vec.data = (1.0 - relax) * gf.vec + relax * gnew.vec
        cur = np.array(gf.vec)
        d = np.linalg.norm(cur - prev) / (np.linalg.norm(cur) + 1e-30); prev = cur.copy()
        nit = it + 1
        if d < tol:
            break
    B = grad(gf)
    Bg = float(sqrt(B[0] * B[0] + B[1] * B[1])(mesh(0.0, AY - T / 2)))
    return {"B_gap": Bg, "iters": nit}


def validate_vs_fem(T_throat=0.006, NI_list=(2000.0, 8000.0, 30000.0, 80000.0)):
    """Compare the circuit B_gap(NI) to the nonlinear FEM at sample drives (one mesh,
    reused across NI).  Returns the agreement + the Picard iteration counts."""
    mesh = _build_mesh(T_throat)
    rows = []
    for NI in NI_list:
        fem = nonlinear_fem(float(NI), T_throat, mesh=mesh)
        cir = circuit_oneshot(float(NI), T_throat)
        rows.append({"NI": float(NI), "B_gap_fem": fem["B_gap"], "B_gap_circuit": cir,
                     "rel_err": abs(fem["B_gap"] - cir) / (abs(fem["B_gap"]) + 1e-30),
                     "iters": fem["iters"]})
    return {"T_throat": float(T_throat), "rows": rows,
            "max_rel_err": float(max(r["rel_err"] for r in rows)),
            "mean_iters": float(np.mean([r["iters"] for r in rows]))}


def run(with_fem=False):
    mp = design_map()
    out = {"design_map": mp}
    if with_fem:
        val = validate_vs_fem()
        out["validation"] = val
        out["cost"] = {
            "map_points": mp["n_points"], "map_ms": mp["map_seconds"] * 1e3,
            "one_shot_solves": mp["n_points"], "mean_picard_iters": val["mean_iters"],
            "fem_equiv_linear_solves": int(round(mp["n_points"] * val["mean_iters"])),
        }
    return out


def _plot(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(13.8, 4.2), dpi=150)
    mp = out["design_map"]
    colors = ["C0", "C1", "C2", "C3"]

    # Panel A: B_gap(NI) operating curves per throat width, with FEM validation points
    for i, c in enumerate(mp["curves"]):
        axA.plot(np.array(c["NI"]) * 1e-3, np.array(c["B_gap"]), "-", color=colors[i % 4],
                 label=f"throat {c['T_throat']*1e3:.0f} mm")
        axA.plot(c["knee_NI"] * 1e-3, c["knee_Bgap"], "o", color=colors[i % 4], ms=6)
    val = out.get("validation")
    if val is not None:
        ni = [r["NI"] * 1e-3 for r in val["rows"]]
        bf = [r["B_gap_fem"] for r in val["rows"]]
        axA.plot(ni, bf, "kx", ms=9, mew=2, label=f"nonlinear FEM (throat {val['T_throat']*1e3:.0f} mm)")
    axA.set_xlabel("drive  NI  [kA-turns]"); axA.set_ylabel("gap field  $B_{gap}$  [T]")
    axA.set_title("operating curve $B_{gap}$(NI), per throat width\n(each circuit point = ONE root-find; o = knee)")
    axA.legend(fontsize=8, loc="lower right")

    # Panel B: the saturation knob -- the knee drive vs throat width
    Tt = [c["T_throat"] * 1e3 for c in mp["curves"]]
    knee_NI = [c["knee_NI"] * 1e-3 for c in mp["curves"]]
    knee_B = [c["knee_Bgap"] for c in mp["curves"]]
    axB.plot(Tt, knee_NI, "C0o-", label="knee drive [kA-t]")
    axB2 = axB.twinx()
    axB2.plot(Tt, knee_B, "C3s--", label="knee $B_{gap}$ [T]")
    axB2.plot(Tt, [c["Bgap_knee_pred"] for c in mp["curves"]], "C3:", lw=1,
              label="$B_k\\,A_{throat}/A_{gap}$")
    axB.set_xlabel("iron throat width  [mm]")
    axB.set_ylabel("knee drive NI [kA-t]", color="C0")
    axB2.set_ylabel("knee $B_{gap}$ [T]", color="C3")
    axB.set_title("the saturation design knob:\na thinner throat saturates EARLIER")
    h1, l1 = axB.get_legend_handles_labels(); h2, l2 = axB2.get_legend_handles_labels()
    axB.legend(h1 + h2, l1 + l2, fontsize=7, loc="upper left")

    # Panel C: the linear-cost win
    cost = out.get("cost")
    if cost is not None:
        bars = axC.bar(["circuit\n(1 root-find each)", "nonlinear FEM\n(Picard each)"],
                       [cost["one_shot_solves"], cost["fem_equiv_linear_solves"]], color=["C0", "C3"])
        axC.set_yscale("log"); axC.set_ylabel("linear solves for the whole map")
        for b, v in zip(bars, [cost["one_shot_solves"], cost["fem_equiv_linear_solves"]]):
            axC.text(b.get_x() + b.get_width() / 2, v, f"{v}", ha="center", va="bottom", fontsize=9)
        axC.set_title(f"{cost['map_points']} operating points in {cost['map_ms']:.0f} ms\n"
                      f"(FEM: ~{cost['mean_picard_iters']:.0f} Picard iters EACH, "
                      f"agree {out['validation']['max_rel_err']*100:.0f}%)")
    else:
        axC.text(0.5, 0.5, "FEM validation + cost:\nrun with --fem", ha="center", va="center",
                 transform=axC.transAxes, fontsize=11)
        axC.set_xticks([]); axC.set_yticks([]); axC.set_title("the linear-cost win")

    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout(); fig.savefig(png, bbox_inches="tight"); plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    with_fem = "--fem" in sys.argv
    print("Saturation in the dipole: the iron flux path is a Chaplygin guide, "
          "B_gap(NI) at linear cost\n")
    out = run(with_fem=with_fem)
    mp = out["design_map"]
    print(f"  operating-curve design map: {mp['n_points']} points "
          f"({len(mp['curves'])} throat widths x {len(mp['NI'])} drives) in "
          f"{mp['map_seconds']*1e3:.1f} ms (circuit root-find each):")
    for c in mp["curves"]:
        print(f"    throat {c['T_throat']*1e3:4.0f} mm:  saturation knee at NI = "
              f"{c['knee_NI']*1e-3:6.1f} kA-t, B_gap = {c['knee_Bgap']:.2f} T "
              f"(predicted B_k*A_throat/A_gap = {c['Bgap_knee_pred']:.2f} T)")
    print(f"    -> a thinner iron throat saturates the magnet EARLIER (lower knee).")
    if with_fem:
        val, cost = out["validation"], out["cost"]
        print(f"  nonlinear FEM validation (throat {val['T_throat']*1e3:.0f} mm, A-formulation Picard):")
        for r in val["rows"]:
            print(f"    NI = {r['NI']*1e-3:5.1f} kA-t:  FEM B_gap = {r['B_gap_fem']:.3f} T  vs  "
                  f"circuit {r['B_gap_circuit']:.3f} T  (rel.err {r['rel_err']*100:.1f}%, {r['iters']} Picard iters)")
        print(f"    max rel.err = {val['max_rel_err']*100:.1f}% (the lumped-circuit error: "
              f"fringing + corner crowding)")
        print(f"  THE LINEAR-COST WIN: {cost['map_points']} operating points in {cost['map_ms']:.0f} ms "
              f"= {cost['one_shot_solves']} root-finds; the equivalent")
        print(f"    nonlinear FEM = {cost['map_points']} x ~{cost['mean_picard_iters']:.0f} Picard iters "
              f"= ~{cost['fem_equiv_linear_solves']} curved-mesh linear solves.")
    else:
        print("  nonlinear FEM validation + cost: run with  --fem  (slow)")
    print("\n  => the iron flux RETURN path is a saturable Chaplygin guide, so the dipole's")
    print("     gap-field operating curve B_gap(NI) -- including iron saturation -- is a")
    print("     magnetic-circuit 1-shot: a whole saturation design space at LINEAR cost.")
    _plot(out)


if __name__ == "__main__":
    main()
