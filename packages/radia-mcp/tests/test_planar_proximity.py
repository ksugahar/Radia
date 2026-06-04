"""Multi-conductor circuit (skin + PROXIMITY effect) -- validates
radia_ngsolve.solve.solve_planar_eddy_multi against the exact round-wire
Kelvin (ber/bei) AC resistance.

Three checks against the same q = a*sqrt(w*mu0*sigma) reference used by
test_planar_eddy:

  1. ONE conductor (series or parallel path) reproduces the single-wire
     Kelvin Rac/Rdc  -- the new multi-conductor assembly reduces correctly.
  2. TWO well-separated wires in SERIES (same current I in each) dissipate
     2x the isolated-wire loss: Rac_circuit = 2*Rac_single  (negligible
     coupling at large spacing).
  3. Bringing the two wires CLOSE raises the circuit resistance above
     2*Rac_single -- the proximity effect (each wire sits in the other's AC
     field). We require a clear, monotone increase.
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from scipy.special import kelvin
from ngsolve import Mesh, CoefficientFunction, TaskManager
from netgen.occ import OCCGeometry, MoveTo, WorkPlane, Glue, X, Y
from radia_mcp.radia_ngsolve.solve import solve_planar_eddy, solve_planar_eddy_multi
from radia_mcp.radia_ngsolve.force import ohmic_loss_2d

MU0 = 4e-7 * math.pi
A, SIGMA, I_TOT, Q = 1e-3, 5.8e7, 1.0, 4.0
OMEGA = Q * Q / (A * A * MU0 * SIGMA)
RDC = 1.0 / (SIGMA * math.pi * A * A)
DELTA = math.sqrt(2.0 / (OMEGA * MU0 * SIGMA))


def rac_over_rdc_exact(q):
    be, ke, bep, kep = kelvin(q)
    ber, bei, berp, beip = be.real, be.imag, bep.real, bep.imag
    return (q / 2.0) * (ber * beip - bei * berp) / (berp ** 2 + beip ** 2)


def build_wire_mesh(centers, pad=12 * A):
    """Air box containing round wires named wire0, wire1, ... at ``centers``."""
    xs = [c[0] for c in centers]
    ys = [c[1] for c in centers]
    x0, x1 = min(xs) - pad, max(xs) + pad
    y0, y1 = min(ys) - pad, max(ys) + pad
    wires = []
    for k, (cx, cy) in enumerate(centers):
        w = WorkPlane().Circle(cx, cy, A).Face()
        w.faces.name = f"wire{k}"
        w.maxh = DELTA / 4.0
        wires.append(w)
    box = MoveTo(x0, y0).Rectangle(x1 - x0, y1 - y0).Face()
    air = box
    for w in wires:
        air = air - w
    air.faces.name = "air"
    for sel in (air.edges.Max(X), air.edges.Min(X), air.edges.Max(Y), air.edges.Min(Y)):
        sel.name = "outer"
    return Mesh(OCCGeometry(Glue([air] + wires), dim=2).GenerateMesh(maxh=pad / 6))


def circuit_rac(mesh, conductors, connection):
    sigma = mesh.MaterialCF({c: SIGMA for c in conductors}, default=0.0)
    nu = CoefficientFunction(1.0 / MU0)
    gfu = solve_planar_eddy_multi(mesh, nu, sigma, OMEGA, conductors,
                                  connection=connection, total_current=I_TOT,
                                  order=3)
    comps = gfu.components
    Az = comps[0]
    P = 0.0
    for k, c in enumerate(conductors):
        Vc = comps[1] if connection == "parallel" else comps[1 + k]
        P += ohmic_loss_2d(-1j * OMEGA * Az + Vc, mesh, sigma, region=c)
    return 2.0 * P / (I_TOT ** 2)


def main():
    rac_single_exact = RDC * rac_over_rdc_exact(Q)
    print(f"  q={Q}  delta/a={DELTA/A:.3f}  Rac_single(Kelvin)={rac_single_exact:.6e} ohm/m")

    with TaskManager():
        # 1. single conductor via the multi path == Kelvin single wire
        m1 = build_wire_mesh([(0.0, 0.0)])
        rac1 = circuit_rac(m1, ["wire0"], "series")
        err1 = (rac1 - rac_single_exact) / rac_single_exact
        print(f"  [1] single (series path)  Rac={rac1:.6e}  err={100*err1:+.2f}%")
        assert abs(err1) < 0.02, f"single-conductor multi off {100*err1:.2f}%"

        # 2. two FAR wires in series -> 2 x single
        D_far = 16 * A
        mfar = build_wire_mesh([(-D_far / 2, 0.0), (D_far / 2, 0.0)])
        rac_far = circuit_rac(mfar, ["wire0", "wire1"], "series")
        err2 = (rac_far - 2 * rac_single_exact) / (2 * rac_single_exact)
        print(f"  [2] two far (D={D_far/A:.0f}a) series Rac={rac_far:.6e}  "
              f"vs 2*single err={100*err2:+.2f}%")
        assert abs(err2) < 0.02, f"two far wires off from 2x single {100*err2:.2f}%"

        # 3. two CLOSE wires in series -> proximity raises Rac
        D_near = 2.4 * A
        mnear = build_wire_mesh([(-D_near / 2, 0.0), (D_near / 2, 0.0)])
        rac_near = circuit_rac(mnear, ["wire0", "wire1"], "series")
        prox = rac_near / rac_far
        print(f"  [3] two near (D={D_near/A:.1f}a) series Rac={rac_near:.6e}  "
              f"proximity factor Rac_near/Rac_far = {prox:.4f}")
        assert prox > 1.02, f"expected proximity increase, got factor {prox:.4f}"

    print(f"\n[OK] multi-conductor circuit validated: single->Kelvin ({100*abs(err1):.2f}%), "
          f"2 far->2x single ({100*abs(err2):.2f}%), proximity +{100*(prox-1):.1f}%.")


if __name__ == "__main__":
    main()
