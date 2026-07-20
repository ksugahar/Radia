"""Litz AC-loss demonstration -- tie the build123d Litz archetype to round-wire eddy-current physics.

WHY Litz wire exists: subdividing a solid conductor into ``n`` insulated, transposed strands of radius
``R0/sqrt(n)`` (same total copper) cuts the high-frequency eddy loss.  Two mechanisms, two scalings:

  * **skin** (self-field) excess  ``F(q)-1`` with ``q = r*sqrt(w*mu*sigma)`` propto strand radius, so at
    low ``q`` the bundle/solid excess ratio -> ``1/n^2``;
  * **proximity** (external-field) excess: per strand propto ``r^4`` (a cylinder's induced eddy dipole in
    a uniform AC field), times ``n`` strands -> ``1/n`` of the solid -- usually the DOMINANT term, and the
    one helical transposition (the twist :func:`litz_wire` / :func:`hierarchical_litz` build) averages out.

The single-strand skin ratio is anchored to the FE solver ``solve_planar_eddy`` (it reproduces the exact
Kelvin ber/bei AC resistance to < 1 %), and the proximity penalty that transposition cancels is shown with
``solve_planar_eddy_multi`` (two close strands: ``Rac_near/Rac_far > 1``).  The strand count / radius / fill
come straight from :mod:`radia_mcp.build123d.archetypes`.
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
from radia_mcp.build123d.archetypes import litz_wire, litz_packing_radius, litz_fill_factor

MU0 = 4e-7 * math.pi
A_W, SIGMA, I_TOT, Q = 1e-3, 5.8e7, 1.0, 4.0          # 1 mm copper strand, q=4 working point
OMEGA = Q * Q / (A_W * A_W * MU0 * SIGMA)
DELTA = math.sqrt(2.0 / (OMEGA * MU0 * SIGMA))
RDC = 1.0 / (SIGMA * math.pi * A_W * A_W)


def F_skin(q):
    """Exact round-wire skin AC resistance ratio Rac/Rdc via Kelvin functions (q = a*sqrt(w*mu0*sigma))."""
    be, ke, bep, kep = kelvin(q)
    ber, bei, berp, beip = be.real, be.imag, bep.real, bep.imag
    return (q / 2.0) * (ber * beip - bei * berp) / (berp ** 2 + beip ** 2)


def prox_lowf(r, omega, H, sigma=SIGMA, mu=MU0):
    """Ferreira low-frequency proximity loss per unit length [W/m] of a round wire of radius r in a
    uniform transverse AC field of amplitude H -- the induced-eddy-dipole result, propto r^4."""
    return math.pi / 8.0 * sigma * omega ** 2 * mu ** 2 * r ** 4 * H ** 2


def build_wire_mesh(centers, a=A_W, pad=12 * A_W):
    """Air box containing round wires named wire0, wire1, ... at ``centers`` (after test_planar_*)."""
    xs = [c[0] for c in centers]; ys = [c[1] for c in centers]
    x0, x1 = min(xs) - pad, max(xs) + pad
    y0, y1 = min(ys) - pad, max(ys) + pad
    wires = []
    for k, (cx, cy) in enumerate(centers):
        w = WorkPlane().Circle(cx, cy, a).Face()
        w.faces.name = f"wire{k}"
        w.maxh = DELTA / 4.0
        wires.append(w)
    air = MoveTo(x0, y0).Rectangle(x1 - x0, y1 - y0).Face()
    for w in wires:
        air = air - w
    air.faces.name = "air"
    for sel in (air.edges.Max(X), air.edges.Min(X), air.edges.Max(Y), air.edges.Min(Y)):
        sel.name = "outer"
    return Mesh(OCCGeometry(Glue([air] + wires), dim=2).GenerateMesh(maxh=pad / 6))


def circuit_rac(mesh, conductors, connection="series"):
    sigma = mesh.MaterialCF({c: SIGMA for c in conductors}, default=0.0)
    nu = CoefficientFunction(1.0 / MU0)
    gfu = solve_planar_eddy_multi(mesh, nu, sigma, OMEGA, conductors,
                                  connection=connection, total_current=I_TOT, order=3)
    comps = gfu.components
    Az = comps[0]
    P = 0.0
    for k, c in enumerate(conductors):
        Vc = comps[1] if connection == "parallel" else comps[1 + k]
        P += ohmic_loss_2d(-1j * OMEGA * Az + Vc, mesh, sigma, region=c)
    return 2.0 * P / (I_TOT ** 2)


def test_litz_subdivision_scalings():
    """Analytic: splitting a solid wire into n strands of radius R0/sqrt(n) (same copper) cuts the skin
    self-field excess as 1/n^2 (low q) and the proximity excess as 1/n (exact)."""
    R0 = 1e-3
    # skin: evaluate at low q0 where F-1 ~ q^4 so the ratio is the clean (1/sqrt(n))^4 = 1/n^2
    q0 = 0.6
    for n in (3, 7, 19):
        qs = q0 / math.sqrt(n)
        ratio = (F_skin(qs) - 1.0) / (F_skin(q0) - 1.0)
        assert abs(ratio - 1.0 / n ** 2) < 0.01 / n ** 2, f"skin excess ~1/n^2 (n={n}, got {ratio:.4g})"
    # proximity: n strands of r=R0/sqrt(n) in the same field -> total = 1/n of the solid (algebraic, exact)
    omega = 1.0e5
    for n in (3, 7, 19):
        r = R0 / math.sqrt(n)
        ratio = n * prox_lowf(r, omega, 1.0) / prox_lowf(R0, omega, 1.0)
        assert abs(ratio - 1.0 / n) < 1e-12, f"proximity excess is exactly 1/n (n={n})"


def test_litz_strand_skin_fe_vs_kelvin():
    """FE anchor: one Litz strand (current-driven round wire) reproduces the exact Kelvin Rac/Rdc < 1 %,
    so the analytic F(q) used for the Litz scaling laws is the FE-validated one."""
    with TaskManager():
        mesh = build_wire_mesh([(0.0, 0.0)])
        sigma = mesh.MaterialCF({"wire0": SIGMA}, default=0.0)
        gfu = solve_planar_eddy(mesh, CoefficientFunction(1.0 / MU0), sigma, OMEGA,
                                driven_region="wire0", total_current=I_TOT, order=3)
        Az, Vc = gfu.components
        P = ohmic_loss_2d(-1j * OMEGA * Az + Vc, mesh, sigma, region="wire0")
    ratio_fem = (2.0 * P / I_TOT ** 2) / RDC
    ratio_exact = F_skin(Q)
    assert abs(ratio_fem - ratio_exact) / ratio_exact < 0.01, (ratio_fem, ratio_exact)


def test_litz_proximity_motivates_transposition():
    """FE: two strands carrying the same current dissipate ~2x a single strand when far apart, but MORE
    when brought close -- the proximity effect each strand sees in its neighbour's AC field.  This is the
    penalty helical transposition (every strand visits every position along the twist) averages away."""
    with TaskManager():
        rac_far = circuit_rac(build_wire_mesh([(-8 * A_W, 0.0), (8 * A_W, 0.0)]), ["wire0", "wire1"])
        rac_near = circuit_rac(build_wire_mesh([(-1.2 * A_W, 0.0), (1.2 * A_W, 0.0)]), ["wire0", "wire1"])
    assert rac_near / rac_far > 1.02, f"proximity should raise Rac (got {rac_near / rac_far:.4f})"


def test_litz_geometry_to_physics_handoff():
    """The strand count / radius / fill from the build123d archetype feed the loss laws directly: an
    n-strand bundle of strand radius rs is copper-equivalent to a solid wire of radius rs*sqrt(n), and its
    proximity loss is 1/n of that solid's."""
    n, rs = 7, 0.378e-3
    Rb = litz_packing_radius(n, rs)
    litz = litz_wire(n, rs, Rb, 6.0e-3, 4.0e-3, name="cu")
    assert len(litz.children) == n, "archetype builds n strands"
    ff = litz_fill_factor(n, rs, Rb + rs)                 # physical single-layer envelope
    assert 0.0 < ff < 1.0
    R_eq = rs * math.sqrt(n)                               # copper-equivalent solid radius
    assert abs(n * math.pi * rs ** 2 - math.pi * R_eq ** 2) < 1e-18, "equal copper area"
    prox_reduction = n * prox_lowf(rs, 1e5, 1.0) / prox_lowf(R_eq, 1e5, 1.0)
    assert abs(prox_reduction - 1.0 / n) < 1e-12, "n-strand bundle proximity loss is 1/n of the solid"


def main():
    print(f"  working point: f={OMEGA/2/math.pi:.0f} Hz, copper, solid radius {A_W*1e3:.1f} mm "
          f"(q={Q}, delta/a={DELTA/A_W:.3f})")
    print(f"  solid wire             Rac/Rdc = {F_skin(Q):.3f}   ({100*(F_skin(Q)-1):.0f}% eddy excess)")
    for n in (7, 19, 49):
        qs = Q / math.sqrt(n)
        print(f"  litz n={n:2d} (r={A_W/math.sqrt(n)*1e3:.3f} mm)  strand Rac/Rdc = {F_skin(qs):.4f}   "
              f"skin excess x1/{n**2}, proximity x1/{n}")
    test_litz_subdivision_scalings()
    test_litz_strand_skin_fe_vs_kelvin()
    test_litz_proximity_motivates_transposition()
    test_litz_geometry_to_physics_handoff()
    print("[OK] Litz AC-loss: subdivision scalings (skin 1/n^2, proximity 1/n), single-strand skin FE vs "
          "Kelvin < 1%, proximity penalty (transposition target) shown, archetype geometry -> loss laws.")


if __name__ == "__main__":
    main()
