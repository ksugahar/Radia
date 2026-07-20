"""ANTI-PERIODIC BC on a real MOTOR SECTOR (annular wedge, ROTATIONAL boundary
identification) -- the rotational counterpart of the translational strip check in
test_planar_periodic.py. This is the geometry used to solve ONE pole of a machine
and recover the whole cross-section by cyclic (anti-)symmetry.

Rigorous self-consistency (no closed form): the source Jz = cos(2*theta) =
(x^2-y^2)/r^2 is a 4-pole pattern that FLIPS SIGN under a +90deg rotation, i.e. it
is ANTI-PERIODIC over a quarter (pi/2) sector. So the full-annulus solution
(A_z=0 on both arcs) restricted to one quarter MUST equal the solution on a single
quarter annular wedge whose two radial cuts are ANTI-PERIODICALLY identified
(rotation by pi/2, phase = -1) with A_z=0 on the arcs. Agreement to ~2e-3 rel-RMS
(cross-mesh, full vs single sector at this resolution) confirms the rotational
periodic_h1 identification.

netgen gotcha (encoded below): a spline3 CIRCULAR arc takes the TANGENT-
INTERSECTION control point as its middle point ((R,R) for a 90deg (R,0)->(0,R)
arc), NOT a point on the arc -- and the boundary loop must be CCW or GenerateMesh
stalls.
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import (Mesh, H1, GridFunction, BilinearForm, LinearForm, grad, dx,
                     x, y, CoefficientFunction, Integrate, TaskManager)
from netgen.geom2d import SplineGeometry
from radia_mcp.radia_ngsolve.solve import periodic_h1

Ri, Ro = 1.0, 2.0
Jz = (x * x - y * y) / (x * x + y * y)            # = cos(2 theta) (4-pole, anti-periodic over pi/2)


def _solve(fes):
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += grad(u) * grad(v) * dx
    f = LinearForm(fes)
    f += Jz * v * dx
    a.Assemble()
    f.Assemble()
    g = GridFunction(fes)
    g.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    return g


def full_annulus_mesh(maxh=0.06):
    geo = SplineGeometry()
    geo.AddCircle((0, 0), Ro, leftdomain=1, rightdomain=0, bc="arc")
    geo.AddCircle((0, 0), Ri, leftdomain=0, rightdomain=1, bc="arc")
    return Mesh(geo.GenerateMesh(maxh=maxh))


def quarter_wedge_mesh(maxh=0.06):
    """Quarter annular wedge [0, pi/2] with the two RADIAL cuts identified for a
    +90deg rotation (copy=). CCW loop; spline3 arc control pt = tangent
    intersection (R,R)."""
    geo = SplineGeometry()
    a_ = geo.AppendPoint(Ri, 0.0)
    b_ = geo.AppendPoint(Ro, 0.0)
    cc = geo.AppendPoint(0.0, Ro)
    d_ = geo.AppendPoint(0.0, Ri)
    rad0 = geo.Append(["line", a_, b_], bc="cut", leftdomain=1, rightdomain=0)
    geo.Append(["spline3", b_, geo.AppendPoint(Ro, Ro), cc], bc="arc",
               leftdomain=1, rightdomain=0)
    geo.Append(["line", cc, d_], bc="cut", leftdomain=1, rightdomain=0, copy=rad0)
    geo.Append(["spline3", d_, geo.AppendPoint(Ri, Ri), a_], bc="arc",
               leftdomain=1, rightdomain=0)
    return Mesh(geo.GenerateMesh(maxh=maxh))


def main():
    mesh_full = full_annulus_mesh()
    mesh_sec = quarter_wedge_mesh()
    area = Integrate(CoefficientFunction(1.0), mesh_sec)
    area_expect = (math.pi / 4) * (Ro**2 - Ri**2)
    assert abs(area - area_expect) < 0.02, f"wedge area {area:.4f} != {area_expect:.4f}"

    with TaskManager():
        A_full = _solve(H1(mesh_full, order=3, dirichlet="arc"))
        A_sec = _solve(periodic_h1(mesh_sec, order=3, dirichlet="arc", antiperiodic=True))

    # Compare at interior sample points; normalise the worst by the FIELD SCALE
    # (max|A_full| over the samples) so the cos(2t)=0 zero-crossing at theta=45deg
    # does not create a spurious relative blow-up.
    diffs, refs = [], []
    for fr in (0.3, 0.5, 0.7):
        r = Ri + fr * (Ro - Ri)
        for th in (math.pi / 8, math.pi / 4, 3 * math.pi / 8):
            px, py = r * math.cos(th), r * math.sin(th)
            af = A_full(mesh_full(px, py))
            as_ = A_sec(mesh_sec(px, py)).real
            diffs.append(as_ - af)
            refs.append(af)
    rms = math.sqrt(sum(d * d for d in diffs) / sum(rr * rr for rr in refs))
    scale = max(abs(rr) for rr in refs)
    worst = max(abs(d) for d in diffs) / scale
    print(f"  wedge area = {area:.4f} (expect {area_expect:.4f})")
    print(f"  A_sec vs A_full: rel-RMS = {rms:.2e}, worst |d|/max|A| = {worst:.2e} "
          f"(over {len(diffs)} interior pts)")

    assert rms < 5e-3, f"anti-periodic wedge rel-RMS {rms:.2e} too large"
    assert worst < 1e-2, f"anti-periodic wedge worst {worst:.2e} too large"
    print("\n[OK] anti-periodic BC on a rotational motor-sector wedge validated "
          f"(full-vs-sector self-consistency, rel-RMS {rms:.1e}).")


if __name__ == "__main__":
    main()
