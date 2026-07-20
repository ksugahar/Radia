"""Cross-check the axisymmetric H1Henrotte solver against a stored open-boundary
(Kelvin) regression reference -- not an analytical formula, but an independent
open-region FE solution kept as an opaque stored reference data file.

Reference data (stored regression reference):
    An axisymmetric permanent magnet -- a cylinder r in [0, 10] mm,
    z in [-10, 10] mm, material Hc = 3e5 A/m, mu_r = 1, magnetized axially
    (+z, theta = 90 deg) -- inside a circular air domain closed by a Kelvin
    (open-boundary) transformation.  The reference samples B_z along r = 7 mm,
    z in [-50, 50] mm (100 points) on a converged 100 mm domain and stores
    {Z, B} in a .mat data file.

The permanent-magnet field of a mu_r = 1 body is SCALE-INVARIANT (depends only
on shape ratios and M, not absolute size), so we reproduce the identical
geometry numerically in SI (lengths used as-is, mu0 = 4 pi e-7, Hc = 3e5) and
the B values match the reference's millimeter problem in tesla directly.

We solve the same magnet with solve_axi_magnetostatic on a large truncated
domain (open-boundary approximation, R_out = 500 >> 50 so truncation < 0.2 %),
sample B_z at the same (r=7, z) points, and require agreement with the stored
reference to within 3 % of the peak field.  This validates the axihenrotte
open-region magnetostatics end-to-end against an independent FE reference.

The reference .mat is loaded from RADIA_REGRESSION_DATA (or tests/data/) and is
treated as opaque; the test SKIPS cleanly if the file is not present.
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

# Stored open-boundary regression reference (opaque data file). Location is
# configurable; defaults to a local tests/data/ directory. The test skips if
# the file is absent, so it is portable and carries no internal path.
REF_DATA_DIR = os.environ.get(
    "RADIA_REGRESSION_DATA",
    os.path.join(os.path.dirname(__file__), "data"))
REF_MAT = os.path.join(REF_DATA_DIR, "axi_kelvin_magnet_ref.mat")

MU0 = 4e-7 * math.pi
A_MAG = 10.0     # magnet radius
H_MAG = 10.0     # magnet half-height
HC = 3.0e5       # coercivity [A/m]; B_rem = mu0*Hc = 0.377 T (mu_r = 1)
R_OUT = 500.0    # truncated open-domain radius
R_SAMPLE = 7.0   # reference sampling radius


def build_mesh(maxh_mag=0.4, maxh_near=2.0, maxh_far=50.0,
               r_near=70.0, r_out=R_OUT):
    from netgen.occ import OCCGeometry, MoveTo, WorkPlane, Glue, X, Y

    magnet = MoveTo(0, -H_MAG).Rectangle(A_MAG, 2 * H_MAG).Face()
    magnet.faces.name = "magnet"
    magnet.maxh = maxh_mag
    magnet.edges.Min(X).name = "axis"

    near_disk = WorkPlane().Circle(0, 0, r_near).Face()
    far_disk = WorkPlane().Circle(0, 0, r_out).Face()
    half = MoveTo(0, -r_out).Rectangle(r_out, 2 * r_out).Face()
    near = (near_disk * half) - magnet
    near.faces.name = "air"
    near.maxh = maxh_near
    near.edges.Min(X).name = "axis"
    far = (far_disk * half) - (near_disk * half)
    far.faces.name = "air"
    far.maxh = maxh_far
    far.edges.Min(X).name = "axis"
    far.edges.Max(X).name = "outer"
    far.edges.Max(Y).name = "outer"
    far.edges.Min(Y).name = "outer"

    from ngsolve import Mesh
    geo = OCCGeometry(Glue([far, near, magnet]), dim=2)
    return Mesh(geo.GenerateMesh(maxh=maxh_far))


def main():
    import numpy as np

    if not os.path.exists(REF_MAT):
        print(f"[SKIP] stored regression reference not found: {REF_MAT}")
        return
    import scipy.io as sio
    d = sio.loadmat(REF_MAT)
    Z = np.asarray(d["Z"], dtype=float).ravel()
    B_ref = np.asarray(d["B"], dtype=float).ravel()

    from ngsolve import (Mesh, CoefficientFunction, grad, TaskManager)
    from ngsolve import x as r_cf
    from radia_mcp.radia_ngsolve.solve import solve_axi_magnetostatic

    mesh = build_mesh()
    nu = CoefficientFunction(1.0 / MU0)   # mu_r = 1 everywhere

    with TaskManager():
        gfu = solve_axi_magnetostatic(
            mesh, nu, magnets={"magnet": (HC, 90.0)}, order=2)

    B_z = grad(gfu)[0] + gfu / r_cf

    # Sample B_z at the reference points (r = 7, z in [-50, 50]); off-axis so the
    # H1Henrotte gradient is reliable.
    B_ng = np.empty_like(B_ref)
    for i, z in enumerate(Z):
        B_ng[i] = B_z(mesh(R_SAMPLE, float(z)))

    peak = np.max(np.abs(B_ref))
    dev = np.abs(B_ng - B_ref)
    max_rel = np.max(dev) / peak
    rms_rel = np.sqrt(np.mean(dev**2)) / peak
    i0 = int(np.argmin(np.abs(Z)))

    print(f"  B_z@(r=7,z=0):  ref {B_ref[i0]:.6f} T   NGSolve {B_ng[i0]:.6f} T")
    print(f"  peak|B_ref| = {peak:.6f} T")
    print(f"  max dev / peak = {100*max_rel:.3f} %   rms dev / peak = {100*rms_rel:.3f} %")
    assert max_rel < 3e-2, f"axihenrotte vs stored ref off by {100*max_rel:.2f}% of peak (>3%)"
    print(f"\n[OK] axihenrotte open-region magnet cross-checked against the stored "
          f"open-boundary reference to {100*max_rel:.2f}% of peak.")


if __name__ == "__main__":
    main()
