r"""AGE step 3 -- COMPLEX (eddy-current) gap coupling: the eddy-current gap-coupling axis.

The air gap is current-free (sigma=0), so even in the TIME-HARMONIC (eddy) regime it
stays Laplace -- the SAME real harmonic DtN couples a CONDUCTING (eddy-current) FE region
across the UN-MESHED gap.  That is exactly the capability a magnetostatic air-gap element
cannot give, and that mainstream machine codes only get via a meshed sliding band.

Validation (independent, fully-meshed reference -- no Bessel needed): a conducting annulus
ro<r<Rext (sigma!=0, jw*mu*sigma eddy term, ~3.6 skin depths thick) driven by a rotor-ring
phasor at ri.  Solve it TWO ways and compare the field over ro<r<Rext:
  (A) AGE       : mesh ONLY ro<r<Rext; the gap is its real DtN Robin operator (un-meshed);
  (B) reference : mesh the WHOLE ri<r<Rext (gap sigma=0 + conductor), Dirichlet at ri.
Agreement proves the complex AGE represents the un-meshed gap EXACTLY with eddy currents
present.  Single harmonic n; complex NGSolve solve.
"""
import math
import os
import sys

from ngsolve import (Mesh, H1, GridFunction, BilinearForm, LinearForm, grad, dx, ds,
                     x, y, atan2, cos, BND, TaskManager)
from netgen.geom2d import SplineGeometry

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.radia_ngsolve.airgap_element import annular_dtn_matrix

RI, RO, REXT = 0.10, 0.11, 0.20
N = 2
A_IN, A_EXT = 1.0 + 0j, 0.3 + 0j
MU0 = 4e-7 * math.pi
SIGMA, FREQ = 2.0e6, 200.0          # ~2.5 cm skin depth in a 9 cm conductor -> strong eddy
OMEGA = 2.0 * math.pi * FREQ


def _inv(mat, freedofs):
    try:
        return mat.Inverse(freedofs, inverse="umfpack")
    except Exception:
        return mat.Inverse(freedofs, inverse="pardiso")


def age_solve():
    geo = SplineGeometry()
    geo.AddCircle((0, 0), REXT, leftdomain=1, rightdomain=0, bc="outer")
    geo.AddCircle((0, 0), RO, leftdomain=0, rightdomain=1, bc="inner")
    mesh = Mesh(geo.GenerateMesh(maxh=0.006)); mesh.Curve(5)
    cosn = cos(N * atan2(y, x))
    (_, _), (M21, M22) = annular_dtn_matrix(N, RI, RO)
    fes = H1(mesh, order=5, complex=True, dirichlet="outer")
    A, v = fes.TnT()
    dsi = ds(definedon=mesh.Boundaries("inner"))
    a = BilinearForm(fes)
    a += grad(A) * grad(v) * dx + 1j * OMEGA * MU0 * SIGMA * A * v * dx + M22 * A * v * dsi
    f = LinearForm(fes)
    f += (-M21 * A_IN) * cosn * v * dsi
    with TaskManager():
        a.Assemble(); f.Assemble()
        gfu = GridFunction(fes)
        gfu.Set(A_EXT * cosn, definedon=mesh.Boundaries("outer"))
        gfu.vec.data += _inv(a.mat, fes.FreeDofs()) * (f.vec - a.mat * gfu.vec)
    return mesh, gfu


def ref_solve():
    geo = SplineGeometry()
    geo.AddCircle((0, 0), REXT, leftdomain=1, rightdomain=0, bc="outer")
    geo.AddCircle((0, 0), RO, leftdomain=2, rightdomain=1)
    geo.AddCircle((0, 0), RI, leftdomain=0, rightdomain=2, bc="inner")
    geo.SetMaterial(1, "fe"); geo.SetMaterial(2, "gap")
    mesh = Mesh(geo.GenerateMesh(maxh=0.006)); mesh.Curve(5)
    cosn = cos(N * atan2(y, x))
    sigma = mesh.MaterialCF({"fe": SIGMA, "gap": 0.0}, default=0.0)
    fes = H1(mesh, order=5, complex=True, dirichlet="inner|outer")
    A, v = fes.TnT()
    a = BilinearForm(fes)
    a += grad(A) * grad(v) * dx + 1j * OMEGA * MU0 * sigma * A * v * dx
    f = LinearForm(fes)
    with TaskManager():
        a.Assemble(); f.Assemble()
        gfu = GridFunction(fes)
        gfu.Set(mesh.BoundaryCF({"inner": A_IN * cosn, "outer": A_EXT * cosn}), BND)
        gfu.vec.data += _inv(a.mat, fes.FreeDofs()) * (f.vec - a.mat * gfu.vec)
    return mesh, gfu


def test_complex_age_matches_fully_meshed_with_eddy():
    mesh_a, gfu_a = age_solve()
    mesh_r, gfu_r = ref_solve()
    diffs, refs = [], []
    for frac in (0.1, 0.5, 0.9):
        rr = RO + frac * (REXT - RO)
        for th in (0.3, 1.2, 2.5, -0.8):
            px, py = rr * math.cos(th), rr * math.sin(th)
            diffs.append(abs(gfu_a(mesh_a(px, py)) - gfu_r(mesh_r(px, py))))
            refs.append(abs(gfu_r(mesh_r(px, py))))
    rel = max(diffs) / max(refs)
    pmid = gfu_r(mesh_r((RO + REXT) / 2.0, 0.0))
    print(f"complex AGE vs fully-meshed (eddy, {SIGMA:.0e} S/m, {FREQ:.0f} Hz): "
          f"max rel err = {rel:.2e}; sample A = {pmid:.4f} "
          f"(|Im/Re| = {abs(pmid.imag / pmid.real):.2f} -> genuine eddy phase)")
    assert abs(pmid.imag) > 0.05 * abs(pmid.real), "eddy should make A complex (else trivial test)"
    assert rel < 1e-2, f"complex AGE off the fully-meshed reference by {rel:.2e}"


def main():
    test_complex_age_matches_fully_meshed_with_eddy()
    print("[OK] AGE step 3: the COMPLEX gap DtN couples a conducting (eddy-current) FE region "
          "across the UN-MESHED gap, matching the fully-meshed reference. Eddy x mesh-free gap "
          "-- the eddy x mesh-free-gap axis that magnetostatic air-gap elements and "
          "sliding-band methods do not offer.")


if __name__ == "__main__":
    main()
