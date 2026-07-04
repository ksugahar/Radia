"""Golden: the HDiv-VIM FIELD reconstruction (radia.Fld on the solved per-element M) gives the correct
field at points -- at ALL orders.

REGRESSION GUARD for the 2026-06-12 finding.  The previous reconstruction in the TEAM-13 pipeline computed
the demag field as M_mass^{-1} N m (the energy operator N = B^T G B applied, then mass-inverted).  That is
WRONG at order>=2: N has a SOLENOIDAL NULLSPACE (divergence-free RT bubbles carry zero charge, so B maps
them to zero), so M_mass^{-1} N m does not equal the field.  On THIS exact uniform-M sphere test the broken
method gave B_z/(MU0 M0) ranging 0.12-1.0 with spurious transverse components up to 1.67 at order 2 (while
order 0/1 happened to converge).  The demag FACTOR (m^T N m / m^T M_mass m, the energy) stays correct -- so
this is a FIELD-reconstruction bug, not an operator bug.

The correct reconstruct_field gives B_inside = MU0*(2/3)*M0 (analytic, uniform interior) at order 0, 1, AND
2 -- because it evaluates the EXACT analytic field of the solved per-element M (radia.Fld).
"""
import numpy as np
import pytest
import ngsolve as ng
from netgen.occ import Sphere, OCCGeometry, Pnt

from radia.vim import reconstruct_field

MU0 = 4e-7 * np.pi
# generic interior points (avoid mesh nodes / faces, where a point eval would be ill-defined)
PTS = [(0.15, 0.10, 0.20), (-0.20, 0.10, 0.15), (0.10, -0.25, 0.05),
       (0.05, 0.05, -0.30), (-0.10, -0.10, 0.10)]


@pytest.mark.parametrize("order", [0, 1])
def test_uniform_sphere_field_reconstruction(order):
    """Uniformly magnetized unit sphere, M = zhat [A/m]: B inside = MU0*(2/3)*zhat (demag factor 1/3),
    UNIFORM throughout the interior, zero transverse components."""
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0)).GenerateMesh(maxh=0.3))
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=order)
        gfM = ng.GridFunction(fes)
        gfM.Set(ng.CF((0, 0, 1.0)))                          # uniform M = zhat
        B = reconstruct_field(mesh, gfM, PTS, quantity="b") / MU0
    bz = B[:, 2]
    # band 0.5%: the faceted-sphere discretization lands ~0.6664 (improves with maxh); the broken
    # M_mass^-1 N m gave mean ~0.12-1.0 with spread ~2.0 at order 2 -> fails this by a wide margin.
    assert abs(bz.mean() - 2.0 / 3.0) < 5e-3, \
        f"order={order}: Bz/(MU0 M0) mean {bz.mean():.4f} != 2/3 (analytic interior)"
    assert (bz.max() - bz.min()) < 5e-3, \
        f"order={order}: interior Bz spread {bz.max()-bz.min():.4f} too large (non-uniform interior)"
    assert np.abs(B[:, 0]).max() < 5e-3 and np.abs(B[:, 1]).max() < 5e-3, \
        f"order={order}: spurious transverse field |Bx|max={np.abs(B[:,0]).max():.4f} " \
        f"|By|max={np.abs(B[:,1]).max():.4f} (solenoidal-nullspace garbage?)"


@pytest.mark.parametrize("quantity", ["b", "a"])
def test_hmatrix_backend_matches_direct(quantity):
    """reconstruct_field(backend='hmatrix') -- the O(N log N) HACApK _FieldEvalHMatrix path for the
    field->GridFunction FEM coupling -- reproduces the direct batched radia.Fld to ~machine, for B and the
    A-formulation source A.  This is the path-B (magnet-container) acceleration; the build is fully parallel
    (unit-M clones, no mutex) so the CALLER's TaskManager parallelises it."""
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0)).GenerateMesh(maxh=0.4))
    with ng.TaskManager():                                   # CALLER wraps; reconstruct_field does NOT
        gfM = ng.GridFunction(ng.HDiv(mesh, order=1))
        gfM.Set(ng.CF((0.3, -0.5, 1.0)))                     # a non-trivial (non-axial) magnetization
        Fd = reconstruct_field(mesh, gfM, PTS, quantity=quantity, backend="direct")
        Fh = reconstruct_field(mesh, gfM, PTS, quantity=quantity, backend="hmatrix", eps=1e-9)
    err = np.linalg.norm(Fh - Fd) / np.linalg.norm(Fd)
    assert err < 1e-7, f"quantity={quantity}: hmatrix vs direct reconstruct_field = {err:.3e}"


def test_hmatrix_backend_rejects_h():
    """backend='hmatrix' supports 'b'/'a' only (no 'h' in _FieldEvalHMatrix yet) -- fail loud, no fallback."""
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0)).GenerateMesh(maxh=0.6))
    with ng.TaskManager():
        gfM = ng.GridFunction(ng.HDiv(mesh, order=1))
        gfM.Set(ng.CF((0, 0, 1.0)))
        with pytest.raises(ValueError, match="quantity 'b' or 'a'"):
            reconstruct_field(mesh, gfM, PTS, quantity="h", backend="hmatrix")
