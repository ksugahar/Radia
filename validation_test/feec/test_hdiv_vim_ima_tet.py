"""Golden: RT1 pure-TET HDiv-VIM supports IMA (image method / mirror symmetry), so HDiv-VIM can serve as
radia's MMM (2026-07-04).  A reduced (1/2, 1/4, 1/8) symmetry model WITH the mirror image reproduces the
FULL model's demag factor (a z-magnetised sphere = 1/3), while the reduced model WITHOUT the image is wrong.

The C++ highorder charge Gram folds the mirror-image charge interactions into every entry
  G_IMA(a,b) = G(a,b) + sum_i sign_i * 0.5*(QuadDotRefl(a,b,mask_i) + QuadDotRefl(b,a,mask_i))
with QuadDotRefl(tgt,src,mask) = the source's PhiInner potential at tgt's outer points reflected on the mask
axes (mirror isometry).  Physics + reflection/sign convention validated against the existing RT0 analytic
IMA (memory hdiv-tet-hex-coupling-pyramid-gated).  IMA is wired for the FLAT pure-TET RT1 path; hex/wedge
and curved tet fail loud toward collocation MMMM (locked below).
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")
import ngsolve as ng                                       # noqa: E402
from netgen.csg import CSGeometry, Sphere, OrthoBrick, Pnt  # noqa: E402
from radia.vim._solve import hdiv_demag_solve              # noqa: E402

_H0 = 1000.0


def _sphere(box_lo, maxh=0.45):
    g = CSGeometry()
    g.Add(Sphere(Pnt(0, 0, 0), 1.0) * OrthoBrick(Pnt(*box_lo), Pnt(2, 2, 2)))
    with ng.TaskManager():
        return ng.Mesh(g.GenerateMesh(maxh=maxh))


def _demag(box_lo, image):
    with ng.TaskManager():
        out = hdiv_demag_solve(_sphere(box_lo), mu_r=100.0,
                               H_ext=ng.CoefficientFunction((0, 0, _H0)), image=image)
    return out["demag"], out["n_el"]


@pytest.mark.parametrize("tag,box_lo,image", [
    ("half", (0, -2, -2), "+x"),           # z-field parallel to the x=0 mirror -> +x
    ("quarter", (0, 0, -2), "+x+y"),
    ("eighth", (0, 0, 0), "+x+y-z"),       # z-field perpendicular to the z=0 mirror -> -z
])
def test_tet_ima_reduced_reproduces_full_demag(tag, box_lo, image):
    """The image-augmented reduced tet model reproduces the full-sphere demag 1/3 (the IMA payoff)."""
    D, ne = _demag(box_lo, image)
    assert abs(D - 1.0 / 3.0) < 0.02, f"{tag}: image demag {D:.4f} off 1/3 (n_el={ne})"


def test_tet_reduced_without_image_is_wrong():
    """Sanity: the SAME half mesh WITHOUT the image is NOT a full sphere -> demag far from 1/3 (so the
    passing image cases above are genuinely the IMA folding, not a coincidence)."""
    D, _ = _demag((0, -2, -2), None)
    assert D < 0.30, f"half-model no-image demag {D:.4f} should be well below 1/3"


def test_hex_wedge_curved_ima_fail_loud():
    """IMA is wired for the FLAT pure-TET path only; hex / wedge / curved-tet fail loud toward collocation
    MMMM (No-Fallbacks -- never silently drop the image)."""
    from ngsolve.meshes import MakeStructured3DMesh
    mp = lambda x, y, z: (0.01 * (x - 0.5), 0.01 * (y - 0.5), 0.01 * (z - 0.5))  # noqa: E731
    with ng.TaskManager():
        hexm = MakeStructured3DMesh(hexes=True, nx=2, ny=2, nz=2, mapping=mp)
    with pytest.raises((ValueError, NotImplementedError), match="TET"):
        with ng.TaskManager():
            hdiv_demag_solve(hexm, mu_r=100.0, H_ext=ng.CoefficientFunction((0, 0, _H0)), image="+x")
    # curved tet + image also fails loud
    tetm = _sphere((0, -2, -2))
    with pytest.raises((ValueError, NotImplementedError), match="CURVED|curve"):
        with ng.TaskManager():
            hdiv_demag_solve(tetm, mu_r=100.0, H_ext=ng.CoefficientFunction((0, 0, _H0)),
                             image="+x", curve_order=2)
