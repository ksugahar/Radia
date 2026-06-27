"""Golden test for rad.ObjPyramid -- the 5-face pyramid MSC element solved by the multipole-moment MMM path.

A square-base pyramid has 5 faces (1 quad base + 4 triangles) -> 5 surface-charge DOF, the same count as
the wedge.  Its single quadrupole moment row is the per-element RESIDUAL EIGENMODE (the in-plane dx^2-dy^2
type mode for a symmetric pyramid -- DISTINCT from the wedge's axial mode; both derived in
examples/vim/eigenmode_quadrupole_derivation.wls).  So the pyramid solves through the SAME moment kernel as
hex (6) and wedge (5) with no element-specific code -- the eigenmode closure is what makes that possible.

Locks:
  (1) ObjPyramid creates a valid element and its permanent-magnet field is finite/non-zero;
  (2) a soft-iron pyramid in a uniform applied field magnetizes (moment path, 5 DOF), axially by symmetry;
  (3) a cube DECOMPOSED into 6 pyramids (apex at the cube center) reproduces the single-hexahedron cube
      demagnetization (both ~ 3*H0 = N=1/3) to within a few percent -- the cross-element-type consistency
      check that certifies the pyramid faces / normals / eigenmode are correct.
"""
from pathlib import Path
import sys
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
import radia as rad

MU0 = 4e-7 * np.pi
H0 = 1000.0


def _pyramid(base4, apex, M=(0, 0, 0)):
    """base4 = the 4 base-quad corners (v0..v3), apex = v4."""
    return rad.ObjPyramid([list(b) for b in base4] + [list(apex)], list(M))


def test_pyramid_exists_and_pm_field_finite():
    rad.UtiDelAll()
    L = 0.01
    base = [(-L, -L, 0), (L, -L, 0), (L, L, 0), (-L, L, 0)]
    p = _pyramid(base, (0, 0, 0.02), M=(0, 0, 954930.0))   # permanent magnet
    B = rad.Fld(p, 'b', [0.05, 0.0, 0.05])
    rad.UtiDelAll()
    assert len(B) == 3 and all(np.isfinite(B)) and np.linalg.norm(B) > 0


def test_pyramid_soft_iron_magnetizes_axially():
    """A symmetric soft-iron pyramid in a uniform Hz magnetizes in +z (moment path), transverse ~ 0."""
    rad.UtiDelAll(); rad.set_demag_backend("collocation_mmmm")
    L = 0.01
    base = [(-L, -L, -0.005), (L, -L, -0.005), (L, L, -0.005), (-L, L, -0.005)]
    p = _pyramid(base, (0, 0, 0.015))
    rad.MatApl(p, rad.MatLin(1000.0))
    rad.Solve(rad.ObjCnt([p, rad.ObjBckg(lambda q: [0, 0, MU0 * H0])]), 1e-6, 1000, 0)
    M = np.asarray(rad.ObjM(p)["magnetization"], float)
    rad.UtiDelAll()
    assert M[2] > 0, f"pyramid soft iron should magnetize in +z: M={M}"
    assert abs(M[0]) < 1e-3 * abs(M[2]) and abs(M[1]) < 1e-3 * abs(M[2]), f"transverse M not ~0: {M}"


def _cube_Mz(build):
    rad.UtiDelAll(); rad.set_demag_backend("collocation_mmmm")
    objs = build()
    for o in objs:
        rad.MatApl(o, rad.MatLin(1000.0))
    rad.Solve(rad.ObjCnt(objs + [rad.ObjBckg(lambda q: [0, 0, MU0 * H0])]), 1e-6, 1000, 0)
    Mz = float(np.mean([rad.ObjM(o)["magnetization"][2] for o in objs]))
    rad.UtiDelAll()
    return Mz


def test_cube_as_six_pyramids_matches_hex():
    """A cube meshed as 6 pyramids (each cube face -> apex at the cube center) gives the same
    demagnetization as the single-hexahedron cube (both ~ 3*H0, N=1/3) -- cross-element-type consistency."""
    L = 0.01
    c = [(-L, -L, -L), (L, -L, -L), (L, L, -L), (-L, L, -L),
         (-L, -L, L), (L, -L, L), (L, L, L), (-L, L, L)]
    faces = [[0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]]

    Mz_hex = _cube_Mz(lambda: [rad.ObjHexahedron([list(v) for v in c], [0, 0, 0])])
    Mz_pyr = _cube_Mz(lambda: [_pyramid([c[i] for i in f], (0, 0, 0)) for f in faces])

    assert 2.0 * H0 < Mz_hex < 4.0 * H0, f"hex cube M_z={Mz_hex:.1f} not ~3*H0 (N=1/3)"
    assert 2.0 * H0 < Mz_pyr < 4.0 * H0, f"6-pyramid cube M_z={Mz_pyr:.1f} not ~3*H0 (N=1/3)"
    rel = abs(Mz_pyr - Mz_hex) / abs(Mz_hex)
    assert rel < 0.12, f"6-pyramid cube M_z={Mz_pyr:.1f} != hex M_z={Mz_hex:.1f} (rel {rel:.2%})"


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
