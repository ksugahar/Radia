"""Golden test for the EIEM2-retirement fail-loud guards (Phase 3b).

The parameter-free moment-yano formulation (BuildMomentSystemCore) is the SOLE surface-charge (MSC)
demag solver; the EIEM2 collocation kernel has been retired.  Two element-composition cases the EIEM2
kernel used to cover are NOT representable by moment and are now rejected fail-loud (No Fallbacks -- a
silent wrong number is worse than an error), raised from radTApplication::MakeAutoRelax:

  (1) Radia::Error204 -- a single demag Solve that MIXES MMM elements (tetrahedron / RecMag, 3 DOF) with
      MSC surface-charge elements (hexahedron / wedge, 5/6 DOF).  Solve them as separate containers.
  (2) Radia::Error205 -- B-input hysteresis (b_input_newton / b_input_hantila) on MSC elements; the
      B-input Newton/Hantila Jacobian is block-3x3 (MMM) only.  Hysteresis stays on tetrahedron / RecMag.

This test also locks that the KEPT paths are unaffected: a pure-hex (moment) and a pure-tet (MMM)
soft-iron solve still magnetize in an applied field, and permanent-magnet field evaluation of a mixed
container (no Solve) is untouched.
"""
import sys
import numpy as np
import pytest

sys.path.insert(0, r"S:\Radia\01_GitHub\src\radia")
import radia as rad

MU0 = 4e-7 * np.pi
H0 = 1000.0


def _hex(cx, L=0.01):
    return [[cx - L, -L, -L], [cx + L, -L, -L], [cx + L, L, -L], [cx - L, L, -L],
            [cx - L, -L, L], [cx + L, -L, L], [cx + L, L, L], [cx - L, L, L]]


def _tet(cx, L=0.01):
    return [[cx, 0, 0], [cx + L, 0, 0], [cx + L * 0.5, L * 0.87, 0], [cx + L * 0.5, L * 0.29, L * 0.82]]


def test_mixed_mmm_msc_solve_raises():
    """A single demag Solve mixing a tet (MMM, 3 DOF) and a hex (MSC, 6 DOF) soft iron -> Error204."""
    rad.UtiDelAll(); rad.set_demag_backend("yano")
    t = rad.ObjTetrahedron(_tet(0.0), [0, 0, 0]); rad.MatApl(t, rad.MatLin(1000.0))
    h = rad.ObjHexahedron(_hex(0.1), [0, 0, 0]); rad.MatApl(h, rad.MatLin(1000.0))
    cont = rad.ObjCnt([t, h, rad.ObjBckg(lambda p: [0, 0, MU0 * H0])])
    with pytest.raises(RuntimeError, match="mixes MMM"):
        rad.Solve(cont, 1e-6, 100, 0)
    rad.UtiDelAll()


def test_binput_hysteresis_on_msc_raises():
    """B-input hysteresis (b_input_newton) on a hex (MSC) soft iron -> Error205 (3-DOF MMM only)."""
    rad.UtiDelAll(); rad.set_demag_backend("yano")
    K = 3
    eta = np.array([0.0, 0.5, 1.0])
    r = np.linspace(0, 2.0, 20)
    f_k = [(r.tolist(), (1000.0 * (k + 1) * r).tolist()) for k in range(K)]
    mat = rad.MatPlayHysteresis(K, eta, f_k)
    h = rad.ObjHexahedron(_hex(0.0), [0, 0, 0]); rad.MatApl(h, mat)
    cont = rad.ObjCnt([h, rad.ObjBckg(lambda p: [0, 0, MU0 * H0])])
    rad.SolverConfig(b_input_newton=True)
    try:
        with pytest.raises(RuntimeError, match="B-input hysteresis"):
            rad.Solve(cont, 1e-6, 100, 0)
    finally:
        rad.SolverConfig(b_input_newton=False)
        rad.UtiDelAll()


def test_pure_hex_moment_still_solves():
    """KEPT: a pure-hex soft iron magnetizes in an applied field via the moment-yano solver."""
    rad.UtiDelAll(); rad.set_demag_backend("yano")
    h = rad.ObjHexahedron(_hex(0.0), [0, 0, 0]); rad.MatApl(h, rad.MatLin(1000.0))
    rad.Solve(rad.ObjCnt([h, rad.ObjBckg(lambda p: [0, 0, MU0 * H0])]), 1e-6, 1000, 0)
    Mz = rad.ObjM(h)["magnetization"][2]
    rad.UtiDelAll()
    assert 2.0 * H0 < Mz < 4.0 * H0, f"pure-hex moment M_z={Mz:.1f} (expected ~3*H0)"


def test_pure_tet_mmm_still_solves():
    """KEPT: a pure-tet soft iron magnetizes in an applied field via the MMM solver (unaffected)."""
    rad.UtiDelAll(); rad.set_demag_backend("auto")
    t = rad.ObjTetrahedron(_tet(0.0), [0, 0, 0]); rad.MatApl(t, rad.MatLin(1000.0))
    rad.Solve(rad.ObjCnt([t, rad.ObjBckg(lambda p: [0, 0, MU0 * H0])]), 1e-6, 1000, 0)
    Mz = rad.ObjM(t)["magnetization"][2]
    rad.UtiDelAll()
    assert Mz > 0, f"pure-tet MMM M_z={Mz:.1f} (should magnetize)"


def test_mixed_container_permanent_magnet_field_untouched():
    """KEPT: a mixed tet+hex container of PERMANENT magnets (no Solve) evaluates fields normally --
    the guard only fires on a demag Solve, not on field evaluation."""
    rad.UtiDelAll()
    t = rad.ObjTetrahedron(_tet(0.0), [0, 0, 954930])
    h = rad.ObjHexahedron(_hex(0.1), [0, 0, 954930])
    cont = rad.ObjCnt([t, h])
    B = rad.Fld(cont, 'b', [0.05, 0, 0.05])
    rad.UtiDelAll()
    assert len(B) == 3 and all(np.isfinite(B))


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
