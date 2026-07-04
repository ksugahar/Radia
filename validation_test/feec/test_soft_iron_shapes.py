"""radia.vim mesh-less-SHAPE intent constructors: soft_iron_box / soft_iron_hex -> auto -> HDiv-VIM.

Item (b) of the HDiv-only gate (CLAUDE.md DIRECTION 2026-07-04): the mesh-less soft-iron capability is KEPT
via intent constructors that STRUCTURE-MESH a simple shape (subdivided structured hex mesh) and register it
through soft_iron_from_mesh, so rad.Solve auto-routes it to the FEEC HDiv-VIM (RT1) -- the API-compatible
replacement for the legacy ObjHexahedron + MatApl(MatLin) mesh-less collocation-MMMM route.  The returned
container IS the sub-mesh, so rad.Fld reflects the resolved per-sub-element M (no write-back plumbing).

VALIDATED 2026-07-05: box demag 0.33284 (== the earlier full-hex HDiv ground truth), box==hex(box corners)
to 0.0 (trilinear CHEXA order correct), nonlinear bh_table path works.
"""
import math

import numpy as np
import pytest

pytest.importorskip("ngsolve")
import radia as rad  # noqa: E402
import radia.vim as vim  # noqa: E402
import ngsolve as ng  # noqa: E402

A = 0.01
MU_R = 1000.0
MU0 = 4.0e-7 * math.pi
H0 = 1.0e4
PROBE = [0.0, 0.0, 3 * A]


@pytest.mark.flaky(reruns=3, reruns_delay=1)   # NGSolve GetTrafo lattice first-touch flake under contention
def test_soft_iron_box_auto_routes_to_hdiv():
    """soft_iron_box -> rad.Solve auto -> HDiv-VIM (dict return with demag), cube demag ~1/3, Mz_avg
    matches a direct full-hex HDiv solve; rad.Fld reflects the resolved field."""
    rad.UtiDelAll()
    with ng.TaskManager():
        cont = vim.soft_iron_box(center=(0, 0, 0), size=(2 * A, 2 * A, 2 * A), mu_r=MU_R, nsub=4)
        from radia.vim import _radsolve
        assert _radsolve.is_hdiv_eligible(cont)       # mesh-backed hex sub-mesh -> HDiv-eligible
        res = rad.Solve(rad.ObjCnt([cont, rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])]), 1e-6, 2000, 0)
        B = np.array(rad.Fld(cont, 'b', PROBE), float)
    assert isinstance(res, dict) and "demag" in res, "soft_iron_box did not auto-route to the HDiv-VIM"
    assert 0.30 < res["demag"] < 0.36, f"box demag {res['demag']:.4f} not ~1/3 (energy demag)"
    # Mz_avg matches the direct full-hex HDiv ground truth (~36111) to a few %
    assert 3.4e4 < res["M_avg"][2] < 3.8e4, f"box Mz_avg {res['M_avg'][2]:.1f} out of band"
    assert abs(B[2]) > 1e-4, f"box external Bz {B[2]:.2e} implausibly small (Fld not reflecting resolved M)"
    rad.UtiDelAll()


@pytest.mark.flaky(reruns=3, reruns_delay=1)
def test_soft_iron_hex_box_corners_equals_box():
    """soft_iron_hex given the 8 corners of the SAME box == soft_iron_box -- locks the trilinear CHEXA
    corner ordering (a twisted order would change the mesh -> a different demag)."""
    rad.UtiDelAll()
    with ng.TaskManager():
        box = vim.soft_iron_box(center=(0, 0, 0), size=(2 * A, 2 * A, 2 * A), mu_r=MU_R, nsub=4)
        rb = rad.Solve(rad.ObjCnt([box, rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])]), 1e-6, 2000, 0)
        rad.UtiDelAll()
        verts = [[-A, -A, -A], [A, -A, -A], [A, A, -A], [-A, A, -A],
                 [-A, -A, A], [A, -A, A], [A, A, A], [-A, A, A]]
        hexi = vim.soft_iron_hex(verts, mu_r=MU_R, nsub=4)
        rh = rad.Solve(rad.ObjCnt([hexi, rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])]), 1e-6, 2000, 0)
    assert abs(rh["demag"] - rb["demag"]) < 1e-3, \
        f"soft_iron_hex(box corners) demag {rh['demag']:.5f} != soft_iron_box {rb['demag']:.5f} (trilinear order?)"
    rad.UtiDelAll()


@pytest.mark.flaky(reruns=3, reruns_delay=1)
def test_soft_iron_box_nonlinear():
    """soft_iron_box with a bh_table (nonlinear) auto-routes to the HDiv-VIM energy-Newton path."""
    rad.UtiDelAll()
    BH = [[0, 0], [100, 0.5], [1000, 1.2], [10000, 1.8], [100000, 2.0]]
    with ng.TaskManager():
        cont = vim.soft_iron_box(center=(0, 0, 0), size=(2 * A, 2 * A, 2 * A), bh_table=BH, nsub=3)
        res = rad.Solve(rad.ObjCnt([cont, rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * 5 * H0])]), 1e-6, 2000, 0)
    assert isinstance(res, dict) and res.get("nonlinear") is True, "box bh_table did not run the nonlinear HDiv path"
    assert res["M_avg"][2] > 1.0e5, f"nonlinear box Mz_avg {res['M_avg'][2]:.1f} too small (expect saturated)"
    rad.UtiDelAll()


def test_soft_iron_box_input_validation():
    """Fail-loud (No-Fallbacks) on bad input: both/neither of mu_r/bh_table, non-positive size."""
    with pytest.raises(ValueError):
        vim.soft_iron_box(center=(0, 0, 0), size=(2 * A, 2 * A, 2 * A))                 # neither
    with pytest.raises(ValueError):
        vim.soft_iron_box(center=(0, 0, 0), size=(2 * A, 2 * A, 2 * A), mu_r=1000, bh_table=[[0, 0], [1, 1]])
    with pytest.raises(ValueError):
        vim.soft_iron_box(center=(0, 0, 0), size=(0.0, 2 * A, 2 * A), mu_r=1000)          # non-positive size
    with pytest.raises(ValueError):
        vim.soft_iron_hex(np.zeros((4, 3)), mu_r=1000)                                    # not 8x3
