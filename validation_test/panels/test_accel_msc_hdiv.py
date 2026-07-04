"""calc_accel_msc.solve_msc with the FEEC HDiv-VIM backend (the only backend exposed by this panel).

The accelerator-magnet MSC panel drives the HDiv-VIM (radia.vim.MeshSoftIron + rad.Solve, which
routes a mesh-backed TET soft iron to HDiv-VIM RT1).  This panel intentionally does not expose the
mesh-less multipole-moment MMM backend.  The HDiv-VIM is KELVIN-less / iron-only / TET-only, so the .vol
must contain only the 'yoke' volume material, and IMA symmetry is not supported there.

Locks:
  (1) on an iron-only yoke + coil, the hdiv backend magnetizes the iron (converged, M_avg large) -- the
      panel routes the coil field into the HDiv-VIM solve correctly;
  (2) demag_backend='hdiv' on a MULTI-material .vol (yoke + air) returns a clean error (iron-only);
  (3) demag_backend='hdiv' with an IMA symmetry string returns a clean error (not supported yet);
  (4) demag_backend='collocation_mmmm' returns a clean error (not a backend for this panel).

Also guards the 2026-06-17 coil-bug fix: solve_msc used to rad.UtiDelAll() AFTER building the coil
(destroying it); now the iron responds to the coil, so M_avg is large (~1e5 A/m), not ~0.

Slow: each solve runs the real HDiv-VIM solver on a small tet yoke.
"""
import math
import os
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "src", "radia"))
sys.path.insert(0, os.path.join(REPO, "src", "radia", "panels"))

pytest.importorskip("ngsolve")
pytest.importorskip("radia")

pytestmark = [
    pytest.mark.filterwarnings("ignore:Gimbal lock:UserWarning"),
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
]

import ngsolve as ng  # noqa: E402

L = 0.01          # 10 mm cube yoke
MU_R = 1000.0


@pytest.fixture(scope="module")
def coil_script(tmp_path_factory):
    """A tiny 360-deg current loop (radius 30 mm, 2000 A) -> axial field through the yoke at origin."""
    p = tmp_path_factory.mktemp("accel_msc_hdiv") / "loop_coil.py"
    p.write_text(
        "import sys, os\n"
        f"sys.path.insert(0, r'{os.path.join(REPO, 'src', 'radia')}')\n"
        "from coil_builder import CoilBuilder\n"
        "def build_coil():\n"
        "    return (CoilBuilder(current=2000).set_start([0.03, 0, 0])\n"
        "            .set_cross_section(0.006, 0.006).add_arc(0.03, 360))\n"
    )
    return str(p)


def _iron_only_yoke_vol(path):
    """A small tet iron cube centred at origin, single material 'yoke', saved as a Netgen .vol."""
    from netgen.csg import CSGeometry, OrthoBrick, Pnt
    geo = CSGeometry()
    geo.Add(OrthoBrick(Pnt(-L / 2, -L / 2, -L / 2), Pnt(L / 2, L / 2, L / 2)).mat("yoke"))
    with ng.TaskManager():
        mesh = ng.Mesh(geo.GenerateMesh(maxh=L / 3))
    mesh.ngmesh.Save(str(path))
    return str(path)


def _yoke_plus_air_vol(path):
    """A two-material tet mesh (yoke | air) -> NOT iron-only, to exercise the hdiv guard."""
    from netgen.csg import CSGeometry, OrthoBrick, Pnt
    geo = CSGeometry()
    geo.Add(OrthoBrick(Pnt(-L, -L, -L), Pnt(0, L, L)).mat("yoke"))
    geo.Add(OrthoBrick(Pnt(0, -L, -L), Pnt(L, L, L)).mat("air"))
    with ng.TaskManager():
        mesh = ng.Mesh(geo.GenerateMesh(maxh=L))
    mesh.ngmesh.Save(str(path))
    return str(path)


def _linear_mat():
    from em_material import EMMaterial
    return EMMaterial(name="linear", sigma=2e6, mu_r=MU_R, bh_curve=None)


def test_hdiv_panel_magnetizes(tmp_path, coil_script):
    """The hdiv backend solves the iron-only yoke + coil: converged, HDiv-VIM, and the iron is actually
    magnetised by the coil (M_avg large -- the coil-bug-fix guard)."""
    from calc_accel_msc import solve_msc
    vol = _iron_only_yoke_vol(tmp_path / "yoke.vol")
    rh = solve_msc(coil_script=coil_script, vol_file=vol, mat=_linear_mat(),
                   demag_backend="hdiv", solver=0, tol=1e-7, max_iter=800)
    assert "error" not in rh, rh
    assert rh["converged"]
    assert rh["demag_backend"] == "hdiv" and rh["solver"] == "HDiv-VIM"
    mz_h = rh["M_avg"][2]
    # coil-bug-fix guard: the iron must actually be magnetised by the coil (not ~0)
    assert abs(mz_h) > 1e4, f"hdiv M_avg_z={mz_h:.1f} too small -- coil field missing?"


def test_collocation_mmmm_backend_rejected(tmp_path, coil_script):
    """demag_backend='collocation_mmmm' returns a clean error -- this panel is HDiv-VIM only."""
    from calc_accel_msc import solve_msc
    vol = _iron_only_yoke_vol(tmp_path / "yoke.vol")
    r = solve_msc(coil_script=coil_script, vol_file=vol, mat=_linear_mat(),
                  demag_backend="collocation_mmmm", solver=0, tol=1e-6, max_iter=400)
    assert "error" in r, r
    assert "collocation_mmmm" in r["error"].lower() and "panel" in r["error"].lower()


def test_hdiv_rejects_multimaterial(tmp_path, coil_script):
    """demag_backend='hdiv' on a yoke+air .vol -> clean error (the VIM is iron-only / KELVIN-less)."""
    from calc_accel_msc import solve_msc
    vol = _yoke_plus_air_vol(tmp_path / "yoke_air.vol")
    r = solve_msc(coil_script=coil_script, vol_file=vol, mat=_linear_mat(),
                  demag_backend="hdiv", solver=0, tol=1e-6, max_iter=400)
    assert "error" in r, r
    assert "iron-only" in r["error"].lower()


def test_hdiv_ima_rejected(tmp_path, coil_script):
    """demag_backend='hdiv' with an IMA image string returns a clean error; IMA is collocation MMMM scope."""
    from calc_accel_msc import solve_msc
    p = _iron_only_yoke_vol(tmp_path / "half_yoke.vol")
    r = solve_msc(coil_script=coil_script, vol_file=str(p), mat=_linear_mat(),
                  demag_backend="hdiv", ima="-z", solver=0, tol=1e-7, max_iter=800)
    assert "error" in r, r
    assert "ima" in r["error"].lower() or "image" in r["error"].lower()


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
