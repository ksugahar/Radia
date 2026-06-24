"""calc_accel_msc.solve_msc with the FEEC HDiv-VIM backend (the only backend exposed by this panel).

The accelerator-magnet MSC panel drives the HDiv-VIM (radia.vim.soft_iron_from_mesh + rad.Solve, which
auto-routes a mesh-backed soft iron to the HDiv-VIM).  This panel intentionally does not expose the
mesh-less moment-yano backend.  The HDiv-VIM is KELVIN-less / iron-only, so the .vol must contain only the 'yoke' volume
material, and IMA symmetry is not supported there yet.

Locks:
  (1) on an iron-only yoke + coil, the hdiv backend magnetizes the iron (converged, M_avg large) -- the
      panel routes the coil field into the HDiv-VIM solve correctly;
  (2) demag_backend='hdiv' on a MULTI-material .vol (yoke + air) returns a clean error (iron-only);
  (3) demag_backend='hdiv' with an IMA symmetry string returns a clean error (not supported yet);
  (4) demag_backend='yano' returns a clean error (not a backend for this panel).

Also guards the 2026-06-17 coil-bug fix: solve_msc used to rad.UtiDelAll() AFTER building the coil
(destroying it); now the iron responds to the coil, so M_avg is large (~1e5 A/m), not ~0.

Slow: each solve runs the real HDiv-VIM solver on a small (64-hex) yoke.
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
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402

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
    """A 4^3-hex iron cube centred at origin, single material 'yoke', saved as a Netgen .vol."""
    mp = lambda x, y, z: (L * (x - 0.5), L * (y - 0.5), L * (z - 0.5))  # noqa: E731
    with ng.TaskManager():
        mesh = MakeStructured3DMesh(hexes=True, nx=4, ny=4, nz=4, mapping=mp)
    mesh.ngmesh.SetMaterial(1, "yoke")
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


def test_yano_backend_removed(tmp_path, coil_script):
    """demag_backend='yano' returns a clean error -- this panel is HDiv-VIM only."""
    from calc_accel_msc import solve_msc
    vol = _iron_only_yoke_vol(tmp_path / "yoke.vol")
    r = solve_msc(coil_script=coil_script, vol_file=vol, mat=_linear_mat(),
                  demag_backend="yano", solver=0, tol=1e-6, max_iter=400)
    assert "error" in r, r
    assert "yano" in r["error"].lower() and "panel" in r["error"].lower()


def test_hdiv_rejects_multimaterial(tmp_path, coil_script):
    """demag_backend='hdiv' on a yoke+air .vol -> clean error (the VIM is iron-only / KELVIN-less)."""
    from calc_accel_msc import solve_msc
    vol = _yoke_plus_air_vol(tmp_path / "yoke_air.vol")
    r = solve_msc(coil_script=coil_script, vol_file=vol, mat=_linear_mat(),
                  demag_backend="hdiv", solver=0, tol=1e-6, max_iter=400)
    assert "error" in r, r
    assert "iron-only" in r["error"].lower()


def test_hdiv_ima_runs(tmp_path, coil_script):
    """demag_backend='hdiv' with an IMA image string now RUNS via the image-charge Gram (the panel passes
    it through to the HDiv-VIM).  Uses a HALF yoke (z>0) + image='-z' -- a valid reduced model.  The IMA
    physics (reduced+IMA == full) itself is locked in tests/feec/test_hdiv_vim_ima.py."""
    from calc_accel_msc import solve_msc
    p = tmp_path / "half_yoke.vol"
    mp = lambda x, y, z: (L * (x - 0.5), L * (y - 0.5), 0.5 * L * z)   # [-L/2,L/2]^2 x [0,L/2] half  # noqa: E731
    with ng.TaskManager():
        m = MakeStructured3DMesh(hexes=True, nx=3, ny=3, nz=2, mapping=mp)
    m.ngmesh.SetMaterial(1, "yoke")
    m.ngmesh.Save(str(p))
    r = solve_msc(coil_script=coil_script, vol_file=str(p), mat=_linear_mat(),
                  demag_backend="hdiv", ima="-z", solver=0, tol=1e-7, max_iter=800)
    assert "error" not in r, r
    assert r["converged"] and r["ima"] == "-z"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
