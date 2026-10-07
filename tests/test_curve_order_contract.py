"""Pin the NGSolve geometry-order behaviour behind radia.mesh_curve.

Each scenario runs in a fresh interpreter because Netgen keeps a
process-global geometry: a CAD-less mesh loaded after another
geometry-bearing mesh is curved against that other geometry.
See docs/ngsolve_integration/curve_order.md.
"""
import json
import math
import pathlib
import subprocess
import sys
import textwrap

import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")

R = 0.01
AREA = 4.0 * math.pi * R * R
# Load the pure-Python helper from this checkout without importing the
# native radia package (the checkout may be unbuilt).
HELPER = pathlib.Path(__file__).resolve().parents[1] / "src" / "radia" / "mesh_curve.py"

PRELUDE = textwrap.dedent(f"""
    import json, math, sys
    import ngsolve as ng
    from netgen.occ import Sphere, Pnt, OCCGeometry
    R = {R!r}
    AREA = 4*math.pi*R*R
    def rel(m):
        return ng.Integrate(ng.CoefficientFunction(1), m, ng.BND)/AREA - 1
    def curve(m, p):
        with ng.TaskManager():
            m.Curve(p)
    def helper():
        import importlib.util
        spec = importlib.util.spec_from_file_location("mesh_curve", {str(HELPER)!r})
        mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
        return mod
""")


def _run(tmp_path, body):
    code = PRELUDE + textwrap.dedent(body)
    proc = subprocess.run([sys.executable, "-c", code], cwd=tmp_path,
                          capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, proc.stderr[-2000:]
    return json.loads(proc.stdout.strip().splitlines()[-1])


def _write_meshes(tmp_path):
    return _run(tmp_path, """
        geo = OCCGeometry(Sphere(Pnt(0, 0, 0), R))
        m = ng.Mesh(geo.GenerateMesh(maxh=R/2)); curve(m, 3); m.ngmesh.Save("c3.vol")
        flat = geo.GenerateMesh(maxh=R/2); flat.Save("flat.vol")
        mm = OCCGeometry(Sphere(Pnt(0, 0, 0), R*1e3)).GenerateMesh(maxh=R*1e3/2)
        mm.Scale(1e-3); mm.Save("mm_scaled.vol")
        lines = open("c3.vol").read().splitlines()
        open("c3_nocad.vol", "w").write("\\n".join(lines[:lines.index("endmesh")+1]) + "\\n")
        print(json.dumps({"ok": True}))
    """)


@pytest.fixture(scope="module")
def meshes(tmp_path_factory):
    path = tmp_path_factory.mktemp("curve_order")
    _write_meshes(path)
    return path


def test_loading_keeps_stored_order(meshes):
    out = _run(meshes, """
        m = ng.Mesh("c3.vol")
        print(json.dumps({"order": m.GetCurveOrder(), "rel": rel(m)}))
    """)
    assert out["order"] == 3
    assert abs(out["rel"]) < 1e-3


def test_embedded_cad_recurves_correctly(meshes):
    out = _run(meshes, """
        m = ng.Mesh("flat.vol"); curve(m, 2)
        print(json.dumps({"rel": rel(m)}))
    """)
    assert abs(out["rel"]) < 5e-3


def test_cadless_curve_flattens_silently(meshes):
    out = _run(meshes, """
        m = ng.Mesh("c3_nocad.vol"); curve(m, 3)
        print(json.dumps({"order": m.GetCurveOrder(), "rel": rel(m)}))
    """)
    assert out["order"] == 3          # reported order ...
    assert out["rel"] < -0.02         # ... but the surface is the flat polytope


def test_cadless_curve_borrows_process_geometry(meshes):
    out = _run(meshes, """
        ng.Mesh("mm_scaled.vol")
        m = ng.Mesh("c3_nocad.vol"); curve(m, 2)
        print(json.dumps({"rel": rel(m)}))
    """)
    assert out["rel"] > 1e3


def test_rescaled_mesh_curve_explodes(meshes):
    out = _run(meshes, """
        m = ng.Mesh("mm_scaled.vol"); curve(m, 2)
        print(json.dumps({"rel": rel(m)}))
    """)
    assert out["rel"] > 1e3


def test_ensure_curve_order_guards(meshes):
    out = _run(meshes, """
        mc = helper(); CurveOrderError, ensure_curve_order = mc.CurveOrderError, mc.ensure_curve_order
        res = {}
        res["kept"] = ensure_curve_order(ng.Mesh("c3.vol"), 2)["action"]
        r = ensure_curve_order(ng.Mesh("flat.vol"), 2)
        res["curved"] = r["action"]; res["d_area"] = r["area_rel_change"]
        for name in ("mm_scaled.vol",):
            try:
                ensure_curve_order(ng.Mesh(name), 2); res[name] = "no error"
            except CurveOrderError as e:
                res[name] = "extent" if "extent" in str(e) else str(e)
        print(json.dumps(res))
    """)
    assert out["kept"] == "kept"
    assert out["curved"] == "curved" and 0 < out["d_area"] < 0.5
    assert out["mm_scaled.vol"] == "extent"


def test_ensure_curve_order_rejects_cadless(meshes):
    out = _run(meshes, """
        mc = helper(); CurveOrderError, ensure_curve_order = mc.CurveOrderError, mc.ensure_curve_order
        lines = open("flat.vol").read().splitlines()
        open("flat_nocad.vol", "w").write("\\n".join(lines[:lines.index("endmesh")+1]) + "\\n")
        try:
            ensure_curve_order(ng.Mesh("flat_nocad.vol"), 2); res = "no error"
        except CurveOrderError as e:
            res = "no CAD" if "no CAD" in str(e) else str(e)
        print(json.dumps({"res": res}))
    """)
    assert out["res"] == "no CAD"
