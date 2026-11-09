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
        mc = helper(); CurveOrderError, ensure = mc.CurveOrderError, mc.ensure_curve_order
        res = {}
        with ng.TaskManager():
            res["kept"] = ensure(ng.Mesh("c3.vol"), 2, vol_path="c3.vol")["action"]
            r = ensure(ng.Mesh("flat.vol"), 2, vol_path="flat.vol")
            res["curved"] = r["action"]; res["ratio"] = r["boundary_ratio"]
            try:
                ensure(ng.Mesh("mm_scaled.vol"), 2, vol_path="mm_scaled.vol"); res["mm"] = "no error"
            except CurveOrderError as e:
                res["mm"] = ("box" if "box" in str(e) else str(e), e.mesh_modified)
        print(json.dumps(res))
    """)
    assert out["kept"] == "kept"
    assert out["curved"] == "curved" and 1.0 < out["ratio"] < 1.1
    assert out["mm"] == ["box", False]


def test_ensure_curve_order_rejects_cadless_even_with_borrowed_geometry(meshes):
    # A CAD-less file loaded after a geometry-bearing one sees that geometry
    # through GetGeometry(); provenance must come from the file itself.
    out = _run(meshes, """
        mc = helper(); CurveOrderError, ensure = mc.CurveOrderError, mc.ensure_curve_order
        lines = open("flat.vol").read().splitlines()
        open("flat_nocad.vol", "w").write("\\n".join(lines[:lines.index("endmesh")+1]) + "\\n")
        ng.Mesh("c3.vol")
        m = ng.Mesh("flat_nocad.vol")
        borrowed = getattr(m.ngmesh.GetGeometry(), "shape", None) is not None
        try:
            with ng.TaskManager():
                ensure(m, 2, vol_path="flat_nocad.vol")
            res = "no error"
        except CurveOrderError as e:
            res = ("no CAD" if "no CAD" in str(e) else str(e), e.mesh_modified)
        print(json.dumps({"borrowed": borrowed, "res": res, "order": m.GetCurveOrder()}))
    """)
    assert out["borrowed"]
    assert out["res"] == ["no CAD", False]
    assert out["order"] == 1


def test_embedded_cad_stays_with_its_mesh(meshes):
    out = _run(meshes, """
        a = ng.Mesh("flat.vol")
        ng.Mesh("mm_scaled.vol")          # load a different CAD afterwards
        mc = helper()
        with ng.TaskManager():
            mc.ensure_curve_order(a, 2, vol_path="flat.vol")
        print(json.dumps({"rel": rel(a)}))
    """)
    assert abs(out["rel"]) < 5e-3


def test_ensure_curve_order_accepts_coarse_correct_mesh(tmp_path):
    out = _run(tmp_path, """
        geo = OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0))
        geo.GenerateMesh(maxh=2.0, curvaturesafety=1).Save("coarse.vol")
        m = ng.Mesh("coarse.vol")
        with ng.TaskManager():
            r = helper().ensure_curve_order(m, 2, vol_path="coarse.vol")
        print(json.dumps({"action": r["action"], "domain_ratio": r["domain_ratio"]}))
    """)
    assert out["action"] == "curved"
    assert 1.0 < out["domain_ratio"] < 3.0


def test_trailing_text_is_not_a_geometry_archive(meshes):
    out = _run(meshes, """
        mc = helper(); CurveOrderError, ensure = mc.CurveOrderError, mc.ensure_curve_order
        lines = open("flat.vol").read().splitlines()
        body = lines[:lines.index("endmesh")+1]
        open("flat_comment.vol", "w").write("\\n".join(body + ["# harmless trailing comment"]) + "\\n")
        ng.Mesh("c3.vol")
        m = ng.Mesh("flat_comment.vol")
        try:
            with ng.TaskManager():
                ensure(m, 2, vol_path="flat_comment.vol")
            res = "no error"
        except CurveOrderError as e:
            res = ("not supported" if "not supported" in str(e) else str(e), e.mesh_modified)
        print(json.dumps({"kind": mc.vol_geometry_kind("flat_comment.vol"), "res": res,
                          "order": m.GetCurveOrder()}))
    """)
    assert out["kind"] == "other"
    assert out["res"] == ["not supported", False]
    assert out["order"] == 1


def test_load_vol_keeps_only_a_restored_occ_geometry(tmp_path):
    """Guards: radia.mesh_curve.load_vol keeps an attached geometry only for a restored OCC archive.

    Assumes: Netgen binds a .vol whose tail it does not restore to the last geometry of the
    process, here a CSG one, also when the tail looks like a CSG archive; a CSG restore cannot
    be told apart from that fallback, so a CSG geometry is never kept.
    """
    out = _run(tmp_path, """
        import gzip, shutil
        from netgen.csg import CSGeometry, Sphere as CSphere, Pnt as CPnt
        from ngsolve.meshes import MakeStructured3DMesh
        g = CSGeometry(); g.Add(CSphere(CPnt(0, 0, 0), 1.0)); g.GenerateMesh(maxh=0.5).Save("sphere_csg.vol")
        MakeStructured3DMesh(hexes=False, nx=2, ny=2, nz=2).ngmesh.Save("cube.vol")
        cube = open("cube.vol").read()
        open("cube_csgtail.vol", "w").write(cube + "csgsurfaces-extra\\n")
        open("cube_csgbare.vol", "w").write(cube + "csgsurfaces\\n")
        OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=R/2).Save("occ.vol")
        for plain in ("cube.vol", "sphere_csg.vol"):
            with open(plain, "rb") as src, gzip.open(plain + ".gz", "wb") as dst:
                shutil.copyfileobj(src, dst)
        mc = helper(); res = {}
        seed = ng.Mesh("sphere_csg.vol")
        plain = ng.Mesh("cube_csgtail.vol").ngmesh.GetGeometry()
        res["premise"] = [mc.vol_geometry_kind("cube_csgtail.vol"), type(plain).__name__,
                          plain is seed.ngmesh.GetGeometry()]
        one = ng.CoefficientFunction(1)
        for name in ("cube.vol", "cube.vol.gz", "cube_csgtail.vol", "cube_csgbare.vol",
                     "sphere_csg.vol", "sphere_csg.vol.gz"):
            raw = ng.Mesh(name); m = mc.load_vol(name)
            res[name] = [type(m.ngmesh.GetGeometry()).__name__, m.GetCurveOrder(),
                         ng.Integrate(one, m) == ng.Integrate(one, raw)]
        m = mc.load_vol("occ.vol")
        res["occ"] = [type(m.ngmesh.GetGeometry()).__name__, m.GetCurveOrder()]
        print(json.dumps(res))
    """)
    assert out["premise"] == ["csg", "CSGeometry", True], out["premise"]
    for name in ("cube.vol", "cube.vol.gz", "cube_csgtail.vol", "cube_csgbare.vol",
                 "sphere_csg.vol", "sphere_csg.vol.gz"):
        assert out[name] == ["NetgenGeometry", 1, True], (name, out[name])
    assert out["occ"] == ["OCCGeometry", 1], out["occ"]


def test_radia_routes_never_curve_and_straight_order_is_explicit(tmp_path):
    """Guards: radia.vim._vim._require_curve_order; straight_curve_order detaches a borrowed CAD.

    Assumes: Curve of a mesh without CAD keeps the straight map and only raises the reported
    order, and Curve onto a borrowed geometry is wrong or crashes, so no Radia route curves.
    """
    pytest.importorskip("radia._radia_pybind")
    out = _run(tmp_path, """
        from netgen.csg import CSGeometry, OrthoBrick, Pnt as CPnt
        from ngsolve.meshes import MakeStructured3DMesh
        from radia.mesh_curve import CurveOrderError, load_vol, straight_curve_order
        from radia.vim._vim import _require_curve_order
        res = {}
        MakeStructured3DMesh(hexes=False, nx=2, ny=2, nz=2).ngmesh.Save("cube.vol")
        q3 = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=R/2)); curve(q3, 3)
        for name, m in (("structured", MakeStructured3DMesh(hexes=True, nx=1, ny=1, nz=1)),
                        ("cadfree_vol", load_vol("cube.vol")),
                        ("occ_uncurved", ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=R/2)))):
            try:
                _require_curve_order(m, 2, "probe"); res[name] = ["accepted", m.GetCurveOrder()]
            except CurveOrderError as error:
                res[name] = ["refused", m.GetCurveOrder(), "straight_curve_order" in str(error)]
        res["q3"] = [_require_curve_order(q3, 2, "probe") is q3, q3.GetCurveOrder()]
        g = CSGeometry(); g.Add(OrthoBrick(CPnt(-1, -1, -1), CPnt(2, 2, 2))); g.GenerateMesh(maxh=1.0).Save("brick.vol")
        brick = ng.Mesh("brick.vol")
        m = ng.Mesh("cube.vol")
        res["borrowed"] = type(m.ngmesh.GetGeometry()).__name__
        with ng.TaskManager():
            straight_curve_order(m, 2)
        res["straight"] = [m.GetCurveOrder(), type(m.ngmesh.GetGeometry()).__name__,
                           ng.Integrate(ng.CoefficientFunction(1), m)]
        print(json.dumps(res))
    """)
    for name in ("structured", "cadfree_vol", "occ_uncurved"):
        assert out[name] == ["refused", 1, True], (name, out[name])
    assert out["q3"] == [True, 3], out["q3"]
    assert out["borrowed"] == "CSGeometry", out
    order, geometry, volume = out["straight"]
    assert order == 2 and geometry == "NetgenGeometry" and abs(volume - 1.0) < 1e-12, out["straight"]


def test_public_entries_never_curve_a_borrowed_geometry(tmp_path):
    """Guards: vim.Solve / MagnetizationSource refuse an uncurved mesh; the path seam refuses too.

    Assumes: after a CSG GenerateMesh, a plain Mesh() of a CAD-free .vol (with or without
    trailing text or a fake CSG tail, gzipped or not) is bound to that CSG geometry, and Curve
    onto it is wrong or crashes; a stored order and an explicitly curved CAD are kept.
    """
    pytest.importorskip("radia._radia_pybind")
    out = _run(tmp_path, """
        import gzip, shutil
        from netgen.csg import CSGeometry, Sphere as CSphere, Pnt as CPnt
        from ngsolve.meshes import MakeStructured3DMesh
        from radia import vim
        from radia.mesh_curve import CurveOrderError, ensure_curve_order, load_vol, vol_geometry_kind
        from radia.vim._vim import _require_curve_order
        occ = OCCGeometry(Sphere(Pnt(0, 0, 0), R))
        occ.GenerateMesh(maxh=R/2).Save("occ_flat.vol")
        q2 = ng.Mesh(occ.GenerateMesh(maxh=R/2)); curve(q2, 2); q2.ngmesh.Save("q2.vol")
        lines = open("q2.vol").read().splitlines()
        open("q2_nocad.vol", "w").write("\\n".join(lines[:lines.index("endmesh")+1]) + "\\n")
        MakeStructured3DMesh(hexes=False, nx=2, ny=2, nz=2).ngmesh.Save("cube.vol")
        cube = open("cube.vol").read()
        open("cube_comment.vol", "w").write(cube + "# written by another tool\\n")
        open("cube_csgtail.vol", "w").write(cube + "csgsurfaces-extra\\n")
        for plain in ("cube.vol", "cube_comment.vol"):
            with open(plain, "rb") as src, gzip.open(plain + ".gz", "wb") as dst:
                shutil.copyfileobj(src, dst)
        g = CSGeometry(); g.Add(CSphere(CPnt(0, 0, 0), 1.0)); keep = g.GenerateMesh(maxh=0.5)
        H = ng.CoefficientFunction((0.0, 0.0, 1.0e3))
        res = {"kind_comment": vol_geometry_kind("cube_comment.vol")}
        for entry in ("solve", "source"):
            m = ng.Mesh("cube.vol")
            res["borrowed"] = type(m.ngmesh.GetGeometry()).__name__
            try:
                with ng.TaskManager():
                    if entry == "solve":
                        vim.Solve(m, mu_r=100.0, H_ext=H, order=1, curve_order=2)
                    else:
                        vim.MagnetizationSource(m, (0.0, 0.0, 1.0e5), order=1, curve_order=2)
                res[entry] = ["ran", m.GetCurveOrder()]
            except CurveOrderError as error:
                res[entry] = ["refused", m.GetCurveOrder(), "straight_curve_order" in str(error)]
        for name in ("cube.vol", "cube_comment.vol", "cube.vol.gz", "cube_comment.vol.gz",
                     "cube_csgtail.vol", "occ_flat.vol"):
            m = load_vol(name)
            res[name] = [type(m.ngmesh.GetGeometry()).__name__]
            try:
                _require_curve_order(m, 2, "HDiv device prepare"); res[name].append("accepted")
            except CurveOrderError:
                res[name].append("refused")
        one = ng.CoefficientFunction(1.0)
        stored = load_vol("q2_nocad.vol")
        before = [ng.Integrate(one, stored, ng.BND), ng.Integrate(one, stored)]
        with ng.TaskManager():
            vim.Solve(stored, mu_r=100.0, H_ext=H, order=1, curve_order=2)
        after = [ng.Integrate(one, stored, ng.BND), ng.Integrate(one, stored)]
        res["stored"] = [stored.GetCurveOrder(), max(abs(a/b - 1) for a, b in zip(after, before)), rel(stored)]
        cad = load_vol("occ_flat.vol")
        with ng.TaskManager():
            action = ensure_curve_order(cad, 2, vol_path="occ_flat.vol")["action"]
            vim.Solve(cad, mu_r=100.0, H_ext=H, order=1, curve_order=2)
        res["cad"] = [action, cad.GetCurveOrder(), rel(cad)]
        print(json.dumps(res))
    """)
    assert out["borrowed"] == "CSGeometry" and out["kind_comment"] == "other", out
    assert out["solve"] == ["refused", 1, True] and out["source"] == ["refused", 1, True], out
    for name in ("cube.vol", "cube_comment.vol", "cube.vol.gz", "cube_comment.vol.gz", "cube_csgtail.vol"):
        assert out[name] == ["NetgenGeometry", "refused"], (name, out[name])
    assert out["occ_flat.vol"] == ["OCCGeometry", "refused"], out["occ_flat.vol"]
    order, change, area = out["stored"]
    assert order == 2 and change == 0.0 and abs(area) < 5e-3, out["stored"]
    action, order, area = out["cad"]
    assert action == "curved" and order == 2 and abs(area) < 5e-3, out["cad"]


def test_production_vol_entries_load_through_load_vol():
    """Guards: the HDiv device prepare and soft_iron_from_vol path loads.

    Assumes: a direct Mesh(path) in these modules would reattach a borrowed geometry.
    """
    import ast
    root = pathlib.Path(__file__).resolve().parents[1] / "src" / "radia"
    for relative in ("hdiv_gpu.py", "vim/_radsolve.py"):
        tree = ast.parse((root / relative).read_text(encoding="utf-8"))
        direct = [node.lineno for node in ast.walk(tree)
                  if isinstance(node, ast.Call) and getattr(node.func, "attr", getattr(node.func, "id", None)) == "Mesh"
                  and node.args and isinstance(node.args[0], ast.Call) and getattr(node.args[0].func, "id", None) == "str"]
        loads = [node for node in ast.walk(tree)
                 if isinstance(node, ast.Call) and getattr(node.func, "id", None) == "load_vol"]
        assert direct == [] and loads, (relative, direct, len(loads))


def test_production_curved_entries_require_the_stored_order():
    """Guards: no curved Radia entry calls Curve; each requires the stored geometry order.

    Assumes: a mesh object or a loaded .vol cannot prove where its geometry came from.
    """
    import ast
    root = pathlib.Path(__file__).resolve().parents[1] / "src" / "radia"
    for relative in ("hdiv_gpu.py", "vim/_solve.py", "vim/_solver.py", "vim/_magnetization_source.py"):
        tree = ast.parse((root / relative).read_text(encoding="utf-8"))
        calls = [node for node in ast.walk(tree) if isinstance(node, ast.Call)]
        names = [getattr(node.func, "attr", getattr(node.func, "id", None)) for node in calls]
        curving = [node.lineno for node, name in zip(calls, names) if name in ("Curve", "_curve_mesh")]
        assert curving == [] and "_require_curve_order" in names, (relative, curving)


@pytest.mark.parametrize("factor", [0.8, 1.2, 1e-3])
def test_mesh_transformed_after_loading_is_rejected(meshes, factor):
    out = _run(meshes, f"""
        mc = helper(); CurveOrderError, ensure = mc.CurveOrderError, mc.ensure_curve_order
        m = ng.Mesh("flat.vol"); m.ngmesh.Scale({factor!r}); m = ng.Mesh(m.ngmesh)
        try:
            with ng.TaskManager():
                ensure(m, 2, vol_path="flat.vol")
            res = "no error"
        except CurveOrderError as e:
            res = ("moved" if "moved or rescaled" in str(e) else str(e), e.mesh_modified)
        print(json.dumps({{"res": res, "order": m.GetCurveOrder()}}))
    """)
    assert out["res"] == ["moved", False]
    assert out["order"] == 1

