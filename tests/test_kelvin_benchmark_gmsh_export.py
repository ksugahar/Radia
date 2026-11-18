"""Gmsh B/H export of the Kelvin magnetic-sphere benchmark.

Guards: the exported B is the total physical field. In air that is the
reduced potential corrected by the source lift, plus the applied field.
The Kelvin (mapped) region is excluded, each element carries its own
interface trace (ElementNodeData), and export errors propagate.
Assumes: the analytical sphere (mu_r=100, H0=1 A/m along z, a=0.05 m) on the
coarse Kelvin mesh used by validation_test/panels/test_kelvin_benchmark_direct.
The bounds are loose code-breakage guards, not accuracy criteria: the
reduced-field-only export gives a far-field Bz near 0 and the unlifted air
field a 46 % median normal-B jump. Neither is near the bounds.
"""
import math

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
gmsh = pytest.importorskip("gmsh")

MU0 = 4e-7 * math.pi
A = 0.05


def _sphere_vol(tmp_path):
    from netgen.occ import (Glue, IdentificationType, OCCGeometry, Pnt,
                            Sphere, Vertex)
    iron = Sphere(Pnt(0, 0, 0), A)
    iron.mat("magnetic")
    iron.maxh = .015
    iron.faces.name = "sphere"
    inner = Sphere(Pnt(0, 0, 0), .2)
    inner.faces.name = "kelvin_int"
    air = inner - iron
    air.mat("air")
    outer = Sphere(Pnt(0, 0, .6), .2)
    outer.mat("kelvin")
    outer.faces.name = "kelvin_ext"
    ground = Vertex(Pnt(0, 0, .6))
    ground.name = "GND"
    shape = Glue([air, iron, outer, ground])
    inside = [f for f in shape.faces if f.name == "kelvin_int"][0]
    outside = [f for f in shape.faces if f.name == "kelvin_ext"][0]
    inside.Identify(outside, "kelvin_periodic", IdentificationType.PERIODIC)
    with ng.TaskManager():
        mesh = OCCGeometry(shape).GenerateMesh(maxh=.035)
    path = tmp_path / "sphere.vol"
    mesh.Save(str(path))
    return path


def _read_b(path):
    gmsh.initialize(["gmsh", "-noconfig"])
    try:
        gmsh.open(str(path))
        material_of = {}
        for dim, ptag in gmsh.model.getPhysicalGroups(3):
            name = gmsh.model.getPhysicalName(dim, ptag)
            for ent in gmsh.model.getEntitiesForPhysicalGroup(dim, ptag):
                for tags in gmsh.model.mesh.getElements(dim, ent)[1]:
                    material_of.update({int(t): name for t in tags})
        node_tags, coords, _ = gmsh.model.mesh.getNodes()
        xyz = dict(zip(node_tags, np.asarray(coords).reshape(-1, 3)))
        views = {gmsh.view.option.getString(t, "Name"): t
                 for t in gmsh.view.getTags()}
        kind, elems, data, _t, ncomp = gmsh.view.getModelData(views["B"], 0)
        samples = []
        for elem, values in zip(elems, data):
            conn = gmsh.model.mesh.getElement(int(elem))[1]
            for node, b in zip(conn, np.asarray(values).reshape(-1, 3)):
                samples.append((material_of[int(elem)], int(node),
                                xyz[node], b))
        return kind, ncomp, sorted(views), samples
    finally:
        gmsh.finalize()


def test_kelvin_benchmark_exports_total_b_with_interface_traces(tmp_path):
    from radia.panels.calc_kelvin_benchmark import solve_kelvin_benchmark

    out = tmp_path / "sphere.msh"
    result = solve_kelvin_benchmark(str(_sphere_vol(tmp_path)), fes_order=2,
                                    output_msh=str(out))
    assert result["msh_file"] == str(out)

    kind, ncomp, names, samples = _read_b(out)
    assert (kind, ncomp, names) == ("ElementNodeData", 3, ["B", "H"])
    mats = {m for m, *_ in samples}
    assert mats == {"magnetic", "air"}, f"exported materials {mats}"

    far = [b[2] for m, _n, p, b in samples
           if m == "air" and np.linalg.norm(p) > 0.15]
    far_bz = float(np.mean(far))
    assert far_bz == pytest.approx(MU0, rel=0.05), (
        f"far-field Bz {far_bz:.4e} T, expected about mu0*H0 = {MU0:.4e} T")

    sides = {}
    for m, n, p, b in samples:
        if abs(np.linalg.norm(p) - A) < 1e-6:
            sides.setdefault(n, {}).setdefault(m, []).append(b)
    jumps = []
    for n, side in sides.items():
        if len(side) != 2:
            continue
        p = next(pp for m, nn, pp, _b in samples if nn == n)
        r_hat = p / np.linalg.norm(p)
        b_mag = np.mean(side["magnetic"], axis=0)
        b_air = np.mean(side["air"], axis=0)
        jumps.append(abs((b_mag - b_air) @ r_hat) / np.linalg.norm(b_mag))
    assert len(jumps) > 50, f"only {len(jumps)} two-sided interface nodes"
    median_jump = float(np.median(jumps))
    assert median_jump < 0.1, (
        f"median normal-B jump {median_jump:.3f} across the sphere, "
        "expected near 0 (normal B is continuous)")


def test_kelvin_benchmark_export_failure_propagates(tmp_path, monkeypatch):
    import radia.gmsh_post_export as gpe
    from radia.panels.calc_kelvin_benchmark import solve_kelvin_benchmark

    def boom(self, filename, time=0.0, timestep=0):
        raise RuntimeError("export failed")

    monkeypatch.setattr(gpe.GmshPostExport, "write", boom)
    with pytest.raises(RuntimeError, match="export failed"):
        solve_kelvin_benchmark(str(_sphere_vol(tmp_path)), fes_order=2,
                               output_msh=str(tmp_path / "sphere.msh"))
