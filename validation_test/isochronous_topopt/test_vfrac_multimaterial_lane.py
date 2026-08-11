"""Multi-material Sculpt x HDiv-MMM lane: partitioned volume fractions ->
conformal all-hex interface -> named `.vol` regions -> per-region solve.

Body: a concentric sphere (core r<=0.3 inside shell 0.3<r<=0.6) whose
UNION is the same r<=0.6 sphere as the single-material lane.  Two gates
follow from that construction:

* the interface is CONFORMAL -- the two regions share interface nodes, and
  the exported `.vol` carries the outer skin (exactly the topological
  skin) plus the material-interface faces;
* giving both regions the SAME permeability must reproduce the sphere's
  closed-form demagnetizing factor 1/3, which checks the partition, the
  interface, the label plumbing and the per-region assembly at once.

A contrasting pair (weak core inside a permeable shell) is then locked on
the per-region mean magnetization -- the demag factor is a GEOMETRIC
quantity of the whole body and cannot show that the per-region mu_r took
effect, so it is the wrong discriminator there (measured 2026-08-10:
identical to five digits across a 200x permeability contrast).

Requires a Cubit license (sculpt.exe + headless import/export); ~5 min.
"""
import json
from pathlib import Path

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netCDF4")
pytest.importorskip("scipy")

from netgen.occ import OCCGeometry, Pnt, Sphere  # noqa: E402
from ngsolve import Integrate, Mesh, TaskManager, VOL  # noqa: E402
import ngsolve as ngs  # noqa: E402

from radia.topopt_cad import (  # noqa: E402
    nodal_from_element_density, write_vfrac_exodus)

pytest.importorskip("radia_mcp.cubit.server")

R_CORE, R_BODY = 0.30, 0.60


@pytest.fixture(scope="module")
def artifacts(tmp_path_factory):
    root = tmp_path_factory.mktemp("vfrac_multimat")
    with TaskManager():
        fine = Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0))
                    .GenerateMesh(maxh=0.15))
        vols = np.asarray(Integrate(ngs.CoefficientFunction(1.0), fine,
                                    VOL, element_wise=True), float)
        centroid = [
            np.asarray(Integrate(axis, fine, VOL, element_wise=True),
                       float) / vols
            for axis in (ngs.x, ngs.y, ngs.z)]
        radius = np.sqrt(sum(c ** 2 for c in centroid))
        core = nodal_from_element_density(
            fine, (radius <= R_CORE).astype(float))
        shell = nodal_from_element_density(
            fine, ((radius > R_CORE) & (radius <= R_BODY)).astype(float))
    vf = write_vfrac_exodus(fine, {"core": core, "shell": shell},
                            root / "two", level=0.5, cells=24, supersample=4)

    from radia_mcp.cubit.server import cubit_vfrac_to_vol
    report = json.loads(cubit_vfrac_to_vol(
        vf["path"], out_vol=str(root / "two.vol"),
        out_msh=str(root / "two.msh"), material_names="core,shell"))
    return {"vf": vf, "mesh_report": report, "root": root}


def test_multimaterial_bridge_gates_and_labels(artifacts):
    report = artifacts["mesh_report"]
    assert report.get("status") == "ok", report.get("gates", report)
    assert report["gates"] == {"closure_ok": True,
                               "no_inverted_elements": True,
                               "boundary_faces_ok": True,
                               "per_material_closure_ok": True}
    assert [row["name"] for row in report["materials"]] == ["core", "shell"]
    # a multi-material mesh carries interface faces ON TOP of the skin;
    # the OUTER count is what must match the topological skin
    assert report["vol_outer_faces"] == report["msh_skin_faces"]
    assert report["vol_interface_faces"] > 0, (
        "no material-interface faces exported -- the regions are not "
        "sharing an interface")


def test_each_region_closes_against_its_own_volume_fraction(artifacts):
    """The TOTAL closure conserves a partition bias -- core over-meshed at
    the shell's expense sums to the same number -- yet the per-region
    ``mu_r`` solve is only as trustworthy as the per-region geometry.  So
    every block is gated against its OWN fraction integral.
    """
    report = artifacts["mesh_report"]
    per_material = report["per_material_closure"]
    assert set(per_material) == {"core", "shell"}
    for name, entry in per_material.items():
        assert entry["ok"] is True, (name, entry)
        assert entry["closure"] <= report["closure_tolerance"], (name, entry)
    # the per-region volumes must reconstitute the total the shared gate
    # used, otherwise the two closures are measuring different meshes
    total = sum(entry["mesh_volume"] for entry in per_material.values())
    assert total == pytest.approx(report["mesh_volume"], rel=1e-6)


def test_multimaterial_interface_is_conformal(artifacts):
    mesh = Mesh(str(artifacts["mesh_report"]["vol"]))
    assert sorted(set(mesh.GetMaterials())) == ["core", "shell"]
    names = mesh.GetMaterials()
    nodes = {"core": set(), "shell": set()}
    for el in mesh.Elements(VOL):
        target = nodes[names[el.index]]
        for vertex in el.vertices:
            target.add(vertex.nr)
    shared = nodes["core"] & nodes["shell"]
    assert len(shared) > 0, (
        "core and shell share no nodes -- the Sculpt interface is NOT "
        "conformal, and normal flux would not be continuous across it")


def test_multimaterial_same_mu_recovers_sphere_demag(artifacts):
    from radia_mcp.radia_ngsolve.server import hdiv_vim_demag_eval
    vol = Path(artifacts["mesh_report"]["vol"])
    result = json.loads(hdiv_vim_demag_eval(
        str(vol), mu_r="core=1000,shell=1000", h_ext="0,0,1e5"))
    assert result.get("status") == "ok", result
    assert sorted(result["material_volumes"]) == ["core", "shell"]
    # same body as the single-material lane -> same closed-form 1/3 band
    assert 0.30 < result["demag_factor"] < 0.37, result["demag_factor"]
    magnetization = result["mean_magnetization"]
    ratio = (magnetization["core"][2] / magnetization["shell"][2])
    # one material in two blocks: the regions magnetize alike (the shell
    # sits nearer the surface, so a modest spread is physical)
    assert 0.7 < ratio < 1.4, ratio


def test_per_region_mu_changes_the_regions(artifacts):
    from radia_mcp.radia_ngsolve.server import hdiv_vim_demag_eval
    vol = Path(artifacts["mesh_report"]["vol"])
    contrast = json.loads(hdiv_vim_demag_eval(
        str(vol), mu_r="core=5,shell=1000", h_ext="0,0,1e5"))
    assert contrast.get("status") == "ok", contrast
    magnetization = contrast["mean_magnetization"]
    ratio = magnetization["core"][2] / magnetization["shell"][2]
    # a 200x permeability contrast must leave the weak core far less
    # magnetized than the shell around it
    assert ratio < 0.5, ratio


def test_per_region_mu_must_cover_every_material(artifacts):
    from radia_mcp.radia_ngsolve.server import hdiv_vim_demag_eval
    vol = Path(artifacts["mesh_report"]["vol"])
    partial = json.loads(hdiv_vim_demag_eval(str(vol), mu_r="core=1000"))
    assert partial["status"] == "error"
    assert "shell" in partial["error"]
    unknown = json.loads(hdiv_vim_demag_eval(
        str(vol), mu_r="core=1000,shell=50,ghost=7"))
    assert unknown["status"] == "error"
    assert "ghost" in unknown["error"]
