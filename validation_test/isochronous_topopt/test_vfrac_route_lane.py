"""Sculpt volume-fraction route lane: the full Sculpt x HDiv-MMM combination
on a body with a closed-form answer.

density -> ``radia.topopt_cad.write_vfrac_exodus`` (Sculpt-native volume
fractions; no STL) -> ``cubit_vfrac_to_vol`` (standalone ``sculpt.exe
--input_vfrac`` -> Exodus -> headless Cubit -> gated ``.vol``) ->
``hdiv_vim_demag_eval`` (air-mesh-free charge-Gram demag solve).

The body is a sphere, whose demagnetizing factor is exactly 1/3 -- the one
number that verifies the WHOLE chain at once: a broken export reads 0, a
volume-biased regeneration shifts it, and a broken solver misses the band.
Measured on this lane's configuration (2026-08-10): closure 2e-3 band,
demag 0.333 +- band.  Requires a Cubit license (sculpt.exe + headless
import/export); runs in ~2 minutes.
"""
import json
from pathlib import Path

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netCDF4")
pytest.importorskip("scipy")

from netgen.occ import OCCGeometry, Pnt, Sphere  # noqa: E402
from ngsolve import Integrate, Mesh, TaskManager  # noqa: E402

from radia.topopt_cad import (  # noqa: E402
    nodal_from_element_density, write_vfrac_exodus)

pytest.importorskip("radia_mcp.cubit.server")


@pytest.fixture(scope="module")
def artifacts(tmp_path_factory):
    root = tmp_path_factory.mktemp("vfrac_lane")
    import ngsolve as ngs
    with TaskManager():
        fine = Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0))
                    .GenerateMesh(maxh=0.15))
        vols = np.asarray(Integrate(ngs.CoefficientFunction(1.0), fine,
                                    ngs.VOL, element_wise=True), float)
        cx = np.asarray(Integrate(ngs.x, fine, ngs.VOL,
                                  element_wise=True), float) / vols
        cy = np.asarray(Integrate(ngs.y, fine, ngs.VOL,
                                  element_wise=True), float) / vols
        cz = np.asarray(Integrate(ngs.z, fine, ngs.VOL,
                                  element_wise=True), float) / vols
        rho = (np.sqrt(cx ** 2 + cy ** 2 + cz ** 2) <= 0.6).astype(float)
        nodal = nodal_from_element_density(fine, rho)
    vf = write_vfrac_exodus(fine, nodal, root / "sphere_vf", level=0.5,
                            cells=24, supersample=4)

    from radia_mcp.cubit.server import cubit_vfrac_to_vol
    report = json.loads(cubit_vfrac_to_vol(
        vf["path"], out_vol=str(root / "sphere.hex.vol"),
        out_msh=str(root / "sphere.hex.msh")))
    return {"vf": vf, "mesh_report": report, "root": root}


def test_vfrac_bridge_gates(artifacts):
    report = artifacts["mesh_report"]
    assert report.get("status") == "ok", report.get("gates", report)
    assert report["gates"] == {"closure_ok": True,
                               "no_inverted_elements": True,
                               "boundary_faces_ok": True}
    # interface reconstruction honors the fraction integral (measured
    # 2e-3 on this configuration; the STL route's gate is 3e-2)
    assert report["closure"] < 1e-2, report["closure"]


def test_vfrac_sphere_demag_third(artifacts):
    from radia_mcp.radia_ngsolve.server import hdiv_vim_demag_eval
    vol = Path(artifacts["mesh_report"]["vol"])
    assert vol.is_file()
    result = json.loads(hdiv_vim_demag_eval(str(vol), mu_r=1000.0,
                                            h_ext="0,0,1e5"))
    assert result.get("status") == "ok", result
    assert result["ne"] > 100
    demag = result["demag_factor"]
    # sphere: exactly 1/3; band covers the P0->P1 body blur + the hex
    # staircase-free reconstruction at 24 cells across the bbox
    assert 0.30 < demag < 0.37, demag
