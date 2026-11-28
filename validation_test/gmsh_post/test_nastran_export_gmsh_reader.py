"""radia.nastran_export read back by Gmsh's own Nastran reader.

Guards: the second-order BDF written from a curved .vol keeps the curved
geometry, with positive Jacobians, the right edge-node placement and the
Nastran node order of CTETRA/CHEXA/CPENTA/CTRIA6/CQUAD8.
Assumes: Gmsh 5 reads BDF with its own node tables, independent of the writer.
NGSolve's integral of 1 over the curved mesh is the reference volume.

Two limits of the oracle, handled explicitly:
- Gmsh reads a CPYRAM as a 5-node pyramid, so pyramid mid-edge nodes are
  compared through a 5-node truncation; their card order follows the
  Nastran spec only.
- Gmsh expects a continuation after a full 8-field CPENTA/CTRIA6 line.  The
  BDF needs none; the copy handed to Gmsh gets an empty one.
"""
import gzip
from pathlib import Path

import gmsh
import numpy as np
import ngsolve as ng
import pytest
from netgen.occ import Box, Cylinder, Glue, OCCGeometry, Pnt, WorkPlane, Z

from radia.gmsh_post_export import GmshPostExport
from radia.nastran_export import write_nastran_bdf

pytestmark = pytest.mark.usefixtures("ngsolve_taskmanager")

ROOT = Path(__file__).resolve().parents[2]
MIXED_O2 = ROOT / "tests" / "fixtures" / "cubit" / "mixed_cylinder_o2.vol.gz"
TOL = 1e-9  # GRID* prints 10 significant digits


def _gmsh_copy(path: Path) -> Path:
    out = []
    for line in path.read_text().splitlines():
        if line.startswith(("CPENTA ", "CTRIA6 ")) and len(line) == 72:
            out += [line + "+", "+"]
        else:
            out.append(line)
    copy = path.with_name(path.stem + ".gmsh.bdf")
    copy.write_text("\n".join(out) + "\n")
    return copy


def _gmsh_measure(path: Path, dim: int) -> dict:
    """Total measure from Gmsh's Jacobians, element counts, smallest det."""
    gmsh.initialize(["-nt", "1"])
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        gmsh.open(str(_gmsh_copy(path) if path.suffix == ".bdf" else path))
        total, types, min_det = 0.0, {}, np.inf
        for t in gmsh.model.mesh.getElementTypes(dim=dim):
            loc, w = gmsh.model.mesh.getIntegrationPoints(t, "Gauss8")
            _, det, _ = gmsh.model.mesh.getJacobians(t, loc)
            det = np.asarray(det).reshape(-1, len(w))
            types[gmsh.model.mesh.getElementProperties(t)[0]] = det.shape[0]
            total += float((det * np.asarray(w)).sum())
            min_det = min(min_det, float(det.min()))
        return {"measure": total, "types": types, "min_det": min_det}
    finally:
        gmsh.finalize()


def _netgen_two_materials():
    geo = OCCGeometry(Glue([
        Cylinder(Pnt(0, 0, 0), Z, r=0.5, h=1).mat("iron"),
        Box(Pnt(0.6, -0.3, 0), Pnt(1.2, 0.3, 1)).mat("air")]))
    mesh = ng.Mesh(geo.GenerateMesh(maxh=0.25))
    mesh.Curve(2)
    return mesh


def test_curved_tets_keep_volume_and_boundary_area(tmp_path):
    mesh = _netgen_two_materials()
    out = tmp_path / "cyl_box.bdf"
    info = write_nastran_bdf(mesh, out, order=2)
    n_bnd = len(list(mesh.Elements(ng.BND)))
    assert info["elements"] == {"CTETRA": mesh.ne, "CTRIA6": n_bnd}, info
    # Netgen lists tet corners negatively: every element must be flipped.
    assert info["flipped_to_positive"] == mesh.ne, info
    vol = _gmsh_measure(out, 3)
    area = _gmsh_measure(out, 2)
    ref_vol = ng.Integrate(1, mesh, order=8)
    ref_area = ng.Integrate(1, mesh, ng.BND, order=8)
    assert vol["types"] == {"Tetrahedron 10": mesh.ne}, vol
    assert vol["min_det"] > 0, vol
    assert abs(vol["measure"] / ref_vol - 1) < TOL, (vol["measure"], ref_vol)
    assert abs(area["measure"] / ref_area - 1) < TOL, (area["measure"], ref_area)
    # The straight order-1 file is measurably different: the check can fail.
    write_nastran_bdf(mesh, tmp_path / "linear.bdf", order=1)
    lin = _gmsh_measure(tmp_path / "linear.bdf", 3)["measure"]
    assert abs(lin / ref_vol - 1) > 1e-3, (lin, ref_vol)


def test_mixed_cubit_mesh_matches_the_gmsh_writer_nodes(tmp_path):
    """Hex/pyramid/tet .vol from Cubit: the BDF equals the existing Gmsh
    writer's high-order nodes truncated to the Nastran node counts."""
    vol = tmp_path / "mixed.vol"
    vol.write_bytes(gzip.decompress(MIXED_O2.read_bytes()))
    mesh = ng.Mesh(str(vol))
    assert mesh.GetCurveOrder() == 2
    bdf = tmp_path / "mixed.bdf"
    info = write_nastran_bdf(mesh, bdf, order=2, boundaries=False)
    assert set(info["elements"]) == {"CHEXA", "CTETRA", "CPYRAM"}, info

    msh = tmp_path / "mixed.msh"
    GmshPostExport(mesh).write_mesh(str(msh))
    keep = {12: (17, 20), 14: (7, 5), 11: (11, 10)}  # Hex27->20, Pyr14->5, Tet10
    trunc = tmp_path / "mixed_trunc.msh"
    gmsh.initialize(["-nt", "1"])
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        gmsh.open(str(msh))
        tags, coords, _ = gmsh.model.mesh.getNodes()
        blocks = []
        for t in gmsh.model.mesh.getElementTypes(3):
            _, nodes = gmsh.model.mesh.getElementsByType(t)
            n_all = gmsh.model.mesh.getElementProperties(t)[3]
            new_type, n_keep = keep[t]
            blocks.append((new_type, np.asarray(nodes).reshape(-1, n_all)[:, :n_keep]))
        gmsh.model.add("truncated")
        entity = gmsh.model.addDiscreteEntity(3)
        gmsh.model.mesh.addNodes(3, entity, tags, coords)
        for new_type, nodes in blocks:
            gmsh.model.mesh.addElementsByType(entity, new_type, [], nodes.reshape(-1).tolist())
        gmsh.write(str(trunc))
    finally:
        gmsh.finalize()

    got = _gmsh_measure(bdf, 3)
    want = _gmsh_measure(trunc, 3)
    assert got["types"] == want["types"], (got, want)
    assert got["min_det"] > 0, got
    assert abs(got["measure"] / want["measure"] - 1) < TOL, (got["measure"], want["measure"])


def test_curved_2d_mesh_keeps_area(tmp_path):
    face = WorkPlane().Circle(0.5).Face()
    mesh = ng.Mesh(OCCGeometry(face, dim=2).GenerateMesh(maxh=0.2, quad_dominated=True))
    mesh.Curve(2)
    out = tmp_path / "disk.bdf"
    info = write_nastran_bdf(mesh, out, order=2)
    assert set(info["elements"]) <= {"CQUAD8", "CTRIA6"} and "CQUAD8" in info["elements"], info
    got = _gmsh_measure(out, 2)
    ref = ng.Integrate(1, mesh, order=8)
    assert got["min_det"] > 0, got
    assert abs(got["measure"] / ref - 1) < TOL, (got["measure"], ref)
