"""Periodic RVE lane: Sculpt `--periodic` on the volume-fraction route.

A periodic unit cell is the mesh building block for arrays -- magnet
arrays, lamination RVEs, motor unit cells.  What makes a cell tileable is
NOT a flat boundary (Sculpt deliberately leaves it ragged) but that the
node set is INVARIANT under the period: the nodes of one period face are
the nodes of the opposite face, translated, so two copies laid side by
side share their seam exactly.

Body: a simple-cubic array of overlapping spheres (R = 0.55 L, so the
material crosses every face -- a body that does not touch the boundary
would make the test vacuous), sampled with the minimum-image distance so
the density is periodic by construction.

Gates, with `periodic=False` as the negative control on the SAME input:
periodic gives matched +period / -period node counts that are equal and
nonzero on all three axes; plain gives none.

Requires a Cubit license (sculpt.exe + headless import/export); ~4 min.
"""
import json
from pathlib import Path

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netCDF4")
pytest.importorskip("scipy")

from ngsolve import TaskManager  # noqa: E402
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402

from radia.topopt_cad import write_vfrac_exodus  # noqa: E402

pytest.importorskip("radia_mcp.cubit.server")

PERIOD = 0.02
RADIUS = 0.55 * PERIOD


def _minimum_image(delta):
    return delta - PERIOD * np.round(delta / PERIOD)


@pytest.fixture(scope="module")
def artifacts(tmp_path_factory):
    root = tmp_path_factory.mktemp("vfrac_periodic")
    with TaskManager():
        cell = MakeStructured3DMesh(
            nx=12, ny=12, nz=12,
            mapping=lambda x, y, z: (PERIOD * x, PERIOD * y, PERIOD * z))
    pts = np.array([list(v.point) for v in cell.vertices])
    radius = np.sqrt(sum(_minimum_image(pts[:, i] - 0.5 * PERIOD) ** 2
                         for i in range(3)))
    nodal = (radius <= RADIUS).astype(float)
    vf = write_vfrac_exodus(
        cell, nodal, root / "rve", level=0.5, cells=24, supersample=4,
        bounds=((0.0, 0.0, 0.0), (PERIOD, PERIOD, PERIOD)))

    from radia_mcp.cubit.server import cubit_vfrac_to_vol
    reports = {}
    for tag, flag in (("plain", False), ("periodic", True)):
        reports[tag] = json.loads(cubit_vfrac_to_vol(
            vf["path"], out_vol=str(root / f"{tag}.vol"),
            out_msh=str(root / f"{tag}.msh"), periodic=flag))
    return {"vf": vf, "reports": reports, "root": root}


def test_bounds_override_spans_exactly_one_period(artifacts):
    vf = artifacts["vf"]
    # no guard ring, and the box is the requested period exactly -- a pad
    # would put the grid off the period and break Sculpt's check
    np.testing.assert_allclose(vf["bounds_min"], [0.0, 0.0, 0.0], atol=0)
    np.testing.assert_allclose(vf["bounds_max"],
                               [PERIOD, PERIOD, PERIOD], rtol=1e-12)
    assert vf["nel"] == [24, 24, 24]
    np.testing.assert_allclose(vf["cell_size"], PERIOD / 24, rtol=1e-12)


def test_periodic_mesh_node_set_is_translation_invariant(artifacts):
    report = artifacts["reports"]["periodic"]
    assert report.get("status") == "ok", report.get("gates", report)
    assert report["gates"]["periodic_nodes_match"] is True
    match = report["periodic_node_match"]
    assert match["ok"] is True
    for entry in match["axes"]:
        # equal counts = a bijection between the two period faces
        assert entry["matched_plus"] == entry["matched_minus"], entry
        # the body crosses every face, so every axis must carry a seam
        assert entry["matched_plus"] > 50, entry
        # ... and that seam is COMPLETE (measured 1.000 on all three axes)
        assert entry["seam_fill_rate"] >= 0.98, entry


def test_plain_mesh_is_the_negative_control(artifacts):
    report = artifacts["reports"]["plain"]
    assert report.get("status") == "ok", report.get("gates", report)
    # the plain route does not claim periodicity and must not report it
    assert "periodic_nodes_match" not in report["gates"]
    assert report["periodic_node_match"] is None
    # ... and the same input without --periodic must FAIL the seam test,
    # otherwise the gate proves nothing.  Note the discriminator is the
    # fill RATE, not the raw count: a plain Sculpt mesh keeps a handful of
    # untouched lattice corners that match by accident (measured 9).
    from radia_mcp.cubit.server import _periodic_node_match
    cell = artifacts["vf"]["cell_size"]
    measured = _periodic_node_match(Path(report["vol"]),
                                    [PERIOD, PERIOD, PERIOD],
                                    [cell, cell, cell])
    assert measured["axes"], measured
    worst = max(entry["seam_fill_rate"] for entry in measured["axes"])
    assert worst < 0.5, (
        f"the non-periodic mesh filled {worst:.3f} of its seam -- the "
        "gate would not discriminate")


def test_periodic_cell_still_passes_the_solver_gates(artifacts):
    """Periodicity must not cost the solver contract: closure, no
    inverted elements, and the outer skin still exported."""
    report = artifacts["reports"]["periodic"]
    assert report["gates"]["closure_ok"] is True
    assert report["gates"]["no_inverted_elements"] is True
    assert report["gates"]["boundary_faces_ok"] is True
    # the ragged periodic boundary costs nothing in volume fidelity
    assert report["closure"] < 5e-3, report["closure"]


def test_period_faces_arrive_as_named_boundaries(artifacts):
    """The six period faces must reach the `.vol` as NAMED boundaries.

    Sculpt leaves the periodic boundary deliberately ragged, so a
    downstream coordinate filter cannot separate a period face from the
    interior -- carrying Sculpt's own sideset knowledge through as
    ``rve_xmin`` .. ``rve_zmax`` is the only reliable handle for periodic
    BCs and for a strict-label contract on an RVE mesh.
    """
    report = artifacts["reports"]["periodic"]
    assert report["gates"]["rve_labels_ok"] is True
    names = set((report.get("rve_sideset_names") or {}).values())
    expected = {"rve_xmin", "rve_xmax", "rve_ymin", "rve_ymax",
                "rve_zmin", "rve_zmax"}
    assert expected <= names, names
    assert "skin" in names, "the interior material surface must stay named"

    # the labels must survive INTO the .vol, not just the report
    text = Path(report["vol"]).read_text(encoding="utf-8", errors="replace")
    marker = text.find("bcnames")
    assert marker >= 0, "the exported .vol carries no bcnames section"
    lines = text[marker:].splitlines()
    count = int(lines[1])
    in_vol = {parts[1] for parts in
              (line.split() for line in lines[2:2 + count])
              if len(parts) >= 2}
    assert expected <= in_vol, in_vol

    # ... and each named face must actually HUG its plane (the report's
    # own gate re-verified here from the geometry, independently of the
    # Exodus-side classification that produced the names)
    labels = report["rve_boundary_labels"]
    cell = artifacts["vf"]["cell_size"]
    for axis, low, high in ((0, "rve_xmin", "rve_xmax"),
                            (1, "rve_ymin", "rve_ymax"),
                            (2, "rve_zmin", "rve_zmax")):
        assert labels[low]["hi"][axis] <= 1.5 * cell, labels[low]
        assert labels[high]["lo"][axis] >= PERIOD - 1.5 * cell, labels[high]
        # opposite faces, not the same one twice
        assert labels[low]["n_faces"] > 0 and labels[high]["n_faces"] > 0


def test_plain_route_reports_no_rve_labels(artifacts):
    """The negative control: without --periodic there are no period faces
    to name, and the tool must not invent them."""
    report = artifacts["reports"]["plain"]
    assert report.get("rve_sideset_names") is None
    assert report.get("rve_boundary_labels") is None
    assert "rve_labels_ok" not in report["gates"]
