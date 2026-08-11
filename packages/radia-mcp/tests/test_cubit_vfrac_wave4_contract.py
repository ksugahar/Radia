"""License-free contracts of the ``cubit_vfrac_to_vol`` wave-4 additions:
the numeric/path hardening, the Sculpt CLI construction, and the two report
helpers (``_classify_rve_sidesets`` / ``_vol_bcname_face_stats``).

Everything here runs without Cubit, Sculpt, or a license.  The real-mesher
evidence lives in ``validation_test/isochronous_topopt/test_vfrac_*_lane.py``;
this file locks the parts that only need the input contract and the argv.
"""
import json
import math
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")

from radia_mcp.cubit import server                        # noqa: E402
from radia_mcp.cubit.server import cubit_vfrac_to_vol     # noqa: E402


# ----------------------------------------------------------------------
# a minimal but VALID vfrac Exodus, so the tool proceeds past validation
# ----------------------------------------------------------------------
def _write_vfrac(path: Path, n_mat: int = 1, nel=(2, 2, 2),
                 lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0)) -> None:
    nc = pytest.importorskip("netCDF4")
    nx, ny, nz = nel
    n_cells = nx * ny * nz
    ds = nc.Dataset(str(path), "w", format="NETCDF3_64BIT_OFFSET")
    try:
        ds.createDimension("len_name", 256)
        ds.createDimension("time_step", None)
        ds.createDimension("num_elem_var", 1 + n_mat)
        ds.createDimension("num_glo_var", 21 + n_mat)
        ds.createDimension("num_el_in_blk1", n_cells)

        names = ["VOID"] + [f"MAT_{i + 1}" for i in range(n_mat)]
        var = ds.createVariable("name_elem_var", "S1",
                                ("num_elem_var", "len_name"))
        arr = np.zeros((1 + n_mat, 256), dtype="S1")
        for row, name in enumerate(names):
            for col, ch in enumerate(name):
                arr[row, col] = ch.encode()
        var[:] = arr

        gnames = ["gxint", "gyint", "gzint", "xmin", "ymin", "zmin",
                  "xmax", "ymax", "zmax"]
        gnames += [f"pad{i}" for i in range(21 + n_mat - len(gnames))]
        gvar = ds.createVariable("name_glo_var", "S1",
                                 ("num_glo_var", "len_name"))
        garr = np.zeros((21 + n_mat, 256), dtype="S1")
        for row, name in enumerate(gnames):
            for col, ch in enumerate(name):
                garr[row, col] = ch.encode()
        gvar[:] = garr
        gvals = ds.createVariable("vals_glo_var", "f8",
                                  ("time_step", "num_glo_var"))
        row = [float(nx), float(ny), float(nz),
               lo[0], lo[1], lo[2], hi[0], hi[1], hi[2]]
        row += [0.0] * (21 + n_mat - len(row))
        gvals[0, :] = row

        # VOID is index 1; materials are 2..n_mat+1, each half-filled
        for index in range(1 + n_mat):
            v = ds.createVariable(f"vals_elem_var{index + 1}eb1", "f8",
                                  ("time_step", "num_el_in_blk1"))
            v[0, :] = (0.0 if index == 0 else 1.0 / n_mat)
    finally:
        ds.close()


def _sculpt_stub(tmp_path: Path, monkeypatch, recorder: dict,
                 returncode: int = 1):
    """Point the tool at a fake sculpt.exe and capture its argv."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir(exist_ok=True)
    (bin_dir / "sculpt.exe").write_text("stub", encoding="utf-8")
    monkeypatch.setattr("radia_mcp.cubit.session.get_cubit_bin_dir",
                        lambda: str(bin_dir))

    class _Run:
        def __init__(self):
            self.returncode = returncode
            self.stdout = "stub sculpt"
            self.stderr = ""

    import subprocess as sp

    def fake_run(cli, **kwargs):
        recorder["cli"] = list(cli)
        recorder["kwargs"] = kwargs
        return _Run()

    monkeypatch.setattr(sp, "run", fake_run)


# ----------------------------------------------------------------------
# numeric + path contract (the bf8ab4c0b hardening its STL sibling had)
# ----------------------------------------------------------------------
@pytest.mark.parametrize("kwargs, needle", [
    ({"timeout_s": float("nan")}, "timeout_s"),
    ({"timeout_s": 0.0}, "timeout_s"),
    ({"closure_tolerance": float("inf")}, "closure_tolerance"),
    ({"closure_tolerance": -0.1}, "closure_tolerance"),
    ({"closure_tolerance": "abc"}, "numeric argument"),
    ({"smooth_method": "eight"}, "numeric argument"),
    ({"gq_threshold": float("nan")}, "gq_threshold"),
])
def test_numeric_contract_fails_before_any_subprocess(tmp_path, kwargs,
                                                      needle):
    vf = tmp_path / "design.e.1.0"
    _write_vfrac(vf)
    report = json.loads(cubit_vfrac_to_vol(str(vf), **kwargs))
    assert report["status"] == "error"
    assert report["kind"] == "input"
    assert needle in report["error"]


def test_coinciding_output_paths_are_rejected(tmp_path):
    vf = tmp_path / "design.e.1.0"
    _write_vfrac(vf)
    same = str(tmp_path / "both.out")
    report = json.loads(cubit_vfrac_to_vol(str(vf), out_vol=same,
                                           out_msh=same))
    assert report["status"] == "error"
    assert report["kind"] == "input"
    assert "distinct" in report["error"]


def test_material_names_count_mismatch_is_rejected(tmp_path):
    vf = tmp_path / "design.e.1.0"
    _write_vfrac(vf, n_mat=2)
    report = json.loads(cubit_vfrac_to_vol(str(vf),
                                           material_names="only_one"))
    assert report["status"] == "error"
    assert report["kind"] == "input"
    assert "2 material" in report["error"]


def test_duplicate_material_names_are_rejected(tmp_path):
    vf = tmp_path / "design.e.1.0"
    _write_vfrac(vf, n_mat=2)
    report = json.loads(cubit_vfrac_to_vol(str(vf),
                                           material_names="iron,iron"))
    assert report["status"] == "error"
    assert "distinct" in report["error"]


def test_missing_sculpt_is_an_environment_error(tmp_path, monkeypatch):
    """kind -- not just stage -- must say 'environment': the server
    instructions route retries on kind, and the default needle scan
    classifies on the MESSAGE, which would call a missing install an
    input error and invite an agent to retry it forever."""
    vf = tmp_path / "design.e.1.0"
    _write_vfrac(vf)
    monkeypatch.setattr("radia_mcp.cubit.session.get_cubit_bin_dir",
                        lambda: None)
    report = json.loads(cubit_vfrac_to_vol(str(vf)))
    assert report["status"] == "error"
    assert report["stage"] == "environment"
    assert report["kind"] == "environment"
    assert "sculpt.exe" in report["error"]


def test_missing_netcdf4_is_an_environment_error(tmp_path, monkeypatch):
    vf = tmp_path / "design.e.1.0"
    _write_vfrac(vf)
    import builtins
    real_import = builtins.__import__

    def no_netcdf(name, *args, **kwargs):
        if name == "netCDF4":
            raise ImportError("stubbed out")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", no_netcdf)
    report = json.loads(cubit_vfrac_to_vol(str(vf)))
    assert report["status"] == "error"
    assert report["kind"] == "environment"
    assert "netCDF4" in report["error"]


# ----------------------------------------------------------------------
# Sculpt CLI construction -- every flag must actually reach sculpt.exe
# ----------------------------------------------------------------------
def test_plain_route_cli_flags(tmp_path, monkeypatch):
    vf = tmp_path / "design.e.1.0"
    _write_vfrac(vf)
    rec = {}
    _sculpt_stub(tmp_path, monkeypatch, rec)
    report = json.loads(cubit_vfrac_to_vol(str(vf)))
    assert report["status"] == "error" and report["stage"] == "cubit"
    cli = rec["cli"]
    assert "-ivf" in cli and "-cv" in cli
    assert cli[cli.index("-SS") + 1] == "2"          # variable, not rve
    assert "--periodic" not in cli
    assert "-mn" in cli and "iron" in cli            # default single label


def test_periodic_route_switches_sideset_mode(tmp_path, monkeypatch):
    vf = tmp_path / "design.e.1.0"
    _write_vfrac(vf)
    rec = {}
    _sculpt_stub(tmp_path, monkeypatch, rec)
    json.loads(cubit_vfrac_to_vol(str(vf), periodic=True))
    cli = rec["cli"]
    assert cli[cli.index("-SS") + 1] == "5"          # rve
    assert "--periodic" in cli


def test_quality_and_smooth_knobs_reach_the_command(tmp_path, monkeypatch):
    vf = tmp_path / "design.e.1.0"
    _write_vfrac(vf)
    rec = {}
    _sculpt_stub(tmp_path, monkeypatch, rec)
    json.loads(cubit_vfrac_to_vol(
        str(vf), smooth_method=8, gq_iters=3, gq_threshold=0.35,
        pillow_surfaces=True, defeature=1, min_vol_cells=7))
    cli = " ".join(rec["cli"])
    assert "--smooth 8" in cli
    for token in ("3", "0.35", "1", "7"):
        assert token in rec["cli"]
    # gq_threshold must not be emitted when the iteration count is zero
    rec2 = {}
    _sculpt_stub(tmp_path, monkeypatch, rec2)
    json.loads(cubit_vfrac_to_vol(str(vf), gq_iters=0, gq_threshold=0.44))
    assert "0.44" not in rec2["cli"]


def test_multi_material_names_become_mn_pairs(tmp_path, monkeypatch):
    vf = tmp_path / "design.e.1.0"
    _write_vfrac(vf, n_mat=2)
    rec = {}
    _sculpt_stub(tmp_path, monkeypatch, rec)
    json.loads(cubit_vfrac_to_vol(str(vf), material_names="core,shell"))
    cli = rec["cli"]
    assert cli[cli.index("core") - 1] == "1"
    assert cli[cli.index("shell") - 1] == "2"


def test_sculpt_timeout_is_reported_as_a_cubit_stage_error(tmp_path,
                                                           monkeypatch):
    vf = tmp_path / "design.e.1.0"
    _write_vfrac(vf)
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir(exist_ok=True)
    (bin_dir / "sculpt.exe").write_text("stub", encoding="utf-8")
    monkeypatch.setattr("radia_mcp.cubit.session.get_cubit_bin_dir",
                        lambda: str(bin_dir))
    import subprocess as sp

    def fake_run(cli, **kwargs):
        raise sp.TimeoutExpired(cli, kwargs.get("timeout", 1))

    monkeypatch.setattr(sp, "run", fake_run)
    report = json.loads(cubit_vfrac_to_vol(str(vf), timeout_s=5))
    assert report["status"] == "error"
    assert report["stage"] == "cubit"
    assert "timed out" in report["error"]


# ----------------------------------------------------------------------
# report helpers
# ----------------------------------------------------------------------
def test_closure_miss_diagnosis_names_both_knobs_and_prescribes_neither():
    """A closure miss must NOT be reported as 'raise cells and retry'.

    Measured 2026-08-11 (6-pole ring, analytic volume known, design mesh
    nz=10): cells 56 -> 92 turned a 1.27 % gate FAILURE into a 0.37 %
    pass while moving the mesh from 0.16 % to 1.8 % away from the true
    volume -- the fraction field was the inaccurate quantity, and the
    finer lattice aliased an under-resolved design field.  The report
    therefore has to name the design mesh as a candidate too.
    """
    report = server._closure_miss_diagnosis([56, 56, 20], 0.0127, 0.01)
    assert report["lattice_cells"] == [56, 56, 20]
    # a miss by 1.27x asks for a proportionally finer lattice IF the
    # mesher is what is limiting -- strictly larger, and capped
    assert 56 < report["cells_if_mesher_limited"] <= 512
    note = report["note"]
    assert "DESIGN MESH" in note, "the second knob must be named"
    assert "WORSE" in note, (
        "the note must state that a finer lattice can make the mesh "
        "worse, otherwise it reads as a prescription")
    assert "independent reference" in note

    # the suggestion must saturate at the writer's own cap rather than
    # emitting a lattice write_vfrac_exodus would reject
    huge = server._closure_miss_diagnosis([400, 400, 400], 0.9, 0.01)
    assert huge["cells_if_mesher_limited"] == 512


def test_bcname_face_stats_groups_faces_by_name(tmp_path):
    """A hand-written .vol: two named boundaries on opposite z planes."""
    vol = tmp_path / "toy.vol"
    vol.write_text(
        "mesh3d\ndimension\n3\ngeomtype\n0\n\n"
        "surfaceelements\n2\n"
        "  1 1 1 0 3 1 2 3\n"
        "  2 2 1 0 3 4 5 6\n\n"
        "bcnames\n2\n1\tbottom\n2\ttop\n\n"
        "points\n6\n"
        "0 0 0\n1 0 0\n0 1 0\n"
        "0 0 2\n1 0 2\n0 1 2\n",
        encoding="utf-8")
    stats = server._vol_bcname_face_stats(vol)
    assert set(stats) == {"bottom", "top"}
    assert stats["bottom"]["n_faces"] == 1
    assert stats["bottom"]["hi"][2] == pytest.approx(0.0)
    assert stats["top"]["lo"][2] == pytest.approx(2.0)


def test_bcname_face_stats_returns_empty_without_points(tmp_path):
    vol = tmp_path / "nopoints.vol"
    vol.write_text("mesh3d\nsurfaceelements\n0\n", encoding="utf-8")
    assert server._vol_bcname_face_stats(vol) == {}


def test_rve_classification_rejects_a_non_rve_exodus(tmp_path):
    """Fewer than six sidesets means Sculpt's rve mode did not fire --
    fail loudly instead of naming whatever is there."""
    nc = pytest.importorskip("netCDF4")
    exo = tmp_path / "plain_sculpt.e.1.0"
    ds = nc.Dataset(str(exo), "w", format="NETCDF3_64BIT_OFFSET")
    try:
        ds.createDimension("num_nodes", 8)
        ds.createDimension("num_el_blk", 1)
        ds.createDimension("num_el_in_blk1", 1)
        ds.createDimension("num_nod_per_el1", 8)
        ds.createDimension("num_side_sets", 1)
        for name, values in (("coordx", [0, 1, 1, 0, 0, 1, 1, 0]),
                             ("coordy", [0, 0, 1, 1, 0, 0, 1, 1]),
                             ("coordz", [0, 0, 0, 0, 1, 1, 1, 1])):
            v = ds.createVariable(name, "f8", ("num_nodes",))
            v[:] = np.asarray(values, dtype=float)
        blk = ds.createVariable("eb_prop1", "i4", ("num_el_blk",))
        blk[:] = [1]
        conn = ds.createVariable("connect1", "i4",
                                 ("num_el_in_blk1", "num_nod_per_el1"))
        conn[:] = np.arange(1, 9).reshape(1, 8)
        ss = ds.createVariable("ss_prop1", "i4", ("num_side_sets",))
        ss[:] = [1]
        ds.createDimension("num_side_ss1", 1)
        e1 = ds.createVariable("elem_ss1", "i4", ("num_side_ss1",))
        e1[:] = [1]
    finally:
        ds.close()
    names, error = server._classify_rve_sidesets(
        exo, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [0.5, 0.5, 0.5])
    assert names == {}
    assert "rve sideset mode did not fire" in error
