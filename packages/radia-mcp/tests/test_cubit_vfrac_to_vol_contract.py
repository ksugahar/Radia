"""License-free contracts of ``cubit_vfrac_to_vol``: input validation must
fail loudly BEFORE any Sculpt / Cubit process is spawned."""
import json

import pytest

np = pytest.importorskip("numpy")

from radia_mcp.cubit.server import cubit_vfrac_to_vol  # noqa: E402


def test_missing_file_is_an_input_error(tmp_path):
    report = json.loads(cubit_vfrac_to_vol(str(tmp_path / "nope")))
    assert report["status"] == "error"
    assert report["kind"] == "input"
    assert "write_vfrac_exodus" in report["error"]


def test_wrong_element_variables_are_rejected(tmp_path):
    nc = pytest.importorskip("netCDF4")
    out = tmp_path / "bad.e.1.0"
    ds = nc.Dataset(str(out), "w", format="NETCDF3_64BIT_OFFSET")
    try:
        ds.createDimension("len_name", 256)
        ds.createDimension("time_step", None)
        ds.createDimension("num_elem_var", 1)
        var = ds.createVariable("name_elem_var", "S1",
                                ("num_elem_var", "len_name"))
        arr = np.zeros((1, 256), dtype="S1")
        for col, ch in enumerate("LSD"):
            arr[0, col] = ch.encode()
        var[:] = arr
    finally:
        ds.close()
    report = json.loads(cubit_vfrac_to_vol(str(out)))
    assert report["status"] == "error"
    assert report["kind"] == "input"
    assert "MAT_1" in report["error"]


def test_base_name_resolves_to_e_1_0(tmp_path):
    # base name without the suffix must resolve to <base>.e.1.0 and then
    # report ITS absence (not a directory or half-suffixed path)
    report = json.loads(cubit_vfrac_to_vol(str(tmp_path / "design_vf")))
    assert report["status"] == "error"
    assert "design_vf.e.1.0" in report["error"]
