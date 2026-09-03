from __future__ import annotations

import gzip
import importlib.util
import json
import shutil
from pathlib import Path

import numpy as np
import pytest

from radia.simulink import ih_operator_assembly as assembly


@pytest.fixture(autouse=True)
def _taskmanager():
    if importlib.util.find_spec("ngsolve") is None:
        yield
        return

    import ngsolve as ng

    with ng.TaskManager():
        yield


def _workpiece(path: Path, *, material: str = "workpiece", boundary: str = "sibc") -> Path:
    pytest.importorskip("ngsolve")
    from netgen.occ import Box, OCCGeometry, Pnt

    shape = Box(Pnt(0, 0, 0), Pnt(0.01, 0.01, 0.01))
    shape.mat(material)
    for face in shape.faces:
        face.name = boundary
    mesh = OCCGeometry(shape).GenerateMesh(maxh=0.006)
    mesh.Save(str(path))
    return path


def _bema_coil(path: Path) -> Path:
    pytest.importorskip("ngsolve")
    from netgen.occ import Box, OCCGeometry, Pnt, X

    shape = Box(Pnt(0, 0, 0), Pnt(0.01, 0.002, 0.002))
    shape.mat("coil")
    shape.faces.name = "body"
    shape.faces.Min(X).name = "source"
    shape.faces.Max(X).name = "sink"
    mesh = OCCGeometry(shape).GenerateMesh(maxh=0.003)
    mesh.Save(str(path))
    return path


def _fake_unit_current(workpiece, coil, backend, options, run_dir):
    from ngsolve import BND, CF, H1, GridFunction, Integrate, Mesh

    from radia.gmsh_post_export import GmshPostExport

    mesh = Mesh(str(workpiece))
    fes = H1(mesh, order=1)
    qsurf = GridFunction(fes)
    qsurf.vec[:] = 125.0
    qsurf_path = run_dir / "fake_qsurf.sol"
    qsurf.Save(str(qsurf_path))
    field_path = run_dir / "fake_electromagnetic_fields.msh"
    fields = GmshPostExport(mesh)
    fields.add_scalar_field("fake_heat_flux", qsurf)
    fields.write(str(field_path))
    power = float(
        Integrate(
            CF(125.0),
            mesh,
            BND,
            definedon=mesh.Boundaries(options.workpiece_label),
        )
    )
    return assembly.UnitCurrentResult(
        np.full(fes.ndof, 125.0),
        {
            "P_wp_W": power,
            "method": f"fake-{backend}",
            "msh_file": str(field_path),
        },
        qsurf_path,
        field_path,
    )


def test_default_output_path_uses_both_geometry_names(tmp_path):
    output = assembly.default_output_path(tmp_path / "workpiece.vol.gz", tmp_path / "coil.step")
    assert output.name == "workpiece_coil_ih_native.json"


def test_options_reject_nonphysical_values():
    with pytest.raises(ValueError, match="frequency_hz"):
        assembly.IHOperatorAssemblyOptions(frequency_hz=0.0).checked()
    with pytest.raises(ValueError, match="thermal_order"):
        assembly.IHOperatorAssemblyOptions(thermal_order=2).checked()
    with pytest.raises(ValueError, match="versioned IH label contract"):
        assembly.IHOperatorAssemblyOptions(workpiece_label="surface").checked()


def test_geometry_to_native_config_preserves_power(monkeypatch, tmp_path):
    pytest.importorskip("cubit_mesh_export")
    workpiece = _workpiece(tmp_path / "workpiece.vol")
    coil = tmp_path / "coil.step"
    coil.write_text("STEP fixture is not parsed by this focused test\n", encoding="ascii")
    monkeypatch.setattr(assembly, "_solve_unit_current", _fake_unit_current)

    output = tmp_path / "native.json"
    run_dir = tmp_path / "run"
    config = assembly.assemble_ih_operators(
        workpiece,
        coil,
        output=output,
        run_dir=run_dir,
        options=assembly.IHOperatorAssemblyOptions(peec_proximity=False),
    )

    saved = json.loads(output.read_text(encoding="utf-8"))
    assert saved == config
    assert config["schema"] == assembly.CONFIG_SCHEMA
    assert config["operator_basis"] == "exact-single-current-linear-response"
    assert config["surrogate"] is False
    assert config["eddy_solver"] == "peec"
    assert config["n_eddy_unknown"] == 1
    assert config["n_heat"] > 0
    assert config["n_temperature"] >= config["n_heat"]
    assert len(config["heat_projection"]) == config["n_heat"]
    assert len(config["heat_to_temperature_projection"]) == (
        config["n_temperature"] * config["n_heat"]
    )
    assert config["unit_current"]["relative_power_error"] < 1.0e-12
    assert all(Path(path).is_file() for path in config["vol_check_reports"])
    gmsh = [Path(path) for path in config["artifacts"]["gmsh"]]
    assert len(gmsh) == 3
    assert len({path.name for path in gmsh}) == 3
    assert all(path.is_file() for path in gmsh)
    assert (run_dir / "run.log").is_file()
    result = json.loads((run_dir / "result.json").read_text(encoding="utf-8"))
    assert result["radia_result"]["status"] == "passed"


def test_crlf_vol_gz_is_materialized_for_ngsolve(monkeypatch, tmp_path):
    pytest.importorskip("cubit_mesh_export")
    workpiece = _workpiece(tmp_path / "workpiece.vol")
    compressed = tmp_path / "workpiece.vol.gz"
    with workpiece.open("rb") as source, gzip.open(compressed, "wb") as target:
        shutil.copyfileobj(source, target)
    coil = tmp_path / "coil.step"
    coil.write_text("STEP fixture is not parsed by this focused test\n", encoding="ascii")
    seen_workpiece = None

    def capture_materialized(solver_workpiece, *args, **kwargs):
        nonlocal seen_workpiece
        seen_workpiece = Path(solver_workpiece)
        return _fake_unit_current(solver_workpiece, *args, **kwargs)

    monkeypatch.setattr(assembly, "_solve_unit_current", capture_materialized)
    run_dir = tmp_path / "run"
    config = assembly.assemble_ih_operators(
        compressed,
        coil,
        output=tmp_path / "native.json",
        run_dir=run_dir,
        options=assembly.IHOperatorAssemblyOptions(peec_proximity=False),
    )
    assert seen_workpiece == run_dir / "workpiece.solver.vol"
    assert seen_workpiece.is_file()
    assert b"\r" not in seen_workpiece.read_bytes()
    assert config["geometry"]["workpiece_vol"] == str(compressed.resolve())
    assert config["geometry"]["solver_workpiece_vol"] == str(seen_workpiece)
    assert config["geometry"]["gzip_materialized"] is True


def test_bema_coil_vol_gz_is_materialized_after_contract_check(monkeypatch, tmp_path):
    pytest.importorskip("cubit_mesh_export")
    workpiece = _workpiece(tmp_path / "workpiece.vol")
    coil = _bema_coil(tmp_path / "coil.vol")
    compressed = tmp_path / "coil.vol.gz"
    with coil.open("rb") as source, gzip.open(compressed, "wb") as target:
        shutil.copyfileobj(source, target)
    seen_coil = None

    def capture_materialized(solver_workpiece, solver_coil, *args, **kwargs):
        nonlocal seen_coil
        seen_coil = Path(solver_coil)
        return _fake_unit_current(solver_workpiece, solver_coil, *args, **kwargs)

    monkeypatch.setattr(assembly, "_solve_unit_current", capture_materialized)
    run_dir = tmp_path / "run"
    config = assembly.assemble_ih_operators(
        workpiece,
        compressed,
        output=tmp_path / "native.json",
        run_dir=run_dir,
    )
    assert seen_coil == run_dir / "coil.solver.vol"
    assert seen_coil.is_file()
    assert b"\r" not in seen_coil.read_bytes()
    assert config["eddy_solver"] == "bem-a"
    assert config["geometry"]["coil_file"] == str(compressed.resolve())
    assert config["geometry"]["solver_coil_file"] == str(seen_coil)
    assert config["geometry"]["gzip_materialized"] is True


def test_strict_workpiece_labels_fail_before_solver(monkeypatch, tmp_path):
    pytest.importorskip("cubit_mesh_export")
    workpiece = _workpiece(tmp_path / "bad.vol", material="steel", boundary="surface")
    coil = tmp_path / "coil.step"
    coil.write_text("not reached\n", encoding="ascii")
    called = False

    def should_not_run(*args, **kwargs):
        nonlocal called
        called = True
        raise AssertionError("electromagnetic solve must follow check-vol")

    monkeypatch.setattr(assembly, "_solve_unit_current", should_not_run)
    with pytest.raises(RuntimeError, match="check-vol failed"):
        assembly.assemble_ih_operators(
            workpiece,
            coil,
            output=tmp_path / "bad.json",
            run_dir=tmp_path / "bad-run",
        )
    assert called is False
    report = json.loads(
        (tmp_path / "bad-run" / "workpiece.vol-check.json").read_text(encoding="utf-8")
    )
    assert report["passed"] is False
    assert report["labels"]["missing"]["materials"] == ["workpiece"]


def test_hole_requires_solved_cohomology_mode(monkeypatch, tmp_path):
    pytest.importorskip("cubit_mesh_export")
    workpiece = _workpiece(tmp_path / "workpiece.vol")
    coil = tmp_path / "coil.step"
    coil.write_text("STEP fixture is not parsed by this focused test\n", encoding="ascii")

    def missing_loop(*args, **kwargs):
        result = _fake_unit_current(*args, **kwargs)
        payload = dict(result.solver_payload)
        payload.update(
            wp_genus=1,
            wp_loop_dof=False,
            wp_loop_dof_skip_reason="focused missing-loop fixture",
        )
        return assembly.UnitCurrentResult(
            result.heat_flux_W_per_m2,
            payload,
            result.qsurf_solution,
            result.field_mesh,
        )

    monkeypatch.setattr(assembly, "_solve_unit_current", missing_loop)
    output = tmp_path / "native.json"
    run_dir = tmp_path / "run"
    with pytest.raises(RuntimeError, match="requires a solved cohomology loop DOF"):
        assembly.assemble_ih_operators(
            workpiece,
            coil,
            output=output,
            run_dir=run_dir,
            options=assembly.IHOperatorAssemblyOptions(peec_proximity=False),
        )
    assert not output.exists()
    result = json.loads((run_dir / "result.json").read_text(encoding="utf-8"))
    assert result["radia_result"]["status"] == "failed"
    assert result["radia_result"]["config"] is None


def test_bema_coil_requires_versioned_source_sink_contract(tmp_path):
    workpiece = tmp_path / "workpiece.vol"
    coil = tmp_path / "coil.vol"
    workpiece.write_text("placeholder", encoding="ascii")
    coil.write_text("placeholder", encoding="ascii")
    _, _, backend = assembly._checked_geometry_paths(workpiece, coil)
    assert backend == "bem-a"
    contract = json.loads(
        assembly._contract_path("ih_coil_bema_v1.json").read_text(encoding="utf-8")
    )
    assert contract["schema"] == "radia.vol-label-contract.v1"
    assert contract["required"]["boundaries"] == ["body", "source", "sink"]
    assert contract["allowed"]["boundaries"] == ["body", "source", "sink"]
