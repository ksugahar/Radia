"""IH notebook comparison and application-interface manifest contracts.

Only IH keeps a Jupyter workbench while its Simulink operation is evaluated.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
import time

from radia.ih_design import IHDesignSpec, METHOD_PEEC_IND
from radia.ih_notebook import IHWorkbench
from radia.notebook_workbench import CommandWorkbench


class EchoSpec:
    def build_command(self) -> list[str]:
        return [sys.executable, "-c", "print('notebook-runner-ok')"]

    def missing_required_inputs(self) -> list[str]:
        return []


class SleepSpec:
    def build_command(self) -> list[str]:
        return [
            sys.executable,
            "-c",
            "import time; print('sleeping', flush=True); time.sleep(10)",
        ]

    def missing_required_inputs(self) -> list[str]:
        return []


class TimedOutputSpec:
    def __init__(self, output_path: Path):
        self.output_path = output_path

    def build_command(self) -> list[str]:
        path_literal = json.dumps(str(self.output_path))
        code = (
            "import json, pathlib; "
            f"p=pathlib.Path({path_literal}); "
            "p.parent.mkdir(parents=True, exist_ok=True); "
            "p.write_text(json.dumps({"
            "'t_mesh_s': 1.5, "
            "'t_solve_s': 9.0, "
            "'t_post_s': 0.2, "
            "'t_assembly_s': 4.0, "
            "'t_io_s': 0.1, "
            "'nested': {'t_factor_s': 6.0}"
            "}), encoding='utf-8')"
        )
        return [sys.executable, "-c", code, "--output", str(self.output_path)]

    def missing_required_inputs(self) -> list[str]:
        return []


def test_run_local_writes_radia_result_artifact(tmp_path: Path):
    wb = CommandWorkbench(EchoSpec(), run_root=tmp_path / "runs", timeout_s=5)

    record = wb.run_local()

    assert record.status == "passed"
    assert record.returncode == 0
    assert record.log_path.is_file()
    assert record.result_path.is_file()
    assert "notebook-runner-ok" in record.log_path.read_text(encoding="utf-8")
    payload = json.loads(record.result_path.read_text(encoding="utf-8"))
    assert list(payload) == ["radia_result"]
    result = payload["radia_result"]
    assert result["schema"] == "radia.notebook_panel_run.v2"
    assert result["status"] == "passed"
    assert result["runtime_radia_version"]
    assert result["executed_at_utc"]
    assert result["completed_at_utc"]
    assert result["timing"]["wall_elapsed_s"] >= 0
    assert result["timing"]["top_stages"] == []


def test_run_local_collects_top_four_cli_timing_stages(tmp_path: Path):
    output_path = tmp_path / "solver_output.json"
    wb = CommandWorkbench(
        TimedOutputSpec(output_path),
        run_root=tmp_path / "runs",
        timeout_s=5,
    )

    record = wb.run_local()

    payload = json.loads(record.result_path.read_text(encoding="utf-8"))
    stages = payload["radia_result"]["timing"]["top_stages"]
    assert [stage["name"] for stage in stages] == [
        "t_solve_s",
        "nested.t_factor_s",
        "t_assembly_s",
        "t_mesh_s",
    ]
    assert [stage["elapsed_s"] for stage in stages] == [9.0, 6.0, 4.0, 1.5]
    assert all(stage["source"] == str(output_path) for stage in stages)


def test_run_local_timeout_is_recorded(tmp_path: Path):
    wb = CommandWorkbench(SleepSpec(), run_root=tmp_path / "runs", timeout_s=1)

    record = wb.run_local(timeout_s=1)

    assert record.status == "timeout"
    payload = json.loads(record.result_path.read_text(encoding="utf-8"))
    assert payload["radia_result"]["status"] == "timeout"
    assert payload["radia_result"]["timeout_s"] == 1


def test_background_run_can_be_cancelled(tmp_path: Path):
    wb = CommandWorkbench(SleepSpec(), run_root=tmp_path / "runs", timeout_s=30)

    thread = wb.start_background_run(timeout_s=30)
    time.sleep(0.3)
    wb.cancel()
    thread.join(timeout=5)

    assert not thread.is_alive()
    assert wb.last_run is not None
    assert wb.last_run.status == "cancelled"
    payload = json.loads(wb.last_run.result_path.read_text(encoding="utf-8"))
    assert payload["radia_result"]["status"] == "cancelled"


def test_ih_comparison_workbench_has_app_specific_run_root():
    assert str(IHWorkbench().run_root).endswith("radia_ih")


def test_ih_comparison_workbench_builds_headless_command():
    workbench = IHWorkbench(
        IHDesignSpec(method=METHOD_PEEC_IND, peec_step="coil.step")
    )
    command = workbench.build_command()
    assert any("calc_inductance.py" in part for part in command), command


def test_spec_cell_source_makes_notebook_initial_values_canonical():
    workbench = IHWorkbench(
        IHDesignSpec(method=METHOD_PEEC_IND, peec_step="coil.step")
    )

    source = workbench.spec_cell_source()

    assert "from radia.ih_design import IHDesignSpec" in source
    assert "from radia.ih_notebook import IHWorkbench" in source
    assert "IHDesignSpec(**" in source
    assert "'peec_step': 'coil.step'" in source
    assert "json" not in source.lower()


def test_application_manifest_records_simulink_first_and_ih_dual_state():
    manifest_path = (
        Path(__file__).resolve().parents[2]
        / "src"
        / "radia"
        / "panels"
        / "application_interface_manifest.json"
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    states = {item["id"]: item["state"] for item in manifest["applications"]}

    for application_id in (
        "radia-em",
        "radia-pcb",
        "radia-motor",
        "radia-streamfunction",
    ):
        assert states[application_id] == "active-simulink-block"
    assert states["radia-ih"] == "active-dual-comparison"
    assert states["radia-export-menu"] == "active-cubit-toolbar"
    assert manifest["library"]["artifact"] == "matlab/radia_simulink_library.slx"
    assert manifest["library"]["backend"] == "python-headless-cli"


def test_panel_notebooks_do_not_import_pyside():
    notebook_dir = (
        Path(__file__).resolve().parents[2]
        / "src"
        / "radia"
        / "panels"
        / "notebooks"
    )
    for notebook_path in notebook_dir.glob("*.ipynb"):
        text = notebook_path.read_text(encoding="utf-8")
        assert "active-ipywidgets" not in text
        assert "PySide6" not in text
        assert "PyQt" not in text


def test_only_ih_analysis_notebook_remains():
    notebook_dir = (
        Path(__file__).resolve().parents[2]
        / "src"
        / "radia"
        / "panels"
        / "notebooks"
    )
    analysis_notebooks = {
        path.name for path in notebook_dir.glob("radia_*.ipynb")
        if path.name != "radia_export_menu.ipynb"
    }
    assert analysis_notebooks == {"radia_ih.ipynb"}
    text = (notebook_dir / "radia_ih.ipynb").read_text(encoding="utf-8")
    assert "from radia.ih_design import IHDesignSpec" in text
    assert "spec = IHDesignSpec()" in text
    assert "JSON files are run artifacts, not preset storage" in text


def test_ih_notebook_includes_comparison_workbench_notes():
    notebook_dir = (
        Path(__file__).resolve().parents[2]
        / "src"
        / "radia"
        / "panels"
        / "notebooks"
    )
    text = (notebook_dir / "radia_ih.ipynb").read_text(encoding="utf-8")
    assert "## Notebook Panel Notes" in text
    assert "IH comparison workbench" in text
    assert "Simulink block" in text
    assert "app-specific `calc_*.py` CLI" in text
    assert "DesignSpec" in text
    assert "Run local" in text
    assert "netgen.webgui" in text
    assert "GMSH `.msh v4.1`" in text


def test_ih_notebook_carries_esim_and_previous_result_notes():
    notebook_path = (
        Path(__file__).resolve().parents[2]
        / "src"
        / "radia"
        / "panels"
        / "notebooks"
        / "radia_ih.ipynb"
    )
    text = notebook_path.read_text(encoding="utf-8")

    assert "## IH Workpiece Notes" in text
    assert "Linear SIBC" in text
    assert "Nonlinear ESIM" in text
    assert "1-D cell problem" in text
    assert "`Z_s`" in text
    assert "esim_converged" in text
    assert "esim_iterations" in text
    assert "docs/ih_esim_benchmark/sweep_data_dense/" in text
    assert "I100_f50k_scalar.json" in text
    assert "I100_f50k_per_panel.json" in text
