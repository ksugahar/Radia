"""Notebook panel workbench contracts.

These tests keep the promoted ipynb panels usable from Jupyter without
requiring PySide6 or a live notebook server.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
import time

from radia.em_notebook import EMWorkbench
from radia.em_design import EMDesignSpec
from radia.ih_design import IHDesignSpec, METHOD_PEEC_IND
from radia.ih_notebook import IHWorkbench
from radia.motor_design import ANALYSIS_LAMINATION, MotorDesignSpec
from radia.motor_notebook import MotorWorkbench
from radia.notebook_workbench import CommandWorkbench
from radia.pcb_design import PCBDesignSpec
from radia.pcb_notebook import PCBWorkbench
from radia.streamfunction_design import StreamFunctionDesignSpec
from radia.streamfunction_notebook import StreamFunctionWorkbench


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


def test_promoted_workbenches_have_app_specific_run_roots():
    assert str(IHWorkbench().run_root).endswith("radia_ih")
    assert str(EMWorkbench().run_root).endswith("radia_em")
    assert str(PCBWorkbench().run_root).endswith("radia_pcb")
    assert str(MotorWorkbench().run_root).endswith("radia_motor")
    assert str(StreamFunctionWorkbench().run_root).endswith("radia_streamfunction")


def test_promoted_workbenches_build_headless_commands():
    cases = [
        (EMWorkbench(EMDesignSpec(coil_script="coil.py")), "calc_accel_magnet.py"),
        (
            IHWorkbench(IHDesignSpec(method=METHOD_PEEC_IND, peec_step="coil.step")),
            "calc_inductance.py",
        ),
        (PCBWorkbench(PCBDesignSpec(inp="board.inp")), "calc_pcb_peec.py"),
        (
            MotorWorkbench(MotorDesignSpec(analysis=ANALYSIS_LAMINATION)),
            "calc_motor_lamination.py",
        ),
        (
            StreamFunctionWorkbench(
                StreamFunctionDesignSpec(coil_vol="coil.vol", eval_vol="eval.vol")
            ),
            "calc_streamfunction.py",
        ),
    ]

    for workbench, script_name in cases:
        command = workbench.build_command()
        assert any(script_name in part for part in command), command


def test_spec_cell_source_makes_notebook_initial_values_canonical():
    workbench = EMWorkbench(EMDesignSpec(coil_script="coil.py", n_steps=3))

    source = workbench.spec_cell_source()

    assert "from radia.em_design import EMDesignSpec" in source
    assert "from radia.em_notebook import EMWorkbench" in source
    assert "EMDesignSpec(**" in source
    assert "'coil_script': 'coil.py'" in source
    assert "'n_steps': 3" in source
    assert "json" not in source.lower()


def test_panel_notebooks_are_marked_as_local_runner():
    manifest_path = (
        Path(__file__).resolve().parents[2]
        / "src"
        / "radia"
        / "panels"
        / "notebooks"
        / "panel_notebook_manifest.json"
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    states = {panel["id"]: panel["state"] for panel in manifest["panels"]}

    for panel_id in (
        "radia-ih",
        "radia-em",
        "radia-pcb",
        "radia-motor",
        "radia-streamfunction",
    ):
        assert states[panel_id] == "active-local-runner"
    assert states["radia-export-menu"] == "active-cubit-toolbar"
    assert "calc_*.py CLI arguments" in manifest["policy"]
    assert "DesignSpec settings" in manifest["policy"]
    assert "previous in-repo result artifacts" in manifest["policy"]


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


def test_active_panel_notebooks_use_designspec_cells_for_initial_values():
    notebook_dir = (
        Path(__file__).resolve().parents[2]
        / "src"
        / "radia"
        / "panels"
        / "notebooks"
    )
    active = {
        "radia_em.ipynb": "EMDesignSpec",
        "radia_ih.ipynb": "IHDesignSpec",
        "radia_pcb.ipynb": "PCBDesignSpec",
        "radia_motor.ipynb": "MotorDesignSpec",
        "radia_streamfunction.ipynb": "StreamFunctionDesignSpec",
    }

    for notebook_name, spec_name in active.items():
        text = (notebook_dir / notebook_name).read_text(encoding="utf-8")
        assert f"from radia." in text
        assert f"import {spec_name}" in text
        assert f"spec = {spec_name}()" in text
        assert "JSON files are run artifacts, not preset storage" in text


def test_active_panel_notebooks_include_panel_notes():
    notebook_dir = (
        Path(__file__).resolve().parents[2]
        / "src"
        / "radia"
        / "panels"
        / "notebooks"
    )
    active_names = (
        "radia_em.ipynb",
        "radia_ih.ipynb",
        "radia_pcb.ipynb",
        "radia_motor.ipynb",
        "radia_streamfunction.ipynb",
    )

    for notebook_name in active_names:
        text = (notebook_dir / notebook_name).read_text(encoding="utf-8")
        assert "## Notebook Panel Notes" in text
        assert "This notebook is the panel" in text
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
