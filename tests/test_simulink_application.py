from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import sys

import pytest

from radia.simulink import application as app


@dataclass
class EchoDesignSpec:
    input_file: str = ""
    value: float = 12.5

    def missing_required_inputs(self) -> list[str]:
        return [] if self.input_file else ["Input file"]

    def build_command(self, *, python: str | None = None) -> list[str]:
        code = (
            "import json,sys; "
            "p=sys.argv[sys.argv.index('--output')+1]; "
            "m=sys.argv[sys.argv.index('--msh-output')+1]; "
            f"open(p,'w',encoding='utf-8').write(json.dumps({{'metric': {self.value}}})); "
            "open(m,'w',encoding='ascii',newline='\\n').write("
            "'$MeshFormat\\n4.1 0 8\\n$EndMeshFormat\\n')"
        )
        return [
            python or sys.executable,
            "-c",
            code,
            "--msh-output",
            "ignored.msh",
            "--output",
            "ignored.json",
        ]


@pytest.fixture
def echo_application(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setitem(
        app.APPLICATIONS,
        "echo",
        app.ApplicationDefinition(__name__, "EchoDesignSpec", ("metric",)),
    )


def test_application_runner_writes_stable_artifacts(tmp_path: Path, echo_application):
    config = tmp_path / "config.json"
    config.write_text(
        json.dumps({
            "schema": app.CONFIG_SCHEMA,
            "settings": {"input_file": "model.vol", "value": 7.25},
        }),
        encoding="utf-8",
    )

    payload = app.run_application("echo", config, tmp_path / "run", timeout_s=30)

    result = payload["radia_result"]
    assert result["schema"] == app.RESULT_SCHEMA
    assert result["backend"] == "python-headless-cli"
    assert result["status"] == "passed"
    assert result["primary"] == {"key": "metric", "value": 7.25}
    assert (tmp_path / "run" / "run.log").is_file()
    assert (tmp_path / "run" / "command.txt").is_file()
    assert (tmp_path / "run" / "solver_result.json").is_file()
    gmsh_path = (tmp_path / "run" / "echo_fields.msh").resolve()
    assert gmsh_path.is_file()
    assert result["artifacts"] == {
        "gmsh_policy": "required",
        "gmsh_format": "msh-v4.1",
        "gmsh_primary": str(gmsh_path),
        "gmsh": [str(gmsh_path)],
    }
    command = (tmp_path / "run" / "command.txt").read_text(encoding="utf-8")
    assert str(gmsh_path) in command
    on_disk = json.loads((tmp_path / "run" / "result.json").read_text(encoding="utf-8"))
    assert on_disk == payload


def test_application_runner_records_missing_inputs(tmp_path: Path, echo_application):
    config = tmp_path / "config.json"
    config.write_text("{}", encoding="utf-8")
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "command.txt").write_text("stale command\n", encoding="utf-8")
    (run_dir / "solver_result.json").write_text(
        '{"metric": 99}\n', encoding="utf-8"
    )
    (run_dir / "echo_fields.msh").write_text(
        "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n", encoding="ascii"
    )
    (run_dir / "old_auxiliary.msh").write_text(
        "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n", encoding="ascii"
    )

    payload = app.run_application("echo", config, run_dir, timeout_s=5)

    result = payload["radia_result"]
    assert result["status"] == "failed"
    assert "Missing required inputs: Input file" in result["error"]
    assert result["artifacts"]["gmsh"] == []
    assert not (run_dir / "command.txt").exists()
    assert not (run_dir / "solver_result.json").exists()
    assert not (run_dir / "echo_fields.msh").exists()
    assert not (run_dir / "old_auxiliary.msh").exists()
    assert json.loads((run_dir / "result.json").read_text())["radia_result"][
        "status"
    ] == "failed"


def test_application_cli_returns_nonzero_but_keeps_result(tmp_path: Path):
    config = tmp_path / "missing.json"
    run_dir = tmp_path / "run"

    returncode = app.main([
        "--application",
        "em",
        "--config",
        str(config),
        "--run-dir",
        str(run_dir),
    ])

    assert returncode == 1
    result = json.loads((run_dir / "result.json").read_text(encoding="utf-8"))
    assert result["radia_result"]["status"] == "failed"


@pytest.mark.parametrize(
    ("application", "settings"),
    [
        ("em", {"coil_script": "coil.py"}),
        ("pcb", {"inp": "board.inp"}),
        ("motor", {"vol": "motor.vol"}),
        ("streamfunction", {"coil_vol": "coil.vol", "eval_vol": "eval.vol"}),
        ("ih", {"peec_step": "coil.step"}),
        ("field", {"request_file": "field-study.json"}),
    ],
)
def test_production_design_specs_expose_runner_output_contract(
    tmp_path: Path,
    application: str,
    settings: dict[str, str],
):
    spec, _ = app._load_spec(application, settings)
    assert spec.missing_required_inputs() == []

    output = tmp_path / application / "solver_result.json"
    command = app._replace_output_path(
        spec.build_command(python=sys.executable),
        output,
    )

    index = command.index("--output")
    assert command[index + 1] == str(output)
