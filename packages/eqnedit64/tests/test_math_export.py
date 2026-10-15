"""Subprocess contract; these tests never import or run the native core."""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import subprocess
import pytest

SPEC = spec_from_file_location("equation_api_contract", Path(__file__).parents[1] / "src/eqnedit64/api.py")
api = module_from_spec(SPEC)
SPEC.loader.exec_module(api)
M = "http://schemas.openxmlformats.org/officeDocument/2006/math"
A = "http://schemas.openxmlformats.org/drawingml/2006/main"
VALID = f'<m:oMathPara xmlns:m="{M}" xmlns:a="{A}"><m:oMath><m:r><a:rPr/><m:t>x</m:t></m:r></m:oMath></m:oMathPara>'


def test_export_checked_atomic_output_and_source_cleanup(tmp_path, monkeypatch):
    calls = []
    def invoke(arguments, timeout):
        source, output = map(Path, arguments)
        calls.append((source, output))
        assert source.read_text("utf-8") == "x"
        assert output.suffix == ".pptxmath" and not output.exists()
        output.write_text(VALID, "utf-8")
        return subprocess.CompletedProcess(arguments, 0, "", "")
    monkeypatch.setattr(api, "_invoke", invoke)
    target = tmp_path / "math.pptxmath"
    target.write_text("previous", "utf-8")
    assert api.export_powerpoint_math("x", target) == target
    assert target.read_text("utf-8") == VALID
    assert all(not path.exists() for call in calls for path in call)


@pytest.mark.parametrize("payload,status", [(None, 1), (None, 0), ("<bad/>", 0),
    (VALID.replace("m:t", "a:t"), 0), (VALID.replace(A, "urn:unsupported"), 0)])
def test_rejected_export_preserves_destination(tmp_path, monkeypatch, payload, status):
    def invoke(arguments, timeout):
        if payload:
            Path(arguments[1]).write_text(payload, "utf-8")
        return subprocess.CompletedProcess(arguments, status, "", "unsupported TeX")
    monkeypatch.setattr(api, "_invoke", invoke)
    target = tmp_path / "math.pptxmath"
    target.write_text("previous", "utf-8")
    with pytest.raises(RuntimeError):
        api.export_powerpoint_math(r"\unsupported{x}", target)
    assert target.read_text("utf-8") == "previous"


def test_export_rejects_empty_tex_before_backend(tmp_path):
    with pytest.raises(ValueError, match="non-empty"):
        api.export_powerpoint_math(" ", tmp_path / "math.pptxmath")
