"""Public runtime provenance cannot expose machines; exact routing stays usable."""
import importlib.util
from pathlib import Path
import pytest

spec = importlib.util.spec_from_file_location("runtime_privacy", Path(__file__).parents[1] / "tools/policy_lint.py")
lint = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lint)

@pytest.mark.parametrize("value", ["INTEL11", "mdx1", "mdx2", "hibino", "192.168.121.102", "10.2.3.4", r"C:\Users\someone\data", r"W:\00_CAE\result"])
def test_identifying_result_provenance_rejected(value):
    assert lint.public_runtime_identity_hits("validation_test/results/result.json", value)

@pytest.mark.parametrize("value", ["Windows", "Linux", "threads=8 Python=3.12 NGSolve=6.2.2607", "laboratory", "Hibino2024"])
def test_non_identifying_metadata_and_citations_allowed(value):
    assert not lint.public_runtime_identity_hits("docs/runtime.md", value)

def test_only_exact_route_exception_is_allowed():
    assert not lint.public_runtime_identity_hits("tools/ci_preflight_mdx.py", "mdx1")
    assert lint.public_runtime_identity_hits("tools/new_result_writer.py", "mdx1")
    assert lint.public_runtime_identity_hits("docs/ci_preflight_mdx.py", "mdx1")

def test_host_filename_rejected_without_reading_runtime_field():
    assert lint.public_runtime_identity_hits("results/intel11-result.json", "{}")

@pytest.mark.parametrize("content", [
    r'{"host": "lab"}',
    r'{"extra": "\u0069ntel11"}',
    r'{"runtime_path": "\\\\storage\\share\\job"}',
])
def test_tracked_json_structural_scan_has_no_raw_grep_bypass(tmp_path, monkeypatch, content):
    target = tmp_path / "result.json"
    target.write_text(content, encoding="utf-8")
    monkeypatch.setattr(lint, "REPO", str(tmp_path))
    monkeypatch.setattr(lint, "_git_ls_files", lambda *args: ["result.json"])
    monkeypatch.setattr(lint, "_git_grep", lambda *args, **kwargs: [])
    assert lint.check_public_runtime_privacy()

@pytest.mark.parametrize("content", [
    r'{"platform_class": "Windows", "threads": 8}',
    r'{"description": "laboratory lab_style labels collaborate"}',
    r'{"formula": "\\rho\\beta"}',
])
def test_structural_scan_keeps_non_identifying_metadata_and_math(tmp_path, monkeypatch, content):
    (tmp_path / "result.json").write_text(content, encoding="utf-8")
    monkeypatch.setattr(lint, "REPO", str(tmp_path))
    monkeypatch.setattr(lint, "_git_ls_files", lambda *args: ["result.json"])
    monkeypatch.setattr(lint, "_git_grep", lambda *args, **kwargs: [])
    assert not lint.check_public_runtime_privacy()

def test_every_exception_is_an_exact_file_with_a_reviewable_reason():
    assert lint.PUBLIC_RUNTIME_OPERATIONAL_ALLOW == frozenset(lint.PUBLIC_RUNTIME_OPERATIONAL_REASONS)
    for path, reason in lint.PUBLIC_RUNTIME_OPERATIONAL_REASONS.items():
        assert not any(character in path for character in "*?[]")
        assert not path.endswith("/")
        assert len(reason) > 20


@pytest.mark.parametrize("path", ["results/lab-record.json", "validation_test/results/mdx_timing.json"])
def test_host_role_filename_without_metadata_is_rejected(path):
    assert lint.public_runtime_identity_hits(path, "{}")

def test_scientific_lab_identifier_filename_is_preserved():
    assert not lint.public_runtime_identity_hits("src/figure/lab_style.py", "# laboratory")
