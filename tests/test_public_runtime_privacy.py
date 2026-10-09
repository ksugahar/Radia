"""Public runtime provenance cannot expose machines; exact routing stays usable."""
import importlib.util
from pathlib import Path
import pytest

spec = importlib.util.spec_from_file_location("runtime_privacy", Path(__file__).parents[1] / "tools/policy_lint.py")
lint = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lint)

@pytest.mark.parametrize("value", ["INTEL11", "mdx1", "mdx2", "hibino", "192.168.254.254", "10.2.3.4", r"C:\Users\someone\data", r"W:\synthetic-project\result"])
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


def test_bare_private_scratch_and_user_roots_are_identifying():
    for value in (r"C:\temp", r"C:\Users"):
        assert lint.public_runtime_identity_hits("result.json", __import__("json").dumps({"note": value}))


def test_citation_exception_is_exact_and_does_not_hide_other_fields():
    path = "packages/radia-mcp/src/radia_mcp/bibliography/data/references.bib"
    citation = "@misc{hibino2026aca,\n author = {Hibino, Yoshihiko and Other, Author},\n title = {Example},\n year = {2026}\n}"
    assert not lint.public_runtime_identity_hits(path, citation)
    assert lint.public_runtime_identity_hits(path, citation.replace("Example", "192.168.254.254"))
    assert lint.public_runtime_identity_hits(path, citation.replace("hibino2026aca", "other-key"))

def test_numbered_host_name_rejected():
    assert lint.public_runtime_identity_hits("result.txt", "100\u53f7\u6a5f")


def test_lowercase_runtime_machine_field_is_rejected():
    assert lint.public_runtime_identity_hits("result.json", '{"machine":"lab"}')
    assert not lint.public_runtime_identity_hits("result.json", '{"machine":"8-pole PMSM benchmark"}')


@pytest.mark.parametrize("host", ["lab", "worker-one", "&#108;ab"])
def test_xml_host_attribute_is_checked_by_full_pipeline(tmp_path, monkeypatch, host):
    (tmp_path / "focused.xml").write_text(f'<testsuite hostname="{host}" tests="1"/>', encoding="utf-8")
    monkeypatch.setattr(lint, "REPO", str(tmp_path))
    monkeypatch.setattr(lint, "_git_ls_files", lambda *args: ["focused.xml"])
    monkeypatch.setattr(lint, "_git_grep", lambda *args, **kwargs: [])
    assert lint.check_public_runtime_privacy()

def test_xml_without_hostname_preserves_numeric_attributes():
    assert not lint.public_runtime_identity_hits("focused.xml", '<testsuite tests="7" time="9.195"/>')
