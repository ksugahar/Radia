"""Guards: external-policy integrity and scan surfaces; assumes reviewed private criteria."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import zipfile

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture
def guard(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location("public_policy_contract", ROOT / "tools/check_public_names.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    raw = json.dumps({"schema": "radia.private-public-name-rules.v1",
                      "flags": ["IGNORECASE"], "patterns": [r"blocked[ _-]*token"]}).encode()
    path = tmp_path / "criteria.json"
    path.write_bytes(raw)
    monkeypatch.setattr(module, "APPROVED_RULES_SHA256", hashlib.sha256(raw).hexdigest())
    monkeypatch.setenv("RADIA_PUBLIC_NAME_RULES_FILE", str(path))
    monkeypatch.delenv("RADIA_PUBLIC_NAME_RULES_JSON", raising=False)
    module.FORBIDDEN = module.load_rules()
    return module, path


def test_missing_ambiguous_or_changed_criteria_fail_closed(guard, monkeypatch):
    module, path = guard
    monkeypatch.delenv("RADIA_PUBLIC_NAME_RULES_FILE")
    with pytest.raises(module.PolicyInputError, match="Exactly one"):
        module.load_rules()
    monkeypatch.setenv("RADIA_PUBLIC_NAME_RULES_FILE", str(path))
    monkeypatch.setenv("RADIA_PUBLIC_NAME_RULES_JSON", path.read_text())
    with pytest.raises(module.PolicyInputError, match="Exactly one"):
        module.load_rules()
    monkeypatch.delenv("RADIA_PUBLIC_NAME_RULES_JSON")
    path.write_bytes(b"{}")
    with pytest.raises(module.PolicyInputError, match="fingerprint"):
        module.load_rules()


@pytest.mark.parametrize("raw", [
    b'{"schema":"wrong","flags":["IGNORECASE"],"patterns":["x"]}',
    b'{"schema":"radia.private-public-name-rules.v1","flags":[],"patterns":["x"]}',
    b'{"schema":"radia.private-public-name-rules.v1","flags":["IGNORECASE"],"patterns":[]}',
    b'{"schema":"radia.private-public-name-rules.v1","flags":["IGNORECASE"],"patterns":[""]}',
    b'{"schema":"radia.private-public-name-rules.v1","flags":["IGNORECASE"],"patterns":["["]}',
    b'{"schema":"a","schema":"b","flags":["IGNORECASE"],"patterns":["x"]}',
])
def test_malformed_schema_flags_and_patterns_rejected(guard, monkeypatch, raw):
    module, path = guard
    path.write_bytes(raw)
    monkeypatch.setattr(module, "APPROVED_RULES_SHA256", hashlib.sha256(raw).hexdigest())
    with pytest.raises(module.PolicyInputError):
        module.load_rules()


def test_semantic_and_encoded_binary_and_wheel_members(guard, tmp_path):
    module, _ = guard
    for text in ("BLOCKED_token", "blocked token", "blocked---token"):
        assert module.FORBIDDEN.search(text)
        for encoding in ("utf-8", "utf-16-le", "utf-16-be"):
            assert module.binary_hits(b"x" + text.encode(encoding) + b"\0\0")
    assert not module.binary_hits(b"allowed material")
    for name in ("blocked_token.exe", "package/tool.pyd", "package.dist-info/METADATA"):
        wheel = tmp_path / "artifact.whl"
        with zipfile.ZipFile(wheel, "w") as archive:
            archive.writestr(name, "blocked token".encode("utf-16-le"))
        assert module.artifact_hits(wheel)


def test_notebook_scan_keeps_semantic_text_and_excludes_image_bytes(guard):
    module, _ = guard
    doc = {"cells": [{"source": ["allowed"], "outputs": [{"text": ["blocked token"],
            "data": {"image/png": "binary payload", "text/plain": ["blocked_token"]}}]}]}
    assert list(module.notebook_text(doc)) == ["allowed", "blocked token", "blocked_token"]


def test_source_scan_and_report_are_executed(guard, tmp_path, monkeypatch):
    module, _ = guard
    (tmp_path / "source.txt").write_text("blocked token", encoding="utf-8")
    monkeypatch.setattr(module, "ROOT", tmp_path)
    def git(_command, **kwargs):
        return "a" * 40 if kwargs.get("text") else b"source.txt\0"
    monkeypatch.setattr(module.subprocess, "check_output", git)
    report = tmp_path / "report.json"
    monkeypatch.setattr(sys, "argv", ["check_public_names.py", "--report", str(report)])
    with pytest.raises(SystemExit, match="source.txt: text"):
        module.main()
    result = json.loads(report.read_text())
    assert result["status"] == "FAIL" and result["source_semantic_files"] == 1
    assert result["criteria_sha256"] == module.APPROVED_RULES_SHA256


def test_policy_workflow_does_not_expose_rules_to_forks_or_target_events():
    path = ROOT / ".github/workflows/policy-lint.yml"
    data = yaml.load(path.read_text(encoding="utf-8"), Loader=yaml.BaseLoader)
    assert set(data["on"]) == {"push", "pull_request"}
    step = next(s for s in data["jobs"]["policy-checks"]["steps"]
                if "python tools/policy_lint.py" in s.get("run", ""))
    env = step["env"]["RADIA_PUBLIC_NAME_RULES_JSON"]
    assert "github.event_name != 'pull_request'" in env
    assert "github.event.pull_request.head.repo.full_name == github.repository" in env
    assert "secrets.RADIA_PUBLIC_NAME_RULES_JSON" in env
    assert "|| ''" in env
