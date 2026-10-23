"""Local operations data stays external and retrieval never runs commands."""
import json

import pytest

from radia_mcp.meta import operations


def test_local_registration_lists_and_reads_provenance(tmp_path, monkeypatch):
    path = tmp_path / "knowledge.json"
    monkeypatch.setenv("RADIA_OPERATIONS_KNOWLEDGE", str(path))
    assert operations.lookup()["status"] == "not_configured"
    entry = {"id": "sample-pbs", "title": "Sample scheduler", "summary": "Submission notes",
             "content": "qsub job.pbs (reference only)",
             "sources": [{"path": "private-manual.pdf", "pdf_pages": [23, 42], "sha256": "a" * 64}]}
    path.write_text(json.dumps({"schema": "radia.operations-knowledge.v1", "topics": [entry]}), encoding="utf-8")
    listed = operations.lookup()
    assert listed["topics"] == [{k: entry[k] for k in ("id", "title", "summary")}]
    assert operations.lookup("sample-pbs")["topic"] == entry
    assert operations.lookup("../private-manual.pdf")["status"] == "unknown_topic"
    path.write_text(json.dumps({"schema": "radia.operations-knowledge.v1", "topics": [entry, entry]}), encoding="utf-8")
    with pytest.raises(ValueError, match="Duplicate"):
        operations.lookup()


def test_corrupt_registration_fails_loudly(tmp_path, monkeypatch):
    path = tmp_path / "knowledge.json"
    monkeypatch.setenv("RADIA_OPERATIONS_KNOWLEDGE", str(path))
    path.write_text('{"schema":"wrong","topics":[]}', encoding="utf-8")
    with pytest.raises(ValueError, match="schema"):
        operations.lookup()
