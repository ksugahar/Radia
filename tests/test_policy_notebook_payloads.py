"""Policy lint must treat encoded notebook payloads as bytes, never as prose."""
import base64
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("policy_lint", ROOT / "tools" / "policy_lint.py")
LINT = importlib.util.module_from_spec(spec)
spec.loader.exec_module(LINT)

# Valid base64 whose text happens to contain flagged words: PNG signature plus
# one byte (12 characters, no padding), then 4-character groups.
# Words are split so this file does not itself trip the lint it tests.
MDX_WORD, LAB_WORD, INTEL_WORD = "md" + "x", "LA" + "B", "intel" + "11"
PNG_TEXT = base64.b64encode(b"\x89PNG\r\n\x1a\n\x00").decode() + "/" + MDX_WORD + "/" + LAB_WORD + INTEL_WORD + "A"


def notebook(png_text, points=("AAAA/" + MDX_WORD,), source="plain prose"):
    return json.dumps({
        "cells": [{"cell_type": "code", "source": [source], "outputs": [
            {"output_type": "display_data", "data": {"image/png": png_text}}]}],
        "metadata": {"widgets": {"application/vnd.jupyter.widget-state+json": {"state": {
            "m1": {"state": {"value": {"points3d": list(points)}}}}}}},
        "nbformat": 4, "nbformat_minor": 5})


def test_accidental_flagged_words_in_payloads_are_not_privacy_hits():
    content = notebook(PNG_TEXT)
    assert base64.b64decode(PNG_TEXT, validate=True).startswith(b"\x89PNG")
    assert LINT.public_runtime_identity_hits("docs/demo.ipynb", content) == []
    assert LINT.check_notebook_binary_integrity("docs/demo.ipynb", content) == []


def test_flagged_words_in_notebook_prose_are_still_hits():
    content = notebook(PNG_TEXT, source="measured on " + LAB_WORD)
    assert LINT.public_runtime_identity_hits("docs/demo.ipynb", content)


def test_text_substitution_inside_payloads_fails_integrity():
    broken = notebook(PNG_TEXT.replace(LAB_WORD, "worker_a"), points=("AAAA/worker_a",))
    hits = LINT.check_notebook_binary_integrity("docs/demo.ipynb", broken)
    assert any("image/png: invalid base64" in hit for hit in hits)
    assert any("points3d/0: invalid base64" in hit for hit in hits)


def test_image_payload_must_keep_its_signature():
    not_png = base64.b64encode(b"GIF89a").decode()
    hits = LINT.check_notebook_binary_integrity("docs/demo.ipynb", notebook(not_png))
    assert any("lost its signature" in hit for hit in hits)


def test_bare_mdx_count_ignores_payload_lines(tmp_path, monkeypatch):
    (tmp_path / "nb.ipynb").write_text(json.dumps(json.loads(notebook(PNG_TEXT)), indent=1))
    monkeypatch.setattr(LINT, "REPO", str(tmp_path))
    assert LINT._bare_mdx_notebook_lines("nb.ipynb") == 0
    (tmp_path / "nb.ipynb").write_text(json.dumps(json.loads(notebook(PNG_TEXT, source="ran on " + MDX_WORD)), indent=1))
    assert LINT._bare_mdx_notebook_lines("nb.ipynb") == 1


def test_public_service_allowlist_is_exact_path_and_context():
    path = "packages/radia-mcp/src/radia_mcp/grant_writing/tools.py"
    lines = (ROOT / path).read_text(encoding="utf-8").splitlines()
    line = next(line for line in lines if '"shared_compute":' in line)
    assert LINT.bare_private_host_lines(path, line) == []
    assert LINT.bare_private_host_lines("other.py", line)
    assert LINT.bare_private_host_lines(path, line + " # host " + MDX_WORD)


def test_strict_private_host_scan_does_not_need_a_baseline():
    assert LINT.bare_private_host_lines("ordinary.md", "run on " + MDX_WORD)
    assert LINT.bare_private_host_lines("docs/demo.ipynb", notebook(PNG_TEXT)) == []
    assert LINT.bare_private_host_lines("docs/demo.ipynb", notebook(PNG_TEXT, source="run on " + MDX_WORD))
