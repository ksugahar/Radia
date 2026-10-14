"""Error-kind serialization boundary for the build123d MCP server."""

import json

from radia_mcp.build123d.server import _dumps, _with_kind


def test_error_kind_marks_any_missing_python_module_as_environment():
    payload = {"status": "error",
               "error": "ModuleNotFoundError: No module named 'meshio'"}

    classified = _with_kind(payload)

    assert classified["kind"] == "environment"
    assert "kind" not in payload


def test_error_kind_preserves_explicit_and_non_error_statuses():
    explicit = {"status": "error", "kind": "internal", "error": "boom"}
    invalid = {"status": "invalid_input", "error": "bad JSON"}

    assert _with_kind(explicit) is explicit
    assert _with_kind(invalid) is invalid


def test_dumps_adds_input_kind_to_an_ordinary_error():
    encoded = _dumps({"status": "error", "error": "file not found"})

    assert json.loads(encoded)["kind"] == "input"
