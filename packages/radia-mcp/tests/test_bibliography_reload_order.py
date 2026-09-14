"""New parent helpers must not be imported eagerly during child-first reload."""
from pathlib import Path
import importlib
import pytest


@pytest.mark.parametrize("name", ["T13_check_surname_braces", "T9_self_citation_ratio"])
def test_child_can_reload_before_parent_gains_new_name_helpers(monkeypatch, name):
    from radia_mcp.bibliography import _bibparse
    module = importlib.import_module("radia_mcp.bibliography.plans." + name)
    namespace = {"__name__": module.__name__, "__package__": module.__package__}
    with monkeypatch.context() as patch:
        patch.delattr(_bibparse, "_split_name_parts")
        patch.delattr(_bibparse, "first_author_family")
        exec(compile(Path(module.__file__).read_bytes(), module.__file__, "exec"), namespace)
    if name.startswith("T13"):
        assert namespace["_wrap_surname"]("Jane Doe") == "Jane {Doe}"
    else:
        assert namespace["_entry_authors"]("Doe, Jane AND Roe, John") == ["doe", "roe"]
