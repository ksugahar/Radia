"""Generated Cubit artifacts must not escape or overwrite user files."""

import json
import threading

from cubit_mesh_export.mcp import server


def test_scaffold_rejects_path_components_and_existing_destination(tmp_path):
    buttons = json.dumps([{"label": "Mesh", "script_basename": "mesh"}])
    result = server.cubit_scaffold_toolbar("../outside", buttons, str(tmp_path))
    assert result.startswith("ERROR:")
    assert not (tmp_path.parent / "outside").exists()

    bad_button = json.dumps([{"label": "Mesh", "script_basename": "../outside"}])
    result = server.cubit_scaffold_toolbar("safe", bad_button, str(tmp_path))
    assert result.startswith("ERROR:")
    assert not (tmp_path / "safe").exists()

    result = server.cubit_scaffold_toolbar("safe", buttons, str(tmp_path))
    assert result.startswith("Toolbar scaffold written")
    xml = tmp_path / "safe" / "toolbar.xml"
    before = xml.read_bytes()
    result = server.cubit_scaffold_toolbar("safe", buttons, str(tmp_path))
    assert "refusing to overwrite" in result
    assert xml.read_bytes() == before


def test_journal_requires_new_jou_destination(tmp_path, monkeypatch):
    class FakeSession:
        _lock = threading.RLock()
        _command_history = []

        def native_journal_snapshot(self):
            return {"journal": "reset\n", "paths": [],
                    "generation_count": 1, "errors": []}

    monkeypatch.setattr(server._cs, "_SINGLETON", FakeSession())
    wrong = json.loads(server.cubit_session_journal(str(tmp_path / "journal.txt")))
    assert wrong["status"] == "error"
    assert not (tmp_path / "journal.txt").exists()

    target = tmp_path / "journal.jou"
    target.write_text("user-owned\n", encoding="utf-8")
    existing = json.loads(server.cubit_session_journal(str(target)))
    assert existing["status"] == "error"
    assert target.read_text(encoding="utf-8") == "user-owned\n"

    fresh = tmp_path / "fresh.jou"
    created = json.loads(server.cubit_session_journal(str(fresh)))
    assert created["status"] == "ok"
    assert fresh.read_text(encoding="utf-8") == "reset\n"
