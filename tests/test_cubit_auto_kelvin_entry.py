import json
import runpy
import sys
import types
from pathlib import Path


AUTO_KELVIN_ENTRY = (
    Path(__file__).parents[1]
    / "packages"
    / "cubit-mesh-export"
    / "src"
    / "cubit_mesh_export"
    / "cubit_helpers"
    / "auto_kelvin_entry.py"
)


def test_auto_kelvin_loads_deployed_sibling_despite_stale_module(tmp_path, monkeypatch):
    entry_dir = tmp_path / "stale_entry_location"
    entry_dir.mkdir()
    entry = entry_dir / "auto_kelvin_entry.py"
    entry.write_text(AUTO_KELVIN_ENTRY.read_text(encoding="utf-8"), encoding="utf-8")
    helper_dir = tmp_path / "deployed_helpers"
    helper_dir.mkdir()
    marker = helper_dir / "called.json"
    (helper_dir / "add_kelvin.py").write_text(
        "import json, os\n"
        "def auto_add_kelvin_from_current_model(**kwargs):\n"
        "    with open(os.environ['AUTO_KELVIN_TEST_MARKER'], 'w', encoding='utf-8') as f:\n"
        "        json.dump(kwargs, f, sort_keys=True)\n"
        "    return {'center': (1.0, 2.0, 3.0), 'R': 4.0, 'symmetry': 'none'}\n",
        encoding="utf-8",
    )
    config = tmp_path / "config.json"
    config.write_text(
        json.dumps({"kelvin_air_block": "air", "kelvin_block_name": "kelvin"}),
        encoding="utf-8",
    )

    stale = types.ModuleType("add_kelvin")
    stale.auto_add_kelvin_from_current_model = lambda **_: (_ for _ in ()).throw(
        AssertionError("stale add_kelvin module was reused"))
    monkeypatch.setitem(sys.modules, "add_kelvin", stale)
    monkeypatch.setenv("CUBIT_HELPERS_DIR", str(helper_dir))
    monkeypatch.setenv("RADIA_LAUNCHER_CONFIG", str(config))
    monkeypatch.setenv("AUTO_KELVIN_TEST_MARKER", str(marker))

    runpy.run_path(str(entry), run_name="__main__")

    assert json.loads(marker.read_text(encoding="utf-8")) == {
        "air_block": "air",
        "kelvin_block": "kelvin",
        "mesh_size": None,
        "reduction": None,
    }
