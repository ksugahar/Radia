"""Cubit probe preserves the shared cross-CAD entity vocabulary."""
from radia_mcp.common.server_hardening import PROBE_FACE_CORE_KEYS, PROBE_SOLID_CORE_KEYS


def test_cubit_probe_entities_carries_core_keys():
    from radia_mcp.cubit import probe_ops
    import importlib.util
    from pathlib import Path
    spec = importlib.util.spec_from_file_location(
        "_cubit_probe_ops_helpers", Path(__file__).with_name("test_cubit_probe_ops.py"))
    helpers = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helpers)
    _MockCubit = helpers._MockCubit

    ent = probe_ops.op_probe(_MockCubit(), ["entities"])
    assert PROBE_SOLID_CORE_KEYS <= set(ent["volumes"][0])
    assert PROBE_FACE_CORE_KEYS <= set(ent["surfaces"][0])
