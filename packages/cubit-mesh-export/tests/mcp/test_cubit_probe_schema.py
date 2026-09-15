"""Cubit probe preserves the shared cross-CAD entity vocabulary."""
from cubit_mesh_export.mcp._support.server_hardening import PROBE_FACE_CORE_KEYS, PROBE_SOLID_CORE_KEYS


def test_cubit_probe_entities_carries_core_keys():
    from cubit_mesh_export.mcp import probe_ops
    from test_cubit_probe_ops import _MockCubit

    ent = probe_ops.op_probe(_MockCubit(), ["entities"])
    assert PROBE_SOLID_CORE_KEYS <= set(ent["volumes"][0])
    assert PROBE_FACE_CORE_KEYS <= set(ent["surfaces"][0])
