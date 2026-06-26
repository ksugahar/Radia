"""Basic Coreform Cubit 2025.12 plugin-shape tests."""

from tests.cubit.cubit_202512_helpers import add_cubit_mesh_curver_to_path


def test_cubit_mesh_curver_low_level_module_shape():
    """The .pyd exposes the internal builder used by the .ccm plugin."""
    import netgen  # noqa: F401 - import first so its DLLs are registered

    add_cubit_mesh_curver_to_path()
    import cubit_mesh_curver

    assert callable(cubit_mesh_curver.build_curved_mesh)
    assert cubit_mesh_curver.has_netgen is True
    assert not hasattr(cubit_mesh_curver, "extract_curved_mesh")
    assert not hasattr(cubit_mesh_curver, "extract_mesh_data")
