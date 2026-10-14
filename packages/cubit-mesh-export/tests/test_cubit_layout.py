from cubit_mesh_export._layout import has_python_binding
import pytest


@pytest.mark.parametrize("init_present,native_present,expected", [
    (False, False, False), (True, False, False),
    (False, True, False), (True, True, True),
])
def test_packaged_binding_requires_both_files(tmp_path, init_present, native_present, expected):
    package = tmp_path / "cubit"
    package.mkdir()
    if init_present:
        (package / "__init__.py").touch()
    if native_present:
        (package / "_cubit3.pyd").touch()
    assert has_python_binding(tmp_path) is expected
