"""Source contract for the independently owned native Cubit plugin."""
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
PLUGIN = ROOT / "src" / "cubit_plugin"


def test_native_plugin_has_one_bdf_command_and_product_owned_identity():
    registry = (PLUGIN / "CubitMeshExportPlugin.cpp").read_text(encoding="utf-8")
    header = (PLUGIN / "ExportNastranCommand.hpp").read_text(encoding="utf-8")
    implementation = (PLUGIN / "ExportNastranCommand.cpp").read_text(
        encoding="utf-8"
    )

    assert registry.count('keys.push_back("ExportNastranCommand")') == 1
    assert header.count("class ExportNastranCommand") == 1
    assert implementation.count('"export nastran_bdf ') == 2
    assert "CubitMeshExportPlugin" in registry
    assert "CoilCommand" not in registry
    assert not (PLUGIN / "CoilCommand.cpp").exists()
    assert not (PLUGIN / "CoilCommand.hpp").exists()


def test_native_sources_use_product_owned_helpers():
    cmake = (PLUGIN / "CMakeLists.txt").read_text(encoding="utf-8")
    message_filter = (PLUGIN / "LearnEditionMessageFilter.hpp").read_text(
        encoding="utf-8"
    )

    assert "project(CubitMeshExportPlugin)" in cmake
    assert "namespace cubit_mesh_export" in message_filter
    assert "RADIA_NETGEN_" not in cmake

    native_source = "\n".join(
        path.read_text(encoding="utf-8", errors="replace")
        for pattern in ("*.cpp", "*.hpp")
        for path in PLUGIN.glob(pattern)
    )
    assert "RADIA_" not in native_source
    assert "radia_kelvin_config" not in native_source
