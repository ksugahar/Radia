"""Guards: Gmsh 4/native mismatches and old mesh formats cannot pass quality gates.

Assumes: official Gmsh 5 is installed in the caller's isolated environment;
its pip install report and retained original wheel bind the loaded bytes. The real
Gmsh API roundtrip is interoperability evidence, not numerical accuracy proof.
"""
import gmsh
import pytest

from cubit_mesh_export import _gmsh_runtime as runtime
from cubit_mesh_export.mcp._support._gmsh_subprocess import run_gmsh_json_subprocess
from cubit_mesh_export.mcp.gmsh_v42 import ELEMENT_TYPES


def test_loaded_gmsh5_roundtrips_msh42_without_leaving_a_model():
    assert not gmsh.isInitialized()
    result = runtime.probe_gmsh5()
    assert result["distribution_version"].startswith("5.")
    assert result["native_version"].startswith("5.")
    assert result["msh42_roundtrip"] is True
    assert result["msh_format"] == "4.2"
    assert len(result["native_sha256"]) == len(result["wrapper_sha256"]) == 64
    assert len(result["wheel_sha256"]) == 64
    assert not gmsh.isInitialized()


def test_major4_is_rejected_before_initializing_native_state(monkeypatch):
    monkeypatch.setattr(runtime.metadata, "version", lambda name: "4.99.0")
    with pytest.raises(RuntimeError, match="major 5 required"):
        runtime.probe_gmsh5()
    assert not gmsh.isInitialized()


def test_missing_wheel_identity_is_not_silently_reconstructed(monkeypatch):
    monkeypatch.delenv("CUBIT_GMSH_INSTALL_REPORT", raising=False)
    with pytest.raises(RuntimeError, match="wheel identity missing"):
        runtime.probe_gmsh5()
    assert not gmsh.isInitialized()


def test_original_wheel_is_required_instead_of_trusting_record(monkeypatch):
    monkeypatch.delenv("CUBIT_GMSH_WHEEL", raising=False)
    with pytest.raises(RuntimeError, match="original wheel missing"):
        runtime.probe_gmsh5()
    assert not gmsh.isInitialized()


def test_changed_original_wheel_is_rejected_before_native_initialization(monkeypatch, tmp_path):
    wheel = tmp_path / "changed.whl"
    wheel.write_bytes(b"changed archive")
    monkeypatch.setenv("CUBIT_GMSH_WHEEL", str(wheel))
    with pytest.raises(RuntimeError, match="SHA256 differs"):
        runtime.probe_gmsh5()
    assert not gmsh.isInitialized()


@pytest.mark.parametrize("header", ["4.3 0 8", "4.2 1 8", "4.2 0 4", "invalid"])
def test_unknown_binary_or_malformed_mesh_never_reaches_the_operation(tmp_path, header):
    path = tmp_path / "input.msh"
    marker = tmp_path / "operation-executed"
    path.write_text(f"$MeshFormat\n{header}\n$EndMeshFormat\n", encoding="ascii")
    script = f"""import json,sys
from pathlib import Path
Path({str(marker)!r}).write_text('executed')
Path(sys.argv[-1]).write_text(json.dumps({{'ok': True}}))
"""
    result = run_gmsh_json_subprocess(script, [str(path)], timeout_s=60)
    assert result["ok"] is False
    assert "MSH 4.2" in result["error"]
    assert not marker.exists()


def test_child_result_keeps_runtime_identity():
    script = "import json,sys;open(sys.argv[-1], 'w').write(json.dumps({'ok': True}))"
    result = run_gmsh_json_subprocess(script, [], timeout_s=60)
    assert result["ok"] is True
    assert result["gmsh_runtime"]["native_version"].startswith("5.")
    assert result["gmsh_runtime"]["msh42_roundtrip"] is True


def test_inventory_connectivity_contract_matches_gmsh5_element_api():
    runtime.probe_gmsh5()
    gmsh.initialize(["-noconfig"])
    try:
        for kind, (family, order, nodes) in ELEMENT_TYPES.items():
            name, dimension, actual_order, actual_nodes, coordinates, primary = (
                gmsh.model.mesh.getElementProperties(kind))
            assert (actual_order, actual_nodes) == (order, nodes), (kind, name, family)
    finally:
        gmsh.finalize()
