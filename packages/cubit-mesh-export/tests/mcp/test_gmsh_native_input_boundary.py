"""Guards: oversized field tags and indirect formats cannot reach native Gmsh.

Assumes: mesh-only quality APIs take direct ASCII MSH 4.2. Native access is
replaced by a process marker so large tags never allocate native field arrays;
the independent inventory reader and mesh-only size_t tags stay unrestricted.
"""
import json
from pathlib import Path

import pytest

from cubit_mesh_export import _gmsh_runtime as runtime
from cubit_mesh_export.mcp.gmsh_v42 import summarize_gmsh_v42_ascii
from cubit_mesh_export.mcp._support import _gmsh_subprocess as child


def _mesh(tag):
    return (f"$MeshFormat\n4.2 0 8\n$EndMeshFormat\n"
            f"$Nodes\n1 1 {tag} {tag}\n0 1 0 1\n{tag}\n0 0 0\n$EndNodes\n"
            f"$Elements\n1 1 {tag} {tag}\n0 1 15 1\n{tag} {tag}\n$EndElements\n")


def _field(kind, tag, maximum=None):
    # MSH 4.2 stores field tags before block values, not interleaved rows.
    block = "15 1 1" if kind == "ElementNodeData" else "1"
    return (f"${kind}\n1\n\"value\"\n1\n0\n3\n0\n1\n1\n"
            f"1 1 {tag} {tag if maximum is None else maximum}\n{block}\n{tag}\n1\n$End{kind}\n")


def _operation(monkeypatch, tmp_path, path):
    native = tmp_path / "native-entered"
    operation = tmp_path / "operation-entered"
    original = child.subprocess.run

    def marked_run(command, **kwargs):
        script_path = Path(command[1])
        text = script_path.read_text(encoding="utf-8")
        # Instrument only readiness/native access; execute the real boundary in
        # the real child process, without ingesting a large field into Gmsh.
        text = text.replace("    runtime = probe_gmsh5()", (
            f"    Path({str(native)!r}).write_text('entered', encoding='utf-8')\n"
            "    runtime = {'native_access': 'marker only'}"))
        script_path.write_text(text, encoding="utf-8")
        return original(command, **kwargs)

    monkeypatch.setattr(child.subprocess, "run", marked_run)
    script = (f"Path({str(operation)!r}).write_text('entered', encoding='utf-8')\n"
              "Path(sys.argv[-1]).write_text(json.dumps({'ok': True}), encoding='utf-8')")
    result = child.run_gmsh_json_subprocess(script, [str(path)], timeout_s=30)
    return result, native.exists(), operation.exists()


@pytest.mark.parametrize("kind", ["NodeData", "ElementData", "ElementNodeData"])
@pytest.mark.parametrize("tag", [2**28, 2147483647, 2147483648, 2147483665])
def test_native_boundary_checks_each_field_before_any_native_call(monkeypatch, tmp_path, kind, tag):
    path = tmp_path / "input.msh"
    text = _mesh(tag) + _field(kind, tag)
    path.write_text(text, encoding="utf-8")
    assert summarize_gmsh_v42_ascii(text)["node_count"] == 1
    result, native, operation = _operation(monkeypatch, tmp_path, path)
    allowed = tag <= 2147483647
    assert result["ok"] is allowed, result
    assert native is allowed and operation is allowed, (native, operation, result)
    if not allowed:
        assert str(tag) in result["error"] and "2147483647" in result["error"], result
        assert "field tag" in result["error"], result


@pytest.mark.parametrize("kind", ["NodeData", "ElementData", "ElementNodeData"])
def test_actual_field_tags_are_checked_even_if_maximum_header_lies(tmp_path, kind):
    path = tmp_path / "input.msh"
    path.write_text(_mesh(1) + _field(kind, 2147483648, maximum=1), encoding="utf-8")
    with pytest.raises(ValueError, match="field tag 2147483648 exceeds 2147483647"):
        runtime.validate_native_msh_fields(path)


def test_mesh_only_size_t_tags_are_not_capped(monkeypatch, tmp_path):
    path = tmp_path / "input.msh"
    text = _mesh(2**40)
    path.write_text(text, encoding="utf-8")
    assert summarize_gmsh_v42_ascii(text)["node_count"] == 1
    result, native, operation = _operation(monkeypatch, tmp_path, path)
    assert result["ok"] and native and operation, result


@pytest.mark.parametrize("suffix,contents", [
    (".geo", 'Merge "legacy.msh";'),
    (".geo", 'Include "nested.geo";'),
    (".geo", 'Merge StrCat("prefix", name);'),
    (".msh", 'Merge "legacy.msh";'),
    (".txt", '$MeshFormat\n4.2 0 8\n$EndMeshFormat\n'),
])
def test_indirect_or_disguised_input_is_refused_before_native_access(monkeypatch, tmp_path, suffix, contents):
    path = tmp_path / ("input" + suffix)
    path.write_text(contents, encoding="utf-8")
    result, native, operation = _operation(monkeypatch, tmp_path, path)
    assert result["ok"] is False and not native and not operation, result
    assert "MSH 4.2" in result["error"], result


@pytest.mark.parametrize("suffix,header", [(".bdf", "$ producer reference\nSOL 101\nCEND\nBEGIN BULK\n"), (".vtk", "# vtk DataFile Version 3.0\nreference\nASCII\n")])
def test_bdf_vtk_are_explicit_independent_reference_formats_only(tmp_path, suffix, header):
    path = tmp_path / ("reference" + suffix)
    path.write_text(header, encoding="utf-8")
    with pytest.raises(ValueError, match="MSH 4.2"):
        runtime.preflight_native_mesh_input(path)
    runtime.preflight_native_mesh_input(path, allow_bdf_vtk=True)


@pytest.mark.parametrize("suffix", [".bdf", ".vtk"])
def test_reference_suffix_does_not_hide_an_oversized_msh_field(tmp_path, suffix):
    path = tmp_path / ("reference" + suffix)
    path.write_text(_mesh(2147483648) + _field("NodeData", 2147483648), encoding="utf-8")
    with pytest.raises(ValueError, match="field tag 2147483648 exceeds 2147483647"):
        runtime.preflight_native_mesh_input(path, allow_bdf_vtk=True)


@pytest.mark.parametrize("suffix", [".geo", ".bdf", ".vtk"])
def test_reference_mode_does_not_execute_a_geo_program(tmp_path, suffix):
    path = tmp_path / ("reference" + suffix)
    path.write_text('Merge "mesh.msh";', encoding="utf-8")
    with pytest.raises(ValueError, match="GEO"):
        runtime.preflight_native_mesh_input(path, allow_bdf_vtk=True)
