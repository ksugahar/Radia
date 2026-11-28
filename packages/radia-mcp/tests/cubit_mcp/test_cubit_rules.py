"""Exporter-owned tests for Cubit lint rules.

Each test provides intentionally bad code and verifies the rule fires,
then provides correct code and verifies no findings.
"""

import pytest
from radia_mcp.cubit.rules import (
    ALL_RULES,
    check_missing_block_registration,
    check_missing_mesh_command,
    check_geometry_block_for_2nd_order,
    check_missing_cubit_init,
    check_wrong_connectivity_for_2nd_order,
    check_element_type_before_add,
    check_hardcoded_absolute_paths,
    check_missing_boundary_block,
    check_ambiguous_face_block,
    check_manual_curve,
    check_missing_block_names,
    check_no_pyqt5_imports,
)


# ============================================================
# Helper
# ============================================================

def _run(rule_func, code: str, filepath: str = "test_script.py"):
    lines = code.splitlines()
    return rule_func(filepath, lines)


# ============================================================
# Cubit rules
# ============================================================

class TestMissingBlockRegistration:
    def test_detects_export_without_block(self):
        code = (
            'import cubit\n'
            'cubit.cmd("export vol output.vol order 2")\n'
        )
        findings = _run(check_missing_block_registration, code)
        assert len(findings) == 1
        assert findings[0]['severity'] == 'CRITICAL'

    def test_allows_with_block(self):
        code = (
            'import cubit\n'
            'cubit.cmd("block 1 add tet all")\n'
            'cubit.cmd("export vol output.vol order 2")\n'
        )
        assert _run(check_missing_block_registration, code) == []

    def test_no_export_no_finding(self):
        code = (
            'import cubit\n'
            'cubit.cmd("mesh volume all")\n'
        )
        assert _run(check_missing_block_registration, code) == []


class TestMissingMeshCommand:
    def test_detects_export_without_mesh(self):
        code = (
            'cubit.cmd("block 1 add tet all")\n'
            'cubit.cmd("export vol output.vol order 2")\n'
        )
        findings = _run(check_missing_mesh_command, code)
        assert len(findings) == 1
        assert findings[0]['severity'] == 'CRITICAL'

    def test_allows_with_mesh(self):
        code = (
            'cubit.cmd("mesh volume all")\n'
            'cubit.cmd("block 1 add tet all")\n'
            'cubit.cmd("export vol output.vol order 2")\n'
        )
        assert _run(check_missing_mesh_command, code) == []


class TestGeometryBlockFor2ndOrder:
    def test_detects_geometry_block_type_change(self):
        code = (
            'cubit.cmd("block 1 add volume all")\n'
            'cubit.cmd("block 1 element type tetra10")\n'
        )
        findings = _run(check_geometry_block_for_2nd_order, code)
        assert len(findings) == 1
        assert findings[0]['severity'] == 'HIGH'

    def test_allows_element_block(self):
        code = (
            'cubit.cmd("block 1 add tet all")\n'
            'cubit.cmd("block 1 element type tetra10")\n'
        )
        assert _run(check_geometry_block_for_2nd_order, code) == []


class TestMissingCubitInit:
    def test_detects_missing_init(self):
        code = (
            'import cubit\n'
            'cubit.cmd("create brick")\n'
        )
        findings = _run(check_missing_cubit_init, code)
        assert len(findings) == 1
        assert findings[0]['severity'] == 'HIGH'

    def test_allows_with_init(self):
        code = (
            'import cubit\n'
            "cubit.init(['cubit', '-nojournal'])\n"
            'cubit.cmd("create brick")\n'
        )
        assert _run(check_missing_cubit_init, code) == []

    def test_skips_no_cubit(self):
        code = 'import numpy\n'
        assert _run(check_missing_cubit_init, code) == []


class TestWrongConnectivity2ndOrder:
    def test_detects_wrong_connectivity(self):
        code = (
            'cubit.cmd("block 1 element type tetra10")\n'
            'conn = cubit.get_connectivity("tet", eid)\n'
        )
        findings = _run(check_wrong_connectivity_for_2nd_order, code)
        assert len(findings) == 1

    def test_allows_expanded(self):
        code = (
            'cubit.cmd("block 1 element type tetra10")\n'
            'conn = cubit.get_expanded_connectivity("tet", eid)\n'
        )
        assert _run(check_wrong_connectivity_for_2nd_order, code) == []


class TestElementTypeBeforeAdd:
    def test_detects_type_before_add(self):
        code = (
            'cubit.cmd("block 1 element type tetra10")\n'
            'cubit.cmd("block 1 add tet all")\n'
        )
        findings = _run(check_element_type_before_add, code)
        assert len(findings) == 1

    def test_allows_add_before_type(self):
        code = (
            'cubit.cmd("block 1 add tet all")\n'
            'cubit.cmd("block 1 element type tetra10")\n'
        )
        assert _run(check_element_type_before_add, code) == []


class TestCubitHardcodedAbsolutePaths:
    def test_detects_windows_path(self):
        code = "sys.path.insert(0, r'D:\\myproject\\lib')"
        findings = _run(check_hardcoded_absolute_paths, code)
        assert len(findings) == 1

    def test_allows_cubit_path(self):
        code = "sys.path.insert(0, r'C:\\Program Files\\Coreform Cubit 2025.12\\bin')"
        assert _run(check_hardcoded_absolute_paths, code) == []

    def test_allows_ngsolve_path(self):
        code = "sys.path.insert(0, r'C:\\NGSolve\\lib')"
        assert _run(check_hardcoded_absolute_paths, code) == []


class TestMissingBoundaryBlock:
    def test_detects_missing_surface_block(self):
        code = (
            'cubit.cmd("block 1 add tet all")\n'
            'cubit.cmd("export vol output.vol order 2")\n'
        )
        findings = _run(check_missing_boundary_block, code)
        assert len(findings) == 1

    def test_allows_with_surface_block(self):
        code = (
            'cubit.cmd("block 1 add tet all")\n'
            'cubit.cmd("block 2 add tri all")\n'
            'cubit.cmd("export vol output.vol order 2")\n'
        )
        assert _run(check_missing_boundary_block, code) == []


class TestAmbiguousFaceBlock:
    def test_detects_generic_face_block(self):
        code = (
            'cubit.cmd("block 1 add hex all")\n'
            'cubit.cmd("block 2 add face all in surface all")\n'
        )
        findings = _run(check_ambiguous_face_block, code)
        assert len(findings) == 1
        assert findings[0]['rule'] == 'ambiguous-face-block'

    def test_allows_explicit_quad_block(self):
        code = 'cubit.cmd("block 2 add quad all in surface all")\n'
        assert _run(check_ambiguous_face_block, code) == []

    def test_allows_explicit_tri_block(self):
        code = 'cubit.cmd("block 2 add tri all in surface all")\n'
        assert _run(check_ambiguous_face_block, code) == []

    def test_allows_documented_mixed_face_block(self):
        code = (
            'cubit.cmd("block 3 add face all in surface all")  '
            '# mixed tri/quad boundary\n'
        )
        assert _run(check_ambiguous_face_block, code) == []


class TestManualCurve:
    def test_detects_manual_curve_without_native_export_vol(self):
        code = (
            'mesh = Mesh("output.vol")\n'
            'mesh.Curve(3)\n'
        )
        findings = _run(check_manual_curve, code)
        assert [f["rule"] for f in findings] == ["curve-without-export-vol"], findings

    def test_allows_native_export_vol(self):
        code = (
            'cubit.cmd("export vol output.vol order 3")\n'
        )
        assert _run(check_manual_curve, code) == []

    def test_allows_occ_generate_mesh(self):
        code = (
            'mesh = geo.GenerateMesh(maxh=0.1)\n'
            'mesh.Curve(3)\n'
        )
        assert _run(check_manual_curve, code) == []

    def test_flags_mesh_curve_on_live_cubit_mesh(self):
        code = (
            'cubit.cmd("export vol output.vol order 2")\n'
            'mesh = Mesh(vol_mesh.mesh_from_cubit(2))\n'
            'mesh.Curve(3)\n'
        )
        findings = _run(check_manual_curve, code)
        assert [(f["rule"], f["severity"], f["line"]) for f in findings] == [
            ("curve-on-live-cubit-mesh", "HIGH", 3)], findings

    def test_flags_curve_through_a_mesh_wrapper(self):
        code = (
            'live = vol_mesh.mesh_from_cubit(2)\n'
            'mesh = ngsolve.Mesh(live)\n'
            'mesh.Curve(3)\n'
        )
        findings = _run(check_manual_curve, code)
        assert [(f["rule"], f["line"]) for f in findings] == [
            ("curve-on-live-cubit-mesh", 3)], findings

    def test_another_mesh_in_the_same_file_is_not_live(self):
        code = (
            'live = vol_mesh.mesh_from_cubit(2)\n'
            'occ_mesh = Mesh(geo.GenerateMesh(maxh=0.1))\n'
            'occ_mesh.Curve(3)\n'
        )
        assert _run(check_manual_curve, code) == []

    def test_allows_vol_mesh_curve(self):
        code = (
            'mesh = Mesh(vol_mesh.mesh_from_cubit(2))\n'
            'vol_mesh.curve(mesh, 3)\n'
        )
        assert _run(check_manual_curve, code) == []


class TestMissingBlockNames:
    def test_detects_unnamed_blocks(self):
        code = (
            'cubit.cmd("block 1 add tet all")\n'
            'cubit.cmd("export vol output.vol order 2")\n'
        )
        findings = _run(check_missing_block_names, code)
        assert len(findings) == 1
        assert findings[0]['severity'] == 'LOW'

    def test_allows_named_blocks(self):
        code = (
            'cubit.cmd("block 1 add tet all")\n'
            'cubit.cmd("block 1 name \\"domain\\"")\n'
            'cubit.cmd("export vol output.vol order 2")\n'
        )
        assert _run(check_missing_block_names, code) == []


class TestNoPyQt5Imports:
    def test_detects_pyqt5_import(self):
        code = (
            'from PyQt5.QtWidgets import QDialog\n'
            'dlg = QDialog()\n'
        )
        findings = _run(check_no_pyqt5_imports, code)
        assert len(findings) == 1
        assert findings[0]['rule'] == 'pyqt5-import-forbidden'

    def test_allows_pyside6_import(self):
        code = (
            'from PySide6.QtWidgets import QDialog\n'
            'dlg = QDialog()\n'
        )
        assert _run(check_no_pyqt5_imports, code) == []


# ============================================================
# Meta-tests
# ============================================================

class TestAllRulesList:
    def test_all_rules_count(self):
        assert len(ALL_RULES) == 17, len(ALL_RULES)

    def test_all_rules_callable(self):
        for rule in ALL_RULES:
            assert callable(rule)

    def test_clean_code_no_findings(self):
        """Clean code should produce zero findings across all rules."""
        clean = (
            'import numpy as np\n'
            'x = np.linspace(0, 1, 100)\n'
        )
        lines = clean.splitlines()
        for rule in ALL_RULES:
            findings = rule("clean.py", lines)
            assert findings == [], f"{rule.__name__} fired on clean code"
