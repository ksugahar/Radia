"""Tests for Cubit lint rules.

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
    check_setgeominfo_without_geometry,
    check_nodeset_sideset_usage,
    check_missing_name_occ_faces,
    check_missing_step_reimport,
    check_noheal_for_named_workflow,
    check_hardcoded_absolute_paths,
    check_missing_boundary_block,
    check_export_file_extension,
    check_curve_without_setgeominfo,
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
            'cme.export_Gmsh_ver4(cubit, "output.msh")\n'
        )
        findings = _run(check_missing_block_registration, code)
        assert len(findings) == 1
        assert findings[0]['severity'] == 'CRITICAL'

    def test_allows_with_block(self):
        code = (
            'import cubit\n'
            'cubit.cmd("block 1 add tet all")\n'
            'cme.export_Gmsh_ver4(cubit, "output.msh")\n'
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
            'cme.export_Gmsh_ver4(cubit, "output.msh")\n'
        )
        findings = _run(check_missing_mesh_command, code)
        assert len(findings) == 1
        assert findings[0]['severity'] == 'CRITICAL'

    def test_allows_with_mesh(self):
        code = (
            'cubit.cmd("mesh volume all")\n'
            'cubit.cmd("block 1 add tet all")\n'
            'cme.export_Gmsh_ver4(cubit, "output.msh")\n'
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


class TestSetGeomInfoWithoutGeometry:
    def test_stub_returns_empty(self):
        code = 'set_cylinder_geominfo(mesh, surface_id)\n'
        # Rule is now a stub (old API detection handled by check_missing_step_reimport)
        assert _run(check_setgeominfo_without_geometry, code) == []


class TestNodesetSidesetUsage:
    def test_detects_nodeset_with_non_exodus(self):
        code = (
            'cubit.cmd("nodeset 1 add node all")\n'
            'cme.export_Gmsh_ver4(cubit, "output.msh")\n'
        )
        findings = _run(check_nodeset_sideset_usage, code)
        assert len(findings) == 1
        assert findings[0]['severity'] == 'HIGH'

    def test_no_finding_without_export(self):
        code = 'cubit.cmd("nodeset 1 add node all")\n'
        assert _run(check_nodeset_sideset_usage, code) == []


class TestMissingNameOccFaces:
    def test_stub_returns_empty(self):
        code = 'name_occ_faces(shape)\n'
        # Rule is now a stub (old API removed)
        assert _run(check_missing_name_occ_faces, code) == []


class TestMissingStepReimport:
    def test_detects_deleted_api(self):
        code = (
            'cme.export_netgen(cubit, geometry=True)\n'
        )
        findings = _run(check_missing_step_reimport, code)
        assert len(findings) >= 1

    def test_allows_radia_export_netgen(self):
        code = (
            'mesh = cubit.cmd("export netgen output.vol order 2")\n'
        )
        assert _run(check_missing_step_reimport, code) == []


class TestNohealForNamedWorkflow:
    def test_stub_returns_empty(self):
        code = 'name_occ_faces(shape)\n'
        # Rule is now a stub (old name-based workflow removed)
        assert _run(check_noheal_for_named_workflow, code) == []


class TestCubitHardcodedAbsolutePaths:
    def test_detects_windows_path(self):
        code = "sys.path.insert(0, r'D:\\myproject\\lib')"
        findings = _run(check_hardcoded_absolute_paths, code)
        assert len(findings) == 1

    def test_allows_cubit_path(self):
        code = "sys.path.insert(0, r'C:\\Program Files\\Coreform Cubit 2024.3\\bin')"
        assert _run(check_hardcoded_absolute_paths, code) == []

    def test_allows_ngsolve_path(self):
        code = "sys.path.insert(0, r'C:\\NGSolve\\lib')"
        assert _run(check_hardcoded_absolute_paths, code) == []


class TestMissingBoundaryBlock:
    def test_detects_missing_surface_block(self):
        code = (
            'cubit.cmd("block 1 add tet all")\n'
            'cubit.cmd("export netgen output.vol order 2")\n'
        )
        findings = _run(check_missing_boundary_block, code)
        assert len(findings) == 1

    def test_allows_with_surface_block(self):
        code = (
            'cubit.cmd("block 1 add tet all")\n'
            'cubit.cmd("block 2 add tri all")\n'
            'cubit.cmd("export netgen output.vol order 2")\n'
        )
        assert _run(check_missing_boundary_block, code) == []


class TestExportFileExtension:
    def test_detects_wrong_extension(self):
        code = 'cme.export_Gmsh_ver4(cubit, "output.txt")\n'
        findings = _run(check_export_file_extension, code)
        assert len(findings) == 1
        assert '.txt' in findings[0]['message']

    def test_allows_correct_extension(self):
        code = 'cme.export_Gmsh_ver4(cubit, "output.msh")\n'
        assert _run(check_export_file_extension, code) == []


class TestCurveWithoutSetGeomInfo:
    def test_detects_manual_curve_without_radia_export_netgen(self):
        code = (
            'mesh = Mesh("output.vol")\n'
            'mesh.Curve(3)\n'
        )
        findings = _run(check_curve_without_setgeominfo, code)
        assert len(findings) == 1

    def test_allows_radia_export_netgen(self):
        code = (
            'cubit.cmd("export netgen output.vol order 3")\n'
        )
        assert _run(check_curve_without_setgeominfo, code) == []

    def test_allows_occ_generate_mesh(self):
        code = (
            'mesh = geo.GenerateMesh(maxh=0.1)\n'
            'mesh.Curve(3)\n'
        )
        assert _run(check_curve_without_setgeominfo, code) == []


class TestMissingBlockNames:
    def test_detects_unnamed_blocks(self):
        code = (
            'cubit.cmd("block 1 add tet all")\n'
            'cme.export_Gmsh_ver4(cubit, "output.msh")\n'
        )
        findings = _run(check_missing_block_names, code)
        assert len(findings) == 1
        assert findings[0]['severity'] == 'LOW'

    def test_allows_named_blocks(self):
        code = (
            'cubit.cmd("block 1 add tet all")\n'
            'cubit.cmd("block 1 name \\"domain\\"")\n'
            'cme.export_Gmsh_ver4(cubit, "output.msh")\n'
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
        # 19 original + 3 Cubit in-process toolbar rules + PySide6-only gate.
        assert len(ALL_RULES) == 23

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
