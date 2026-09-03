"""Unit test for src/radia/panels/radia_export_menu.py (PySide6).

Exercises the new PySide6 Radia Export toolbar that replaced the legacy
C++ Qt5 .ccl plugin in 2026-05.  Tests are headless (QT_QPA_PLATFORM=
offscreen) and do NOT require Cubit or NGSolve -- only PySide6.

Coverage:
    A. Module imports + symbol surface
    B. Export actions are registered through Cubit's Claro API
    C. ExportDialog.cubit_command() round-trip for all 6 formats
       (default options yield the expected `export ...` APREPRO
        commands so the .ccm side is unchanged)
    D. Netgen Kelvin options (add_kelvin, kelvin_mesh, kelvin_sym_*)
    E. Nastran nopyramid flag
    F. FEMEEM scale flag
    G. MEG labels (block table mocking)
    H. install_menu() idempotency on a stub QMainWindow
    I. Current journal is an optional output-name hint, never a prerequisite
    J. Settings persistence round-trip
    K. Mesh Evaluation is not part of the Cubit menu; the
       p-convergence demonstration lives under docs/cubit_mesh_export.

Run::

    set QT_QPA_PLATFORM=offscreen
    # Direct python (simplest -- bypasses pytest's import hooks):
    python validation_test/panels/test_radia_export_menu.py

    # Or pytest with --confcutdir to isolate this Cubit-embedded PySide test
    # from unrelated validation fixtures:
    python -m pytest validation_test/panels/test_radia_export_menu.py -xvs \
        --confcutdir=validation_test/panels
"""

from __future__ import annotations

import os
import sys
import tempfile
import types
import unittest
from unittest.mock import MagicMock, patch

# Force offscreen Qt platform BEFORE PySide6 is imported.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

# Ensure the editable radia package is importable as a side-load.
_REPO_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_REPO_ROOT, "src", "radia", "panels"))

try:
    from PySide6.QtWidgets import (  # noqa: E402
        QApplication, QMainWindow,
    )
except ImportError as exc:
    message = (
        "PySide6 is required only in the Cubit panel runtime; "
        "skip radia_export_menu validation on normal Radia Python."
    )
    if __name__ == "__main__":
        print(f"SKIP: {message}")
        raise SystemExit(0) from exc
    raise unittest.SkipTest(message) from exc


# Module under test.  Top-level `import cubit` is forbidden in the
# module (per layer 2 isolation rules) -- this import must succeed
# without Cubit on PATH when PySide6 is available.
import radia_export_menu as rem  # noqa: E402


_QAPP = None


def _ensure_qapp():
    """Create / reuse a single QApplication for all tests."""
    global _QAPP
    if _QAPP is None:
        _QAPP = QApplication.instance() or QApplication(sys.argv)
    return _QAPP


# ----------------------------------------------------------------------
# A. Module imports + symbol surface
# ----------------------------------------------------------------------
class TestModuleSurface(unittest.TestCase):

    def test_format_constants(self):
        self.assertEqual(rem.FMT_NETGEN, "netgen_vol")
        self.assertEqual(rem.FMT_GMSH, "gmsh")
        self.assertEqual(rem.FMT_NASTRAN, "nastran")
        self.assertEqual(rem.FMT_VTK, "vtk")
        self.assertEqual(rem.FMT_FEMEEM, "femeem")
        self.assertEqual(rem.FMT_MEG, "meg")

    def test_format_extensions(self):
        self.assertEqual(rem._FORMAT_EXTS[rem.FMT_NETGEN], ".vol")
        self.assertEqual(rem._FORMAT_EXTS[rem.FMT_GMSH], ".msh")
        self.assertEqual(rem._FORMAT_EXTS[rem.FMT_NASTRAN], ".bdf")
        self.assertEqual(rem._FORMAT_EXTS[rem.FMT_VTK], ".vtk")
        self.assertEqual(rem._FORMAT_EXTS[rem.FMT_FEMEEM], "")  # dir output
        self.assertEqual(rem._FORMAT_EXTS[rem.FMT_MEG], ".meg")

    def test_public_callables_present(self):
        for name in ("install_menu", "find_claro", "_current_journal_hint",
                     "ExportDialog"):
            self.assertTrue(hasattr(rem, name),
                            f"radia_export_menu missing: {name}")
            self.assertTrue(callable(getattr(rem, name)))

    def test_no_top_level_cubit_import(self):
        """Module must not `import cubit` at top level -- breaks pytest /
        any non-Cubit context per Layer-2 architecture (CLAUDE.md)."""
        with open(rem.__file__, "r", encoding="utf-8") as f:
            src = f.read()
        # Allow `import cubit` only INSIDE a function (indented),
        # never at column 0.
        for line in src.splitlines():
            stripped = line.lstrip()
            if stripped.startswith("import cubit") or \
                    stripped.startswith("from cubit "):
                indent = len(line) - len(stripped)
                self.assertGreater(
                    indent, 0,
                    f"top-level cubit import found: {line!r}")


# ----------------------------------------------------------------------
# B. Cubit owns the menu actions through emclaro
# ----------------------------------------------------------------------
class TestClaroOwnership(unittest.TestCase):

    def test_module_does_not_construct_qt_actions(self):
        with open(rem.__file__, "r", encoding="utf-8") as f:
            source = f.read()
        self.assertNotIn("from PySide6.QtGui import QAction", source)
        self.assertNotIn("QMenu(", source)
        self.assertIn("emclaro.add_to_menu", source)


# ----------------------------------------------------------------------
# C-G. ExportDialog.cubit_command() round-trip
# ----------------------------------------------------------------------
class _DialogTestBase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        _ensure_qapp()
        cls._tmpdir = tempfile.mkdtemp(prefix="rem_test_")
        cls._jou = os.path.join(cls._tmpdir, "myjob.jou")
        # Touch a fake .jou so file_path() can derive basename
        with open(cls._jou, "w") as f:
            f.write("# fake journal\n")

    def _make_dialog(self, fmt, cubit_stub=None):
        # All formats except MEG never touch the cubit_mod; MEG calls
        # get_block_id_list -> we stub it as empty so the block table
        # stays empty.
        stub = cubit_stub or MagicMock()
        stub.get_block_id_list = MagicMock(return_value=[])
        d = rem.ExportDialog(fmt, self._jou, stub)
        return d


class TestNetgenCommand(_DialogTestBase):

    def test_default_command(self):
        d = self._make_dialog(rem.FMT_NETGEN)
        cmd = d.cubit_command()
        # Default order is 3 (combo index 2)
        self.assertRegex(
            cmd,
            r'^export netgen "[^"]+myjob\.vol" order 3 overwrite$')
        # No Kelvin by default
        self.assertNotIn("add_kelvin", cmd)

    def test_order_variants(self):
        d = self._make_dialog(rem.FMT_NETGEN)
        for combo_idx, expected_order in [(0, 1), (1, 2), (2, 3),
                                          (3, 4), (4, 5)]:
            d._order.setCurrentIndex(combo_idx)
            cmd = d.cubit_command()
            self.assertIn(f"order {expected_order} ", cmd + " ")

    def test_kelvin_enabled_default(self):
        d = self._make_dialog(rem.FMT_NETGEN)
        d._kelvin_enable.setChecked(True)
        cmd = d.cubit_command()
        self.assertIn("add_kelvin", cmd)
        # Default sym_x/y/z combo is "(none)" -> no flags emitted
        self.assertNotIn("kelvin_sym_x", cmd)
        self.assertNotIn("kelvin_sym_y", cmd)
        self.assertNotIn("kelvin_sym_z", cmd)

    def test_kelvin_with_mesh_size_and_symmetry(self):
        d = self._make_dialog(rem.FMT_NETGEN)
        d._kelvin_enable.setChecked(True)
        d._kelvin_mesh_size.setText("0.05")
        # Find "bn" / "ht" in the combo items (order may vary)
        for combo, val in [(d._kelvin_sym_x, "bn"),
                            (d._kelvin_sym_y, "ht"),
                            (d._kelvin_sym_z, "bn")]:
            for i in range(combo.count()):
                if combo.itemText(i) == val:
                    combo.setCurrentIndex(i)
                    break
        cmd = d.cubit_command()
        self.assertIn("add_kelvin", cmd)
        self.assertIn("kelvin_mesh 0.05", cmd)
        self.assertIn("kelvin_sym_x bn", cmd)
        self.assertIn("kelvin_sym_y ht", cmd)
        self.assertIn("kelvin_sym_z bn", cmd)


class TestGmshCommand(_DialogTestBase):

    def test_default_3d(self):
        d = self._make_dialog(rem.FMT_GMSH)
        cmd = d.cubit_command()
        self.assertRegex(
            cmd,
            r'^export gmsh "[^"]+\.msh" order \d+ dimension 3 '
            r'overwrite$')

    def test_2d_dimension(self):
        d = self._make_dialog(rem.FMT_GMSH)
        for i in range(d._dimension.count()):
            if d._dimension.itemText(i) == "2D":
                d._dimension.setCurrentIndex(i)
                break
        self.assertIn("dimension 2", d.cubit_command())


class TestNastranCommand(_DialogTestBase):

    def test_default_command(self):
        d = self._make_dialog(rem.FMT_NASTRAN)
        cmd = d.cubit_command()
        self.assertRegex(
            cmd,
            r'^export nastran_bdf "[^"]+\.bdf" order \d+ '
            r'dimension 3 overwrite$')
        # nopyramid default = off
        self.assertNotIn("nopyramid", cmd)

    def test_nopyramid_flag(self):
        d = self._make_dialog(rem.FMT_NASTRAN)
        d._nopyramid.setCurrentIndex(1)  # "Yes"
        self.assertIn("nopyramid", d.cubit_command())


class TestVtkCommand(_DialogTestBase):

    def test_default_command(self):
        d = self._make_dialog(rem.FMT_VTK)
        self.assertRegex(
            d.cubit_command(),
            r'^export vtk "[^"]+\.vtk" order \d+ '
            r'dimension 3 overwrite$')


class TestFemeemCommand(_DialogTestBase):

    def test_default_command(self):
        d = self._make_dialog(rem.FMT_FEMEEM)
        cmd = d.cubit_command()
        # FEMEEM uses directory output (no .ext on filename)
        self.assertRegex(
            cmd,
            r'^export femeem "[^"]+" scale [\d.]+ overwrite$')

    def test_custom_scale(self):
        d = self._make_dialog(rem.FMT_FEMEEM)
        d._scale.setText("0.001")
        self.assertIn("scale 0.001", d.cubit_command())


class TestMegCommand(_DialogTestBase):

    def test_default_3d(self):
        d = self._make_dialog(rem.FMT_MEG)
        cmd = d.cubit_command()
        self.assertRegex(
            cmd,
            r'^export meg "[^"]+\.meg" threed overwrite$')

    def test_2d_dimension(self):
        d = self._make_dialog(rem.FMT_MEG)
        # MEG combo items are ["3D (threed)", "2D (twod)", "Axisymmetric"]
        # -- index 1 is the 2D entry (label != "2D" exactly).
        d._dimension.setCurrentIndex(1)
        self.assertIn(" twod ", d.cubit_command())

    def test_axisym_dimension(self):
        d = self._make_dialog(rem.FMT_MEG)
        # 3rd combo entry is "Axisymmetric"
        d._dimension.setCurrentIndex(2)
        self.assertIn(" axisymmetric ", d.cubit_command())

    def test_block_labels_passthrough(self):
        # Stub cubit with a single block "iron" -> default MMB
        stub = MagicMock()
        stub.get_block_id_list = MagicMock(return_value=[1])
        stub.get_block_name = MagicMock(return_value="iron")
        d = rem.ExportDialog(rem.FMT_MEG, self._jou, stub)
        cmd = d.cubit_command()
        self.assertIn('labels "1:MMB"', cmd)

    def test_air_block_skipped(self):
        # AIR / KELVIN should default to "(none) - Skip export"
        stub = MagicMock()
        stub.get_block_id_list = MagicMock(return_value=[7, 8])
        stub.get_block_name = MagicMock(side_effect=["air", "iron"])
        d = rem.ExportDialog(rem.FMT_MEG, self._jou, stub)
        cmd = d.cubit_command()
        # air -> skipped, only iron labeled
        self.assertIn('labels "8:MMB"', cmd)
        self.assertNotIn("7:", cmd)


# ----------------------------------------------------------------------
# H. install_menu() idempotency on stub QMainWindow
# ----------------------------------------------------------------------
class _FakePyAction:
    def __init__(self):
        self.text = ""
        self.menu_text = ""
        self.status_tip = ""
        self.activate_method = ""

    def setText(self, value):
        self.text = value

    def setMenuText(self, value):
        self.menu_text = value

    def setStatusTip(self, value):
        self.status_tip = value

    def setActivateMethod(self, value):
        self.activate_method = value


class _FakePyActionVector(list):
    def push_back(self, value):
        self.append(value)


class TestInstallMenu(unittest.TestCase):

    def setUp(self):
        rem._claro_keepalive.clear()
        self.removed = []
        self.added = []
        self.fake_emclaro = types.SimpleNamespace(
            is_loaded=lambda: True,
            remove_menu_items=self.removed.append,
            PyAction=_FakePyAction,
            PyActionVector=_FakePyActionVector,
            add_to_menu=lambda *args: self.added.append(args),
        )
        sys.modules["emclaro"] = self.fake_emclaro

    def tearDown(self):
        sys.modules.pop("emclaro", None)
        rem._claro_keepalive.clear()

    def test_install_creates_menu_with_6_export_actions(self):
        self.assertTrue(rem.install_menu())
        self.assertEqual(self.removed, [rem._CLARO_COMPONENT])
        self.assertEqual(len(self.added), 1)
        menu_title, actions, component = self.added[0]
        self.assertEqual(menu_title, rem._CLARO_MENU_TITLE)
        self.assertEqual(component, rem._CLARO_COMPONENT)
        self.assertEqual(len(actions), 6)
        self.assertEqual(
            [action.text for action in actions],
            [spec[0] for spec in rem._MENU_SPECS],
        )
        self.assertTrue(all("launch_export" in action.activate_method
                            for action in actions))

    def test_install_is_idempotent(self):
        rem.install_menu()
        rem.install_menu()
        self.assertEqual(self.removed,
                         [rem._CLARO_COMPONENT, rem._CLARO_COMPONENT])
        self.assertEqual(len(self.added), 2)
        self.assertEqual(len(rem._claro_keepalive), 7)


class TestFindClaro(unittest.TestCase):

    def setUp(self):
        _ensure_qapp()
        self._win = QMainWindow()
        self._win.setObjectName("claro")
        self._win.show()

    def tearDown(self):
        self._win.hide()
        self._win.setParent(None)
        self._win.deleteLater()
        self._win = None
        _ensure_qapp().processEvents()

    def test_find_claro_returns_our_stub(self):
        m = rem.find_claro()
        self.assertIs(m, self._win)


# ----------------------------------------------------------------------
# I. Optional current-journal hint
# ----------------------------------------------------------------------
class TestCurrentJournalHint(unittest.TestCase):

    def test_already_loaded_via_cubit(self):
        stub = MagicMock()
        with tempfile.NamedTemporaryFile(suffix=".jou", delete=False) as f:
            jou = f.name
        try:
            stub.get_current_journal_file = MagicMock(return_value=jou)
            got = rem._current_journal_hint(stub)
            self.assertEqual(got.replace("\\", "/"), jou.replace("\\", "/"))
            stub.cmd.assert_not_called()
        finally:
            os.unlink(jou)

    def test_no_journal_is_valid_for_loaded_model(self):
        stub = MagicMock()
        stub.get_current_journal_file = MagicMock(return_value="")
        self.assertEqual(rem._current_journal_hint(stub), "")
        stub.cmd.assert_not_called()

    def test_dialog_without_journal_uses_regular_default(self):
        _ensure_qapp()
        stub = MagicMock()
        stub.get_block_id_list = MagicMock(return_value=[])
        dialog = rem.ExportDialog(rem.FMT_NETGEN, "", stub)
        self.assertTrue(dialog.file_path().endswith("/ExportedMesh.vol"))

    def test_loaded_model_export_never_prompts_to_save_journal(self):
        stub = MagicMock()
        stub.get_volume_count = MagicMock(return_value=1)
        stub.get_current_journal_file = MagicMock(return_value="")
        dialog = MagicMock()
        dialog.exec.return_value = rem.QDialog.Rejected
        with patch.object(rem, "ExportDialog", return_value=dialog) as factory, \
                patch.object(rem.QFileDialog, "getSaveFileName") as save_dialog:
            rem._run_generic_export(rem.FMT_GMSH, stub, parent=None)
        save_dialog.assert_not_called()
        factory.assert_called_once_with(rem.FMT_GMSH, "", stub, parent=None)


# ----------------------------------------------------------------------
# J. Settings persistence round-trip
# ----------------------------------------------------------------------
class TestSettingsRoundTrip(unittest.TestCase):

    def setUp(self):
        # Redirect settings to a temp dir
        self._tmp = tempfile.mkdtemp(prefix="rem_settings_")
        self._orig_settings_dir = rem._settings_dir
        rem._settings_dir = lambda: self._tmp

    def tearDown(self):
        rem._settings_dir = self._orig_settings_dir
        import shutil
        shutil.rmtree(self._tmp, ignore_errors=True)

    def test_load_empty_is_dict(self):
        self.assertEqual(rem._load_settings(), {})

    def test_save_then_load(self):
        rem._save_settings({"netgen_vol": {"order": 4,
                                            "dir": "C:/temp"}})
        got = rem._load_settings()
        self.assertEqual(got["netgen_vol"]["order"], 4)
        self.assertEqual(got["netgen_vol"]["dir"], "C:/temp")

    def test_load_corrupt_returns_empty(self):
        with open(os.path.join(self._tmp, "export_settings.json"),
                  "w") as f:
            f.write("{ not valid json ]")
        self.assertEqual(rem._load_settings(), {})


# ----------------------------------------------------------------------
# K. Mesh evaluation routing
# ----------------------------------------------------------------------
class TestMeshEvaluationDocsRouting(unittest.TestCase):

    def test_mesh_evaluation_runner_removed_from_menu_module(self):
        self.assertFalse(hasattr(rem, "_run_mesh_evaluation"))
        self.assertFalse(hasattr(rem, "_show_mesh_eval_result"))

    def test_docs_demo_is_the_p_convergence_surface(self):
        demo = os.path.join(
            _REPO_ROOT, "docs", "cubit_mesh_export", "netgen",
            "p_convergence_demo.ipynb")
        self.assertTrue(os.path.isfile(demo), demo)
        with open(demo, "r", encoding="utf-8") as f:
            text = f.read()
        self.assertIn("calc_mesh_eval.py", text)
        self.assertIn("p-convergence", text)
        self.assertIn("validation_test", text)
        self.assertIn("p_convergence_demo_results.json", text)


if __name__ == "__main__":
    unittest.main(verbosity=2)
