"""Unit test for src/radia/panels/radia_export_menu.py (PySide6).

Exercises the new PySide6 Radia Export toolbar that replaced the legacy
C++ Qt5 .ccl plugin in 2026-05.  Tests are headless (QT_QPA_PLATFORM=
offscreen) and do NOT require Cubit or NGSolve -- only PySide6.

Coverage:
    A. Module imports + symbol surface
    B. QAction comes from QtGui (PySide6/Qt6 location), not QtWidgets
    C. ExportDialog.cubit_command() round-trip for all 6 formats
       (default options yield the expected `export ...` APREPRO
        commands so the .ccm side is unchanged)
    D. Netgen Kelvin options (add_kelvin, kelvin_mesh, kelvin_sym_*)
    E. Nastran nopyramid flag
    F. FEMEEM scale flag
    G. MEG labels (block table mocking)
    H. install_menu() idempotency on a stub QMainWindow
    I. ensure_jou_path() session-stickiness (_last_jou_path)
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
import unittest
from unittest.mock import MagicMock

# Force offscreen Qt platform BEFORE PySide6 is imported.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

# Ensure the editable radia package is importable as a side-load.
_REPO_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_REPO_ROOT, "src", "radia", "panels"))

try:
    from PySide6.QtCore import Qt  # noqa: E402
    from PySide6.QtGui import QAction  # noqa: E402
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
        for name in ("install_menu", "find_claro", "ensure_jou_path",
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
# B. QAction location (PySide6/Qt6 moved it from QtWidgets to QtGui)
# ----------------------------------------------------------------------
class TestQActionLocation(unittest.TestCase):

    def test_qaction_from_qtgui_not_qtwidgets(self):
        """PySide6 ships QAction in QtGui; importing from QtWidgets fails.
        Agent's bug fix during the port -- this test locks it down."""
        from PySide6 import QtGui, QtWidgets
        self.assertTrue(hasattr(QtGui, "QAction"))
        self.assertFalse(hasattr(QtWidgets, "QAction"))

    def test_module_uses_qtgui_qaction(self):
        # The QAction symbol the module imported must resolve from QtGui.
        self.assertEqual(rem.QAction.__module__, "PySide6.QtGui")


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
class TestInstallMenu(unittest.TestCase):

    def setUp(self):
        _ensure_qapp()
        # Inject a stub cubit module into sys.modules so install_menu's
        # `import cubit` succeeds in this headless context.
        self._fake_cubit = MagicMock()
        sys.modules["cubit"] = self._fake_cubit

        # Create a fake Cubit main window with objectName='claro'.
        self._win = QMainWindow()
        self._win.setObjectName("claro")
        self._win.show()  # offscreen -> still counts as topLevelWidget

    def tearDown(self):
        # Remove the menu and the stub window completely, then process
        # Qt's deleteLater queue so the next test's find_claro() can
        # NOT pick up this window as a stale topLevelWidget.
        bar = self._win.menuBar()
        for act in list(bar.actions()):
            sub = act.menu()
            if sub and sub.objectName() == rem._INSTALLED_OBJECT_NAME:
                bar.removeAction(act)
                sub.deleteLater()
        self._win.hide()
        self._win.setParent(None)
        self._win.deleteLater()
        self._win = None
        _ensure_qapp().processEvents()
        sys.modules.pop("cubit", None)

    def _radia_export_menu(self):
        bar = self._win.menuBar()
        for act in bar.actions():
            sub = act.menu()
            if sub and sub.objectName() == rem._INSTALLED_OBJECT_NAME:
                return sub
        return None

    def test_install_creates_menu_with_6_export_actions(self):
        menu = rem.install_menu()
        self.assertIsNotNone(menu, "install_menu returned None")
        # Radia Export is intentionally export-only.  Mesh evaluation is
        # documented as a docs notebook demo, not a Cubit menu action.
        actions = [a for a in menu.actions() if not a.isSeparator()]
        self.assertEqual(len(actions), 6,
                         f"expected 6 non-separator actions, got "
                         f"{len(actions)}: {[a.text() for a in actions]}")
        labels = [a.text() for a in actions]
        for expected in ("Netgen Vol", "GMSH", "Nastran", "VTK",
                         "FEMEEM", "MEG"):
            self.assertTrue(
                any(expected in lbl for lbl in labels),
                f"missing action containing {expected!r} in {labels}")
        self.assertFalse(
            any("Mesh Evaluation" in lbl for lbl in labels),
            f"Mesh Evaluation must not be exposed in the menu: {labels}")

    def test_install_is_idempotent(self):
        rem.install_menu()
        rem.install_menu()  # twice -> still 1 menu
        # Look in the window install_menu actually picked (might be
        # ours, or another 'claro' if Qt has stale leftovers from a
        # previous test that processEvents could not yet collect).
        target = rem.find_claro()
        self.assertIsNotNone(target, "find_claro() returned None")
        bar = target.menuBar()
        n = sum(1 for a in bar.actions()
                if a.menu() and a.menu().objectName() ==
                rem._INSTALLED_OBJECT_NAME)
        self.assertEqual(n, 1, f"after 2 install_menu calls, "
                              f"{n} Radia Export menus exist (should be 1)")

    def test_find_claro_returns_our_stub(self):
        # find_claro() should find our fake QMainWindow.
        m = rem.find_claro()
        self.assertIs(m, self._win)


# ----------------------------------------------------------------------
# I. ensure_jou_path session-stickiness
# ----------------------------------------------------------------------
class TestEnsureJouPathSession(unittest.TestCase):

    def setUp(self):
        rem._last_jou_path = ""

    def test_already_loaded_via_cubit(self):
        # Branch 1: cubit.get_current_journal_file() returns a valid path
        stub = MagicMock()
        with tempfile.NamedTemporaryFile(suffix=".jou", delete=False) as f:
            jou = f.name
        try:
            stub.get_current_journal_file = MagicMock(return_value=jou)
            got = rem.ensure_jou_path(stub, parent=None)
            self.assertEqual(got.replace("\\", "/"), jou.replace("\\", "/"))
        finally:
            os.unlink(jou)

    def test_session_sticky_after_first_set(self):
        # Branch 2: no cubit-loaded jou, but _last_jou_path was set
        # earlier in the session -> reuse without prompting.
        stub = MagicMock()
        stub.get_current_journal_file = MagicMock(return_value="")
        with tempfile.NamedTemporaryFile(suffix=".jou", delete=False) as f:
            jou = f.name
        try:
            rem._last_jou_path = jou
            got = rem.ensure_jou_path(stub, parent=None)
            self.assertEqual(got, jou)
            # The cubit module's `save` should NOT have been called --
            # the sticky path is reused, no resave.
            stub.cmd.assert_not_called()
        finally:
            os.unlink(jou)


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
        self.assertIn("docs/cubit_mesh_export", text)


if __name__ == "__main__":
    unittest.main(verbosity=2)
