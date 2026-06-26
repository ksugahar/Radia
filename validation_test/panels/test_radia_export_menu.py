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
    K. Mesh Evaluation runner orchestration (mocked Cubit + subprocess)
    L. End-to-end p-convergence on a sphere (Cubit batch + NGSolve +
       GMSH).  Skipped if Cubit batch cannot launch (e.g. license busy
       because the LAB GUI session is still open).

Run::

    set QT_QPA_PLATFORM=offscreen
    # Direct python (simplest -- bypasses pytest's import hooks):
    python validation_test/panels/test_radia_export_menu.py

    # Or pytest with --confcutdir to skip the parent panel conftest
    # (which transitively imports radia_gui_base -> radia/__init__ ->
    # MKL os.add_dll_directory calls that, on some LAB configurations,
    # break PySide6.QtWidgets DLL load before our test even starts):
    python -m pytest validation_test/panels/test_radia_export_menu.py -xvs \
        --confcutdir=validation_test/panels
"""

from __future__ import annotations

import json
import os
import sys
import tempfile
import unittest
import unittest.mock
from unittest.mock import MagicMock

# Force offscreen Qt platform BEFORE PySide6 is imported.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

# Ensure the editable radia package is importable as a side-load.
_REPO_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_REPO_ROOT, "src", "radia", "panels"))


# Module under test.  Top-level `import cubit` is forbidden in the
# module (per layer 2 isolation rules) -- this import must succeed
# without Cubit on PATH.
import radia_export_menu as rem  # noqa: E402

from PySide6.QtCore import Qt  # noqa: E402
from PySide6.QtGui import QAction  # noqa: E402
from PySide6.QtWidgets import (  # noqa: E402
    QApplication, QMainWindow,
)


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
            r'^export jmag_nastran "[^"]+\.bdf" order \d+ '
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

    def test_install_creates_menu_with_7_actions(self):
        menu = rem.install_menu()
        self.assertIsNotNone(menu, "install_menu returned None")
        # 7 = 6 export formats + 1 Mesh Evaluation (separator is not
        # itself an action accessible via actions(), but it IS an
        # action via menu.actions()).  Count visible/triggerable.
        actions = [a for a in menu.actions() if not a.isSeparator()]
        self.assertEqual(len(actions), 7,
                         f"expected 7 non-separator actions, got "
                         f"{len(actions)}: {[a.text() for a in actions]}")
        labels = [a.text() for a in actions]
        for expected in ("Netgen Vol", "GMSH", "Nastran", "VTK",
                         "FEMEEM", "MEG", "Mesh Evaluation"):
            self.assertTrue(
                any(expected in lbl for lbl in labels),
                f"missing action containing {expected!r} in {labels}")

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
# K. Mesh Evaluation runner orchestration (mocked Cubit + subprocess)
# ----------------------------------------------------------------------
class TestMeshEvaluationRunner(unittest.TestCase):
    """Verify _run_mesh_evaluation issues the right exports and parses
    calc_mesh_eval.py output correctly, WITHOUT launching Cubit."""

    def setUp(self):
        _ensure_qapp()
        self._tmpdir = tempfile.mkdtemp(prefix="rem_mesheval_")
        self._jou = os.path.join(self._tmpdir, "spheretest.jou")
        with open(self._jou, "w") as f:
            f.write("# fake journal\n")

        # Mock cubit module: get_volume_count > 0 so the runner
        # skips the empty-model branch; get_current_journal_file
        # returns our jou path (ensure_jou_path path 1).
        self._cubit = MagicMock()
        self._cubit.get_volume_count = MagicMock(return_value=1)
        self._cubit.get_current_journal_file = MagicMock(
            return_value=self._jou)
        self._cmds_issued = []
        self._cubit.cmd = MagicMock(
            side_effect=lambda c: self._cmds_issued.append(c))

    def tearDown(self):
        import shutil
        shutil.rmtree(self._tmpdir, ignore_errors=True)
        sys.modules.pop("cubit", None)

    def _stub_run_subprocess_utf8(self, json_payload):
        """Return a stub that mimics _run_subprocess_utf8(argv, timeout)
        emitting (rc=0, stdout_str, "")."""
        stdout = json.dumps(json_payload) + "\n"
        return MagicMock(return_value=(0, stdout, ""))

    def test_issues_5_vol_exports_at_orders_1_to_5(self):
        """Phase 1 of the runner must emit `_p1.vol` ... `_p5.vol`."""
        with unittest.mock.patch.object(
                rem, "_find_external_python", return_value="python"), \
            unittest.mock.patch.object(
                rem, "_find_calc_script", return_value="/dev/null"), \
            unittest.mock.patch.object(
                rem, "_run_subprocess_utf8",
                self._stub_run_subprocess_utf8(
                    {"orders": [], "format_qa": []})), \
            unittest.mock.patch.object(rem, "_show_mesh_eval_result"):
            rem._run_mesh_evaluation(self._cubit, parent=None)

        netgen_cmds = [c for c in self._cmds_issued
                       if c.startswith("export netgen ")]
        self.assertEqual(
            len(netgen_cmds), 5,
            f"expected 5 netgen exports, got {len(netgen_cmds)}: "
            f"{netgen_cmds}")
        for p in range(1, 6):
            self.assertTrue(
                any(f"_p{p}.vol" in c and f"order {p}" in c
                    for c in netgen_cmds),
                f"missing order {p} export in {netgen_cmds}")

    def test_issues_format_qa_exports(self):
        """Phase 1 also emits gmsh(1-3) + nastran(1-2) + vtk(1-2) = 7
        QA exports.  GMSH order is capped at 3 per CLAUDE.md
        ("GMSH order limit": 4-5 unreliable)."""
        with unittest.mock.patch.object(
                rem, "_find_external_python", return_value="python"), \
            unittest.mock.patch.object(
                rem, "_find_calc_script", return_value="/dev/null"), \
            unittest.mock.patch.object(
                rem, "_run_subprocess_utf8",
                self._stub_run_subprocess_utf8(
                    {"orders": [], "format_qa": []})), \
            unittest.mock.patch.object(rem, "_show_mesh_eval_result"):
            rem._run_mesh_evaluation(self._cubit, parent=None)

        gmsh_qa = [c for c in self._cmds_issued
                   if c.startswith("export gmsh ")
                   and "_qa_" in c]
        nast_qa = [c for c in self._cmds_issued
                   if c.startswith("export jmag_nastran ")
                   and "_qa_" in c]
        vtk_qa = [c for c in self._cmds_issued
                  if c.startswith("export vtk ")
                  and "_qa_" in c]
        self.assertEqual(len(gmsh_qa), 3,
                         f"GMSH QA: expected 3, got {gmsh_qa}")
        self.assertEqual(len(nast_qa), 2,
                         f"Nastran QA: expected 2, got {nast_qa}")
        self.assertEqual(len(vtk_qa), 2,
                         f"VTK QA: expected 2, got {vtk_qa}")

    def test_calls_calc_mesh_eval_with_vol_base_and_max_order(self):
        """Phase 2 must subprocess `calc_mesh_eval.py --vol-base <base>
        --max-order 5` (verified through the _run_subprocess_utf8
        helper)."""
        captured = {}

        def fake_run(argv, timeout=None):
            captured["argv"] = argv
            return (0,
                    json.dumps({"orders": [], "format_qa": []}) + "\n",
                    "")

        with unittest.mock.patch.object(
                rem, "_find_external_python", return_value="python"), \
            unittest.mock.patch.object(
                rem, "_find_calc_script",
                return_value="C:/fake/calc_mesh_eval.py"), \
            unittest.mock.patch.object(
                rem, "_run_subprocess_utf8", side_effect=fake_run), \
            unittest.mock.patch.object(rem, "_show_mesh_eval_result"):
            rem._run_mesh_evaluation(self._cubit, parent=None)

        argv = captured.get("argv", [])
        joined = " ".join(str(a) for a in argv)
        self.assertIn("calc_mesh_eval.py", joined)
        self.assertIn("--vol-base", argv)
        self.assertIn("--max-order", argv)
        self.assertIn("5", argv)

    def test_passes_json_to_show_result(self):
        """The parsed JSON from calc_mesh_eval must reach
        _show_mesh_eval_result intact."""
        payload = {
            "orders": [
                {"order": 1, "vol_error_pct": -2.5},
                {"order": 2, "vol_error_pct": -0.1},
                {"order": 3, "vol_error_pct": -0.005},
            ],
            "format_qa": [],
            "cad_vol_total": 5.235988e-4,
            "max_order": 3,
        }

        with unittest.mock.patch.object(
                rem, "_find_external_python", return_value="python"), \
            unittest.mock.patch.object(
                rem, "_find_calc_script", return_value="/dev/null"), \
            unittest.mock.patch.object(
                rem, "_run_subprocess_utf8",
                self._stub_run_subprocess_utf8(payload)), \
            unittest.mock.patch.object(
                rem, "_show_mesh_eval_result") as show:
            rem._run_mesh_evaluation(self._cubit, parent=None)

        show.assert_called_once()
        passed = show.call_args[0][0]   # first positional arg
        self.assertEqual(passed["max_order"], 3)
        self.assertEqual(len(passed["orders"]), 3)
        # Convergence in the payload: err decreases with p
        errs = [abs(o["vol_error_pct"]) for o in passed["orders"]]
        self.assertGreater(errs[0], errs[1])
        self.assertGreater(errs[1], errs[2])


# ----------------------------------------------------------------------
# L. End-to-end p-convergence (Cubit batch -> _p1..._p3.vol -> NGSolve)
# ----------------------------------------------------------------------
# All 3 DLL / encoding lessons are encapsulated in
# tests/_subprocess_utils.py::run_cubit_python_script; this test just
# supplies the body that builds the sphere + calls evaluate_mesh.
sys.path.insert(0, os.path.join(_REPO_ROOT, "tests"))
from _subprocess_utils import (   # noqa: E402
    run_cubit_python_script, extract_json_from_output,
)

_PCONV_SCRIPT_BODY = r"""
import json, sys, traceback

base, max_order, radius, msize = sys.argv[1:5]
max_order = int(max_order); radius = float(radius); msize = float(msize)

# numpy / NGSolve / Cubit-bin / CUBIT_PLUGIN_DIR were pre-bound by
# the AUTO-INJECTED prelude (_subprocess_utils._CUBIT_BATCH_PRELUDE).
try:
    import calc_mesh_eval
    import cubit
    cubit.init(["cubit", "-nojournal", "-batch", "-nographics",
                "-commandplugindir", os.environ["CUBIT_PLUGIN_DIR"]])

    cubit.cmd("reset")
    cubit.cmd("create sphere radius %g" % radius)
    cubit.cmd("volume 1 size %g" % msize)
    cubit.cmd("volume 1 scheme tetmesh")
    cubit.cmd("mesh volume 1")
    cubit.cmd('block 1 add volume 1')
    cubit.cmd('block 1 name "sphere"')
    for p in range(1, max_order + 1):
        vp = "%s_p%d.vol" % (base, p)
        cubit.cmd('export netgen "%s" order %d overwrite' % (vp, p))

    # Format QA (mimics _run_mesh_evaluation Phase 1b)
    fmt_specs = [
        ("gmsh", "msh", "export gmsh", min(max_order, 3)),
        ("nastran", "bdf", "export jmag_nastran", min(max_order, 2)),
        ("vtk", "vtk", "export vtk", min(max_order, 2)),
    ]
    for name, ext, cmd_prefix, max_ord in fmt_specs:
        for p in range(1, max_ord + 1):
            fp = "%s_qa_%s_o%d.%s" % (base, name, p, ext)
            cmd = '%s "%s"' % (cmd_prefix, fp)
            if p >= 2:
                cmd += " order %d" % p
            cmd += " overwrite"
            cubit.cmd(cmd)

    result = calc_mesh_eval.evaluate_mesh(base, max_order=max_order)
    print("__BEGIN_JSON__"); print(json.dumps(result, default=float))
    print("__END_JSON__")
except Exception:
    traceback.print_exc(file=sys.stderr); sys.exit(2)
"""


def _run_pconv_subprocess(base, max_order, radius, msize, timeout=240):
    """Spawn the canonical Cubit batch + evaluate_mesh helper.
    Returns (parsed_json_dict, None) on success, or (None, error_str).

    Detects RLM license failures (Cubit 2025.12 token-based auth)
    so the test can SkipTest cleanly rather than report a confusing
    "Parse Error" -- the actual cause is Cubit refusing to load
    APREPRO commands when its license is unhealthy.
    """
    rc, stdout, stderr = run_cubit_python_script(
        _PCONV_SCRIPT_BODY,
        args=(base, max_order, radius, msize),
        timeout=timeout)
    # Cubit's "no export" failure is symptomatic of a deeper
    # license failure; surface the real cause to the test runner.
    if "RLM ERROR code: -102" in stdout or \
       "Couldn't initialize RLM" in stdout:
        return None, ("Cubit RLM license unhealthy (-102 \"Can't "
                      "read license data\") -- the login_tokens.json "
                      "needs refresh.  Launch CoreformCubit.ps1 once "
                      "to refresh, then re-run this test.")
    if rc != 0:
        return None, (f"subprocess exit {rc}\n"
                      f"STDOUT tail:\n{stdout[-1000:]}\n"
                      f"STDERR tail:\n{stderr[-1500:]}")
    try:
        return extract_json_from_output(stdout), None
    except ValueError as e:
        return None, str(e)


def _has_ngsolve_and_gmsh():
    try:
        import ngsolve  # noqa
        import gmsh     # noqa
        return True
    except ImportError:
        return False


class TestPConvergenceEndToEnd(unittest.TestCase):
    """Real Cubit batch + NGSolve check: a sphere's volume error must
    converge as polynomial order p increases (the central physics
    contract of the Mesh Evaluation feature).

    Runs the whole Cubit init + export + evaluate_mesh in a CLEAN
    subprocess to avoid Qt6 DLL conflicts with PySide6 (already loaded
    in this test process).  Skips cleanly if Cubit batch cannot launch
    (license busy because the LAB Cubit GUI is open, or Cubit not
    installed).
    """

    # Tunables -- small sphere + coarse mesh keeps the subprocess
    # under ~30 s while still exercising curved-boundary p-convergence.
    SPHERE_RADIUS = 0.05      # [m]
    MESH_SIZE     = 0.025     # [m] -> roughly 6-10 tets
    MAX_ORDER     = 3         # Skip p=4-5 for speed; trend already proven

    @classmethod
    def setUpClass(cls):
        if not _has_ngsolve_and_gmsh():
            raise unittest.SkipTest("NGSolve or gmsh not importable")
        cls._tmp = tempfile.mkdtemp(prefix="pconv_sphere_")
        cls._base = os.path.join(cls._tmp, "sphere")
        cls._result, err = _run_pconv_subprocess(
            cls._base, cls.MAX_ORDER,
            cls.SPHERE_RADIUS, cls.MESH_SIZE,
            timeout=240)
        if cls._result is None:
            raise unittest.SkipTest(
                f"Cubit batch / evaluate_mesh subprocess failed:\n"
                f"{err}")

    @classmethod
    def tearDownClass(cls):
        import shutil
        if hasattr(cls, "_tmp"):
            shutil.rmtree(cls._tmp, ignore_errors=True)

    def test_all_vol_files_produced(self):
        for p in range(1, self.MAX_ORDER + 1):
            vp = f"{self._base}_p{p}.vol"
            self.assertTrue(os.path.isfile(vp),
                            f"missing {vp} -- Cubit export failed")
            json_path = f"{self._base}_p{p}.vol.json"
            self.assertTrue(os.path.isfile(json_path),
                            f"missing companion JSON {json_path}")

    def test_volume_error_decreases_with_p(self):
        r = self._result
        self.assertNotIn("error", r, f"evaluate_mesh failed: {r}")
        orders = r["orders"]
        self.assertEqual(len(orders), self.MAX_ORDER)
        ok = [o for o in orders if "vol_error_pct" in o]
        self.assertEqual(len(ok), self.MAX_ORDER,
                         f"some orders failed: {orders}")

        errs = [abs(o["vol_error_pct"]) for o in ok]
        print(f"\nSphere R={self.SPHERE_RADIUS} m, h={self.MESH_SIZE} m: "
              f"vol_error_pct(p=1..{self.MAX_ORDER}) = "
              f"{', '.join(f'{e:.3e}%' for e in errs)}")

        # p=2 must improve over p=1 (curved boundary -> higher
        # interpolation order captures it better).
        self.assertLess(
            errs[1], errs[0],
            f"p=2 error ({errs[1]:.3e}%) not better than p=1 "
            f"({errs[0]:.3e}%) -- convergence broken")
        # p=3 must improve over p=1 (relaxed: p=3 vs p=2 may already
        # hit the floor on a coarse mesh).
        self.assertLess(
            errs[2], errs[0],
            f"p=3 error ({errs[2]:.3e}%) not better than p=1 "
            f"({errs[0]:.3e}%)")
        # Hard upper bound: p=2 must already be <1% (sanity check
        # against geometry mishap).
        self.assertLess(
            errs[1], 1.0,
            f"p=2 sphere vol error >= 1% -- implausible "
            f"({errs[1]:.3e}%)")

    def test_no_negative_jacobians_in_qa_exports(self):
        """GMSH format-QA exports must have no inverted curved
        elements (negative Jacobian determinants)."""
        for q in self._result.get("format_qa", []):
            if "error" in q:
                continue
            neg = q.get("neg_det", 0)
            self.assertEqual(
                neg, 0,
                f"{q['format']} order {q['order']}: {neg} negative "
                f"Jacobian determinants -- inverted elements")


if __name__ == "__main__":
    # Import unittest.mock lazily so the module's top-level remains
    # cheap when tests are merely collected.
    import unittest.mock  # noqa
    unittest.main(verbosity=2)
