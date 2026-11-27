"""
Netgen Export menu — one-level menu-bar action + PySide6 dialogs.

Clicking **Export Netgen** on Cubit's menu bar opens the Netgen dialog directly.
The dialogs are PySide6 and run inside Cubit's embedded Python
interpreter.  The Netgen action optionally adds the Kelvin open-boundary
sphere (``kelvin.py`` of ``radia_mcp.cubit``, loaded by path), runs the
native plugin's ``export netgen`` command and verifies the result with
``python -m radia.vol_check`` in the external interpreter named by
``RADIA_MCP_CUBIT_PYTHON``.  Dialogs are parented to the Cubit main window.

Layer 2 (Cubit GUI Python).  See ``radia_mcp.cubit.knowledge.custom_toolbar``
(TOOLBAR_PYTHON_SCRIPT_CONVENTIONS, TOOLBAR_CLARO_PARENT_HELPER).  This file
is copied into the toolbar payload, so it imports nothing from ``radia_mcp``.

This module avoids:
  - parenthesized multi-line PySide6 imports (Cubit's importer chokes
    on them) — use backslash continuations
  - top-level `import cubit` — the host already imports it
  - dialogs without `parent=find_claro()` — they get hidden behind the
    Cubit main window on Windows

Qt binding: **PySide6 only**. Target is Cubit 2026.8
which ships PySide6.  No PyQt5 fallback — per CLAUDE.md "No Fallbacks
— Fail Fast, Fail Loud", an old Cubit without PySide6 must raise the
underlying ImportError loudly so the operator can fix the env.
"""

import json
import os
import sys

# Qt binding: PySide6 only. Target is Cubit 2026.8
# which ships PySide6.  Per CLAUDE.md "No Fallbacks — Fail Fast, Fail
# Loud" any ImportError propagates so the operator sees the real cause
# (old Cubit, broken PySide6 install, wrong Python env) instead of a
# silent PyQt5 fallback that masks it.
#
# Use backslash continuations (NOT parens) per Cubit toolbar convention.
# The menu-bar entry is a module-owned QAction; Cubit-owned menu objects
# are never enumerated from PySide6 (that deadlocks Cubit).
from PySide6.QtWidgets import QApplication, QCheckBox, QComboBox, \
    QDialog, QDialogButtonBox, QFileDialog, QFormLayout, QGroupBox, \
    QHBoxLayout, QLabel, QLineEdit, QMainWindow, \
    QMessageBox, QPushButton, QTableWidget, QTableWidgetItem, \
    QVBoxLayout
from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QAction


# ----------------------------------------------------------------------
# Settings persistence
# ----------------------------------------------------------------------

def _settings_dir():
    """Return the exporter-owned per-user settings directory."""
    roaming = os.environ.get("APPDATA") or os.path.join(
        os.path.expanduser("~"), "AppData", "Roaming")
    appdata = os.path.join(roaming, "radia-mcp")
    os.makedirs(appdata, exist_ok=True)
    return appdata


def _settings_path():
    return os.path.join(_settings_dir(), "export_settings.json")


def _load_settings():
    p = _settings_path()
    if not os.path.isfile(p):
        return {}
    try:
        with open(p, "r", encoding="utf-8") as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError):
        return {}


def _save_settings(obj):
    try:
        with open(_settings_path(), "w", encoding="utf-8") as f:
            json.dump(obj, f, indent=2)
    except OSError:
        pass


def _default_start_dir():
    """Read default_dir from export_settings.json (set by Python startup)."""
    s = _load_settings()
    dd = s.get("default_dir")
    if dd and os.path.isdir(dd):
        return dd.replace("\\", "/")
    return ""


# ----------------------------------------------------------------------
# Claro (Cubit main window) helper — see custom_toolbar.py
# ----------------------------------------------------------------------

def find_claro():
    """Return Cubit's main QMainWindow (objectName 'Claro'), or None.

    Filter by objectName, never windowTitle (locale-dependent and
    changes when a file is loaded).  Used as `parent=` for every
    QDialog so the dialog cannot fall behind the main window.

    MEASURED on Coreform Cubit 2025.12: the main window's objectName is
    ``"Claro"`` with a CAPITAL C, and a second QMainWindow named
    ``"QtJournalEditor"`` also exists.  Matching the old lowercase
    ``"claro"`` therefore never hit, and the any-QMainWindow fallback
    silently returned the journal editor instead -- which is why dialogs
    were parented to the wrong window.  Compare case-insensitively and
    only fall back when no window is named claro at all.
    """
    app = QApplication.instance()
    if app is None or not hasattr(app, "topLevelWidgets"):
        return None
    for w in app.topLevelWidgets():
        if isinstance(w, QMainWindow) and w.objectName().lower() == "claro":
            return w
    # Fallback: a QMainWindow that is not one of Cubit's known auxiliary
    # windows (helps in headless / test contexts).
    for w in app.topLevelWidgets():
        if isinstance(w, QMainWindow) and w.objectName() != "QtJournalEditor":
            return w
    return None


def _is_valid_jou_path(p):
    if not p:
        return False
    p_norm = p.replace("\\", "/")
    if p_norm in ("/", ".", ".."):
        return False
    base = os.path.basename(p_norm)
    if not base:
        return False
    # Strip extension if any — basename without ext must be non-empty
    stem, _ = os.path.splitext(base)
    return bool(stem)


def _current_journal_hint(cubit_mod):
    """Return Cubit's current journal as an optional output-name hint."""
    try:
        path = cubit_mod.get_current_journal_file()
    except Exception:
        return ""
    if not path:
        return ""
    normalized = path.replace("\\", "/")
    return normalized if _is_valid_jou_path(normalized) else ""


def _ensure_model(cubit_mod, parent):
    """If no geometry is loaded, prompt for a .jou to play.

    Returns the .jou path that was played (empty string if model was
    already loaded; empty if cancelled).
    """
    if cubit_mod.get_volume_count() > 0:
        return ""

    jou, _ = QFileDialog.getOpenFileName(
        parent, "No model loaded - Select Journal File",
        _default_start_dir(), "Cubit Journal (*.jou);;All Files (*)")
    if not jou:
        return ""
    jou = jou.replace("\\", "/")
    cubit_mod.cmd(f'play "{jou}"')
    if cubit_mod.get_volume_count() > 0:
        return jou
    return ""


# ----------------------------------------------------------------------
# ExportDialog — format-aware QDialog (port of ExportDialog in
# the removed Qt5 component).
# ----------------------------------------------------------------------

# Key of the Netgen options in the exporter-owned settings file.
_SETTINGS_KEY = "netgen_vol"


class ExportDialog(QDialog):
    """Netgen .vol export options dialog.

    Emits the `export netgen` command and the Kelvin options.  Settings are
    persisted to ``%APPDATA%/radia-mcp/export_settings.json`` under
    ``netgen_vol``.
    """

    def __init__(self, source_path_hint, cubit_mod, parent=None):
        super().__init__(parent)
        self._cubit = cubit_mod
        self.setWindowTitle("Export Netgen Vol")
        self.setMinimumWidth(500)

        layout = QVBoxLayout(self)
        form = QFormLayout()

        # A loaded journal may suggest a useful initial directory and
        # basename, but it is never required and never controls the export.
        default_dir = ""
        base_name = "ExportedMesh"

        if source_path_hint:
            jp = source_path_hint.replace("\\", "/")
            slash = jp.rfind("/")
            if slash >= 0:
                default_dir = jp[:slash]
                fname = jp[slash + 1:]
            else:
                fname = jp
            dot = fname.rfind(".")
            if dot > 0:
                base_name = fname[:dot]

        if not default_dir:
            default_dir = _default_start_dir()
        if not default_dir:
            default_dir = os.getcwd().replace("\\", "/")

        # --- Directory + filename row -------------------------------
        row = QHBoxLayout()
        self._dir = QLineEdit(default_dir)
        browse = QPushButton("...")
        browse.setFixedWidth(30)
        browse.clicked.connect(self._browse_dir)
        self._dir.textChanged.connect(self._update_preview)
        row.addWidget(self._dir)
        row.addWidget(browse)
        form.addRow("Directory:", row)
        self._filename = QLineEdit(base_name + ".vol")
        self._filename.textChanged.connect(self._update_preview)
        form.addRow("Filename:", self._filename)

        # --- Order combo ---------------------------------------------
        self._order = QComboBox()
        self._order.addItems(["1", "2", "3", "4", "5"])
        self._order.setCurrentIndex(2)  # default 3
        self._order.currentIndexChanged.connect(self._update_preview)
        form.addRow("Order:", self._order)

        layout.addLayout(form)

        # --- Kelvin / symmetry --------------------------------------
        box = QGroupBox(
            "Kelvin open-boundary  (auto-add exterior sphere; "
            "needs an 'air' block)")
        kv = QVBoxLayout(box)

        row1 = QHBoxLayout()
        self._kelvin_enable = QCheckBox("Add Kelvin (idempotent)")
        self._kelvin_enable.setToolTip(
            "Idempotent: skipped if a 'kelvin' block already exists.\n"
            "Uncheck for Dirichlet-truncation baselines.")
        self._kelvin_enable.toggled.connect(self._update_preview)
        row1.addWidget(self._kelvin_enable)
        row1.addSpacing(20)
        row1.addWidget(QLabel("Mesh size [m]:"))
        self._kelvin_mesh_size = QLineEdit()
        self._kelvin_mesh_size.setPlaceholderText(
            "auto (blank) or e.g. 0.03")
        self._kelvin_mesh_size.setMaximumWidth(140)
        self._kelvin_mesh_size.setToolTip(
            "Tet edge length on the Kelvin shell.  Blank inherits "
            "the size from the air outer surface (copy-mesh).")
        self._kelvin_mesh_size.textChanged.connect(self._update_preview)
        row1.addWidget(self._kelvin_mesh_size)
        row1.addStretch()
        kv.addLayout(row1)

        row2 = QHBoxLayout()
        row2.addWidget(QLabel("Symmetry reduction:"))

        def _make_sym(axis):
            cb = QComboBox()
            cb.addItem("off")
            cb.addItem("bn")
            cb.addItem("ht")
            cb.setMaximumWidth(70)
            cb.setToolTip(
                f"BC label on the {axis}=0 plane when the air block "
                f"is reduced to {axis}>=0.\n"
                "  bn = B.n=0   (flux parallel,   solver '+', "
                "A:Dirichlet, Omega:natural)\n"
                "  ht = HxN=0   (flux perp,       solver '-', "
                "A:natural,    Omega:Dirichlet)\n"
                "  off = no reduction on this axis")
            cb.currentIndexChanged.connect(self._update_preview)
            return cb

        self._kelvin_sym_x = _make_sym("x")
        self._kelvin_sym_y = _make_sym("y")
        self._kelvin_sym_z = _make_sym("z")
        row2.addWidget(QLabel("x:"))
        row2.addWidget(self._kelvin_sym_x)
        row2.addSpacing(8)
        row2.addWidget(QLabel("y:"))
        row2.addWidget(self._kelvin_sym_y)
        row2.addSpacing(8)
        row2.addWidget(QLabel("z:"))
        row2.addWidget(self._kelvin_sym_z)
        row2.addStretch()
        kv.addLayout(row2)

        layout.addWidget(box)

        # --- Command preview row ------------------------------------
        prev_form = QFormLayout()
        self._preview = QLineEdit()
        self._preview.setReadOnly(True)
        self._preview.setAlignment(Qt.AlignLeft)
        prev_form.addRow("Command:", self._preview)
        layout.addLayout(prev_form)

        # --- OK / Cancel buttons ------------------------------------
        buttons = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self._on_accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        # --- Load saved settings ------------------------------------
        # A source-path hint wins for the initial directory. Otherwise reuse
        # the user's previous output directory.
        all_settings = _load_settings()
        s = all_settings.get(_SETTINGS_KEY, {})
        if not source_path_hint and "dir" in s:
            self._dir.setText(s["dir"])
        if "order" in s:
            self._order.setCurrentIndex(int(s["order"]))
        if "kelvin_enable" in s:
            self._kelvin_enable.setChecked(bool(s["kelvin_enable"]))
        if "kelvin_mesh_size" in s:
            self._kelvin_mesh_size.setText(s["kelvin_mesh_size"])
        for cb, key in (
                (self._kelvin_sym_x, "kelvin_sym_x"),
                (self._kelvin_sym_y, "kelvin_sym_y"),
                (self._kelvin_sym_z, "kelvin_sym_z")):
            if key in s:
                idx = cb.findText(s[key])
                if idx >= 0:
                    cb.setCurrentIndex(idx)

        self._update_preview()

    # ------------------------------------------------------------------
    def _browse_dir(self):
        d = QFileDialog.getExistingDirectory(
            self, "Select Output Directory", self._dir.text())
        if d:
            self._dir.setText(d.replace("\\", "/"))

    # ------------------------------------------------------------------
    def file_path(self):
        d = self._dir.text().replace("\\", "/")
        if not d.endswith("/"):
            d += "/"
        return d + self._filename.text()

    # ------------------------------------------------------------------
    def order(self):
        return self._order.currentIndex() + 1

    # ------------------------------------------------------------------
    def cubit_command(self):
        """Build the `export netgen` command string."""
        return export_command(self.file_path(), self.order())

    # ------------------------------------------------------------------
    def kelvin_options(self):
        """The Kelvin options in the form ``kelvin.validate_options`` takes."""
        return kelvin_options(
            self._kelvin_enable.isChecked(),
            self._kelvin_mesh_size.text(),
            self._kelvin_sym_x.currentText(),
            self._kelvin_sym_y.currentText(),
            self._kelvin_sym_z.currentText())

    # ------------------------------------------------------------------
    def _update_preview(self):
        text = self.cubit_command()
        if self._kelvin_enable.isChecked():
            text = "[Auto-Kelvin first] " + text
        self._preview.setText(text)
        self._preview.setCursorPosition(0)

    # ------------------------------------------------------------------
    def _on_accept(self):
        # Persist settings before accepting
        all_settings = _load_settings()
        all_settings[_SETTINGS_KEY] = {
            "dir": self._dir.text(),
            "order": self._order.currentIndex(),
            "kelvin_enable": self._kelvin_enable.isChecked(),
            "kelvin_mesh_size": self._kelvin_mesh_size.text(),
            "kelvin_sym_x": self._kelvin_sym_x.currentText(),
            "kelvin_sym_y": self._kelvin_sym_y.currentText(),
            "kelvin_sym_z": self._kelvin_sym_z.currentText(),
        }
        _save_settings(all_settings)
        self.accept()


def export_command(path, order, kelvin_arguments=""):
    """The native plugin command for one curved ``.vol`` export.

    ``kelvin_arguments`` is ``kelvin.export_arguments(info)`` after Auto-Kelvin;
    the plugin then writes the Kelvin pairs and needs a fresh output path.
    """
    extra = f"{kelvin_arguments} " if kelvin_arguments else ""
    return f'export netgen "{path}" order {int(order)} {extra}overwrite'


def kelvin_options(enabled, mesh_size, sym_x, sym_y, sym_z):
    """Map the dialog fields onto the Kelvin option names of ``kelvin.py``.

    A blank mesh size inherits the air outer-surface size; any other text is
    passed through so ``kelvin.validate_options`` rejects it loudly.
    """
    mesh = str(mesh_size).strip()
    return {"add_kelvin": bool(enabled),
            "kelvin_mesh": mesh or None,
            "kelvin_sym_x": sym_x, "kelvin_sym_y": sym_y, "kelvin_sym_z": sym_z}


def _kelvin_module():
    """Load ``kelvin.py`` by path: beside this file (toolbar payload) or one
    level up (installed ``radia_mcp/cubit`` package).  It is stdlib-only and
    imports ``cubit`` inside its functions."""
    import importlib.util
    here = os.path.dirname(os.path.abspath(__file__))
    for candidate in (os.path.join(here, "kelvin.py"),
                      os.path.join(os.path.dirname(here), "kelvin.py")):
        if os.path.isfile(candidate):
            spec = importlib.util.spec_from_file_location(
                "_radia_mcp_cubit_kelvin", candidate)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            return module
    raise RuntimeError(f"kelvin.py not found beside {here} or its parent; "
                       "run python -m radia_mcp.cubit.install to repair the toolbar")


def _open_in_os(path):
    """Open a file or directory with the OS-associated application."""
    try:
        if sys.platform == "win32":
            os.startfile(path)
        elif sys.platform == "darwin":
            import subprocess
            subprocess.Popen(["open", path])
        else:
            import subprocess
            subprocess.Popen(["xdg-open", path])
    except OSError:
        pass


# ----------------------------------------------------------------------
# Netgen Vol runner — exports + runs the radia.vol_check subprocess and
# shows a consistency-check table dialog.
# ----------------------------------------------------------------------

def _run_subprocess_utf8(argv, timeout):
    """Subprocess.run wrapper hardened for Cubit / NGSolve invocation.

    Three protections rolled into one call:
    1. Strips QT_* and PYSIDE* env vars so the child process does not
       inherit our offscreen Qt platform.
    2. Sets ``PYTHONIOENCODING=utf-8`` so the child's print() emits
       UTF-8 -- Windows console default is cp932 in Japanese locale,
       which crashes our reader thread when Cubit / NGSolve print
       non-ASCII text on stderr (memory-dump headers, file paths
       with Japanese chars, etc.).
    3. Captures bytes (text=False) and decodes ourselves with
       errors="replace" so any unexpected non-UTF-8 byte sequence
       on stderr does not abort the read.

    Returns: (returncode, stdout_str, stderr_str).  On TimeoutExpired,
    re-raises -- the caller decides the policy.
    """
    import subprocess
    env = os.environ.copy()
    for k in list(env.keys()):
        if "QT" in k.upper() or "PYSIDE" in k.upper():
            del env[k]
    env["PYTHONIOENCODING"] = "utf-8"
    proc = subprocess.run(argv, capture_output=True, timeout=timeout,
                            env=env)
    out = (proc.stdout or b"").decode("utf-8", errors="replace")
    err = (proc.stderr or b"").decode("utf-8", errors="replace")
    return proc.returncode, out, err


PYTHON_ENV = "RADIA_MCP_CUBIT_PYTHON"


def _external_python():
    """The external interpreter that has radia and NGSolve (never Cubit's own).

    The startup shim sets ``RADIA_MCP_CUBIT_PYTHON`` to the installer's
    interpreter.  It must be an existing absolute path; PATH is never searched.
    """
    python = os.environ.get(PYTHON_ENV, "")
    if not python or not os.path.isabs(python) or not os.path.isfile(python):
        raise RuntimeError(f"{PYTHON_ENV} must name the absolute path of the Python that "
                           f"has radia installed; got {python!r}")
    return python


def verification_command(python, vol_path, report_path):
    """Isolated ``radia.vol_check`` run that writes the JSON report."""
    return [python, "-I", "-B", "-m", "radia.vol_check", vol_path,
            "--format", "json", "--report-json", report_path]


def _run_netgen_export(cubit_mod, parent):
    """Optional Kelvin, native export, then radia.vol_check verification."""
    source_path_hint = ""
    if cubit_mod.get_volume_count() == 0:
        source_path_hint = _ensure_model(cubit_mod, parent)
        if cubit_mod.get_volume_count() == 0:
            return
    if not source_path_hint:
        source_path_hint = _current_journal_hint(cubit_mod)

    dlg = ExportDialog(source_path_hint, cubit_mod, parent=parent)
    if dlg.exec() != QDialog.Accepted:
        return

    order = dlg.order()
    vol_path = dlg.file_path().replace("\\", "/")
    out_dir = os.path.dirname(vol_path).replace("\\", "/")
    try:
        python = _external_python()
        options = dlg.kelvin_options()
        kelvin = _kelvin_module() if options["add_kelvin"] else None
        if kelvin is not None:
            kelvin.validate_options(options)
    except (RuntimeError, ValueError) as exc:
        print(f"Export not started: {exc}")
        QMessageBox.warning(parent, "Export not started", str(exc))
        return
    if kelvin is not None:
        # The Kelvin pairs publish only to a fresh path; never delete user files.
        existing = [p for p in (vol_path, vol_path + ".json", vol_path + ".partial.vol")
                    if os.path.lexists(p)]
        if existing:
            message = ("A Kelvin export needs a fresh output path. These files already "
                       "exist; choose another name or remove them yourself:\n" + "\n".join(existing))
            print(f"Export not started: {message}")
            QMessageBox.warning(parent, "Export not started", message)
            return
    cubit_mod.cmd(f'cd "{out_dir}"')

    # --- Phase 0: optional Kelvin open-boundary sphere (fails loudly) ---
    command = dlg.cubit_command()
    if kelvin is not None:
        print("Adding the Kelvin open-boundary sphere...")
        try:
            info = kelvin.apply_options(options)
            print(f"Kelvin: {info}")
            if not info:
                raise RuntimeError("the model already has a Kelvin block, so its offset is "
                                   "unknown here; rebuild it from the journal")
            command = export_command(vol_path, order, kelvin.export_arguments(info))
        except Exception as exc:
            print(f"Kelvin FAILED: {type(exc).__name__}: {exc}")
            QMessageBox.warning(parent, "Kelvin failed",
                                f"{type(exc).__name__}: {exc}\nNothing was exported.")
            return

    # --- Phase 1: C++ export .vol + companion JSON ---
    print(f"Exporting Netgen .vol (order {order})...")
    cubit_mod.cmd(command)

    if not os.path.isfile(vol_path):
        QMessageBox.warning(
            parent, "Export failed",
            f"export netgen command did not produce:\n{vol_path}")
        return
    print(f"Exported: {vol_path}")

    # --- Phase 2: radia.vol_check in the external interpreter ---
    print("Verifying mesh with radia.vol_check...")

    import subprocess
    report_path = vol_path + ".check.json"
    if os.path.isfile(report_path):
        os.remove(report_path)

    try:
        rc, stdout, stderr = _run_subprocess_utf8(
            verification_command(python, vol_path, report_path), timeout=300)
    except subprocess.TimeoutExpired:
        print("Verification FAILED: radia.vol_check timed out.")
        QMessageBox.warning(parent, "Verification",
                            "radia.vol_check timed out.")
        return
    except OSError as exc:
        print(f"Verification FAILED: cannot start radia.vol_check: {exc}")
        QMessageBox.warning(parent, "Verification failed",
                            f"Cannot start radia.vol_check:\n{exc}")
        return

    try:
        with open(report_path, encoding="utf-8") as stream:
            r = json.load(stream)
    except (OSError, json.JSONDecodeError):
        print(f"Verification FAILED: no valid radia.vol_check report: "
              f"{stderr[:500]}")
        QMessageBox.warning(parent, "Verification",
                            f"No valid radia.vol_check report:\n{stderr[:2000]}")
        return
    print(_verification_summary(rc, r, report_path))
    if rc != 0 or not r.get("passed", False):
        QMessageBox.warning(parent, "Verification failed",
                            json.dumps(r, indent=2)[:4000])
        return

    _show_netgen_result(r, vol_path, order, parent)
    _open_in_os(vol_path)


def _verification_summary(rc, report, report_path):
    """Console lines for a radia.vol_check report, so the log keeps the verdict."""
    passed = rc == 0 and bool(report.get("passed", False))
    lines = [
        f"Verification {'passed' if passed else 'FAILED'}: "
        f"{report.get('n_elements', '?')} elements, "
        f"{report.get('n_points', '?')} points, "
        f"order {report.get('order', '?')} (report: {report_path})"]
    if report.get("error"):
        lines.append(f"  error: {report['error']}")
    for warning in report.get("warnings", [])[:10]:
        lines.append(f"  warning: {warning}")
    return "\n".join(lines)


def _show_netgen_result(r, vol_path, order, parent):
    """Modeless dialog showing per-label volume/area/length consistency."""
    dlg = QDialog(parent)
    dlg.setAttribute(Qt.WA_DeleteOnClose)
    dlg.setWindowTitle(f"Netgen Vol Export - Order {order}")
    dlg.setMinimumSize(700, 400)
    layout = QVBoxLayout(dlg)

    n_el = r.get("n_elements", 0)
    layout.addWidget(QLabel(
        f"Output: {vol_path}\nElements: {n_el}, Order: {order}"))

    def _make_table(label, headers, data, fill_row):
        layout.addWidget(QLabel(label))
        if not data:
            layout.addWidget(QLabel("  (none)"))
            return
        tbl = QTableWidget(len(data), len(headers), dlg)
        tbl.setHorizontalHeaderLabels(headers)
        tbl.horizontalHeader().setStretchLastSection(True)
        tbl.verticalHeader().setVisible(False)
        tbl.setEditTriggers(QTableWidget.NoEditTriggers)
        for i, item in enumerate(data):
            fill_row(tbl, i, item)
        tbl.resizeColumnsToContents()
        layout.addWidget(tbl)

    def _err_item(err):
        return QTableWidgetItem(f"{err:.2e}")

    # Material table
    mats = r.get("materials", [])

    def _fill_mat(t, i, m):
        t.setItem(i, 0, QTableWidgetItem(m.get("name", "")))
        t.setItem(i, 1, QTableWidgetItem(
            f"{m.get('cad_volume', 0.0):.6e}"))
        t.setItem(i, 2, QTableWidgetItem(
            f"{m['ng_volume']:.6e}" if "ng_volume" in m else "N/A"))
        if "error_pct" in m:
            t.setItem(i, 3, _err_item(m["error_pct"]))
    _make_table(
        "Materials (Volume):",
        ["Name", "CAD Volume", "NGSolve Volume", "Error [%]"],
        mats, _fill_mat)

    # Boundary table
    bnds = r.get("boundaries", [])

    def _fill_bnd(t, i, b):
        t.setItem(i, 0, QTableWidgetItem(b.get("name", "")))
        t.setItem(i, 1, QTableWidgetItem(
            f"{b.get('cad_area', 0.0):.6e}"))
        t.setItem(i, 2, QTableWidgetItem(
            f"{b['ng_area']:.6e}" if "ng_area" in b else "N/A"))
        if "error_pct" in b:
            t.setItem(i, 3, _err_item(b["error_pct"]))
    _make_table(
        "Boundaries (Surface Area):",
        ["Name", "CAD Area", "NGSolve Area", "Error [%]"],
        bnds, _fill_bnd)

    # Edge table
    edges = r.get("edges", [])

    def _fill_edge(t, i, e):
        t.setItem(i, 0, QTableWidgetItem(e.get("name", "")))
        t.setItem(i, 1, QTableWidgetItem(
            f"{e.get('cad_length', 0.0):.6e}"))
        t.setItem(i, 2, QTableWidgetItem(
            f"{e['ng_length']:.6e}" if "ng_length" in e else "N/A"))
        if "error_pct" in e:
            t.setItem(i, 3, _err_item(e["error_pct"]))
    _make_table(
        "Edges (BBND Length):",
        ["Name", "CAD Length", "NGSolve Length", "Error [%]"],
        edges, _fill_edge)

    # Warnings
    warns = r.get("warnings", [])
    if warns:
        html = "<b style='color:red'>Warnings:</b><ul>"
        for w in warns:
            html += f"<li>{w}</li>"
        html += "</ul>"
        wl = QLabel(html)
        wl.setTextFormat(Qt.RichText)
        layout.addWidget(wl)

    buttons = QDialogButtonBox(QDialogButtonBox.Ok, dlg)
    buttons.accepted.connect(dlg.close)
    layout.addWidget(buttons)
    dlg.show()


def launch_netgen_export():
    """Open the Netgen .vol export dialog from the menu or the toolbar.

    ``WorkflowToolbar`` buttons execute small Python scripts.  Those scripts
    call this entry point instead of modifying Cubit's QMenuBar.
    """
    parent = find_claro()
    if parent is None:
        raise RuntimeError("Cubit main window is not available")
    import cubit as cubit_mod
    return _run_netgen_export(cubit_mod, parent)


# ----------------------------------------------------------------------
# Menu installation
# ----------------------------------------------------------------------

# The single top-level menu-bar action.  The toolbar template, its script
# and the toolbar checks carry the same label.
_MENU_LABEL = "Export Netgen (.vol)..."
_MENU_BAR_TEXT = "&Export Netgen"
_MENU_STATUS_TIP = "Export curved mesh with labels (.vol, order 1-5)"
_MENU_OBJECT_NAME = "radia_mcp_export_netgen_action"

# The QAction is created and owned by this module, so removing it never
# touches Cubit-owned menu objects.
_menu_action = None
_shutdown_connected = False
# Bumped by every install/remove.  A deferred re-attach only runs while its
# generation is still current, so remove_menu() cancels a pending one.
_attach_generation = 0


def _on_export_triggered(_checked=False):
    import traceback
    try:
        launch_netgen_export()
    except Exception:
        traceback.print_exc()


def _attach_menu_action():
    """Add the Export action directly to Cubit's menu bar (one level)."""
    global _menu_action
    parent = find_claro()
    if parent is None:
        return None
    _detach_menu_action()
    bar = parent.menuBar()
    action = QAction(_MENU_BAR_TEXT, bar)
    action.setObjectName(_MENU_OBJECT_NAME)
    action.setToolTip(_MENU_LABEL)
    action.setStatusTip(_MENU_STATUS_TIP)
    action.triggered.connect(_on_export_triggered)
    bar.addAction(action)
    _menu_action = action
    return action


def _detach_menu_action():
    global _menu_action
    action = _menu_action
    _menu_action = None
    if action is None:
        return False
    try:
        bar = action.parent()
        if bar is not None:
            bar.removeAction(action)
        action.deleteLater()
    except RuntimeError:
        # Cubit replaced its menu bar during startup; the old bar and its
        # child action are already destroyed on the C++ side.
        pass
    return True


def _deferred_attach(generation):
    if generation == _attach_generation:
        _attach_menu_action()


def install_menu():
    """Install a one-level "Export Netgen" entry on Cubit's menu bar.

    Clicking **Export Netgen** opens the Netgen ``.vol`` dialog directly; there is
    no submenu.  Cubit's Claro API (``emclaro.add_to_menu``) can only add
    actions under a popup menu, so the entry is a module-owned ``QAction``
    added to the main window's ``QMenuBar``.  It is attached once now and
    again after the event loop starts, because Cubit rebuilds its menu bar
    during cold start.  Only this module's own action is ever touched; do
    not enumerate Cubit-owned menus from PySide6 (that deadlocks Cubit).

    Returns True on success, or None when the Cubit main window is
    unavailable (batch mode, or Cubit not present).
    """
    global _shutdown_connected, _attach_generation
    _attach_generation += 1
    generation = _attach_generation
    if _attach_menu_action() is None:
        print("[radia-mcp] Cubit main window not available - "
              "Export menu not installed")
        return None
    QTimer.singleShot(0, lambda: _deferred_attach(generation))
    app = QApplication.instance()
    if app is not None and not _shutdown_connected:
        app.aboutToQuit.connect(remove_menu)
        _shutdown_connected = True
    print("[radia-mcp] Export menu installed (one-level action)")
    return True


def remove_menu():
    """Remove the Export action and cancel any pending re-attach."""
    global _attach_generation
    _attach_generation += 1
    return _detach_menu_action()


# Compatibility entry point for an explicit interactive ``play``.  Normal
# installation uses the official WorkflowToolbar package instead.
if __name__ == "_coreform_cubit":
    install_menu()
