"""
Shared base for Radia standalone analysis windows (PySide6).

Provides:
  - ModePanel: form layout with add_line/add_combo/add_spin/add_browse helpers
  - AnalysisWindow: QMainWindow with Run/Stop/Save/GMSH buttons + QProcess
"""

import json
import os
import sys

from PySide6.QtCore import Qt, QProcess
from PySide6.QtGui import QFont, QIcon
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget,
    QVBoxLayout, QHBoxLayout, QFormLayout,
    QLabel, QLineEdit, QComboBox, QSpinBox,
    QPushButton, QPlainTextEdit, QFileDialog,
    QMessageBox, QSplitter, QGroupBox, QStyle,
)

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PANELS_DIR = os.path.join(_THIS_DIR, "panels")
_RESOURCES_DIR = os.path.join(_THIS_DIR, "resources")
_PYTHON = sys.executable
_SETTINGS_DIR = os.path.join(os.path.expanduser("~"), ".radia")


def _icon_path():
    for ext in (".ico", ".png"):
        p = os.path.join(_RESOURCES_DIR, "radia_icon" + ext)
        if os.path.isfile(p):
            return p
    return ""


def calc_script(name):
    return os.path.join(_PANELS_DIR, name)


def msh_output(vol_path, suffix):
    base = vol_path if vol_path else "output"
    return os.path.splitext(base)[0] + suffix + ".msh"


# ============================================================
# ModePanel: form layout base
# ============================================================

class ModePanel(QWidget):
    """Base class for analysis parameter panels."""

    _RED = "QLineEdit { background-color: #FFD0D0; }"
    _NORMAL = ""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._form = QFormLayout()
        self._form.setLabelAlignment(Qt.AlignRight)
        self.setLayout(self._form)
        self._widgets = {}
        self._row_indices = {}

    def _set_row_visible(self, key, visible):
        row_idx = self._row_indices.get(key)
        if row_idx is None:
            return
        label_item = self._form.itemAt(row_idx, QFormLayout.LabelRole)
        field_item = self._form.itemAt(row_idx, QFormLayout.FieldRole)
        if label_item and label_item.widget():
            label_item.widget().setVisible(visible)
        if field_item and field_item.widget():
            field_item.widget().setVisible(visible)

    def add_line(self, key, label, default="", placeholder=""):
        w = QLineEdit(default)
        if placeholder:
            w.setPlaceholderText(placeholder)
        self._form.addRow(label, w)
        self._widgets[key] = w
        self._row_indices[key] = self._form.rowCount() - 1
        return w

    def add_combo(self, key, label, items, default=0):
        w = QComboBox()
        w.addItems(items)
        w.setCurrentIndex(default)
        self._form.addRow(label, w)
        self._widgets[key] = w
        self._row_indices[key] = self._form.rowCount() - 1
        return w

    def add_spin(self, key, label, value=1, lo=1, hi=999):
        w = QSpinBox()
        w.setRange(lo, hi)
        w.setValue(value)
        self._form.addRow(label, w)
        self._widgets[key] = w
        self._row_indices[key] = self._form.rowCount() - 1
        return w

    def add_browse(self, key, label, default="", filter_str="All files (*.*)"):
        row = QHBoxLayout()
        le = QLineEdit(default)
        btn = QPushButton("...")
        btn.setFixedWidth(30)
        btn.clicked.connect(lambda: self._do_browse(le, filter_str))
        row.addWidget(le)
        row.addWidget(btn)
        container = QWidget()
        container.setLayout(row)
        row.setContentsMargins(0, 0, 0, 0)
        self._form.addRow(label, container)
        self._widgets[key] = le
        self._row_indices[key] = self._form.rowCount() - 1
        return le

    def _do_browse(self, line_edit, filter_str):
        path, _ = QFileDialog.getOpenFileName(self, "Select file", "", filter_str)
        if path:
            line_edit.setText(path)

    def val(self, key):
        w = self._widgets[key]
        if isinstance(w, QLineEdit):
            return w.text().strip()
        elif isinstance(w, QComboBox):
            return w.currentText()
        elif isinstance(w, QSpinBox):
            return str(w.value())
        return ""

    def save_state(self):
        state = {}
        for key, w in self._widgets.items():
            if isinstance(w, QLineEdit):
                state[key] = w.text()
            elif isinstance(w, QComboBox):
                state[key] = w.currentIndex()
            elif isinstance(w, QSpinBox):
                state[key] = w.value()
        return state

    def restore_state(self, state):
        if not state:
            return
        for key, val in state.items():
            w = self._widgets.get(key)
            if w is None:
                continue
            try:
                if isinstance(w, QLineEdit):
                    w.setText(str(val))
                elif isinstance(w, QComboBox):
                    w.setCurrentIndex(int(val))
                elif isinstance(w, QSpinBox):
                    w.setValue(int(val))
            except (ValueError, TypeError):
                pass

    def is_runnable(self):
        return True

    def build_command(self, vol_path):
        raise NotImplementedError


# ============================================================
# AnalysisWindow: main window with Run/Stop/GMSH
# ============================================================

class AnalysisWindow(QMainWindow):
    """Standalone analysis window with .vol input, parameter panel, Run/Stop."""

    def __init__(self, title, vol_path="", settings_key="default"):
        super().__init__()
        self.setWindowTitle(title)
        self.resize(650, 600)
        self.setMinimumSize(500, 400)
        self._settings_key = settings_key

        icon = _icon_path()
        if icon:
            self.setWindowIcon(QIcon(icon))

        self._process = None
        self._last_msh = None
        self._panel = None  # set by subclass via _set_panel()
        self._vol_path = vol_path

        self._build_ui(vol_path)

    def _set_panel(self, panel):
        """Set the parameter panel (call from subclass __init__)."""
        self._panel = panel
        self._panel_area.addWidget(panel)
        panel.validationChanged = self._update_run_state
        self._update_run_state()

    def _build_ui(self, vol_path):
        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(10, 10, 10, 5)

        # Model path
        model_row = QHBoxLayout()
        model_row.addWidget(QLabel("Model (.vol):"))
        self._vol_edit = QLineEdit(vol_path)
        self._vol_edit.setPlaceholderText(
            ".vol file (exported from Cubit or Netgen)")
        model_row.addWidget(self._vol_edit, 1)
        browse_btn = QPushButton("Browse...")
        browse_btn.clicked.connect(self._browse_vol)
        model_row.addWidget(browse_btn)
        root.addLayout(model_row)

        # Splitter: panel | output
        splitter = QSplitter(Qt.Vertical)
        root.addWidget(splitter, 1)

        # Panel area (subclass inserts panel here)
        self._panel_container = QWidget()
        self._panel_area = QVBoxLayout(self._panel_container)
        self._panel_area.setContentsMargins(0, 0, 0, 0)
        splitter.addWidget(self._panel_container)

        # Output
        out_group = QGroupBox("Output")
        out_layout = QVBoxLayout(out_group)
        out_layout.setContentsMargins(5, 5, 5, 5)
        self._output = QPlainTextEdit()
        self._output.setReadOnly(True)
        self._output.setFont(QFont("Consolas", 9))
        out_layout.addWidget(self._output)
        splitter.addWidget(out_group)
        splitter.setSizes([300, 250])

        # Buttons
        style = self.style()
        btn_row = QHBoxLayout()

        self._run_btn = QPushButton(
            style.standardIcon(QStyle.SP_MediaPlay), " Run")
        self._run_btn.setFixedHeight(32)
        self._run_btn.clicked.connect(self._on_run)
        btn_row.addWidget(self._run_btn)

        self._stop_btn = QPushButton(
            style.standardIcon(QStyle.SP_MediaStop), " Stop")
        self._stop_btn.setFixedHeight(32)
        self._stop_btn.setEnabled(False)
        self._stop_btn.clicked.connect(self._on_stop)
        btn_row.addWidget(self._stop_btn)

        self._gmsh_btn = QPushButton(
            style.standardIcon(QStyle.SP_ComputerIcon), " Open GMSH")
        self._gmsh_btn.setFixedHeight(32)
        self._gmsh_btn.setEnabled(False)
        self._gmsh_btn.clicked.connect(self._open_gmsh)
        btn_row.addWidget(self._gmsh_btn)

        btn_row.addStretch()
        root.addLayout(btn_row)

        # Status bar
        self._status = self.statusBar()

    def _browse_vol(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select .vol file",
            os.path.dirname(self._vol_edit.text()),
            "Netgen Vol (*.vol);;All (*)")
        if path:
            self._vol_edit.setText(path)

    def _update_run_state(self):
        if self._process is not None:
            return
        runnable = self._panel.is_runnable() if self._panel else True
        self._run_btn.setEnabled(runnable)

    def _on_run(self):
        if self._process is not None:
            return
        vol = self._vol_edit.text().strip()
        try:
            cmd = self._panel.build_command(vol)
        except ValueError as e:
            QMessageBox.warning(self, "Input Error", str(e))
            return

        self._save_settings()
        self._output.clear()
        self._output.appendPlainText(f"> {' '.join(cmd)}\n")
        self._last_msh = None
        self._gmsh_btn.setEnabled(False)

        work_dir = os.path.dirname(vol) if vol else os.getcwd()
        if not os.path.isdir(work_dir):
            work_dir = os.getcwd()

        self._process = QProcess(self)
        self._process.setWorkingDirectory(work_dir)
        self._process.readyReadStandardOutput.connect(self._read_stdout)
        self._process.readyReadStandardError.connect(self._read_stderr)
        self._process.finished.connect(self._on_finished)
        self._process.start(_PYTHON, cmd[1:])

        self._run_btn.setEnabled(False)
        self._stop_btn.setEnabled(True)
        self._status.showMessage("Running...")

    def _on_stop(self):
        if self._process:
            self._process.kill()
            self._process.waitForFinished(3000)

    def _read_stdout(self):
        data = self._process.readAllStandardOutput().data()
        text = data.decode("utf-8", errors="replace").rstrip("\n")
        if text:
            self._output.appendPlainText(text)

    def _read_stderr(self):
        data = self._process.readAllStandardError().data()
        text = data.decode("utf-8", errors="replace").rstrip("\n")
        if text:
            self._output.appendPlainText(text)

    def _on_finished(self, exit_code, exit_status):
        remaining = self._process.readAllStandardOutput().data()
        if remaining:
            self._output.appendPlainText(
                remaining.decode("utf-8", errors="replace").rstrip("\n"))

        stdout_text = self._output.toPlainText()
        self._process = None
        self._stop_btn.setEnabled(False)
        self._update_run_state()

        if exit_code != 0:
            self._status.showMessage(f"Error (exit code {exit_code})")
            self._output.appendPlainText(
                f"\n*** Process exited with code {exit_code}")
            return

        # Find .msh in output for GMSH button
        result = None
        for line in reversed(stdout_text.split("\n")):
            line = line.strip()
            if line.startswith("{"):
                try:
                    result = json.loads(line)
                except json.JSONDecodeError:
                    pass
                break
        if result and "msh_output" in result:
            msh = result["msh_output"]
            if os.path.isfile(msh):
                self._last_msh = msh
                self._gmsh_btn.setEnabled(True)

        self._status.showMessage("Done." if exit_code == 0 else "Failed.")

    def _open_gmsh(self):
        if not self._last_msh:
            return
        import subprocess
        subprocess.Popen(
            [_PYTHON, "-c",
             "import gmsh; gmsh.initialize(); "
             f"gmsh.merge(r'{self._last_msh}'); "
             "gmsh.fltk.run(); gmsh.finalize()"],
            creationflags=(0x08000000 if sys.platform == "win32" else 0))

    # Settings
    def _settings_path(self):
        return os.path.join(_SETTINGS_DIR, f"radia_{self._settings_key}.json")

    def _save_settings(self):
        os.makedirs(_SETTINGS_DIR, exist_ok=True)
        data = {
            "vol": self._vol_edit.text(),
            "panel": self._panel.save_state() if self._panel else {},
        }
        try:
            with open(self._settings_path(), "w", encoding="utf-8") as f:
                json.dump(data, f, indent=2)
        except OSError:
            pass

    def _restore_settings(self):
        p = self._settings_path()
        if not os.path.isfile(p):
            return
        try:
            with open(p, "r", encoding="utf-8") as f:
                data = json.load(f)
        except (OSError, json.JSONDecodeError):
            return
        if not self._vol_edit.text() and "vol" in data:
            self._vol_edit.setText(data["vol"])
        if self._panel:
            self._panel.restore_state(data.get("panel", {}))

    def closeEvent(self, event):
        try:
            self._save_settings()
        except Exception:
            pass
        if self._process:
            try:
                self._process.kill()
                self._process.waitForFinished(3000)
            except Exception:
                pass
        super().closeEvent(event)


def run_app(window_class, vol_path=""):
    """Entry point helper: create QApplication + window, exec."""
    app = QApplication(sys.argv)
    vol = vol_path or (sys.argv[1] if len(sys.argv) > 1 else "")
    window = window_class(vol)
    window.show()
    sys.exit(app.exec())
