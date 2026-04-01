"""
Register a single "Launch Radia App" button in Coreform Cubit.

This script runs inside Cubit on startup (via ~/.cubit) and adds
a toolbar button that saves the current model as .cub5 and launches
the standalone Radia app (Python 3.12 + tkinter).

Architecture:
  Cubit GUI (Python 3.10 + PySide6):
    1 button -> save .cub5 -> subprocess.run([python3.12, radia_app.py, path.cub5])

  Radia App (Python 3.12 + tkinter):
    All settings, computation, and result display
"""

import sys
import os
import subprocess

# Determine script location (__file__ is not defined when run via Cubit 'play')
try:
    _this_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    import importlib
    try:
        _radia_spec = importlib.util.find_spec("radia")
        if _radia_spec and _radia_spec.origin:
            _this_dir = os.path.join(os.path.dirname(_radia_spec.origin), "panels")
        else:
            raise ImportError
    except (ImportError, AttributeError):
        import glob as _glob
        for _base in [
            os.path.join(sys.prefix, "Lib", "site-packages"),
            os.path.join(sys.prefix, "lib", "python*", "site-packages"),
        ]:
            for _candidate in _glob.glob(os.path.join(_base, "radia", "panels")):
                if os.path.isdir(_candidate):
                    _this_dir = _candidate
                    break
            else:
                continue
            break
        else:
            _this_dir = os.getcwd()

import cubit

# Qt bindings: prefer PySide6, fall back to PyQt5 (Cubit ships PyQt5)
try:
    from PySide6.QtWidgets import (QApplication, QMainWindow, QMessageBox,
                                   QToolBar, QFileDialog, QDialog,
                                   QVBoxLayout, QHBoxLayout, QFormLayout,
                                   QComboBox, QSpinBox, QLineEdit,
                                   QPushButton, QDialogButtonBox)
    from PySide6.QtGui import QAction, QIcon
    from PySide6.QtCore import QSize
except ImportError:
    from PyQt5.QtWidgets import (QApplication, QMainWindow, QMessageBox,
                                 QToolBar, QAction, QFileDialog, QDialog,
                                 QVBoxLayout, QHBoxLayout, QFormLayout,
                                 QComboBox, QSpinBox, QLineEdit,
                                 QPushButton, QDialogButtonBox)
    from PyQt5.QtGui import QIcon
    from PyQt5.QtCore import QSize


def _find_external_python():
    """Find external Python 3.12 (not Cubit's bundled Python 3.10).

    Search order:
      1. RADIA_PYTHON environment variable
      2. py -3 (Windows Python Launcher)
      3. python (from PATH, skip if Cubit's Python)
    """
    radia_py = os.environ.get("RADIA_PYTHON")
    if radia_py and os.path.isfile(radia_py):
        return radia_py

    if sys.platform == "win32":
        try:
            result = subprocess.run(
                ["py", "-3", "-c", "import sys; print(sys.executable)"],
                capture_output=True, text=True, timeout=5
            )
            if result.returncode == 0:
                py_path = result.stdout.strip()
                if os.path.isfile(py_path):
                    return py_path
        except Exception:
            pass

    try:
        result = subprocess.run(
            ["python", "-c", "import sys; print(sys.executable)"],
            capture_output=True, text=True, timeout=5
        )
        if result.returncode == 0:
            py_path = result.stdout.strip()
            if os.path.isfile(py_path) and "Cubit" not in py_path:
                return py_path
    except Exception:
        pass

    return None


def _find_main_window():
    """Find Cubit's main QMainWindow."""
    app = QApplication.instance()
    if app is None:
        return None
    best = None
    best_count = 0
    for widget in app.topLevelWidgets():
        if isinstance(widget, QMainWindow):
            count = len(widget.menuBar().actions())
            if count > best_count:
                best_count = count
                best = widget
    return best


def _find_radia_app():
    """Find radia_app.py path."""
    # Same directory as this file -> ../radia_app.py
    pkg_root = os.path.dirname(_this_dir)  # src/radia/
    app_path = os.path.join(pkg_root, "radia_app.py")
    if os.path.isfile(app_path):
        return app_path
    # Fallback: try importlib to find installed radia package
    try:
        import importlib.util
        spec = importlib.util.find_spec("radia")
        if spec and spec.origin:
            app_path = os.path.join(os.path.dirname(spec.origin), "radia_app.py")
            if os.path.isfile(app_path):
                return app_path
    except Exception:
        pass
    return None


def _has_model():
    """Check if Cubit has geometry loaded."""
    try:
        return cubit.get_volume_count() > 0 or cubit.get_surface_count() > 0
    except Exception:
        return False


def _launch_radia_app():
    """Save current Cubit model as .cub5 and launch Radia app."""
    import tempfile

    ext_python = _find_external_python()
    if not ext_python:
        QMessageBox.critical(
            None, "Error",
            "External Python 3.12 not found.\n\n"
            "Set RADIA_PYTHON environment variable or install Python 3.12."
        )
        return

    # If no model, ask for .jou file
    if not _has_model():
        jou_path, _ = QFileDialog.getOpenFileName(
            None, "Select Journal File",
            os.getcwd(),
            "Cubit Journal (*.jou);;All Files (*)"
        )
        if not jou_path:
            return  # cancelled
        jou_path = jou_path.replace("\\", "/")
        cubit.cmd(f'play "{jou_path}"')

        if not _has_model():
            QMessageBox.warning(
                None, "Warning",
                "Journal file did not create any geometry."
            )
            return

    # Save .cub5 in current working directory
    cub5_path = os.path.join(os.getcwd(), "radia_cubit_model.cub5")
    cub5_path = cub5_path.replace("\\", "/")
    cubit.cmd(f'save cub5 "{cub5_path}" overwrite')

    # Launch standalone app (non-blocking)
    radia_app = _find_radia_app()
    work_dir = os.path.dirname(cub5_path)
    if radia_app:
        cmd = [ext_python, radia_app, cub5_path]
    else:
        cmd = [ext_python, "-m", "radia.radia_app", cub5_path]
    subprocess.Popen(cmd, cwd=work_dir)



class ExportMeshDialog(QDialog):
    """Dialog for mesh export via radia .ccm plugin commands.

    Shows all options per format. Command preview updates live.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Export Mesh")
        self.setMinimumWidth(450)

        layout = QVBoxLayout(self)
        form = QFormLayout()

        # Format
        self._format = QComboBox()
        self._format.addItems(["GMSH", "Nastran BDF", "VTK", "MEG"])
        self._format.currentTextChanged.connect(self._on_format_changed)
        form.addRow("Format:", self._format)

        # Order (GMSH/Nastran/VTK: 1-2, MEG: fixed 1)
        self._order = QSpinBox()
        self._order.setRange(1, 2)
        self._order.setValue(1)
        self._order.valueChanged.connect(self._update_command)
        form.addRow("Order:", self._order)

        # Version (GMSH only)
        self._version = QComboBox()
        self._version.addItems(["2.2", "4.1"])
        self._version.currentTextChanged.connect(self._update_command)
        self._version_row = ("Version:", self._version)
        form.addRow(*self._version_row)

        # Dimension (GMSH, Nastran, VTK)
        self._dimension = QComboBox()
        self._dimension.addItems(["3D", "2D"])
        self._dimension.currentTextChanged.connect(self._update_command)
        self._dimension_row = ("Dimension:", self._dimension)
        form.addRow(*self._dimension_row)

        # NoPyramid (Nastran only)
        self._nopyramid = QComboBox()
        self._nopyramid.addItems(["Keep pyramids", "Convert to degenerate hex (JMAG)"])
        self._nopyramid.currentTextChanged.connect(self._update_command)
        self._nopyramid_row = ("Pyramids:", self._nopyramid)
        form.addRow(*self._nopyramid_row)

        # File
        file_row = QHBoxLayout()
        self._filename = QLineEdit("mesh.msh")
        self._filename.textChanged.connect(self._update_command)
        browse = QPushButton("...")
        browse.setFixedWidth(30)
        browse.clicked.connect(self._browse)
        file_row.addWidget(self._filename)
        file_row.addWidget(browse)
        form.addRow("File:", file_row)

        layout.addLayout(form)

        # Command preview (read-only, shows the Cubit command)
        self._cmd_preview = QLineEdit()
        self._cmd_preview.setReadOnly(True)
        form_cmd = QFormLayout()
        form_cmd.addRow("Command:", self._cmd_preview)
        layout.addLayout(form_cmd)

        buttons = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        self._on_format_changed(self._format.currentText())

    def _on_format_changed(self, fmt):
        is_gmsh = (fmt == "GMSH")
        is_nastran = (fmt == "Nastran BDF")
        is_meg = (fmt == "MEG")

        # Order: MEG is always 1
        self._order.setEnabled(not is_meg)
        if is_meg:
            self._order.setValue(1)

        # Version: GMSH only
        self._version.setVisible(is_gmsh)
        self._version.setEnabled(is_gmsh)

        # Dimension: GMSH, Nastran, VTK (not MEG)
        has_dim = not is_meg
        self._dimension.setVisible(has_dim)
        self._dimension.setEnabled(has_dim)

        # NoPyramid: Nastran only
        self._nopyramid.setVisible(is_nastran)
        self._nopyramid.setEnabled(is_nastran)

        # File extension
        ext_map = {"GMSH": ".msh", "Nastran BDF": ".bdf", "VTK": ".vtk", "MEG": ".meg"}
        base = os.path.splitext(self._filename.text())[0]
        self._filename.setText(base + ext_map.get(fmt, ".msh"))

        self._update_command()

    def _browse(self):
        ext_map = {"GMSH": ".msh", "Nastran BDF": ".bdf", "VTK": ".vtk", "MEG": ".meg"}
        ext = ext_map.get(self._format.currentText(), ".*")
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Mesh",
            os.path.join(os.getcwd(), self._filename.text()),
            f"Mesh (*{ext});;All Files (*)")
        if path:
            self._filename.setText(path)

    def _update_command(self):
        self._cmd_preview.setText(self.get_cubit_command())

    def get_cubit_command(self):
        fmt = self._format.currentText()
        filename = self._filename.text().replace("\\", "/")
        order = self._order.value()

        if fmt == "GMSH":
            ver = "2" if self._version.currentText() == "2.2" else "4"
            dim = "2" if self._dimension.currentText() == "2D" else "3"
            cmd = f'radia export gmsh "{filename}" order {order} version {ver} dimension {dim}'
        elif fmt == "Nastran BDF":
            dim = "2" if self._dimension.currentText() == "2D" else "3"
            cmd = f'radia export nastran "{filename}" order {order} dimension {dim}'
            if self._nopyramid.currentIndex() == 1:
                cmd += " nopyramid"
        elif fmt == "VTK":
            dim = "2" if self._dimension.currentText() == "2D" else "3"
            cmd = f'radia export vtk "{filename}" order {order} dimension {dim}'
        elif fmt == "MEG":
            cmd = f'radia export meg "{filename}"'
        else:
            cmd = ""

        cmd += " overwrite"
        return cmd


def _export_mesh():
    """Show export dialog and run radia export command."""
    if not _has_model():
        QMessageBox.warning(None, "Warning", "No mesh to export.\nMesh the geometry first.")
        return
    dlg = ExportMeshDialog()
    if dlg.exec() == QDialog.Accepted:
        cmd = dlg.get_cubit_command()
        print(f"Running: {cmd}")
        cubit.cmd(cmd)


def register_menu():
    """Register 'Radia-NGSolve' menu in the menu bar."""
    main_window = _find_main_window()
    if main_window is None:
        return

    menu_bar = main_window.menuBar()

    # Remove existing Radia menu (for reload)
    is_reload = False
    for action in list(menu_bar.actions()):
        if action.text().replace("&", "") == "Radia-NGSolve":
            sub = action.menu()
            menu_bar.removeAction(action)
            if sub:
                sub.deleteLater()
            is_reload = True
            break

    app = QApplication.instance()
    if app is not None:
        app.setQuitOnLastWindowClosed(False)

    # Add Radia-NGSolve menu to the menu bar (top level)
    _icon_file = os.path.join(os.path.dirname(_this_dir), "resources", "radia_icon.png")
    radia_menu = menu_bar.addMenu("Radia-NGSolve")

    action_launch = QAction("Launch Radia App...", main_window)
    if os.path.isfile(_icon_file):
        action_launch.setIcon(QIcon(_icon_file))
    action_launch.setStatusTip(
        "Save model as .cub5 and launch Radia standalone app (Python 3.12)"
    )
    action_launch.triggered.connect(_launch_radia_app)
    radia_menu.addAction(action_launch)

    # Export Mesh
    action_export = QAction("Export Mesh...", main_window)
    action_export.setStatusTip("Export mesh (GMSH, Nastran, VTK, MEG) via radia plugin")
    action_export.triggered.connect(_export_mesh)
    radia_menu.addAction(action_export)

    # Separator + Reload (development)
    radia_menu.addSeparator()
    action_reload = QAction("Reload Panels", main_window)
    action_reload.setStatusTip("Re-read register_toolbar.py from disk")
    def _reload_panels():
        startup = os.path.join(_this_dir, "startup.py").replace("\\", "/")
        cubit.cmd(f'play "{startup}"')
    action_reload.triggered.connect(_reload_panels)
    radia_menu.addAction(action_reload)

    # Toolbar button with icon
    # Remove existing toolbar (for reload)
    for tb in main_window.findChildren(QToolBar, "RadiaToolBar"):
        main_window.removeToolBar(tb)
        tb.deleteLater()

    if os.path.isfile(_icon_file):
        toolbar = QToolBar("RadiaToolBar", main_window)
        toolbar.setObjectName("RadiaToolBar")
        toolbar.setIconSize(QSize(32, 32))
        tb_action = toolbar.addAction(QIcon(_icon_file), "Radia")
        tb_action.setToolTip("Launch Radia App (save .cub5 + open standalone app)")
        tb_action.triggered.connect(_launch_radia_app)
        main_window.addToolBar(toolbar)
        toolbar.setVisible(True)
        toolbar.show()

    if is_reload:
        print("Radia-NGSolve menu re-registered.")
    else:
        print("Radia-NGSolve menu registered.")


# Auto-register when this script is executed
register_menu()
