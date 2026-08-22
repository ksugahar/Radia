"""
Radia Export menu — PySide6 port of the legacy Qt5 .ccl plugin.

Replaces the in-tree C++ Claro component (RadiaComp.cpp) with a pure
Python (PySide6) toolbar that runs inside Cubit's embedded Python
interpreter.  The .ccm (APREPRO commands `export gmsh / nastran /
vtk / netgen / femeem / meg`) is unchanged — this module just calls
`cubit.cmd("export ...")` after collecting user input through Qt
dialogs parented to the Cubit main window.

Layer 2 (Cubit GUI Python).  See CLAUDE.md "Cubit Panel Architecture"
and `radia_mcp.cubit.knowledge.custom_toolbar` (TOOLBAR_PYTHON_SCRIPT_
CONVENTIONS, TOOLBAR_CLARO_PARENT_HELPER).

This module avoids:
  - parenthesized multi-line PySide6 imports (Cubit's importer chokes
    on them) — use backslash continuations
  - top-level `import cubit` — the host already imports it
  - dialogs without `parent=find_claro()` — they get hidden behind the
    Cubit main window on Windows

Qt binding: **PySide6 only** (radia 4.80.0).  Target is Cubit 2025.12
which ships PySide6.  No PyQt5 fallback — per CLAUDE.md "No Fallbacks
— Fail Fast, Fail Loud", an old Cubit without PySide6 must raise the
underlying ImportError loudly so the operator can fix the env.
"""

import json
import os
import sys

# Qt binding: PySide6 only (radia 4.80.0).  Target is Cubit 2025.12
# which ships PySide6.  Per CLAUDE.md "No Fallbacks — Fail Fast, Fail
# Loud" any ImportError propagates so the operator sees the real cause
# (old Cubit, broken PySide6 install, wrong Python env) instead of a
# silent PyQt5 fallback that masks it.
#
# Use backslash continuations (NOT parens) per Cubit toolbar convention.
# Note: QAction lives in QtGui in Qt6 / PySide6 (was QtWidgets in Qt5).
from PySide6.QtWidgets import QApplication, QCheckBox, QComboBox, \
    QDialog, QDialogButtonBox, QFileDialog, QFormLayout, QGroupBox, \
    QHBoxLayout, QLabel, QLineEdit, QMainWindow, QMenu, \
    QMessageBox, QPushButton, QTableWidget, QTableWidgetItem, \
    QVBoxLayout
from PySide6.QtCore import Qt
from PySide6.QtGui import QAction


# ----------------------------------------------------------------------
# Settings persistence
# ----------------------------------------------------------------------

def _settings_dir():
    """Per-user settings directory (matches the legacy C++ .ccl path)."""
    appdata = os.path.join(os.path.expanduser("~"),
                           "AppData", "Roaming", "Radia")
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
    """Return Cubit's main QMainWindow (objectName == 'claro'), or None.

    Filter by objectName, never windowTitle (locale-dependent and
    changes when a file is loaded).  Used as `parent=` for every
    QDialog so the dialog cannot fall behind the main window.
    """
    app = QApplication.instance()
    if app is None or not hasattr(app, "topLevelWidgets"):
        return None
    for w in app.topLevelWidgets():
        if isinstance(w, QMainWindow) and w.objectName() == "claro":
            return w
    # Fallback: any QMainWindow (helps in headless / test contexts)
    for w in app.topLevelWidgets():
        if isinstance(w, QMainWindow):
            return w
    return None


# ----------------------------------------------------------------------
# Journal-path helper
# ----------------------------------------------------------------------

# Per-session sticky journal path (matches legacy `s_lastJouPath` in .ccl)
_last_jou_path = ""


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


def ensure_jou_path(cubit_mod, parent):
    """Ensure a .jou path is known.  Returns the path or "" if cancelled.

    Three-step resolution (matches RadiaComp.cpp::ensure_jou_path):
      1. Cubit's own current journal
      2. The session-sticky path from a previous export
      3. Prompt user with QFileDialog::getSaveFileName + `save journal`
    """
    global _last_jou_path

    # 1. Cubit's own current journal file
    try:
        jf = cubit_mod.get_current_journal_file()
    except Exception:
        jf = ""
    if jf:
        p = jf.replace("\\", "/")
        if _is_valid_jou_path(p):
            return p

    # 2. Session-sticky path
    if _last_jou_path and _is_valid_jou_path(_last_jou_path):
        return _last_jou_path

    # 3. Prompt
    default_dir = _default_start_dir() or os.getcwd().replace("\\", "/")
    initial = default_dir.rstrip("/") + "/model.jou"
    jou, _filter = QFileDialog.getSaveFileName(
        parent, "Save Journal (determines output filenames)",
        initial, "Cubit Journal (*.jou);;All Files (*)")
    if not jou:
        return ""
    jou = jou.replace("\\", "/")
    os.makedirs(os.path.dirname(jou), exist_ok=True)

    # Save journal only (no .cub5)
    cubit_mod.cmd(f'save journal "{jou}" overwrite')
    _last_jou_path = jou
    print(f"Saved journal: {jou}")
    return jou


# ----------------------------------------------------------------------
# ExportDialog — format-aware QDialog (port of ExportDialog in
# RadiaComp.cpp).
# ----------------------------------------------------------------------

# Format keys — match the legacy enum order so saved settings remain
# compatible (settings dict keys: netgen_vol / gmsh / nastran / vtk /
# femeem / meg).
FMT_NETGEN = "netgen_vol"
FMT_GMSH = "gmsh"
FMT_NASTRAN = "nastran"
FMT_VTK = "vtk"
FMT_FEMEEM = "femeem"
FMT_MEG = "meg"

_FORMAT_TITLES = {
    FMT_NETGEN: "Export Netgen Vol",
    FMT_GMSH: "Export GMSH",
    FMT_NASTRAN: "Export Nastran BDF",
    FMT_VTK: "Export VTK",
    FMT_FEMEEM: "Export FEMEEM",
    FMT_MEG: "Export MEG (ELF/MAGIC)",
}

_FORMAT_EXTS = {
    FMT_NETGEN: ".vol",
    FMT_GMSH: ".msh",
    FMT_NASTRAN: ".bdf",
    FMT_VTK: ".vtk",
    FMT_FEMEEM: "",  # directory output
    FMT_MEG: ".meg",
}


def _elf_labels_for_dim(dim_idx):
    """ELF label options per Dimension dropdown index.

    dim_idx: 0=3D(T), 1=2D(K), 2=Axisym(R)
    "(none)" = skip the block during MEG export (e.g. air regions).
    """
    if dim_idx == 0:
        return ["(none) - Skip export",
                "MMB - Magnetic body",
                "MMS - Magnetic shell",
                "MMT - Magnetic thin",
                "MMP - Permanent magnet",
                "MWL - Nonlinear magnet (fixed local axis)",
                "MWV - Nonlinear magnet (direction vectors)",
                "MCL - Coil line",
                "MCO - Current conductor",
                "MAB - Armature block",
                "MAT - Armature thin",
                "MBB - Boundary"]
    return ["(none) - Skip export",
            "MMB - Magnetic body",
            "MMT - Magnetic thin",
            "MMP - Permanent magnet",
            "MWL - Nonlinear magnet (fixed local axis)",
            "MWV - Nonlinear magnet (direction vectors)",
            "MCL - Coil line",
            "MCO - Current conductor",
            "MAB - Armature block",
            "MAT - Armature thin",
            "MBB - Boundary"]


def _elf_prefix(text):
    """Extract 3-char prefix from combo text 'MMB - Magnetic body'."""
    return text[:3]


class ExportDialog(QDialog):
    """Format-aware export options dialog.

    Builds the same set of widgets as the legacy C++ ExportDialog so
    the same `export` APREPRO command is emitted.  Per-format
    settings are persisted to ``%APPDATA%/Radia/export_settings.json``
    under the format key (``netgen_vol`` / ``gmsh`` / ``nastran`` /
    ``vtk`` / ``femeem`` / ``meg``).
    """

    def __init__(self, fmt, jou_path, cubit_mod, parent=None):
        super().__init__(parent)
        self._fmt = fmt
        self._cubit = cubit_mod
        self.setWindowTitle(_FORMAT_TITLES[fmt])
        self.setMinimumWidth(500)

        layout = QVBoxLayout(self)
        form = QFormLayout()

        # Default directory + basename derive from .jou path
        default_dir = ""
        if fmt == FMT_FEMEEM:
            base_name = "femeem_output"
        else:
            base_name = "ExportedMesh"

        if jou_path:
            jp = jou_path.replace("\\", "/")
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
        if fmt == FMT_FEMEEM:
            # FEMEEM: single directory path
            row = QHBoxLayout()
            self._dir = QLineEdit(default_dir + "/" + base_name)
            browse = QPushButton("...")
            browse.setFixedWidth(30)
            browse.clicked.connect(self._browse_dir)
            self._dir.textChanged.connect(self._update_preview)
            row.addWidget(self._dir)
            row.addWidget(browse)
            form.addRow("Directory:", row)
            self._filename = QLineEdit("")
            self._filename.setVisible(False)
        else:
            row = QHBoxLayout()
            self._dir = QLineEdit(default_dir)
            browse = QPushButton("...")
            browse.setFixedWidth(30)
            browse.clicked.connect(self._browse_dir)
            self._dir.textChanged.connect(self._update_preview)
            row.addWidget(self._dir)
            row.addWidget(browse)
            form.addRow("Directory:", row)
            self._filename = QLineEdit(base_name + _FORMAT_EXTS[fmt])
            self._filename.textChanged.connect(self._update_preview)
            form.addRow("Filename:", self._filename)

        # --- Order combo ---------------------------------------------
        self._order = QComboBox()
        if fmt in (FMT_FEMEEM, FMT_MEG):
            # Always order 1 — no form row
            self._order.addItem("1 (linear)")
        elif fmt == FMT_NETGEN:
            self._order.addItems(["1", "2", "3", "4", "5"])
            self._order.setCurrentIndex(2)  # default 3
            self._order.currentIndexChanged.connect(self._update_preview)
            form.addRow("Order:", self._order)
        elif fmt == FMT_GMSH:
            # tet/hex up to order 3, wedge up to order 2.
            # order 4-5 unsupported in GMSH writer (see CLAUDE.md
            # "GMSH order limit").
            self._order.addItems(
                ["1", "2", "3 (tet/hex only, wedge stays 2)"])
            self._order.setCurrentIndex(0)
            self._order.currentIndexChanged.connect(self._update_preview)
            form.addRow("Order:", self._order)
        else:
            self._order.addItems(["1 (linear)", "2 (quadratic)"])
            self._order.currentIndexChanged.connect(self._update_preview)
            form.addRow("Order:", self._order)

        # --- Scale (FEMEEM only) -------------------------------------
        self._scale = None
        if fmt == FMT_FEMEEM:
            self._scale = QLineEdit("1.0")
            self._scale.textChanged.connect(self._update_preview)
            form.addRow("Scale:", self._scale)

        # --- Dimension ----------------------------------------------
        self._dimension = QComboBox()
        if fmt == FMT_MEG:
            self._dimension.addItems(
                ["3D (threed)", "2D (twod)", "Axisymmetric"])
        else:
            self._dimension.addItems(["3D", "2D"])
        self._dimension.currentIndexChanged.connect(self._update_preview)
        if fmt not in (FMT_NETGEN, FMT_FEMEEM):
            form.addRow("Dimension:", self._dimension)

        # --- NoPyramid (Nastran only) -------------------------------
        self._nopyramid = None
        if fmt == FMT_NASTRAN:
            self._nopyramid = QComboBox()
            self._nopyramid.addItems(
                ["Keep pyramids",
                 "Convert to degenerate hex (JMAG)"])
            self._nopyramid.currentIndexChanged.connect(self._update_preview)
            form.addRow("Pyramids:", self._nopyramid)

        # --- MEG block label table ----------------------------------
        self._block_table = None
        if fmt == FMT_MEG:
            layout.addLayout(form)
            layout.addWidget(QLabel("Block -> ELF Label:"))
            self._block_table = QTableWidget(0, 3)
            self._block_table.setHorizontalHeaderLabels(
                ["Block", "Name", "ELF Type"])
            self._block_table.horizontalHeader().setStretchLastSection(True)
            self._block_table.setColumnWidth(0, 50)
            self._block_table.setColumnWidth(1, 140)
            self._block_table.setMinimumHeight(120)
            self._block_table.setMaximumHeight(250)
            layout.addWidget(self._block_table)
            self._populate_block_table()
            # Re-populate labels when DIM changes
            self._dimension.currentIndexChanged.connect(self._update_preview)
        else:
            layout.addLayout(form)

        # --- Kelvin / symmetry (Netgen .vol only) -------------------
        self._kelvin_enable = None
        self._kelvin_mesh_size = None
        self._kelvin_sym_x = None
        self._kelvin_sym_y = None
        self._kelvin_sym_z = None
        if fmt == FMT_NETGEN:
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
                    "  bn = B.n=0   (flux parallel,   Radia '+',  "
                    "A:Dirichlet, Omega:natural)\n"
                    "  ht = HxN=0   (flux perp,       Radia '-',  "
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
        # When jou_path was supplied (Cubit Solve menu launch with a
        # current journal), the .jou-derived defaultDir wins over the
        # previous session's directory.
        all_settings = _load_settings()
        s = all_settings.get(fmt, {})
        if not jou_path and "dir" in s:
            self._dir.setText(s["dir"])
        if "order" in s:
            self._order.setCurrentIndex(int(s["order"]))
        if "dimension" in s:
            self._dimension.setCurrentIndex(int(s["dimension"]))
        if self._nopyramid and "nopyramid" in s:
            self._nopyramid.setCurrentIndex(int(s["nopyramid"]))
        if self._scale and "scale" in s:
            self._scale.setText(s["scale"])
        if self._kelvin_enable and "kelvin_enable" in s:
            self._kelvin_enable.setChecked(bool(s["kelvin_enable"]))
        if self._kelvin_mesh_size and "kelvin_mesh_size" in s:
            self._kelvin_mesh_size.setText(s["kelvin_mesh_size"])
        for cb, key in (
                (self._kelvin_sym_x, "kelvin_sym_x"),
                (self._kelvin_sym_y, "kelvin_sym_y"),
                (self._kelvin_sym_z, "kelvin_sym_z")):
            if cb and key in s:
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
        if self._fmt == FMT_FEMEEM:
            return d.rstrip("/")
        if not d.endswith("/"):
            d += "/"
        return d + self._filename.text()

    # ------------------------------------------------------------------
    def order(self):
        return self._order.currentIndex() + 1

    # ------------------------------------------------------------------
    def cubit_command(self):
        """Build the `export ...` APREPRO command string."""
        f = self.file_path()
        order = self.order()
        fmt = self._fmt
        cmd = ""

        if fmt == FMT_NETGEN:
            cmd = f'export netgen "{f}" order {order}'
            if self._kelvin_enable and self._kelvin_enable.isChecked():
                cmd += " add_kelvin"
                if self._kelvin_mesh_size:
                    ms = self._kelvin_mesh_size.text().strip()
                    if ms:
                        try:
                            v = float(ms)
                            if v > 0.0:
                                cmd += f" kelvin_mesh {v:g}"
                        except ValueError:
                            pass

                def _emit_sym(cb, flag):
                    if not cb:
                        return ""
                    s = cb.currentText()
                    if s in ("bn", "ht"):
                        return f" {flag} {s}"
                    return ""

                cmd += _emit_sym(self._kelvin_sym_x, "kelvin_sym_x")
                cmd += _emit_sym(self._kelvin_sym_y, "kelvin_sym_y")
                cmd += _emit_sym(self._kelvin_sym_z, "kelvin_sym_z")

        elif fmt == FMT_GMSH:
            dim = "2" if self._dimension.currentText() == "2D" else "3"
            cmd = (f'export gmsh "{f}" order {order} '
                   f'dimension {dim}')

        elif fmt == FMT_NASTRAN:
            dim = "2" if self._dimension.currentText() == "2D" else "3"
            cmd = (f'export jmag_nastran "{f}" order {order} '
                   f'dimension {dim}')
            if self._nopyramid and self._nopyramid.currentIndex() == 1:
                cmd += " nopyramid"

        elif fmt == FMT_VTK:
            dim = "2" if self._dimension.currentText() == "2D" else "3"
            cmd = (f'export vtk "{f}" order {order} '
                   f'dimension {dim}')

        elif fmt == FMT_FEMEEM:
            sc = self._scale.text() if self._scale else "1.0"
            cmd = f'export femeem "{f}" scale {sc}'

        elif fmt == FMT_MEG:
            dim_idx = self._dimension.currentIndex()
            dim_opt = ("threed", "twod", "axisymmetric")[dim_idx]
            cmd = f'export meg "{f}" {dim_opt}'
            if self._block_table and self._block_table.rowCount() > 0:
                pairs = []
                for r in range(self._block_table.rowCount()):
                    id_item = self._block_table.item(r, 0)
                    combo = self._block_table.cellWidget(r, 2)
                    if id_item and combo:
                        prefix = _elf_prefix(combo.currentText())
                        if prefix == "(no":  # "(none) - Skip export"
                            continue
                        pairs.append(f"{id_item.text()}:{prefix}")
                if pairs:
                    cmd += f' labels "{",".join(pairs)}"'

        return cmd + " overwrite"

    # ------------------------------------------------------------------
    def _populate_block_table(self):
        if self._block_table is None:
            return
        block_ids = self._cubit.get_block_id_list()
        self._block_table.setRowCount(len(block_ids))
        dim_idx = self._dimension.currentIndex() if self._dimension else 0
        labels = _elf_labels_for_dim(dim_idx)
        for r, bid in enumerate(block_ids):
            name = self._cubit.get_block_name(bid) or "(unnamed)"

            id_item = QTableWidgetItem(str(bid))
            id_item.setFlags(id_item.flags() & ~Qt.ItemIsEditable)
            self._block_table.setItem(r, 0, id_item)

            name_item = QTableWidgetItem(name)
            name_item.setFlags(name_item.flags() & ~Qt.ItemIsEditable)
            self._block_table.setItem(r, 1, name_item)

            combo = QComboBox()
            combo.addItems(labels)
            qname = name.upper()
            default_idx = 1  # MMB
            if (qname in ("AIR", "KELVIN")
                    or qname.startswith("AIR_")
                    or qname.startswith("KELVIN_")):
                default_idx = 0  # (none)
            else:
                for i in range(1, len(labels)):
                    if qname.startswith(_elf_prefix(labels[i])):
                        default_idx = i
                        break
            combo.setCurrentIndex(default_idx)
            combo.currentIndexChanged.connect(self._update_preview)
            self._block_table.setCellWidget(r, 2, combo)

    def _update_block_labels(self):
        if self._block_table is None or self._dimension is None:
            return
        labels = _elf_labels_for_dim(self._dimension.currentIndex())
        for r in range(self._block_table.rowCount()):
            combo = self._block_table.cellWidget(r, 2)
            if combo is None:
                continue
            current = _elf_prefix(combo.currentText())
            combo.blockSignals(True)
            combo.clear()
            combo.addItems(labels)
            for i, lbl in enumerate(labels):
                if _elf_prefix(lbl) == current:
                    combo.setCurrentIndex(i)
                    break
            combo.blockSignals(False)

    # ------------------------------------------------------------------
    def _update_preview(self):
        if self._fmt == FMT_MEG:
            self._update_block_labels()
        self._preview.setText(self.cubit_command())
        self._preview.setCursorPosition(0)

    # ------------------------------------------------------------------
    def _on_accept(self):
        # Persist settings before accepting
        all_settings = _load_settings()
        s = {
            "dir": self._dir.text(),
            "order": self._order.currentIndex(),
            "dimension": self._dimension.currentIndex(),
        }
        if self._nopyramid:
            s["nopyramid"] = self._nopyramid.currentIndex()
        if self._scale:
            s["scale"] = self._scale.text()
        if self._kelvin_enable:
            s["kelvin_enable"] = self._kelvin_enable.isChecked()
        if self._kelvin_mesh_size:
            s["kelvin_mesh_size"] = self._kelvin_mesh_size.text()
        if self._kelvin_sym_x:
            s["kelvin_sym_x"] = self._kelvin_sym_x.currentText()
        if self._kelvin_sym_y:
            s["kelvin_sym_y"] = self._kelvin_sym_y.currentText()
        if self._kelvin_sym_z:
            s["kelvin_sym_z"] = self._kelvin_sym_z.currentText()
        all_settings[self._fmt] = s
        _save_settings(all_settings)
        self.accept()


# ----------------------------------------------------------------------
# Generic export runner — applies to GMSH / Nastran / VTK / FEMEEM / MEG
# (Netgen Vol has its own runner with subprocess verification, see below.)
# ----------------------------------------------------------------------

def _run_generic_export(fmt, cubit_mod, parent):
    """Show ExportDialog and execute the chosen `export` command."""
    global _last_jou_path

    if cubit_mod.get_volume_count() == 0:
        jou = _ensure_model(cubit_mod, parent)
        if cubit_mod.get_volume_count() == 0:
            return
        if jou:
            _last_jou_path = jou

    jou_path = ensure_jou_path(cubit_mod, parent)
    if not jou_path:
        return

    dlg = ExportDialog(fmt, jou_path, cubit_mod, parent=parent)
    if dlg.exec() != QDialog.Accepted:
        return

    out_dir = os.path.dirname(dlg.file_path()).replace("\\", "/")
    cubit_mod.cmd(f'cd "{out_dir}"')

    cmd = dlg.cubit_command()
    if not cmd:
        QMessageBox.warning(parent, "Export", "No export command.")
        return

    cubit_mod.cmd(cmd)

    out = dlg.file_path()
    exists = os.path.isdir(out) if fmt == FMT_FEMEEM else os.path.isfile(out)
    if not exists:
        QMessageBox.warning(
            parent, "Export failed",
            f"Export did not produce expected output:\n{out}")
        return
    print(f"Export complete: {out}")

    # Open exported file with OS-associated application
    _open_in_os(out)


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
# Netgen Vol runner — exports + runs calc_verify_vol.py subprocess and
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


def _find_external_python():
    """Find external Python 3.12 (NOT Cubit's bundled Python 3.10)."""
    env = os.environ.get("RADIA_PYTHON")
    if env:
        return env
    if sys.platform == "win32":
        import subprocess
        try:
            r = subprocess.run(
                ["py", "-3", "-c", "print('ok')"],
                capture_output=True, text=True, timeout=5)
            if r.returncode == 0:
                return "py"
        except (FileNotFoundError, subprocess.TimeoutExpired):
            pass
    return "python"


def _find_calc_script(python, module_name):
    """Resolve panels/<module_name>.py via the external Python 3.12."""
    import subprocess
    args = []
    if python == "py":
        args.append("-3")
    args += [
        "-c",
        f"import radia.panels.{module_name} as m; import os; "
        f"print(os.path.abspath(m.__file__))",
    ]
    try:
        r = subprocess.run([python] + args,
                           capture_output=True, text=True, timeout=10)
        if r.returncode == 0:
            return r.stdout.strip()
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass
    return ""


def _run_netgen_export(cubit_mod, parent):
    """Export .vol via APREPRO + run NGSolve subprocess verification."""
    global _last_jou_path

    if cubit_mod.get_volume_count() == 0:
        jou = _ensure_model(cubit_mod, parent)
        if cubit_mod.get_volume_count() == 0:
            return
        if jou:
            _last_jou_path = jou

    jou_path = ensure_jou_path(cubit_mod, parent)
    if not jou_path:
        return

    dlg = ExportDialog(FMT_NETGEN, jou_path, cubit_mod, parent=parent)
    if dlg.exec() != QDialog.Accepted:
        return

    order = dlg.order()
    vol_path = dlg.file_path().replace("\\", "/")
    out_dir = os.path.dirname(vol_path).replace("\\", "/")
    cubit_mod.cmd(f'cd "{out_dir}"')

    # --- Phase 1: C++ export .vol + companion JSON ---
    print(f"Exporting Netgen .vol (order {order})...")
    cubit_mod.cmd(dlg.cubit_command())

    if not os.path.isfile(vol_path):
        QMessageBox.warning(
            parent, "Export failed",
            f"export netgen command did not produce:\n{vol_path}")
        return
    print(f"Exported: {vol_path}")
    _open_in_os(vol_path)

    # --- Phase 2: NGSolve subprocess to verify .vol ---
    python = _find_external_python()
    script = _find_calc_script(python, "calc_verify_vol")
    if not script:
        QMessageBox.information(
            parent, "Verification skipped",
            "calc_verify_vol.py not found. Skipping NGSolve "
            "verification.")
        return

    print("Verifying mesh with NGSolve...")

    import subprocess
    args = []
    if python == "py":
        args.append("-3")
    args += [script, "--vol", vol_path]

    try:
        rc, stdout, stderr = _run_subprocess_utf8(
            [python] + args, timeout=300)
    except subprocess.TimeoutExpired:
        QMessageBox.warning(parent, "Verification",
                            "NGSolve verification timed out.")
        return

    if rc != 0:
        QMessageBox.warning(
            parent, "Verification failed",
            f"NGSolve verification failed:\n\n{stderr[:2000]}")
        return

    lines = [ln for ln in stdout.splitlines() if ln.strip()]
    if not lines:
        QMessageBox.information(parent, "Verification",
                                "No verification output.")
        return
    try:
        r = json.loads(lines[-1])
    except json.JSONDecodeError:
        QMessageBox.warning(parent, "Verification",
                            "Invalid verification JSON.")
        return

    if "error" in r:
        QMessageBox.warning(parent, "Verification",
                            f"Verification error: {r['error']}")
        return

    _show_netgen_result(r, vol_path, order, parent)


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


def launch_export(fmt):
    """Launch one export dialog from Cubit's official custom toolbar.

    ``WorkflowToolbar`` buttons execute small Python scripts.  Those scripts
    call this stable dispatcher instead of modifying Cubit's QMenuBar.
    """
    supported = {
        FMT_NETGEN, FMT_GMSH, FMT_NASTRAN, FMT_VTK, FMT_FEMEEM, FMT_MEG,
    }
    if fmt not in supported:
        raise ValueError(f"unsupported Radia export format: {fmt!r}")
    import cubit as cubit_mod
    parent = find_claro()
    if parent is None:
        raise RuntimeError("Cubit main window is not available")
    if fmt == FMT_NETGEN:
        return _run_netgen_export(cubit_mod, parent)
    return _run_generic_export(fmt, cubit_mod, parent)


# ----------------------------------------------------------------------
# Menu installation
# ----------------------------------------------------------------------

# Sentinel so we don't add duplicate menus when the toolbar script
# re-runs (Cubit replays the startup script if the user reloads the
# Python interpreter from the GUI).
_INSTALLED_OBJECT_NAME = "RadiaExportMenu"

def _find_installed_menu(main):
    """Return the installed Radia menu, or None if Cubit removed it."""
    for action in list(main.menuBar().actions()):
        sub = action.menu()
        if sub is not None and sub.objectName() == _INSTALLED_OBJECT_NAME:
            return sub
    return None


def _install_menu_on_main(cubit_mod, main):
    """Create the menu after Cubit's QMainWindow is known to exist."""
    menu_bar = main.menuBar()

    # Remove any previous instance (idempotent reload)
    for action in list(menu_bar.actions()):
        sub = action.menu()
        if sub is not None and sub.objectName() == _INSTALLED_OBJECT_NAME:
            menu_bar.removeAction(action)
            sub.deleteLater()

    menu = QMenu("&Radia Export", main)
    menu.setObjectName(_INSTALLED_OBJECT_NAME)

    def _action(label, status_tip, callback):
        act = QAction(label, main)
        act.setStatusTip(status_tip)
        act.triggered.connect(callback)
        menu.addAction(act)
        return act

    # Wrap the callbacks with a parent= the Cubit main window so any
    # QDialog they create gets parented to it.
    _action(
        "Netgen Vol (.vol)...",
        "Export curved mesh with labels (.vol, order 1-5)",
        lambda: _run_netgen_export(cubit_mod, main))
    _action(
        "GMSH...",
        "Export mesh to GMSH format (.msh)",
        lambda: _run_generic_export(FMT_GMSH, cubit_mod, main))
    _action(
        "Nastran BDF...",
        "Export mesh to Nastran BDF format (.bdf)",
        lambda: _run_generic_export(FMT_NASTRAN, cubit_mod, main))
    _action(
        "VTK...",
        "Export mesh to VTK format (.vtk)",
        lambda: _run_generic_export(FMT_VTK, cubit_mod, main))
    _action(
        "FEMEEM...",
        "Export mesh to FEMEEM format (1st-order tet, directory output)",
        lambda: _run_generic_export(FMT_FEMEEM, cubit_mod, main))
    _action(
        "MEG (ELF/MAGIC)...",
        "Export mesh to ELF/MAGIC MEG format (1st-order)",
        lambda: _run_generic_export(FMT_MEG, cubit_mod, main))

    menu_bar.addMenu(menu)
    print("[Radia] Radia Export menu installed "
          "(6 export actions, PySide6 toolbar)")
    return menu


def install_menu():
    """Install the deprecated, non-persistent top-level menu manually.

    Radia's supported persistent Cubit surface is the ``WorkflowToolbar``
    package installed through Custom Toolbar Editor.  This compatibility helper
    remains for interactive use, but startup code must not depend on it because
    Cubit may replace its QMenuBar after ``~/.cubit`` has run.
    """
    # Cubit imports `cubit` at host startup; re-importing here is safe
    # (Python caches it).  Doing the import inside the function rather
    # than at module top means the module can still be imported under
    # `pytest` for unit tests that don't need Cubit.
    try:
        import cubit as cubit_mod
    except ImportError:
        print("[Radia] cubit module not available - "
              "Radia Export menu not installed")
        return None

    main = find_claro()
    if main is None:
        print("[Radia] Cubit main window not found - "
              "Radia Export compatibility menu not installed")
        return None

    return _install_menu_on_main(cubit_mod, main)


# Compatibility entry point for an explicit interactive ``play``.  Normal
# installation uses the official WorkflowToolbar package instead.
if __name__ == "_coreform_cubit":
    install_menu()
