"""
Register custom toolbar buttons in Coreform Cubit.

This script runs inside Cubit on startup (via ~/.cubit) and adds
toolbar buttons using PySide6 for mesh export functions.
"""

import sys
import os
import json
import subprocess
import tempfile

# Determine script location (__file__ is not defined when run via Cubit 'play')
try:
	_this_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
	# Fallback: find panels/ directory via radia package or known path
	import importlib
	try:
		_radia_spec = importlib.util.find_spec("radia")
		if _radia_spec and _radia_spec.origin:
			_this_dir = os.path.join(os.path.dirname(_radia_spec.origin), "panels")
		else:
			raise ImportError
	except (ImportError, AttributeError):
		# Last resort: search site-packages for radia/panels/
		import glob as _glob
		for _base in [
			os.path.join(sys.prefix, "Lib", "site-packages"),
			os.path.join(sys.prefix, "lib", "python*", "site-packages"),
			os.path.join(os.path.expanduser("~"), ".local", "lib", "python*", "site-packages"),
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

_pkg_root = os.path.dirname(_this_dir)  # radia package root (src/radia/)
if _pkg_root not in sys.path:
	sys.path.insert(0, _pkg_root)

import cubit

# Qt bindings: prefer PySide6, fall back to PyQt5 (Cubit ships PyQt5)
try:
	from PySide6.QtWidgets import (
		QApplication, QDialog, QVBoxLayout, QHBoxLayout, QGridLayout,
		QLabel, QLineEdit, QComboBox, QPushButton, QSpinBox, QCheckBox,
		QFileDialog, QMainWindow, QMessageBox,
		QTableWidget, QTableWidgetItem, QHeaderView, QGroupBox,
		QPlainTextEdit,
	)
	from PySide6.QtGui import QAction, QFont
	from PySide6.QtCore import Qt, QProcess
except ImportError:
	from PyQt5.QtWidgets import (
		QApplication, QDialog, QVBoxLayout, QHBoxLayout, QGridLayout,
		QLabel, QLineEdit, QComboBox, QPushButton, QSpinBox, QCheckBox,
		QFileDialog, QMainWindow, QMessageBox,
		QTableWidget, QTableWidgetItem, QHeaderView, QGroupBox,
		QPlainTextEdit, QAction,
	)
	from PyQt5.QtGui import QFont
	from PyQt5.QtCore import Qt, QProcess


def _find_main_window():
	"""Find Cubit's main QMainWindow (the one with the most menu items)."""
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


def _find_external_python():
	"""Find external Python with NGSolve (not Cubit's bundled Python).

	Search order:
	  1. RADIA_PYTHON environment variable
	  2. py -3 (Windows Python Launcher)
	  3. python (from PATH, skip if it's Cubit's Python)
	"""
	# 1. Explicit env var
	radia_py = os.environ.get("RADIA_PYTHON")
	if radia_py and os.path.isfile(radia_py):
		return radia_py

	# 2. Windows Python Launcher
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

	# 3. python from PATH
	try:
		result = subprocess.run(
			["python", "-c", "import sys; print(sys.executable)"],
			capture_output=True, text=True, timeout=5
		)
		if result.returncode == 0:
			py_path = result.stdout.strip()
			# Skip if it's Cubit's bundled Python
			if os.path.isfile(py_path) and "Cubit" not in py_path:
				return py_path
	except Exception:
		pass

	return None


def _parse_json_output(process):
	"""Parse JSON from QProcess stdout. Finds the line starting with '{'."""
	stdout = process.readAllStandardOutput().data().decode("utf-8", errors="replace")
	for line in stdout.strip().split("\n"):
		line = line.strip()
		if line.startswith("{"):
			try:
				return json.loads(line)
			except json.JSONDecodeError:
				continue
	QMessageBox.critical(None, "Error", f"No JSON output found:\n{stdout[:500]}")
	return None


def _get_all_volume_ids():
	"""Get all volume IDs in the model."""
	return list(cubit.get_entities("volume"))


# ================================================================
# Export Gmsh Dialog
# ================================================================

class ExportGmshDialog(QDialog):
	"""Dialog for Gmsh export parameters."""

	def __init__(self, parent=None):
		super().__init__(parent)
		self.setWindowTitle("Export Gmsh")
		self.setMinimumWidth(400)
		self._setup_ui()

	def _setup_ui(self):
		layout = QVBoxLayout(self)

		# File name
		row1 = QHBoxLayout()
		row1.addWidget(QLabel("File Name:"))
		self.file_edit = QLineEdit()
		self.file_edit.setPlaceholderText("output.msh")
		row1.addWidget(self.file_edit)
		self.browse_btn = QPushButton("...")
		self.browse_btn.setFixedWidth(30)
		self.browse_btn.clicked.connect(self._browse_file)
		row1.addWidget(self.browse_btn)
		layout.addLayout(row1)

		# Version
		row2 = QHBoxLayout()
		row2.addWidget(QLabel("Version:"))
		self.version_combo = QComboBox()
		self.version_combo.addItems(["2.2", "4.1"])
		row2.addWidget(self.version_combo)
		layout.addLayout(row2)

		# DIM (v4.1 only)
		row3 = QHBoxLayout()
		row3.addWidget(QLabel("DIM (v4.1 only):"))
		self.dim_combo = QComboBox()
		self.dim_combo.addItems(["auto", "2D", "3D"])
		row3.addWidget(self.dim_combo)
		layout.addLayout(row3)

		# Buttons
		btn_row = QHBoxLayout()
		btn_row.addStretch()
		self.export_btn = QPushButton("Export")
		self.export_btn.clicked.connect(self._do_export)
		btn_row.addWidget(self.export_btn)
		self.cancel_btn = QPushButton("Cancel")
		self.cancel_btn.clicked.connect(self.reject)
		btn_row.addWidget(self.cancel_btn)
		layout.addLayout(btn_row)

	def _browse_file(self):
		path, _ = QFileDialog.getSaveFileName(
			self, "Save Gmsh File", "", "Gmsh Files (*.msh);;All Files (*)"
		)
		if path:
			self.file_edit.setText(path)

	def _do_export(self):
		file_name = self.file_edit.text().strip()
		if not file_name:
			QMessageBox.warning(self, "Error", "Please specify an output file name.")
			return

		if not file_name.endswith(".msh"):
			file_name += ".msh"

		version = self.version_combo.currentText()
		dim = self.dim_combo.currentText()

		try:
			import cubit_mesh_export
			if version == "4.1":
				cubit_mesh_export.export_Gmesh(cubit, file_name, version=version, DIM=dim)
			else:
				cubit_mesh_export.export_Gmesh(cubit, file_name, version=version)
			QMessageBox.information(self, "Success", f"Exported: {file_name}")
			self.accept()
		except Exception as e:
			QMessageBox.critical(self, "Export Error", str(e))


# ================================================================
# Volume Calculator Dialog
# ================================================================

class VolumeCalculatorDialog(QDialog):
	"""Dialog for volume calculation with optional NGSolve integration."""

	def __init__(self, parent=None):
		super().__init__(parent)
		self.setWindowTitle("Volume Calculator")
		self.setMinimumWidth(550)
		self.setMinimumHeight(350)
		self._ext_python = _find_external_python()
		self._setup_ui()
		self._calculate_cad_volumes()

	def _setup_ui(self):
		layout = QVBoxLayout(self)

		# Volume table
		self.table = QTableWidget()
		self.table.setColumnCount(4)
		self.table.setHorizontalHeaderLabels(["ID", "Name", "CAD Volume", "NGSolve Volume"])
		self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
		self.table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)
		self.table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeToContents)
		self.table.setEditTriggers(QTableWidget.NoEditTriggers)
		self.table.setSelectionBehavior(QTableWidget.SelectRows)
		layout.addWidget(self.table)

		# NGSolve group
		ngsolve_group = QGroupBox("NGSolve Volume (external Python)")
		ngsolve_layout = QGridLayout()

		ngsolve_layout.addWidget(QLabel("Order:"), 0, 0)
		self.order_spin = QSpinBox()
		self.order_spin.setRange(1, 5)
		self.order_spin.setValue(1)
		ngsolve_layout.addWidget(self.order_spin, 0, 1)

		ngsolve_layout.addWidget(QLabel("Python:"), 1, 0)
		self.python_label = QLabel(self._ext_python or "Not found")
		self.python_label.setStyleSheet(
			"color: green;" if self._ext_python else "color: red;"
		)
		ngsolve_layout.addWidget(self.python_label, 1, 1)

		self.calc_btn = QPushButton("Calculate")
		self.calc_btn.clicked.connect(self._calculate_ngsolve)
		self.calc_btn.setEnabled(self._ext_python is not None)
		ngsolve_layout.addWidget(self.calc_btn, 0, 2)

		ngsolve_group.setLayout(ngsolve_layout)
		layout.addWidget(ngsolve_group)

		# Close button
		btn_row = QHBoxLayout()
		btn_row.addStretch()
		close_btn = QPushButton("Close")
		close_btn.clicked.connect(self.accept)
		btn_row.addWidget(close_btn)
		layout.addLayout(btn_row)

	def _calculate_cad_volumes(self):
		"""Calculate CAD volumes for selected volumes."""
		vol_ids = _get_all_volume_ids()
		self.table.setRowCount(len(vol_ids) + 1)  # +1 for total row

		total = 0.0
		self._vol_ids = vol_ids
		for row, vid in enumerate(vol_ids):
			v = cubit.volume(vid)
			vol = v.volume()
			name = cubit.get_entity_name("volume", vid) or f"Volume {vid}"
			total += vol

			self.table.setItem(row, 0, QTableWidgetItem(str(vid)))
			self.table.setItem(row, 1, QTableWidgetItem(name))
			item = QTableWidgetItem(f"{vol:.6e}")
			item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
			self.table.setItem(row, 2, item)
			self.table.setItem(row, 3, QTableWidgetItem("--"))

		# Total row
		total_row = len(vol_ids)
		font_bold = self.table.font()
		font_bold.setBold(True)

		item_label = QTableWidgetItem("Total")
		item_label.setFont(font_bold)
		self.table.setItem(total_row, 1, item_label)
		item_total = QTableWidgetItem(f"{total:.6e}")
		item_total.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
		item_total.setFont(font_bold)
		self.table.setItem(total_row, 2, item_total)
		self.table.setItem(total_row, 0, QTableWidgetItem(""))
		self.table.setItem(total_row, 3, QTableWidgetItem("--"))

		self._cad_total = total

	def _calculate_ngsolve(self):
		"""Run NGSolve volume calculation via external Python + cub5."""
		if not self._ext_python:
			QMessageBox.warning(self, "Error", "External Python with NGSolve not found.\n"
			                    "Set RADIA_PYTHON environment variable.")
			return

		vol_ids = self._vol_ids
		if not vol_ids:
			return

		order = self.order_spin.value()

		self.calc_btn.setEnabled(False)
		self.calc_btn.setText("Calculating...")

		# Save current model to temp cub5
		tmpdir = tempfile.mkdtemp(prefix="radia_vol_")
		cub5_file = os.path.join(tmpdir, "model.cub5").replace("\\", "/")
		cubit.cmd(f'save cub5 "{cub5_file}" overwrite')

		# Build command
		calc_script = os.path.join(_this_dir, "calc_volume.py")
		args = ["--cub5", cub5_file, "--order", str(order)]

		# Run async via QProcess
		self._process = QProcess(self)
		self._process.finished.connect(self._on_calc_finished)
		self._process.start(self._ext_python, [calc_script] + args)

	def _on_calc_finished(self, exit_code, exit_status):
		"""Handle async calculation result."""
		self.calc_btn.setEnabled(True)
		self.calc_btn.setText("Calculate")

		data = _parse_json_output(self._process)
		self._process = None
		if data is None:
			return

		if "error" in data:
			QMessageBox.critical(self, "NGSolve Error", data["error"])
			return

		if "warning" in data:
			QMessageBox.warning(self, "Warning", data["warning"])

		# Update table with NGSolve results
		ng_total = data["ngsolve_total"]
		ng_volumes = data.get("volumes", [])
		vol_ids = self._vol_ids

		# Update per-volume NGSolve results
		if ng_volumes and len(ng_volumes) == len(vol_ids):
			for row, vol_info in enumerate(ng_volumes):
				ng_vol = vol_info.get("ngsolve_volume")
				if ng_vol is None:
					continue
				cad_vol_text = self.table.item(row, 2).text()
				cad_vol = float(cad_vol_text)
				if cad_vol != 0:
					error_pct = (ng_vol - cad_vol) / cad_vol * 100
					text = f"{ng_vol:.6e} ({error_pct:+.2e}%)"
				else:
					text = f"{ng_vol:.6e}"
				item = QTableWidgetItem(text)
				item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
				self.table.setItem(row, 3, item)

		# Update total row
		total_row = len(vol_ids)
		if self._cad_total != 0:
			error_pct = (ng_total - self._cad_total) / self._cad_total * 100
			text = f"{ng_total:.6e} ({error_pct:+.2e}%)"
		else:
			text = f"{ng_total:.6e}"
		item = QTableWidgetItem(text)
		item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
		font_bold = self.table.font()
		font_bold.setBold(True)
		item.setFont(font_bold)
		self.table.setItem(total_row, 3, item)


# ================================================================
# Surface Area Calculator Dialog
# ================================================================

class SurfaceAreaDialog(QDialog):
	"""Dialog for surface area calculation with optional NGSolve integration."""

	def __init__(self, parent=None):
		super().__init__(parent)
		self.setWindowTitle("Surface Area Calculator")
		self.setMinimumWidth(550)
		self.setMinimumHeight(350)
		self._ext_python = _find_external_python()
		self._setup_ui()
		self._calculate_cad_areas()

	def _setup_ui(self):
		layout = QVBoxLayout(self)

		# Surface table
		self.table = QTableWidget()
		self.table.setColumnCount(4)
		self.table.setHorizontalHeaderLabels(["ID", "Name", "CAD Area", "NGSolve Area"])
		self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
		self.table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)
		self.table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeToContents)
		self.table.setEditTriggers(QTableWidget.NoEditTriggers)
		self.table.setSelectionBehavior(QTableWidget.SelectRows)
		layout.addWidget(self.table)

		# NGSolve group
		ngsolve_group = QGroupBox("NGSolve Surface Area (external Python)")
		ngsolve_layout = QHBoxLayout()
		ngsolve_layout.addWidget(QLabel("Order:"))
		self.order_spin = QSpinBox()
		self.order_spin.setRange(1, 5)
		self.order_spin.setValue(1)
		ngsolve_layout.addWidget(self.order_spin)
		self.calc_btn = QPushButton("Calculate")
		self.calc_btn.clicked.connect(self._calculate_ngsolve)
		self.calc_btn.setEnabled(self._ext_python is not None)
		ngsolve_layout.addWidget(self.calc_btn)
		ngsolve_group.setLayout(ngsolve_layout)
		layout.addWidget(ngsolve_group)

		# Close
		btn_row = QHBoxLayout()
		btn_row.addStretch()
		close_btn = QPushButton("Close")
		close_btn.clicked.connect(self.accept)
		btn_row.addWidget(close_btn)
		layout.addLayout(btn_row)

	def _calculate_cad_areas(self):
		"""Calculate CAD surface areas."""
		vol_ids = _get_all_volume_ids()
		self.table.setRowCount(len(vol_ids) + 1)
		total = 0.0
		self._vol_ids = vol_ids

		for row, vid in enumerate(vol_ids):
			surfaces = cubit.get_relatives("volume", vid, "surface")
			area = sum(cubit.surface(sid).area() for sid in surfaces)
			name = cubit.get_entity_name("volume", vid) or f"Volume {vid}"
			total += area

			self.table.setItem(row, 0, QTableWidgetItem(str(vid)))
			self.table.setItem(row, 1, QTableWidgetItem(name))
			item = QTableWidgetItem(f"{area:.6e}")
			item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
			self.table.setItem(row, 2, item)
			self.table.setItem(row, 3, QTableWidgetItem("--"))

		# Total row
		total_row = len(vol_ids)
		font_bold = self.table.font()
		font_bold.setBold(True)
		item_label = QTableWidgetItem("Total")
		item_label.setFont(font_bold)
		self.table.setItem(total_row, 1, item_label)
		item_total = QTableWidgetItem(f"{total:.6e}")
		item_total.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
		item_total.setFont(font_bold)
		self.table.setItem(total_row, 2, item_total)
		self.table.setItem(total_row, 0, QTableWidgetItem(""))
		self.table.setItem(total_row, 3, QTableWidgetItem("--"))
		self._cad_total = total

	def _calculate_ngsolve(self):
		"""Run NGSolve surface area calculation via QProcess."""
		if not self._ext_python:
			return

		self.calc_btn.setEnabled(False)
		self.calc_btn.setText("Calculating...")

		tmpdir = tempfile.mkdtemp(prefix="radia_surf_")
		cub5_file = os.path.join(tmpdir, "model.cub5").replace("\\", "/")
		cubit.cmd(f'save cub5 "{cub5_file}" overwrite')

		calc_script = os.path.join(_this_dir, "calc_surface.py")
		order = self.order_spin.value()
		args = ["--cub5", cub5_file, "--order", str(order)]

		self._process = QProcess(self)
		self._process.finished.connect(self._on_calc_finished)
		self._process.start(self._ext_python, [calc_script] + args)

	def _on_calc_finished(self, exit_code, exit_status):
		"""Handle async result."""
		self.calc_btn.setEnabled(True)
		self.calc_btn.setText("Calculate")

		data = _parse_json_output(self._process)
		self._process = None
		if data is None or "error" in data:
			if data and "error" in data:
				QMessageBox.critical(self, "Error", data["error"])
			return

		ng_total = data["ngsolve_total"]
		ng_volumes = data.get("volumes", [])
		vol_ids = self._vol_ids
		n_bnd = data.get("n_bnd_elements", 0)

		if ng_volumes and len(ng_volumes) == len(vol_ids):
			for row, vol_info in enumerate(ng_volumes):
				ng_area = vol_info.get("ngsolve_area")
				if ng_area is None:
					continue
				cad_area = float(self.table.item(row, 2).text())
				if cad_area != 0:
					error_pct = (ng_area - cad_area) / cad_area * 100
					text = f"{ng_area:.6e} ({error_pct:+.2e}%)"
				else:
					text = f"{ng_area:.6e}"
				item = QTableWidgetItem(text)
				item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
				self.table.setItem(row, 3, item)

		# Total row
		total_row = len(vol_ids)
		if self._cad_total != 0:
			error_pct = (ng_total - self._cad_total) / self._cad_total * 100
			text = f"{ng_total:.6e} ({error_pct:+.2e}%) [{n_bnd} elems]"
		else:
			text = f"{ng_total:.6e}"
		item = QTableWidgetItem(text)
		item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
		font_bold = self.table.font()
		font_bold.setBold(True)
		item.setFont(font_bold)
		self.table.setItem(total_row, 3, item)


# ================================================================
# Inductance Extractor Dialog
# ================================================================

class InductanceDialog(QDialog):
	"""Dialog for DC inductance extraction via source/sink EFIE."""

	def __init__(self, parent=None):
		super().__init__(parent)
		self.setWindowTitle("IH (BEM): Inductance + SIBC")
		self.setMinimumWidth(600)
		self.setMinimumHeight(700)
		self._ext_python = _find_external_python()
		self._result_file = os.path.join(tempfile.gettempdir(), "radia_inductance_result.json")
		self._solve_start_time = None
		self._setup_ui()
		self._populate_blocks()
		self._try_load_existing_result()

	def _setup_ui(self):
		layout = QVBoxLayout(self)

		# --- Journal editor ---
		jou_label = QLabel("Cubit Journal:")
		jou_label.setStyleSheet("font-weight: bold;")
		layout.addWidget(jou_label)

		self.jou_edit = QPlainTextEdit()
		self.jou_edit.setFont(QFont("Consolas", 9))
		self.jou_edit.setPlainText(self._default_torus_journal())
		self.jou_edit.setMaximumHeight(200)
		layout.addWidget(self.jou_edit)

		jou_btn_row = QHBoxLayout()
		load_jou_btn = QPushButton("Load .jou...")
		load_jou_btn.clicked.connect(self._load_journal)
		jou_btn_row.addWidget(load_jou_btn)
		self.run_jou_btn = QPushButton("Run Journal")
		self.run_jou_btn.clicked.connect(self._run_journal)
		jou_btn_row.addWidget(self.run_jou_btn)
		jou_btn_row.addStretch()
		layout.addLayout(jou_btn_row)

		# --- Port detection (auto from blocks) ---
		port_group = QGridLayout()
		port_group.addWidget(QLabel("Source block:"), 0, 0)
		self.source_label = QLabel("(not found)")
		self.source_label.setStyleSheet("font-weight: bold;")
		port_group.addWidget(self.source_label, 0, 1)

		port_group.addWidget(QLabel("Sink block:"), 1, 0)
		self.sink_label = QLabel("(not found)")
		self.sink_label.setStyleSheet("font-weight: bold;")
		port_group.addWidget(self.sink_label, 1, 1)

		port_group.addWidget(QLabel("Curve order:"), 2, 0)
		self.curve_spin = QSpinBox()
		self.curve_spin.setRange(1, 2)
		self.curve_spin.setValue(2)
		self.curve_spin.setToolTip("Mesh geometry order (2 = curved elements)")
		port_group.addWidget(self.curve_spin, 2, 1)

		port_group.addWidget(QLabel("FES order:"), 3, 0)
		self.fes_spin = QSpinBox()
		self.fes_spin.setRange(0, 2)
		self.fes_spin.setValue(0)
		self.fes_spin.setToolTip("HDivSurface basis order (0 = RWG)")
		port_group.addWidget(self.fes_spin, 3, 1)

		port_group.addWidget(QLabel("Workpiece block:"), 4, 0)
		self.workpiece_label = QLabel("(none)")
		self.workpiece_label.setStyleSheet("font-weight: bold; color: gray;")
		port_group.addWidget(self.workpiece_label, 4, 1)

		port_group.addWidget(QLabel("Air block (Post B):"), 5, 0)
		self.air_label = QLabel("(none)")
		self.air_label.setStyleSheet("font-weight: bold; color: gray;")
		port_group.addWidget(self.air_label, 5, 1)


		layout.addLayout(port_group)

		# --- ESIM / Dowell settings (visible only when workpiece found) ---
		self.esim_group = QGroupBox("Workpiece Surface Impedance")
		esim_layout = QGridLayout(self.esim_group)

		esim_layout.addWidget(QLabel("Impedance:"), 0, 0)
		self.model_combo = QComboBox()
		self.model_combo.addItems(["SIBC", "ESIM", "Dowell"])
		self.model_combo.setToolTip("SIBC: classical Z_s=(1+j)/sigma*delta (default)\n"
		                            "ESIM: 1D cell problem (nonlinear OK)\n"
		                            "Dowell: analytical slab (linear only)")
		esim_layout.addWidget(self.model_combo, 0, 1)

		esim_layout.addWidget(QLabel("Frequency [Hz]:"), 1, 0)
		self.freq_edit = QLineEdit("50000")
		esim_layout.addWidget(self.freq_edit, 1, 1)

		esim_layout.addWidget(QLabel("Sigma [S/m]:"), 2, 0)
		self.sigma_edit = QLineEdit("2.0e6")
		esim_layout.addWidget(self.sigma_edit, 2, 1)

		# mu_r row (SIBC/Dowell only)
		self.mur_label = QLabel("mu_r:")
		esim_layout.addWidget(self.mur_label, 3, 0)
		self.mur_edit = QLineEdit("100")
		esim_layout.addWidget(self.mur_edit, 3, 1)

		# BH curve row (ESIM only)
		self.bh_label = QLabel("BH curve:")
		esim_layout.addWidget(self.bh_label, 4, 0)
		bh_row = QHBoxLayout()
		self.bh_edit = QLineEdit("(built-in Steel)")
		self.bh_edit.setReadOnly(True)
		self.bh_edit.setToolTip("BH curve file (2-column: H[A/m] B[T])")
		bh_row.addWidget(self.bh_edit)
		self.bh_browse = QPushButton("...")
		self.bh_browse.setFixedWidth(30)
		self.bh_browse.clicked.connect(self._browse_bh)
		bh_row.addWidget(self.bh_browse)
		esim_layout.addLayout(bh_row, 4, 1)

		# Curvature row (ESIM only)
		self.curv_label = QLabel("Curvature:")
		esim_layout.addWidget(self.curv_label, 5, 0)
		self.curvature_combo = QComboBox()
		self.curvature_combo.addItems(["Local curvature", "None (flat)"])
		self.curvature_combo.setToolTip(
			"Local curvature: cylindrical cell problem (Bessel I0/I1),\n"
			"None (flat): planar slab cell problem (cosh/sinh).")
		esim_layout.addWidget(self.curvature_combo, 5, 1)

		# Connect impedance model change -> show/hide mu_r vs BH
		self.model_combo.currentTextChanged.connect(self._on_impedance_changed)
		self._on_impedance_changed(self.model_combo.currentText())

		self.esim_group.setVisible(False)
		layout.addWidget(self.esim_group)

		# --- Result table ---
		self.result_table = QTableWidget()
		self.result_table.setColumnCount(2)
		self.result_table.setHorizontalHeaderLabels(["Parameter", "Value"])
		self.result_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
		self.result_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
		self.result_table.setEditTriggers(QTableWidget.NoEditTriggers)
		self.result_table.setRowCount(1)
		self.result_table.setItem(0, 0, QTableWidgetItem("Status"))
		self.result_table.setItem(0, 1, QTableWidgetItem("(not yet computed)"))
		layout.addWidget(self.result_table)

		# Python info
		py_label = QLabel(f"Python: {self._ext_python or 'Not found'}")
		py_label.setStyleSheet("color: green;" if self._ext_python else "color: red;")
		layout.addWidget(py_label)

		# Debug output
		self.debug_text = QLabel("")
		self.debug_text.setWordWrap(True)
		self.debug_text.setStyleSheet("color: gray; font-size: 10px;")
		layout.addWidget(self.debug_text)

		# --- Action buttons ---
		btn_row = QHBoxLayout()
		btn_row.addStretch()
		self.solve_btn = QPushButton("Solve")
		self.solve_btn.clicked.connect(self._extract)
		self.solve_btn.setEnabled(False)
		self.solve_btn.setToolTip("BEM: L extraction + workpiece heating + B-field post")
		btn_row.addWidget(self.solve_btn)
		self.open_gmsh_btn = QPushButton("Open Result")
		self.open_gmsh_btn.clicked.connect(self._open_gmsh_result)
		self.open_gmsh_btn.setEnabled(False)
		btn_row.addWidget(self.open_gmsh_btn)
		close_btn = QPushButton("Close")
		close_btn.clicked.connect(self.accept)
		btn_row.addWidget(close_btn)
		layout.addLayout(btn_row)

	def _default_torus_journal(self):
		"""Default BEM journal: coil surface mesh + workpiece (no air volume)."""
		lines = [
			"# IH (BEM): Torus coil + workpiece + air (B-field post)",
			"# Coil: R=30mm, a=3mm, 355deg (surface mesh for BEM EFIE)",
			"# Workpiece: R=10mm, H=20mm cylinder (SIBC heating)",
			"# Air: R=60mm sphere (B-field distribution output)",
			"reset",
			"",
			"# --- Coil: revolve circle to create gapped torus ---",
			"create surface circle radius 0.003 yplane",
			"move surface 1 x 0.03 include_merged",
			"sweep surface 1 zaxis angle 355",
			"",
			"# Hex sweep mesh (coil surface)",
			"volume 1 scheme sweep source surface 1 target surface 3",
			"surface 1 size 0.003",
			"surface 1 scheme pave",
			"mesh surface 1",
			"mesh volume 1",
			"",
			"# --- Workpiece: cylinder at center ---",
			"create cylinder height 0.020 radius 0.010",
			"volume 2 scheme tetmesh",
			"volume 2 size 0.003",
			"mesh volume 2",
			"",
			"# --- Air: sphere for B-field post-processing ---",
			"create sphere radius 0.060",
			"volume 3 scheme tetmesh",
			"volume 3 size 0.010",
			"mesh volume 3",
			"",
			"# --- Blocks ---",
			"set duplicate block elements on",
			"block 1 add volume 1",
			'block 1 name "conductor"',
			"block 2 add face in surface 1",
			'block 2 name "source"',
			"block 3 add face in surface 3",
			'block 3 name "sink"',
			"block 4 add face in surface all in volume 1",
			'block 4 name "coil"',
			"block 5 add volume 2",
			'block 5 name "workpiece"',
			"block 6 add volume 3",
			'block 6 name "air"',
			"",
			"# Hide air volume (post-processing only)",
			"volume 3 visibility off",
		]
		return "\n".join(lines) + "\n"

	@staticmethod
	def _default_fem_journal():
		"""Default FEM journal: Periodic Kelvin (2-sphere + mesh copy).

		Webcut is required: ACIS spheres have no curves, but
		`copy mesh surface` requires curve/vertex mapping.

		User must select source/sink gap surfaces manually after
		boolean operations change surface IDs.
		"""
		lines = [
			"# IH (FEM): Source/sink + Periodic Kelvin (2-sphere)",
			"# Coil: R=30mm, a=3mm gapped torus (355deg, source/sink on gap)",
			"# Workpiece: R=10mm, H=20mm cylinder (SIBC on interface)",
			"# Kelvin: R=60mm interior + exterior sphere (periodic BC)",
			"reset",
			"",
			"# === Coil: gapped torus ===",
			"create surface circle radius 0.003 yplane",
			"move surface 1 x 0.030 include_merged",
			"sweep surface 1 zaxis angle 355",
			"",
			"# === Workpiece cylinder ===",
			"create cylinder height 0.020 radius 0.010",
			"",
			"# === Interior sphere + webcut for mesh copy topology ===",
			"create sphere radius 0.060",
			"webcut volume 3 with plane zplane",
			"# Webcut creates curves at equator for copy mesh mapping",
			"",
			"# === Copy to create exterior sphere (nomesh) ===",
			"volume 3 4 copy move x 0.300 nomesh",
			"",
			"# === Boolean: air = interior sphere halves - coil - wp ===",
			"volume 1 copy",
			"volume 2 copy",
			"subtract volume 7 8 from volume 3 4",
			"",
			"# === Imprint and merge ===",
			"imprint volume all",
			"merge volume all",
			"",
			"# === Mesh interior domain ===",
			"volume all scheme tetmesh",
			"volume 1 size 0.003",
			"volume 2 size 0.003",
			"volume 3 4 size 0.008",
			"mesh volume 1 2 3 4",
			"# Kelvin: mesh copy + tet done by dialog (_setup_kelvin)",
			"",
			"# === Blocks (volume) ===",
			"set duplicate block elements on",
			"block 1 add volume 1",
			'block 1 name "coil"',
			"block 2 add volume 2",
			'block 2 name "workpiece"',
			"block 3 add volume 3 4",
			'block 3 name "air"',
			"block 4 add volume 5 6",
			'block 4 name "kelvin"',
			"",
			"# === Blocks (surface) ===",
			"# NOTE: Select gap faces in Cubit GUI, surface IDs vary",
			"# after boolean. Use: pick surface on coil gap -> add to block.",
			"# block 5 add face in surface <SOURCE_SID>",
			'# block 5 name "source"',
			"# block 6 add face in surface <SINK_SID>",
			'# block 6 name "sink"',
			"# (Without source/sink, J_theta fallback is used)",
			"",
			"# wp_surface, kelvin_int/ext, outer: auto-created by dialog",
			"",
			"# === Hide non-essential ===",
			"volume 3 4 visibility off",
			"volume 5 6 visibility off",
		]
		return "\n".join(lines) + "\n"

	def _load_journal(self):
		"""Load a .jou file into the editor."""
		path, _ = QFileDialog.getOpenFileName(
			self, "Load Cubit Journal", "",
			"Cubit Journal (*.jou);;All Files (*)")
		if path:
			try:
				with open(path, "r", encoding="utf-8") as f:
					self.jou_edit.setPlainText(f.read())
				self.debug_text.setText(f"Loaded: {path}")
			except Exception as e:
				QMessageBox.warning(self, "Error", f"Failed to load: {e}")

	def _run_journal(self):
		"""Execute the journal text in Cubit, then auto-create source/sink blocks."""
		text = self.jou_edit.toPlainText().strip()
		if not text:
			return

		self.run_jou_btn.setEnabled(False)
		self.run_jou_btn.setText("Running...")
		try:
			for line in text.splitlines():
				line = line.strip()
				if not line or line.startswith("#"):
					continue
				cubit.cmd(line)
			self.debug_text.setText("Journal executed.")
		except Exception as e:
			QMessageBox.warning(self, "Journal Error", str(e))
			self.debug_text.setText(f"Journal error: {e}")
			self.run_jou_btn.setEnabled(True)
			self.run_jou_btn.setText("Run Journal")
			return
		finally:
			self.run_jou_btn.setEnabled(True)
			self.run_jou_btn.setText("Run Journal")

		# Auto-create source/sink blocks from planar gap faces
		self._auto_create_source_sink_blocks()

		# Re-detect blocks
		self._populate_blocks()

	def _auto_create_source_sink_blocks(self):
		"""Find gap faces by y-coordinate and create source/sink blocks."""
		try:
			# Find source/sink surfaces by center y-coordinate
			src_sid, snk_sid = None, None
			for sid in cubit.get_entities("surface"):
				cx, cy, cz = cubit.surface(sid).center_point()
				if cy > 0 and src_sid is None:
					src_sid = sid
				elif cy < 0 and snk_sid is None:
					snk_sid = sid

			if src_sid is None or snk_sid is None:
				self.debug_text.setText(
					"Auto-detect: could not find faces with y>0 and y<0.")
				return

			# Next available block ID
			existing = list(cubit.get_block_id_list())
			next_id = max(existing) + 1 if existing else 1

			# Skip if already defined
			existing_names = set()
			for bid in existing:
				try:
					existing_names.add(cubit.get_exodus_entity_name("block", bid))
				except Exception:
					pass

			has_tri = len(cubit.get_entities("tri")) > 0
			has_quad = len(cubit.get_entities("quad")) > 0
			cubit.cmd("set duplicate block elements on")

			if "source" not in existing_names:
				if has_tri:
					cubit.cmd(f"block {next_id} add tri in surface {src_sid}")
				if has_quad:
					cubit.cmd(f"block {next_id} add face in surface {src_sid}")
				cubit.cmd(f'block {next_id} name "source"')
				next_id += 1

			if "sink" not in existing_names:
				if has_tri:
					cubit.cmd(f"block {next_id} add tri in surface {snk_sid}")
				if has_quad:
					cubit.cmd(f"block {next_id} add face in surface {snk_sid}")
				cubit.cmd(f'block {next_id} name "sink"')
				next_id += 1

			if "boundary" not in existing_names:
				if has_tri:
					cubit.cmd(f"block {next_id} add tri all")
				if has_quad:
					cubit.cmd(f"block {next_id} add face all")
				cubit.cmd(f'block {next_id} name "boundary"')

			self.debug_text.setText(
				f"source=surface {src_sid} (y>0), sink=surface {snk_sid} (y<0)")

		except Exception as e:
			self.debug_text.setText(f"Auto-detect failed: {e}")

	def _set_result(self, param, value):
		"""Set a single-row result in the table."""
		self.result_table.setRowCount(1)
		self.result_table.setItem(0, 0, QTableWidgetItem(param))
		self.result_table.setItem(0, 1, QTableWidgetItem(value))


	def _try_load_existing_result(self):
		"""Load existing result JSON if available (from previous Solve)."""
		if os.path.exists(self._result_file):
			try:
				with open(self._result_file, "r") as f:
					data = json.load(f)
				if "inductance_H" in data:
					self._display_result(data)
					self._enable_gmsh_buttons(data)
					# Restore solve results for Post button
					j_npy = data.get("j_npy", "")
					mesh_vol = data.get("mesh_vol", "")
					if j_npy and os.path.exists(j_npy):
						self._j_npy = j_npy
						self.post_btn.setEnabled(True)
					if mesh_vol and os.path.exists(mesh_vol):
						self._mesh_vol = mesh_vol
					# Restore default post volume
					dl = data.get("default_lxyz")
					dm = data.get("default_maxh")
					if dl:
						self.lxyz_edit.setText(f"{dl}, {dl}, {dl}")
					if dm:
						self.maxh_vol_edit.setText(str(dm))
					self.debug_text.setText(f"Previous result loaded: {self._result_file}")
			except Exception:
				pass

	def _on_impedance_changed(self, text):
		"""Show mu_r for SIBC/Dowell, BH curve for ESIM."""
		is_esim = (text == "ESIM")
		# mu_r: SIBC/Dowell only (linear)
		self.mur_label.setVisible(not is_esim)
		self.mur_edit.setVisible(not is_esim)
		# BH curve + curvature: ESIM only (nonlinear)
		self.bh_label.setVisible(is_esim)
		self.bh_edit.setVisible(is_esim)
		self.bh_browse.setVisible(is_esim)
		self.curv_label.setVisible(is_esim)
		self.curvature_combo.setVisible(is_esim)

	def _browse_bh(self):
		"""Browse for BH curve file (2-column: H[A/m] B[T])."""
		path, _ = QFileDialog.getOpenFileName(
			self, "Load BH Curve", "",
			"Text files (*.txt *.csv *.dat);;All Files (*)")
		if path:
			self.bh_edit.setText(path)

	def _populate_blocks(self):
		"""Auto-detect source/sink/workpiece from Cubit block names."""
		self._source_block = None
		self._sink_block = None
		self._workpiece_block = None
		self._air_block = None
		try:
			for bid in cubit.get_block_id_list():
				try:
					name = cubit.get_exodus_entity_name("block", bid)
					if name == "source":
						self._source_block = name
					elif name == "sink":
						self._sink_block = name
					elif name == "workpiece":
						self._workpiece_block = name
					elif name == "air":
						self._air_block = name
				except Exception:
					pass
		except Exception:
			pass

		if self._source_block:
			self.source_label.setText(self._source_block)
			self.source_label.setStyleSheet("font-weight: bold; color: green;")
		else:
			self.source_label.setText("(not found)")
			self.source_label.setStyleSheet("font-weight: bold; color: red;")

		if self._sink_block:
			self.sink_label.setText(self._sink_block)
			self.sink_label.setStyleSheet("font-weight: bold; color: green;")
		else:
			self.sink_label.setText("(not found)")
			self.sink_label.setStyleSheet("font-weight: bold; color: red;")

		# Air block detection (for B-field post-processing)
		if self._air_block:
			self.air_label.setText(self._air_block)
			self.air_label.setStyleSheet("font-weight: bold; color: green;")
		else:
			self.air_label.setText("(none - no B post)")
			self.air_label.setStyleSheet("font-weight: bold; color: gray;")

		# Workpiece block detection -> show/hide ESIM group
		if self._workpiece_block:
			self.workpiece_label.setText(self._workpiece_block)
			self.workpiece_label.setStyleSheet("font-weight: bold; color: green;")
			self.esim_group.setVisible(True)
		else:
			self.workpiece_label.setText("(none)")
			self.workpiece_label.setStyleSheet("font-weight: bold; color: gray;")
			self.esim_group.setVisible(False)

		# Enable Solve when source+sink found
		can_solve = (self._source_block is not None
		             and self._sink_block is not None
		             and self._ext_python is not None)
		self.solve_btn.setEnabled(can_solve)
		if not can_solve and self._ext_python:
			self.debug_text.setText(
				'Define blocks named "source" and "sink" in your journal.')

	def _extract(self):
		"""Run inductance extraction via external Python."""
		if not self._ext_python:
			QMessageBox.warning(self, "Error", "External Python not found.")
			return

		source = self._source_block
		sink = self._sink_block

		if not source or not sink:
			QMessageBox.warning(self, "Error",
				'Blocks named "source" and "sink" not found.\n'
				"Define them in your Cubit journal:\n"
				'  block N name "source"\n'
				'  block M name "sink"')
			return

		curve_order = self.curve_spin.value()
		fes_order = self.fes_spin.value()

		self.solve_btn.setEnabled(False)
		self.solve_btn.setText("Solving...")
		self._set_result("Status", "Computing...")

		from datetime import datetime
		self._solve_start_time = datetime.now()

		# Save cub5
		tmpdir = tempfile.mkdtemp(prefix="radia_ind_")
		cub5_file = os.path.join(tmpdir, "model.cub5").replace("\\", "/")
		cubit.cmd(f'save cub5 "{cub5_file}" overwrite')

		# Build command
		calc_script = os.path.join(_this_dir, "calc_inductance.py")

		# Output dir = examples/cubit_panels/inductance/results/
		_repo_root = os.path.dirname(os.path.dirname(os.path.dirname(
			os.path.dirname(os.path.abspath(calc_script)))))
		output_dir = os.path.join(_repo_root, "examples", "cubit_panels",
		                          "inductance", "results")
		os.makedirs(output_dir, exist_ok=True)
		msh_output = os.path.join(output_dir, "inductance_J.msh").replace("\\", "/")

		# JSON output file (reliable, avoids stdout parsing issues)
		self._json_output = os.path.join(tmpdir, "result.json").replace("\\", "/")
		args = [
			calc_script,
			"--cub5", cub5_file,
			"--msh-output", msh_output,
			"--source", source,
			"--sink", sink,
			"--order", str(curve_order),
			"--fes-order", str(fes_order),
			"--output", self._json_output,
		]

		# Workpiece parameters (BEM + impedance model)
		if self._workpiece_block:
			_curv = self.curvature_combo.currentText()
			_curv_arg = "local_curvature" if "Local" in _curv else "none"
			_imp = self.model_combo.currentText()  # "SIBC", "ESIM", "Dowell"
			_model_arg = "bem-sibc" if _imp == "SIBC" else _imp.lower()
			# Auto-detect half-thickness from workpiece bounding box
			half_t = 0.005  # fallback
			try:
				for bid in cubit.get_block_id_list():
					name = cubit.get_exodus_entity_name("block", bid)
					if name == self._workpiece_block:
						for vid in cubit.get_block_volumes(bid):
							bb = cubit.get_bounding_box("volume", vid)
							r_wp = max(abs(bb[1]), abs(bb[0]),
							           abs(bb[3]), abs(bb[2]))
							half_t = r_wp
							break
						break
			except Exception:
				pass

			args += [
				"--workpiece", self._workpiece_block,
				"--impedance-model", _model_arg,
				"--frequency", self.freq_edit.text().strip(),
				"--sigma", self.sigma_edit.text().strip(),
				"--half-thickness", str(round(half_t, 6)),
				"--material", "steel",
				"--mu-r", self.mur_edit.text().strip(),
				"--esim-geometry", _curv_arg,
			]
			# BH curve file (if user-specified)
			bh_text = self.bh_edit.text().strip()
			if bh_text and not bh_text.startswith("(") and os.path.isfile(bh_text):
				args += ["--bh-file", bh_text]

		# Run async via QProcess (non-blocking)
		self._process = QProcess(self)
		self._process.readyReadStandardError.connect(self._on_stderr)
		self._process.finished.connect(self._on_extract_finished)
		self._gmsh_launched = False
		self._process.start(self._ext_python, args)

	def _solve_heating(self):
		"""Stage 2: FEM-ESIM workpiece heating (independent from BEM)."""
		if not self._ext_python:
			QMessageBox.warning(self, "Error", "External Python not found.")
			return
		if not self._workpiece_block:
			QMessageBox.warning(self, "Error", 'Block named "workpiece" not found.')
			return

		self.solve_p_btn.setEnabled(False)
		self.solve_p_btn.setText("Solving P...")
		self._set_result("Status", "Computing heating (FEM-ESIM)...")

		# Get workpiece geometry from Cubit block bounding box
		try:
			wp_vids = []
			for bid in cubit.get_block_id_list():
				name = cubit.get_exodus_entity_name("block", bid)
				if name == self._workpiece_block:
					for vid in cubit.get_block_volumes(bid):
						wp_vids.append(vid)
					break
			if not wp_vids:
				# Fallback: use all volumes not in conductor/source/sink blocks
				wp_vids = list(cubit.get_entities("volume"))

			# Bounding box of workpiece
			bb = cubit.get_bounding_box("volume", wp_vids[0])
			# bb = (xmin, xmax, ymin, ymax, zmin, zmax, ...)
			r_wp = max(abs(bb[1]), abs(bb[0]), abs(bb[3]), abs(bb[2]))
			h_wp = abs(bb[5] - bb[4])
		except Exception:
			r_wp = 0.01
			h_wp = 0.02

		# Get coil geometry (from conductor block bounding box)
		try:
			for bid in cubit.get_block_id_list():
				name = cubit.get_exodus_entity_name("block", bid)
				if name == "conductor":
					cond_vids = list(cubit.get_block_volumes(bid))
					if cond_vids:
						bb_c = cubit.get_bounding_box("volume", cond_vids[0])
						r_coil = (abs(bb_c[1]) + abs(bb_c[0])) / 2
						a_coil = min(abs(bb_c[1] - bb_c[0]),
						             abs(bb_c[5] - bb_c[4])) / 4
					break
		except Exception:
			r_coil = 0.03
			a_coil = 0.003

		freq = float(self.freq_edit.text().strip() or "50000")
		sigma_str = self.sigma_edit.text().strip() or "2e6"
		sigma_val = float(sigma_str)
		material = "steel"  # material selection removed; BH/mu_r set directly

		# Build command
		calc_script = os.path.join(_this_dir, "calc_heating.py")
		tmpdir = tempfile.mkdtemp(prefix="radia_heat_")
		self._heat_json = os.path.join(tmpdir, "heat_result.json").replace("\\", "/")

		_curv = self.curvature_combo.currentText()
		_curv_arg = "local_curvature" if "Local" in _curv else "none"
		args = [
			calc_script,
			"--r-coil", str(round(r_coil, 6)),
			"--a-coil", str(round(a_coil, 6)),
			"--r-wp", str(round(r_wp, 6)),
			"--h-wp", str(round(h_wp, 6)),
			"--frequency", str(freq),
			"--sigma", str(sigma_val),
			"--material", material,
			"--esim-geometry", _curv_arg,
			"--output", self._heat_json,
		]

		self._process = QProcess(self)
		self._process.readyReadStandardError.connect(self._on_stderr)
		self._process.finished.connect(self._on_heating_finished)
		self._process.start(self._ext_python, args)

	def _on_heating_finished(self, exit_code, exit_status):
		"""Handle FEM-ESIM heating result."""
		self.solve_p_btn.setEnabled(True)
		self.solve_p_btn.setText("Solve P")
		self._process = None

		data = None
		json_path = getattr(self, '_heat_json', '')
		if json_path and os.path.exists(json_path):
			try:
				with open(json_path, "r") as f:
					data = json.load(f)
			except Exception:
				pass

		if data is None:
			self._set_result("Status", f"Heating error (exit code {exit_code})")
			return
		if "error" in data:
			QMessageBox.critical(self, "Error", data["error"])
			return

		# Display heating results
		P = data.get("P_total", 0)
		Q = data.get("Q_total", 0)
		L = data.get("L_coil", 0)
		delta = data.get("delta", 0)
		freq = data.get("frequency", 0)
		material = data.get("material", "")

		rows = [
			("--- FEM-ESIM Heating ---", f"{material}, {freq:.0f} Hz"),
			("P (workpiece)", f"{P:.4e} W"),
			("Q (workpiece)", f"{Q:.4e} var"),
			("L (coil, FEM)", f"{L*1e9:.2f} nH"),
			("Skin depth", f"{delta*1e3:.3f} mm"),
			("Panels", str(data.get("n_panels", ""))),
			("DOFs (FEM)", str(data.get("ndof", ""))),
			("Time: mesh", f"{data.get('t_mesh', 0):.1f} s"),
			("Time: solve", f"{data.get('t_solve', 0):.1f} s"),
			("Time: ESIM", f"{data.get('t_esim', 0):.1f} s"),
		]

		# P at different currents
		for I in [10, 100, 1000]:
			rows.append((f"P at {I} A", f"{P * I**2:.2f} W"))

		self.result_table.setRowCount(len(rows))
		for i, (param, val) in enumerate(rows):
			self.result_table.setItem(i, 0, QTableWidgetItem(param))
			item = QTableWidgetItem(val)
			item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
			self.result_table.setItem(i, 1, item)

		self.debug_text.setText(
			f"Heating done: P={P:.4e} W (per 1A), "
			f"mesh={data.get('t_mesh',0):.1f}s, solve={data.get('t_solve',0):.1f}s")

	def _post_process(self):
		"""Run post-processing (B-field + GMSH export) via external Python.

		If no Solve result exists, runs Solve first then chains to Post.
		"""
		j_npy = getattr(self, '_j_npy', '')
		if not j_npy or not os.path.exists(j_npy):
			# No solve result -> run Solve first, then chain to Post
			self._chain_post_after_solve = True
			self._extract()
			return

		self._run_post(j_npy)

	def _run_post(self, j_npy):
		"""Launch post-processing subprocess (no Cubit needed)."""
		# Find mesh_vol alongside j_npy
		base_dir = os.path.dirname(j_npy)
		mesh_vol = getattr(self, '_mesh_vol', '')
		if not mesh_vol or not os.path.exists(mesh_vol):
			mesh_vol = os.path.join(base_dir, "surface_mesh.vol").replace("\\", "/")
		if not os.path.exists(mesh_vol):
			QMessageBox.warning(self, "Error",
				f"surface_mesh.vol not found in {base_dir}.\nRun Solve first.")
			return

		self.post_btn.setEnabled(False)
		self.post_btn.setText("Post...")
		self.debug_text.setText("Post-processing (B-field + GMSH export)...")

		fes_order = self.fes_spin.value()

		# Read post volume settings from UI
		try:
			parts = [float(v.strip()) for v in self.lxyz_edit.text().split(",")]
			lx = parts[0] if len(parts) > 0 else 0.07
			ly = parts[1] if len(parts) > 1 else lx
			lz = parts[2] if len(parts) > 2 else lx
		except (ValueError, IndexError):
			lx = ly = lz = 0.07
		try:
			maxh_vol = float(self.maxh_vol_edit.text())
		except ValueError:
			maxh_vol = 0.01

		_repo_root = os.path.dirname(os.path.dirname(os.path.dirname(
			os.path.dirname(os.path.abspath(
				os.path.join(_this_dir, "calc_inductance.py"))))))
		output_dir = os.path.join(_repo_root, "examples", "cubit_panels",
		                          "inductance", "results")
		os.makedirs(output_dir, exist_ok=True)
		msh_output = os.path.join(output_dir, "inductance_J.msh").replace("\\", "/")

		calc_script = os.path.join(_this_dir, "calc_inductance.py")
		self._post_json = os.path.join(base_dir, "post_result.json").replace("\\", "/")
		args = [
			calc_script,
			"--mode", "post",
			"--mesh-vol", mesh_vol,
			"--fes-order", str(fes_order),
			"--msh-output", msh_output,
			"--j-npy", j_npy,
			"--lx", str(lx),
			"--ly", str(ly),
			"--lz", str(lz),
			"--maxh-vol", str(maxh_vol),
			"--output", self._post_json,
		]

		self._process = QProcess(self)
		self._process.readyReadStandardError.connect(self._on_stderr)
		self._process.finished.connect(self._on_post_finished)
		self._process.start(self._ext_python, args)

	def _on_post_finished(self, exit_code, exit_status):
		"""Handle post-processing result."""
		self.post_btn.setEnabled(True)
		self.post_btn.setText("Post")
		self._process = None

		data = None
		json_path = getattr(self, '_post_json', '')
		if json_path and os.path.exists(json_path):
			try:
				with open(json_path, "r") as f:
					data = json.load(f)
			except Exception:
				pass

		if data is None:
			self.debug_text.setText(f"Post error (exit code {exit_code})")
			return

		if "error" in data:
			QMessageBox.critical(self, "Error", data["error"])
			return

		self._enable_gmsh_buttons(data)
		t_post = data.get("t_post", "?")

		# Merge t_post into saved result JSON and re-display
		if os.path.exists(self._result_file):
			try:
				with open(self._result_file, "r") as f:
					saved = json.load(f)
				saved["t_post"] = data.get("t_post")
				saved["gmsh_file"] = data.get("gmsh_file", "")
				with open(self._result_file, "w") as f:
					json.dump(saved, f)
				self._display_result(saved)
			except Exception:
				pass

		self.debug_text.setText(f"Post done ({t_post}s). Click Open Result.")

	def _on_stderr(self):
		"""Handle stderr from subprocess - progress updates."""
		if self._process is None:
			return
		data = self._process.readAllStandardError().data().decode("utf-8", errors="replace")
		for line in data.splitlines():
			if line.startswith("MESH_READY:"):
				self.debug_text.setText("Mesh exported, assembling BEM matrix...")
			elif line.startswith("FIELD_READY:"):
				self.debug_text.setText("J-distribution written, computing B field...")
			elif line.startswith("B_FIELD_READY:"):
				self.debug_text.setText("B-distribution written, finalizing...")
			elif line.startswith("B_FIELD_ERROR:"):
				self.debug_text.setText("B field: " + line.split(":", 1)[1])
			elif line.startswith("BEM-SIBC:") or line.startswith("FEM-ESIM:"):
				self.debug_text.setText(line.strip())
			elif line.startswith("ESIM_START:"):
				self.debug_text.setText("Computing workpiece impedance (ESIM)...")
			elif line.startswith("ESIM_DONE:"):
				self.debug_text.setText("Workpiece impedance done. " + line.split(":", 1)[1])
			elif line.startswith("HEATING_MESH:"):
				self.debug_text.setText("Building FEM-ESIM mesh (auto)...")
			elif line.startswith("HEATING_SOLVE:"):
				self.debug_text.setText("FEM static solve...")
			elif line.startswith("HEATING_ESIM:"):
				self.debug_text.setText("Computing ESIM surface impedance...")

	def _on_extract_finished(self, exit_code, exit_status):
		"""Handle async extraction result."""
		self.solve_btn.setEnabled(True)
		self.solve_btn.setText("Solve")
		self._process = None

		# Read result from JSON file (more reliable than stdout parsing)
		data = None
		json_path = getattr(self, '_json_output', '')
		if json_path and os.path.exists(json_path):
			try:
				with open(json_path, "r") as f:
					data = json.load(f)
			except Exception:
				pass

		if data is None:
			self._set_result("Status", f"Error (exit code {exit_code})")
			return

		if "error" in data:
			QMessageBox.critical(self, "Error", data["error"])
			self._set_result("Status", "Error")
			return

		self._display_result(data)

		# Record to Optuna study
		self._record_to_optuna(data)

		# Save result
		with open(self._result_file, "w") as f:
			json.dump(data, f)

		# Enable Post button (J_coeffs.npy + mesh_vol available)
		j_npy = data.get("j_npy", "")
		mesh_vol = data.get("mesh_vol", "")
		if j_npy and os.path.exists(j_npy):
			self._j_npy = j_npy
			self.post_btn.setEnabled(True)
		if mesh_vol and os.path.exists(mesh_vol):
			self._mesh_vol = mesh_vol

		# Update default post volume from conductor bbox
		dl = data.get("default_lxyz")
		dm = data.get("default_maxh")
		if dl:
			self.lxyz_edit.setText(f"{dl}, {dl}, {dl}")
		if dm:
			self.maxh_vol_edit.setText(str(dm))

		self.open_gmsh_btn.setEnabled(False)  # Post not yet done
		self.debug_text.setText(f"Result saved: {self._result_file}")

		# Chain to Post if requested
		if getattr(self, '_chain_post_after_solve', False):
			self._chain_post_after_solve = False
			j_npy = getattr(self, '_j_npy', '')
			if j_npy and os.path.exists(j_npy):
				self._run_post(j_npy)

	def _enable_gmsh_buttons(self, data):
		"""Enable Open GMSH button if result file exists."""
		gmsh_file = data.get("gmsh_file", "")
		if gmsh_file and os.path.exists(gmsh_file):
			self._gmsh_file = gmsh_file
			self.open_gmsh_btn.setEnabled(True)

	def _open_gmsh_result(self):
		"""Open combined result (.geo with all views) in GMSH."""
		gmsh_file = getattr(self, '_gmsh_file', '')
		if not gmsh_file or not os.path.exists(gmsh_file):
			QMessageBox.warning(self, "Error", "No GMSH result file found. Run Solve first.")
			return
		try:
			import subprocess as _sp
			ext_py = self._ext_python
			if ext_py:
				# Use pythonw.exe to avoid console window
				ext_pyw = ext_py.replace("python.exe", "pythonw.exe")
				if not os.path.exists(ext_pyw):
					ext_pyw = ext_py
				code = (
					"import sys, gmsh; "
					"gmsh.initialize(sys.argv, run=True); "
					"gmsh.finalize()"
				)
				_sp.Popen([ext_pyw, "-c", code, gmsh_file])
				self.debug_text.setText(f"GMSH: {os.path.basename(gmsh_file)}")
			else:
				os.startfile(gmsh_file)
				self.debug_text.setText(f"Opened: {os.path.basename(gmsh_file)}")
		except Exception as e:
			QMessageBox.warning(self, "Error",
				f"Failed to open GMSH: {e}\nFile: {gmsh_file}")

	def _get_cub5_path(self):
		"""Get real cub5 path from Cubit (or fallback to cwd)."""
		try:
			path = cubit.get_file_name()
			if path:
				return path
		except Exception:
			pass
		return os.path.join(os.getcwd(), "untitled.cub5")

	def _record_to_optuna(self, data):
		"""Record result as Optuna trial. No-op if optuna not installed."""
		try:
			from optuna_study_helper import InductanceStudy
		except ImportError:
			return

		cub5_path = self._get_cub5_path()
		study = InductanceStudy(cub5_path)
		if not study.available:
			return

		params = {
			"curve_order": self.curve_spin.value(),
			"fes_order": self.fes_spin.value(),
			"source_block": self._source_block or "",
			"sink_block": self._sink_block or "",
		}

		trial_num = study.record_trial(params, data, self._solve_start_time)
		n = study.n_trials
		self.debug_text.setText(
			f"Trial #{trial_num} recorded ({n} total) -> {study.db_path}"
		)

	def _launch_dashboard(self):
		"""Launch optuna-dashboard for the current model."""
		try:
			from optuna_study_helper import InductanceStudy
		except ImportError:
			QMessageBox.warning(
				self, "Dashboard",
				"optuna not installed.\npip install optuna optuna-dashboard"
			)
			return

		cub5_path = self._get_cub5_path()
		study = InductanceStudy(cub5_path)
		if not study.available or study.n_trials == 0:
			QMessageBox.information(
				self, "Dashboard",
				"No study data yet. Run Solve first."
			)
			return

		import webbrowser
		url = study.launch_dashboard()
		if url:
			webbrowser.open(url)
			self.debug_text.setText(f"Dashboard: {url}")
		else:
			QMessageBox.warning(
				self, "Dashboard",
				"Failed to launch optuna-dashboard."
			)

	def _display_result(self, data):
		"""Display extraction result in the table."""
		L = data.get("inductance_H", 0.0)

		if abs(L) >= 1e-3:
			L_str = f"{L*1e3:.4f} mH"
		elif abs(L) >= 1e-6:
			L_str = f"{L*1e6:.4f} uH"
		elif abs(L) >= 1e-9:
			L_str = f"{L*1e9:.4f} nH"
		else:
			L_str = f"{L:.4e} H"

		rows = [
			("Inductance (coil)", L_str),
			("DOFs (edges)", str(data.get("n_dofs", ""))),
			("Mesh export", f"{data.get('t_export', 0):.1f} s" if data.get('t_export') else "-"),
			("Assembly", f"{data.get('t_assembly', 0):.1f} s" if data.get('t_assembly') else "-"),
			("LU solve", f"{data.get('t_lu', 0):.1f} s" if data.get('t_lu') else "-"),
			("Post", f"{data.get('t_post', 0):.1f} s" if data.get('t_post') else "-"),
			("Faces", str(data.get("n_faces", ""))),
			("Source area", f"{data.get('source_area', 0):.4e} m^2"),
			("Sink area", f"{data.get('sink_area', 0):.4e} m^2"),
			("Surface area", f"{data.get('surface_area', 0):.4e} m^2"),
			("Constraint |D*J-g|", f"{data.get('constraint_residual', 0):.2e}"),
			("Curve order", str(data.get("curve_order", ""))),
			("FES order", str(data.get("fes_order", ""))),
		]

		# Workpiece ESIM/Dowell results
		if "wp_R_effective" in data:
			freq = data.get("wp_frequency", 0)
			R = data.get("wp_R_effective", 0)
			P = data.get("wp_P_total", 0)
			Q = data.get("wp_Q_total", 0)
			delta_min = data.get("wp_delta_min", 0)
			delta_max = data.get("wp_delta_max", 0)
			n_panels = data.get("wp_n_panels", 0)
			model = data.get("wp_model", "")
			material = data.get("wp_material", "")

			rows.append(("", ""))  # separator
			rows.append(("--- Workpiece ---", f"{model.upper()} / {material}"))
			rows.append(("Frequency", f"{freq:.0f} Hz"))
			rows.append(("R (workpiece)", f"{R:.4e} Ohm"))
			rows.append(("P (workpiece)", f"{P:.4e} W"))
			rows.append(("Q (workpiece)", f"{Q:.4e} var"))
			if delta_min > 0:
				rows.append(("Skin depth", f"{delta_min*1e3:.3f} - {delta_max*1e3:.3f} mm"))
			if "wp_P_density" in data and data["wp_P_density"] > 0:
				rows.append(("P density", f"{data['wp_P_density']:.2e} W/m^2"))
			rows.append(("Elements/Panels", str(n_panels)))
			if "wp_bem_ndof" in data:
				rows.append(("BEM DOFs", str(data["wp_bem_ndof"])))
				rows.append(("BEM assembly", f"{data['wp_bem_t_assembly']:.1f} s"))
			if "wp_fem_ndof" in data:
				rows.append(("FEM DOFs", str(data["wp_fem_ndof"])))
				rows.append(("FEM mesh", f"{data['wp_fem_t_mesh']:.1f} s"))
				rows.append(("FEM solve", f"{data['wp_fem_t_solve']:.1f} s"))

			# Total impedance
			omega = 2 * 3.14159265 * freq
			X_coil = omega * L
			Z_total_R = R
			Z_total_X = X_coil + data.get("wp_X_effective", 0)
			rows.append(("", ""))
			rows.append(("--- Total (I=1A) ---", ""))
			rows.append(("R_total", f"{Z_total_R:.4e} Ohm"))
			rows.append(("X_total", f"{Z_total_X:.4e} Ohm"))
			rows.append(("|Z_total|", f"{(Z_total_R**2 + Z_total_X**2)**0.5:.4e} Ohm"))

		self.result_table.setRowCount(len(rows))
		for i, (param, val) in enumerate(rows):
			self.result_table.setItem(i, 0, QTableWidgetItem(param))
			item = QTableWidgetItem(val)
			item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
			self.result_table.setItem(i, 1, item)


# ================================================================
# IH (FEM) Dialog
# ================================================================

class IHFEMDialog(QDialog):
	"""Dialog for FEM-ESIM induction heating with Kelvin transform."""

	def __init__(self, parent=None):
		super().__init__(parent)
		self.setWindowTitle("IH (FEM): Source/Sink + SIBC")
		self.setMinimumWidth(550)
		self.setMinimumHeight(500)
		self._ext_python = _find_external_python()
		self._process = None
		self._setup_ui()

	def _setup_ui(self):
		layout = QVBoxLayout(self)

		# --- Journal editor ---
		jou_label = QLabel("Cubit Journal (FEM: coil + air + source/sink):")
		jou_label.setStyleSheet("font-weight: bold;")
		layout.addWidget(jou_label)
		self.jou_edit = QPlainTextEdit()
		self.jou_edit.setFont(QFont("Consolas", 9))
		self.jou_edit.setPlainText(InductanceDialog._default_fem_journal())
		self.jou_edit.setMaximumHeight(200)
		layout.addWidget(self.jou_edit)

		jou_btn_row = QHBoxLayout()
		load_btn = QPushButton("Load .jou...")
		load_btn.clicked.connect(self._load_journal)
		jou_btn_row.addWidget(load_btn)
		run_btn = QPushButton("Run Journal")
		run_btn.clicked.connect(self._run_journal)
		jou_btn_row.addWidget(run_btn)
		jou_btn_row.addStretch()
		layout.addLayout(jou_btn_row)

		# --- Block detection ---
		block_group = QGridLayout()
		for row, (lbl, attr) in enumerate([
			("coil:", "_coil_block"),
			("workpiece:", "_workpiece_block"),
			("air:", "_air_block"),
			("source:", "_source_block"),
			("sink:", "_sink_block"),
			("wp_surface:", "_wp_surface"),
			("kelvin:", "_kelvin_block"),
			("outer:", "_outer_surface"),
		]):
			block_group.addWidget(QLabel(lbl), row, 0)
			w = QLabel("(not found)")
			w.setStyleSheet("font-weight: bold; color: red;")
			setattr(self, f"label{attr}", w)
			block_group.addWidget(w, row, 1)
		layout.addLayout(block_group)

		# --- Mesh / FE order ---
		order_group = QGridLayout()
		order_group.addWidget(QLabel("Curve order:"), 0, 0)
		self.curve_spin = QSpinBox()
		self.curve_spin.setRange(1, 3)
		self.curve_spin.setValue(2)
		self.curve_spin.setToolTip("Mesh geometry order (2 = curved elements)")
		order_group.addWidget(self.curve_spin, 0, 1)

		order_group.addWidget(QLabel("FES order:"), 0, 2)
		self.fes_spin = QSpinBox()
		self.fes_spin.setRange(1, 3)
		self.fes_spin.setValue(1)
		self.fes_spin.setToolTip("HCurl polynomial order")
		order_group.addWidget(self.fes_spin, 0, 3)
		layout.addLayout(order_group)

		# --- Kelvin (open boundary) ---
		kelvin_group = QGroupBox("Kelvin (Open Boundary)")
		kelvin_layout = QGridLayout(kelvin_group)

		kelvin_layout.addWidget(QLabel("Radius [m]:"), 0, 0)
		self.kelvin_radius_edit = QLineEdit("")
		self.kelvin_radius_edit.setPlaceholderText("auto from bounding box")
		self.kelvin_radius_edit.setToolTip(
			"Kelvin sphere radius [m]. Leave empty to auto-compute "
			"from model bounding box (1.5x enclosing radius).")
		kelvin_layout.addWidget(self.kelvin_radius_edit, 0, 1)

		self.add_kelvin_btn = QPushButton("Add Kelvin")
		self.add_kelvin_btn.setToolTip(
			"Create Kelvin sphere pair in Cubit:\n"
			"1. Sphere at centroid -> webcut air\n"
			"2. Exterior sphere at offset -> kelvin volume\n"
			"3. Mesh copy for periodic identification")
		self.add_kelvin_btn.clicked.connect(self._add_kelvin)
		kelvin_layout.addWidget(self.add_kelvin_btn, 0, 2)

		kelvin_layout.addWidget(QLabel("Symmetry:"), 1, 0)
		self.kelvin_sym_combo = QComboBox()
		self.kelvin_sym_combo.addItems([
			"Full", "Half (Z)", "Quarter (X,Z)", "Eighth (X,Y,Z)"])
		self.kelvin_sym_combo.setToolTip(
			"Full: no symmetry planes (2 hemispheres)\n"
			"Half (Z): z-symmetry plane (1 hemisphere)\n"
			"Quarter (X,Z): x + z symmetry planes (1 quadrant)\n"
			"Eighth (X,Y,Z): x + y + z symmetry planes (1 octant)")
		kelvin_layout.addWidget(self.kelvin_sym_combo, 1, 1)
		self.kelvin_status = QLabel("")
		self.kelvin_status.setWordWrap(True)
		kelvin_layout.addWidget(self.kelvin_status, 2, 0, 1, 3)

		layout.addWidget(kelvin_group)

		# --- Workpiece settings ---
		wp_group = QGroupBox("Workpiece (SIBC)")
		wp_layout = QGridLayout(wp_group)

		wp_layout.addWidget(QLabel("Impedance:"), 0, 0)
		self.model_combo = QComboBox()
		self.model_combo.addItems(["SIBC", "ESIM"])
		wp_layout.addWidget(self.model_combo, 0, 1)

		wp_layout.addWidget(QLabel("Frequency [Hz]:"), 1, 0)
		self.freq_edit = QLineEdit("7000")
		wp_layout.addWidget(self.freq_edit, 1, 1)

		wp_layout.addWidget(QLabel("Sigma [S/m]:"), 2, 0)
		self.sigma_edit = QLineEdit("2.0e6")
		wp_layout.addWidget(self.sigma_edit, 2, 1)

		# mu_r (SIBC only, linear)
		self.mur_label = QLabel("mu_r:")
		wp_layout.addWidget(self.mur_label, 3, 0)
		self.mur_edit = QLineEdit("100")
		wp_layout.addWidget(self.mur_edit, 3, 1)

		# BH curve (ESIM only, nonlinear)
		self.bh_label = QLabel("BH curve:")
		wp_layout.addWidget(self.bh_label, 4, 0)
		bh_row = QHBoxLayout()
		self.bh_edit = QLineEdit("(built-in Steel)")
		self.bh_edit.setReadOnly(True)
		bh_row.addWidget(self.bh_edit)
		self.bh_browse = QPushButton("...")
		self.bh_browse.setFixedWidth(30)
		self.bh_browse.clicked.connect(self._browse_bh)
		bh_row.addWidget(self.bh_browse)
		wp_layout.addLayout(bh_row, 4, 1)

		# Toggle visibility based on impedance model
		self.model_combo.currentTextChanged.connect(self._on_impedance_changed)
		self._on_impedance_changed(self.model_combo.currentText())

		layout.addWidget(wp_group)

		# --- Result ---
		self.result_label = QLabel("")
		self.result_label.setWordWrap(True)
		layout.addWidget(self.result_label)

		# --- Buttons ---
		btn_row = QHBoxLayout()
		btn_row.addStretch()
		self.solve_btn = QPushButton("Solve")
		self.solve_btn.clicked.connect(self._solve)
		self.solve_btn.setToolTip("3D FEM-SIBC: source/sink + HCurl + Karl iteration")
		btn_row.addWidget(self.solve_btn)
		self.open_btn = QPushButton("Open Result")
		self.open_btn.clicked.connect(self._open_result)
		self.open_btn.setEnabled(False)
		btn_row.addWidget(self.open_btn)
		close_btn = QPushButton("Close")
		close_btn.clicked.connect(self.accept)
		btn_row.addWidget(close_btn)
		layout.addLayout(btn_row)

		# Python info
		py_label = QLabel(f"Python: {self._ext_python or 'Not found'}")
		py_label.setStyleSheet("color: green;" if self._ext_python else "color: red;")
		layout.addWidget(py_label)

	def _on_impedance_changed(self, text):
		"""Show mu_r for SIBC (linear), BH for ESIM (nonlinear)."""
		is_esim = (text == "ESIM")
		self.mur_label.setVisible(not is_esim)
		self.mur_edit.setVisible(not is_esim)
		self.bh_label.setVisible(is_esim)
		self.bh_edit.setVisible(is_esim)
		self.bh_browse.setVisible(is_esim)

	def _browse_bh(self):
		path, _ = QFileDialog.getOpenFileName(
			self, "Load BH Curve", "",
			"Text files (*.txt *.csv *.dat);;All Files (*)")
		if path:
			self.bh_edit.setText(path)

	def _load_journal(self):
		path, _ = QFileDialog.getOpenFileName(
			self, "Load Cubit Journal", "",
			"Cubit Journal (*.jou);;All Files (*)")
		if path:
			try:
				with open(path, "r", encoding="utf-8") as f:
					self.jou_edit.setPlainText(f.read())
			except Exception as e:
				QMessageBox.warning(self, "Error", f"Failed to load: {e}")

	def _run_journal(self):
		text = self.jou_edit.toPlainText().strip()
		if not text:
			return
		for line in text.splitlines():
			line = line.strip()
			if not line or line.startswith("#"):
				continue
			cubit.cmd(line)
		self._detect_blocks()

	def _add_kelvin(self):
		"""Create Kelvin sphere pair with symmetry support.

		Symmetry modes:
		  Full:          2 hemispheres (z-webcut), sphere at model centroid
		  Half (Z):      1 hemisphere, sphere centered on z-symmetry plane
		  Quarter (X,Z): 1 quadrant, sphere centered on x,z-symmetry planes
		  Eighth (X,Y,Z):1 octant, sphere centered on x,y,z-symmetry planes
		"""
		import math

		symmetry = self.kelvin_sym_combo.currentText()

		# Check blocks
		found = {}
		try:
			for bid in cubit.get_block_id_list():
				name = cubit.get_exodus_entity_name("block", bid)
				found[name] = bid
		except Exception:
			pass

		if "kelvin" in found:
			self.kelvin_status.setText("Kelvin block already exists.")
			self.kelvin_status.setStyleSheet("color: orange;")
			return

		air_bid = found.get("air")
		if not air_bid:
			QMessageBox.warning(self, "Error",
				'No "air" block found. Run journal first to create blocks.')
			return

		all_vids = set()
		for bid in cubit.get_block_id_list():
			for v in cubit.get_block_volumes(bid):
				all_vids.add(v)

		if not all_vids:
			QMessageBox.warning(self, "Error", "No meshed volumes found.")
			return

		# Bounding box of all model volumes
		bb_min = [1e30, 1e30, 1e30]
		bb_max = [-1e30, -1e30, -1e30]
		for vid in all_vids:
			bb = cubit.get_bounding_box("volume", vid)
			bb_min[0] = min(bb_min[0], bb[0])
			bb_max[0] = max(bb_max[0], bb[1])
			bb_min[1] = min(bb_min[1], bb[3])
			bb_max[1] = max(bb_max[1], bb[4])
			bb_min[2] = min(bb_min[2], bb[6])
			bb_max[2] = max(bb_max[2], bb[7])

		cx = (bb_min[0] + bb_max[0]) / 2
		cy = (bb_min[1] + bb_max[1]) / 2
		cz = (bb_min[2] + bb_max[2]) / 2

		# Sphere center: at symmetry plane intersection
		# Symmetry plane = bounding box face closest to 0
		def _sym_coord(lo, hi):
			return lo if abs(lo) <= abs(hi) else hi

		if "Quarter" in symmetry or "Eighth" in symmetry:
			sx = _sym_coord(bb_min[0], bb_max[0])
		else:
			sx = cx
		if "Eighth" in symmetry:
			sy = _sym_coord(bb_min[1], bb_max[1])
		else:
			sy = cy
		if symmetry != "Full":
			sz = _sym_coord(bb_min[2], bb_max[2])
		else:
			sz = cz

		# R_kelvin: must enclose model from sphere center
		corners = [(x, y, z) for x in (bb_min[0], bb_max[0])
		                      for y in (bb_min[1], bb_max[1])
		                      for z in (bb_min[2], bb_max[2])]
		max_dist = max(math.sqrt((x - sx)**2 + (y - sy)**2 + (z - sz)**2)
		              for x, y, z in corners)

		r_text = self.kelvin_radius_edit.text().strip()
		if r_text:
			try:
				R_kelvin = float(r_text)
			except ValueError:
				QMessageBox.warning(self, "Error", "Invalid radius value.")
				return
		else:
			R_kelvin = max_dist * 1.3
			self.kelvin_radius_edit.setText(f"{R_kelvin:.6f}")

		if R_kelvin <= max_dist * 0.9:
			QMessageBox.warning(self, "Error",
				f"Kelvin radius ({R_kelvin:.4f}) smaller than "
				f"model extent ({max_dist:.4f}) from sphere center.")
			return

		# Offset for exterior sphere (along x, away from model)
		offset = R_kelvin * 10
		ext_sx = sx + offset

		self.kelvin_status.setText(f"Creating Kelvin sphere ({symmetry})...")
		self.kelvin_status.setStyleSheet("color: blue;")
		QApplication.processEvents()

		try:
			cubit.cmd("set journal off")
			air_vids = list(cubit.get_block_volumes(air_bid))

			# 1. Create interior sphere at symmetry center
			cubit.cmd(f"create sphere radius {R_kelvin}")
			sphere_id = cubit.get_last_id("volume")
			if abs(sx) > 1e-10 or abs(sy) > 1e-10 or abs(sz) > 1e-10:
				cubit.cmd(f"move volume {sphere_id} x {sx} y {sy} z {sz}")

			# 2. Webcut air with sphere (trim air to sphere boundary)
			air_bb = cubit.get_bounding_box("volume", air_vids[0])
			air_diag = math.sqrt(
				air_bb[2]**2 + air_bb[5]**2 + air_bb[8]**2) / 2
			if air_diag > R_kelvin * 1.05:
				cubit.cmd(
					f"webcut volume {' '.join(str(v) for v in air_vids)} "
					f"with sheet body {sphere_id}")
				cubit.cmd(f"delete volume {sphere_id}")
				new_all_vids = set(cubit.get_entities("volume"))
				new_vids = new_all_vids - all_vids
				for nv in new_vids:
					bb = cubit.get_bounding_box("volume", nv)
					nc = ((bb[0]+bb[1])/2, (bb[3]+bb[4])/2, (bb[6]+bb[7])/2)
					dist = math.sqrt(
						(nc[0]-sx)**2 + (nc[1]-sy)**2 + (nc[2]-sz)**2)
					if dist > R_kelvin * 0.8:
						cubit.cmd(f"delete volume {nv}")
			else:
				cubit.cmd(f"delete volume {sphere_id}")

			# 3. Webcut air with z-plane (creates topology curves for mesh copy)
			air_vids_new = list(cubit.get_block_volumes(air_bid))
			if not air_vids_new:
				air_vids_new = air_vids
			cubit.cmd(
				f"webcut volume {' '.join(str(v) for v in air_vids_new)} "
				f"with plane zplane offset {sz}")

			# Track all known model volumes (before exterior sphere)
			all_known = set()
			for bid in cubit.get_block_id_list():
				for v in cubit.get_block_volumes(bid):
					all_known.add(v)

			# 4. Create exterior sphere at offset
			cubit.cmd(f"create sphere radius {R_kelvin}")
			ext_sphere = cubit.get_last_id("volume")
			cubit.cmd(
				f"move volume {ext_sphere} x {ext_sx} y {sy} z {sz}")

			# 5. Webcut exterior sphere: z-plane (always)
			cubit.cmd(
				f"webcut volume {ext_sphere} "
				f"with plane zplane offset {sz}")

			# 6. Additional webcuts for symmetry (exterior only)
			if "Quarter" in symmetry or "Eighth" in symmetry:
				ext_now = [v for v in cubit.get_entities("volume")
				           if v not in all_known]
				if ext_now:
					cubit.cmd(
						f"webcut volume "
						f"{' '.join(str(v) for v in ext_now)} "
						f"with plane xplane offset {ext_sx}")

			if "Eighth" in symmetry:
				ext_now = [v for v in cubit.get_entities("volume")
				           if v not in all_known]
				if ext_now:
					cubit.cmd(
						f"webcut volume "
						f"{' '.join(str(v) for v in ext_now)} "
						f"with plane yplane offset {sy}")

			# 7. Find all exterior volumes (near ext_sx)
			kelvin_vids = []
			for vid in cubit.get_entities("volume"):
				if vid not in all_known:
					bb = cubit.get_bounding_box("volume", vid)
					nc_x = (bb[0] + bb[1]) / 2
					if abs(nc_x - ext_sx) < R_kelvin * 1.1:
						kelvin_vids.append(vid)

			if not kelvin_vids:
				all_now = sorted(cubit.get_entities("volume"))
				kelvin_vids = all_now[-2:]

			# 8. For symmetry modes, keep only matching piece(s)
			if symmetry != "Full":
				dir_x = cx - sx  # model direction from sphere center
				dir_y = cy - sy
				dir_z = cz - sz
				eps = R_kelvin * 0.01

				keep = []
				for kv in kelvin_vids:
					bb = cubit.get_bounding_box("volume", kv)
					nc = ((bb[0]+bb[1])/2, (bb[3]+bb[4])/2,
					      (bb[6]+bb[7])/2)
					ok = True
					# z-check (all symmetry modes)
					if dir_z >= 0 and nc[2] < sz - eps:
						ok = False
					elif dir_z < 0 and nc[2] > sz + eps:
						ok = False
					# x-check (Quarter, Eighth)
					if ok and ("Quarter" in symmetry
					           or "Eighth" in symmetry):
						if dir_x >= 0 and nc[0] < ext_sx - eps:
							ok = False
						elif dir_x < 0 and nc[0] > ext_sx + eps:
							ok = False
					# y-check (Eighth)
					if ok and "Eighth" in symmetry:
						if dir_y >= 0 and nc[1] < sy - eps:
							ok = False
						elif dir_y < 0 and nc[1] > sy + eps:
							ok = False

					if ok:
						keep.append(kv)
					else:
						cubit.cmd(f"delete volume {kv}")
				kelvin_vids = keep

			# 9. Set mesh scheme and size
			kelvin_size = R_kelvin / 3
			for kv in kelvin_vids:
				cubit.cmd(f"volume {kv} scheme tetmesh")
				cubit.cmd(f"volume {kv} size {kelvin_size}")

			# 10. Find interior sphere surfaces (free surfaces on air)
			air_vids_final = []
			for bid in cubit.get_block_id_list():
				name = cubit.get_exodus_entity_name("block", bid)
				if name == "air":
					air_vids_final.extend(cubit.get_block_volumes(bid))

			full_area = 4 * math.pi * R_kelvin**2
			if symmetry == "Full":
				min_area = full_area * 0.15  # hemisphere ~2*pi*R^2
			elif "Half" in symmetry:
				min_area = full_area * 0.15
			elif "Quarter" in symmetry:
				min_area = full_area * 0.05  # quadrant ~pi*R^2
			else:  # Eighth
				min_area = full_area * 0.02  # octant ~pi*R^2/2

			int_hemis = []
			for av in air_vids_final:
				for sid in cubit.get_relatives("volume", av, "surface"):
					adj = cubit.get_relatives("surface", sid, "volume")
					if len(adj) == 1:  # free surface (sphere boundary)
						area = cubit.surface(sid).area()
						if area > min_area:
							int_hemis.append(sid)

			# 11. Find exterior sphere surfaces (largest per kelvin vol)
			ext_hemis = []
			for kv in kelvin_vids:
				best_sid, best_area = None, 0
				for sid in cubit.get_relatives("volume", kv, "surface"):
					area = cubit.surface(sid).area()
					if area > best_area:
						best_area = area
						best_sid = sid
				if best_sid is not None:
					ext_hemis.append(best_sid)

			# Sort by (z, x) for consistent matching
			def _surf_key(sid):
				pts = cubit.get_relatives("surface", sid, "vertex")
				if pts:
					cs = [cubit.vertex(v).coordinates() for v in pts]
					z = sum(c[2] for c in cs) / len(cs)
					x = sum(c[0] for c in cs) / len(cs)
					return (z, x)
				return (0, 0)

			int_hemis.sort(key=_surf_key)
			ext_hemis.sort(key=_surf_key)

			# 12. Copy mesh from interior to exterior
			for src_sid, dst_sid in zip(int_hemis, ext_hemis):
				src_curves = cubit.get_relatives(
					"surface", src_sid, "curve")
				dst_curves = cubit.get_relatives(
					"surface", dst_sid, "curve")
				if not src_curves or not dst_curves:
					continue
				src_c = src_curves[0]
				dst_c = dst_curves[0]
				try:
					src_v = cubit.get_relatives(
						"curve", src_c, "vertex")[0]
					dst_v = cubit.get_relatives(
						"curve", dst_c, "vertex")[0]
					cubit.cmd(
						f"copy mesh surface {src_sid} onto surface "
						f"{dst_sid} source curve {src_c} source vertex "
						f"{src_v} target curve {dst_c} target vertex "
						f"{dst_v}")
				except Exception:
					cubit.cmd(
						f"copy mesh surface {src_sid} onto surface "
						f"{dst_sid} source curve {src_c} "
						f"target curve {dst_c}")

			# 13. Tet mesh kelvin volumes
			for kv in kelvin_vids:
				if len(cubit.get_volume_tets(kv)) == 0:
					cubit.cmd(f"mesh volume {kv}")

			# 14. Create blocks: kelvin, kelvin_int, kelvin_ext
			cubit.cmd("set duplicate block elements on")
			bid_next = max(cubit.get_block_id_list()) + 1
			for kv in kelvin_vids:
				cubit.cmd(f"block {bid_next} add volume {kv}")
			cubit.cmd(f'block {bid_next} name "kelvin"')
			bid_next += 1
			for sid in int_hemis:
				tris = cubit.get_surface_tris(sid)
				if tris:
					cubit.cmd(
						f"block {bid_next} add tri in surface {sid}")
			cubit.cmd(f'block {bid_next} name "kelvin_int"')
			bid_next += 1
			for sid in ext_hemis:
				tris = cubit.get_surface_tris(sid)
				if tris:
					cubit.cmd(
						f"block {bid_next} add tri in surface {sid}")
			cubit.cmd(f'block {bid_next} name "kelvin_ext"')

			n_tets = sum(
				len(cubit.get_volume_tets(kv)) for kv in kelvin_vids)
			self.kelvin_status.setText(
				f"Kelvin added ({symmetry}): R={R_kelvin:.4f}m, "
				f"{len(kelvin_vids)} vol, {n_tets} tets, "
				f"{len(int_hemis)} int + {len(ext_hemis)} ext surf")
			self.kelvin_status.setStyleSheet("color: green;")

		except Exception as e:
			import traceback
			traceback.print_exc()
			self.kelvin_status.setText(f"Error: {e}")
			self.kelvin_status.setStyleSheet("color: red;")
		finally:
			cubit.cmd("set journal on")

		self._detect_blocks()

	def _detect_blocks(self):
		"""Detect FEM blocks and auto-create wp_surface/outer from geometry."""
		import math

		found = {}
		try:
			for bid in cubit.get_block_id_list():
				name = cubit.get_exodus_entity_name("block", bid)
				found[name] = bid
		except Exception:
			pass

		# Auto-classify wp_surface and outer if not already defined
		if "wp_surface" not in found and "air" in found:
			self._auto_create_wp_surface(found)
			# Re-scan
			for bid in cubit.get_block_id_list():
				name = cubit.get_exodus_entity_name("block", bid)
				found[name] = bid

		if "kelvin" in found:
			self._setup_kelvin(found)
			for bid in cubit.get_block_id_list():
				name = cubit.get_exodus_entity_name("block", bid)
				found[name] = bid

		if "outer" not in found and "kelvin" in found:
			self._auto_create_outer(found)
			for bid in cubit.get_block_id_list():
				name = cubit.get_exodus_entity_name("block", bid)
				found[name] = bid

		# Auto-create outer from air sphere free surfaces (no kelvin needed)
		if "outer" not in found and "air" in found:
			self._auto_create_outer_from_air(found)
			for bid in cubit.get_block_id_list():
				name = cubit.get_exodus_entity_name("block", bid)
				found[name] = bid

		# Required blocks (red if missing), optional blocks (orange if missing)
		_optional = {"source", "sink", "wp_surface", "kelvin", "outer"}
		for key, attr in [("coil", "_coil_block"),
		                   ("workpiece", "_workpiece_block"),
		                   ("air", "_air_block"),
		                   ("source", "_source_block"),
		                   ("sink", "_sink_block"),
		                   ("wp_surface", "_wp_surface"),
		                   ("kelvin", "_kelvin_block"),
		                   ("outer", "_outer_surface")]:
			label_w = getattr(self, f"label{attr}")
			if key in found:
				label_w.setText(key)
				label_w.setStyleSheet("font-weight: bold; color: green;")
			elif key in _optional:
				label_w.setText("(optional)")
				label_w.setStyleSheet("font-weight: bold; color: orange;")
			else:
				label_w.setText("(not found)")
				label_w.setStyleSheet("font-weight: bold; color: red;")

		# Coil + air required; source/sink optional (fallback to J_theta)
		has_all = all(k in found for k in ("coil", "air"))
		self.solve_btn.setEnabled(has_all and self._ext_python is not None)

		# Warn about J_theta fallback
		if has_all and "source" not in found:
			self.result_label.setText(
				"Note: No source/sink blocks. "
				"J0*e_theta fallback (axisymmetric torus only).")

	def _setup_kelvin(self, found):
		"""Setup Kelvin mesh: mesh copy + kelvin_int/kelvin_ext surface blocks.

		Sphere must be webcut (e.g., zplane) BEFORE meshing to create curves
		for the copy mesh command. The journal should use:
		  create sphere radius R
		  webcut volume X with plane zplane
		  volume X Y copy move x OFFSET nomesh

		Finds hemisphere surface pairs via get_similar_surfaces, copies mesh
		from interior to exterior hemisphere, tets the kelvin volumes,
		then creates kelvin_int/kelvin_ext surface blocks for periodic ID.
		"""
		try:
			kelvin_bid = found.get("kelvin")
			air_bid = found.get("air")
			if not kelvin_bid or not air_bid:
				return

			kelvin_vids = list(cubit.get_block_volumes(kelvin_bid))
			air_vids = list(cubit.get_block_volumes(air_bid))
			if not kelvin_vids or not air_vids:
				return

			# Check if any kelvin volume is already meshed
			all_meshed = all(
				len(cubit.get_volume_tets(v)) + len(cubit.get_volume_hexes(v)) > 0
				for v in kelvin_vids)
			if all_meshed:
				return

			# Find hemisphere surface pairs:
			# Interior = largest free surface per air volume
			# Exterior = corresponding surface per kelvin volume (via similar)
			int_hemis = []
			for av in air_vids:
				best_sid, best_area = None, 0
				for sid in cubit.get_relatives("volume", av, "surface"):
					adj = cubit.get_relatives("surface", sid, "volume")
					if len(adj) == 1:  # free (not shared with coil/wp)
						area = cubit.surface(sid).area()
						if area > best_area:
							best_area = area
							best_sid = sid
				if best_sid is not None:
					int_hemis.append(best_sid)

			ext_hemis = []
			for kv in kelvin_vids:
				best_sid, best_area = None, 0
				for sid in cubit.get_relatives("volume", kv, "surface"):
					area = cubit.surface(sid).area()
					if area > best_area:
						best_area = area
						best_sid = sid
				if best_sid is not None:
					ext_hemis.append(best_sid)

			if not int_hemis or not ext_hemis:
				return

			# Copy mesh from interior to exterior hemispheres
			cubit.cmd("set duplicate block elements on")
			for src_sid, dst_sid in zip(int_hemis, ext_hemis):
				src_curves = cubit.get_relatives("surface", src_sid, "curve")
				dst_curves = cubit.get_relatives("surface", dst_sid, "curve")
				if not src_curves or not dst_curves:
					continue
				src_c = src_curves[0]
				dst_c = dst_curves[0]
				src_v = int(cubit.curve(src_c).vertices()[0]
				            .entity_name().split()[1])
				dst_v = int(cubit.curve(dst_c).vertices()[0]
				            .entity_name().split()[1])
				cubit.cmd(
					f"copy mesh surface {src_sid} onto surface {dst_sid} "
					f"source curve {src_c} source vertex {src_v} "
					f"target curve {dst_c} target vertex {dst_v}")

			# Tet mesh kelvin volumes
			for kv in kelvin_vids:
				if len(cubit.get_volume_tets(kv)) == 0:
					cubit.cmd(f"volume {kv} scheme tetmesh")
					cubit.cmd(f"volume {kv} size 0.015")
					cubit.cmd(f"mesh volume {kv}")

			# Create kelvin_int and kelvin_ext surface blocks
			bid_next = max(cubit.get_block_id_list()) + 1
			for sid in int_hemis:
				cubit.cmd(f"block {bid_next} add tri in surface {sid}")
			cubit.cmd(f'block {bid_next} name "kelvin_int"')
			bid_next += 1
			for sid in ext_hemis:
				cubit.cmd(f"block {bid_next} add tri in surface {sid}")
			cubit.cmd(f'block {bid_next} name "kelvin_ext"')

		except Exception as e:
			import traceback
			traceback.print_exc()

	def _auto_create_wp_surface(self, found):
		"""Auto-detect workpiece surface = shared faces between workpiece and air."""
		try:
			# Find workpiece volume IDs
			wp_vids = set()
			if "workpiece" in found:
				for v in cubit.get_block_volumes(found["workpiece"]):
					wp_vids.add(v)
			if not wp_vids:
				return

			# Find surfaces shared between workpiece and air
			wp_sids = []
			for vid in wp_vids:
				for sid in cubit.get_relatives("volume", vid, "surface"):
					adj = set(cubit.get_relatives("surface", sid, "volume"))
					# Shared with another volume (air) = interface = wp_surface
					if len(adj) > 1 and not adj.issubset(wp_vids):
						wp_sids.append(sid)

			if wp_sids:
				bid_next = max(cubit.get_block_id_list()) + 1
				cubit.cmd("set duplicate block elements on")
				for sid in wp_sids:
					tris = cubit.get_surface_tris(sid)
					if tris:
						cubit.cmd(f"block {bid_next} add tri in surface {sid}")
				cubit.cmd(f'block {bid_next} name "wp_surface"')
		except Exception:
			pass

	def _auto_create_outer(self, found):
		"""Auto-detect outer surface on Kelvin volume."""
		import math
		try:
			kelvin_bid = found["kelvin"]
			kelvin_vids = list(cubit.get_block_volumes(kelvin_bid))
			if not kelvin_vids:
				return
			outer_sids = []
			for vid in kelvin_vids:
				for sid in cubit.get_relatives("volume", vid, "surface"):
					adj = cubit.get_relatives("surface", sid, "volume")
					if len(adj) > 1:
						continue
					# Free surface = outer boundary
					outer_sids.append(sid)

			if outer_sids:
				bid_next = max(cubit.get_block_id_list()) + 1
				cubit.cmd("set duplicate block elements on")
				for sid in outer_sids:
					tris = cubit.get_surface_tris(sid)
					if tris:
						cubit.cmd(f"block {bid_next} add tri in surface {sid}")
				cubit.cmd(f'block {bid_next} name "outer"')
		except Exception:
			pass

	def _auto_create_outer_from_air(self, found):
		"""Auto-detect outer surface on air volume (free surfaces = far boundary)."""
		try:
			air_bid = found.get("air")
			if not air_bid:
				return
			air_vids = list(cubit.get_block_volumes(air_bid))
			if not air_vids:
				return
			# Collect all volume IDs in the model
			all_block_vids = set()
			for bid in cubit.get_block_id_list():
				for v in cubit.get_block_volumes(bid):
					all_block_vids.add(v)
			outer_sids = []
			for vid in air_vids:
				for sid in cubit.get_relatives("volume", vid, "surface"):
					adj = cubit.get_relatives("surface", sid, "volume")
					# Free surface on air = outer boundary
					if len(adj) == 1:
						outer_sids.append(sid)
			if outer_sids:
				bid_next = max(cubit.get_block_id_list()) + 1
				cubit.cmd("set duplicate block elements on")
				for sid in outer_sids:
					tris = cubit.get_surface_tris(sid)
					if tris:
						cubit.cmd(f"block {bid_next} add tri in surface {sid}")
				cubit.cmd(f'block {bid_next} name "outer"')
		except Exception:
			pass

	def _solve(self):
		"""Run 3D FEM-SIBC via external Python (calc_fem_kelvin.py)."""
		if not self._ext_python:
			QMessageBox.warning(self, "Error", "External Python not found.")
			return

		self.solve_btn.setEnabled(False)
		self.solve_btn.setText("Solving...")
		self.result_label.setText("Computing 3D FEM-SIBC...")

		# Save current Cubit model
		tmpdir = tempfile.mkdtemp(prefix="radia_fem_")
		cub5_file = os.path.join(tmpdir, "model.cub5").replace("\\", "/")
		cubit.cmd(f'save cub5 "{cub5_file}" overwrite')

		freq = self.freq_edit.text().strip() or "7000"
		sigma = self.sigma_edit.text().strip() or "2e6"
		impedance = self.model_combo.currentText().lower()  # "sibc" or "esim"

		self._fem_json = os.path.join(tmpdir, "result.json").replace("\\", "/")
		msh_output = os.path.join(tmpdir, "result.msh").replace("\\", "/")
		calc_script = os.path.join(_this_dir, "calc_fem_kelvin.py")

		# Detect coil wire radius from Cubit bounding box
		a_coil = 0.003
		r_wp = 0.010
		try:
			for bid in cubit.get_block_id_list():
				name = cubit.get_exodus_entity_name("block", bid)
				if name == "coil":
					vids = list(cubit.get_block_volumes(bid))
					if vids:
						bb = cubit.get_bounding_box("volume", vids[0])
						a_coil = min(abs(bb[1] - bb[0]),
						             abs(bb[5] - bb[4])) / 4
				elif name == "workpiece":
					vids = list(cubit.get_block_volumes(bid))
					if vids:
						bb = cubit.get_bounding_box("volume", vids[0])
						r_wp = max(abs(bb[1] - bb[0]),
						           abs(bb[3] - bb[2])) / 2
		except Exception:
			pass

		args = [
			calc_script,
			"--cub5", cub5_file,
			"--order", str(self.curve_spin.value()),
			"--fes-order", str(self.fes_spin.value()),
			"--frequency", freq,
			"--sigma", sigma,
			"--impedance", impedance,
			"--a-coil", str(round(a_coil, 6)),
			"--r-wp", str(round(r_wp, 6)),
			"--msh-output", msh_output,
			"--output", self._fem_json,
		]

		# SIBC: pass mu_r; ESIM: pass BH file
		if impedance == "sibc":
			mu_r = self.mur_edit.text().strip() or "100"
			args += ["--mu-r", mu_r]
		else:
			bh = self.bh_edit.text().strip()
			if bh and bh != "(built-in Steel)":
				args += ["--bh-file", bh]
			args += ["--material", "steel"]

		self._process = QProcess(self)
		self._process.finished.connect(self._on_finished)
		self._process.start(self._ext_python, args)

	def _on_finished(self, exit_code, exit_status):
		self.solve_btn.setEnabled(True)
		self.solve_btn.setText("Solve")
		self._process = None

		data = None
		json_path = getattr(self, '_fem_json', '')
		if json_path and os.path.exists(json_path):
			try:
				with open(json_path, "r") as f:
					data = json.load(f)
			except Exception:
				pass

		if data is None:
			self.result_label.setText(f"Error (exit code {exit_code})")
			return
		if "error" in data:
			self.result_label.setText(f"Error: {data['error']}")
			return

		P = data.get("P_total", 0)
		L = data.get("L", 0)
		src = data.get("source_type", "?")
		self.result_label.setText(
			f"P = {P:.4e} W, L = {L*1e9:.2f} nH\n"
			f"H_t = {data.get('H_t_rms', 0):.2f} A/m, "
			f"|Z_s| = {abs(complex(data.get('Z_s', '0'))):.4e}\n"
			f"DOFs={data.get('ndof', 0)}, source={src}, "
			f"t={data.get('t_total', 0):.1f}s")
		self.open_btn.setEnabled(True)

	def _open_result(self):
		"""Open GMSH result file."""
		json_path = getattr(self, '_fem_json', '')
		if not json_path or not os.path.exists(json_path):
			return
		try:
			with open(json_path, "r") as f:
				data = json.load(f)
			msh = data.get("msh_file", "")
			if not msh or not os.path.exists(msh):
				return
			import subprocess as _sp
			ext_py = self._ext_python
			if ext_py:
				ext_pyw = ext_py.replace("python.exe", "pythonw.exe")
				if not os.path.exists(ext_pyw):
					ext_pyw = ext_py
				code = (
					"import sys, gmsh; "
					"gmsh.initialize(sys.argv, run=True); "
					"gmsh.finalize()")
				_sp.Popen([ext_pyw, "-c", code, msh])
			else:
				os.startfile(msh)
		except Exception:
			pass


# ================================================================
# Accelerator Magnet Dialog
# ================================================================

class AccelMagnetDialog(QDialog):
	"""Dialog for accelerator magnet analysis (Omega-reduced / A-formulation).

	Inputs:
	  - Coil Python script: defines build_coil() -> CoilBuilder
	  - Yoke Cubit journal: iron yoke + air mesh
	  - Formulation: omega (scalar H1) or a (vector HCurl)

	Coil is NOT meshed -- Biot-Savart source from wire path.
	"""

	def __init__(self, parent=None):
		super().__init__(parent)
		self.setWindowTitle("Accelerator Magnet")
		self.setMinimumWidth(580)
		self.setMinimumHeight(580)
		self._ext_python = _find_external_python()
		self._process = None
		self._setup_ui()

	def _setup_ui(self):
		layout = QVBoxLayout(self)

		# --- Coil script ---
		coil_group = QGroupBox("Coil (Biot-Savart source)")
		coil_layout = QHBoxLayout(coil_group)
		coil_layout.addWidget(QLabel("Python script:"))
		self.coil_edit = QLineEdit()
		self.coil_edit.setPlaceholderText("path/to/coil.py (must define build_coil())")
		coil_layout.addWidget(self.coil_edit)
		coil_browse = QPushButton("...")
		coil_browse.setFixedWidth(30)
		coil_browse.clicked.connect(self._browse_coil)
		coil_layout.addWidget(coil_browse)
		preview_btn = QPushButton("Preview")
		preview_btn.setToolTip(
			"Run coil script, generate STEP, import into Cubit.\n"
			"Coil geometry is shown but NOT meshed.")
		preview_btn.clicked.connect(self._preview_coil)
		coil_layout.addWidget(preview_btn)
		layout.addWidget(coil_group)

		# --- Yoke journal editor ---
		jou_label = QLabel("Yoke Journal (iron + air blocks, NO coil):")
		jou_label.setStyleSheet("font-weight: bold;")
		layout.addWidget(jou_label)
		self.jou_edit = QPlainTextEdit()
		self.jou_edit.setFont(QFont("Consolas", 9))
		self.jou_edit.setPlainText(self._default_yoke_journal())
		self.jou_edit.setMaximumHeight(180)
		layout.addWidget(self.jou_edit)

		jou_btn_row = QHBoxLayout()
		load_btn = QPushButton("Load .jou...")
		load_btn.clicked.connect(self._load_journal)
		jou_btn_row.addWidget(load_btn)
		run_btn = QPushButton("Run Journal")
		run_btn.clicked.connect(self._run_journal)
		jou_btn_row.addWidget(run_btn)
		jou_btn_row.addStretch()
		layout.addLayout(jou_btn_row)

		# --- Block detection ---
		block_group = QGridLayout()
		for row, (lbl, attr) in enumerate([
			("yoke:", "_yoke_block"),
			("air:", "_air_block"),
			("sym_tangential:", "_sym_tan"),
			("sym_normal:", "_sym_nor"),
			("kelvin:", "_kelvin_block"),
			("outer:", "_outer_surface"),
		]):
			block_group.addWidget(QLabel(lbl), row, 0)
			w = QLabel("(not found)")
			w.setStyleSheet("font-weight: bold; color: red;")
			setattr(self, f"label{attr}", w)
			block_group.addWidget(w, row, 1)
		layout.addLayout(block_group)

		# --- Formulation + order ---
		form_group = QGridLayout()
		form_group.addWidget(QLabel("Formulation:"), 0, 0)
		self.form_combo = QComboBox()
		self.form_combo.addItems([
			"Omega-reduced (H1)", "A-formulation (HCurl)",
			"MMM (Radia tet)", "MSC (Radia hex)"])
		self.form_combo.setToolTip(
			"Omega-reduced: scalar potential FEM (fastest)\n"
			"A-formulation: vector potential FEM (eddy-current ready)\n"
			"MMM: Radia integral (tet, open boundary, no air mesh)\n"
			"MSC: Radia integral (hex, open boundary, no air mesh)")
		self.form_combo.currentTextChanged.connect(self._on_formulation_changed)
		form_group.addWidget(self.form_combo, 0, 1)

		form_group.addWidget(QLabel("Curve order:"), 1, 0)
		self.curve_spin = QSpinBox()
		self.curve_spin.setRange(1, 3)
		self.curve_spin.setValue(2)
		form_group.addWidget(self.curve_spin, 1, 1)

		form_group.addWidget(QLabel("FES order:"), 1, 2)
		self.fes_spin = QSpinBox()
		self.fes_spin.setRange(1, 3)
		self.fes_spin.setValue(1)
		form_group.addWidget(self.fes_spin, 1, 3)
		layout.addLayout(form_group)

		# --- Kelvin (open boundary) ---
		kelvin_group = QGroupBox("Kelvin (Open Boundary)")
		kelvin_layout = QGridLayout(kelvin_group)
		kelvin_layout.addWidget(QLabel("Radius [m]:"), 0, 0)
		self.kelvin_radius_edit = QLineEdit("")
		self.kelvin_radius_edit.setPlaceholderText("auto from bounding box")
		kelvin_layout.addWidget(self.kelvin_radius_edit, 0, 1)
		self.add_kelvin_btn = QPushButton("Add Kelvin")
		self.add_kelvin_btn.setToolTip(
			"Create Kelvin sphere pair in Cubit:\n"
			"1. Sphere at centroid -> webcut air\n"
			"2. Exterior sphere at offset -> kelvin volume\n"
			"3. Mesh copy for periodic identification")
		self.add_kelvin_btn.clicked.connect(self._add_kelvin)
		kelvin_layout.addWidget(self.add_kelvin_btn, 0, 2)
		kelvin_layout.addWidget(QLabel("Symmetry:"), 1, 0)
		self.kelvin_sym_combo = QComboBox()
		self.kelvin_sym_combo.addItems([
			"Full", "Half (Z)", "Quarter (X,Z)", "Eighth (X,Y,Z)"])
		self.kelvin_sym_combo.setToolTip(
			"Full: no symmetry planes (2 hemispheres)\n"
			"Half (Z): z-symmetry plane (1 hemisphere)\n"
			"Quarter (X,Z): x + z symmetry planes (1 quadrant)\n"
			"Eighth (X,Y,Z): x + y + z symmetry planes (1 octant)")
		kelvin_layout.addWidget(self.kelvin_sym_combo, 1, 1)
		self.kelvin_status = QLabel("")
		self.kelvin_status.setWordWrap(True)
		kelvin_layout.addWidget(self.kelvin_status, 2, 0, 1, 3)
		layout.addWidget(kelvin_group)

		# --- Iron material ---
		mat_group = QGroupBox("Iron Yoke Material")
		mat_layout = QGridLayout(mat_group)

		mat_layout.addWidget(QLabel("Model:"), 0, 0)
		self.mat_combo = QComboBox()
		self.mat_combo.addItems(["Nonlinear (BH curve)", "Linear (mu_r)",
		                         "Hysteresis (Play .hys)"])
		self.mat_combo.currentTextChanged.connect(self._on_material_changed)
		mat_layout.addWidget(self.mat_combo, 0, 1)

		self.mur_label = QLabel("mu_r:")
		mat_layout.addWidget(self.mur_label, 1, 0)
		self.mur_edit = QLineEdit("1000")
		mat_layout.addWidget(self.mur_edit, 1, 1)

		self.bh_label = QLabel("BH curve:")
		mat_layout.addWidget(self.bh_label, 2, 0)
		bh_row = QHBoxLayout()
		self.bh_combo = QComboBox()
		self.bh_combo.addItems(["ELF Steel (100 pts)", "Custom file..."])
		self.bh_combo.setToolTip(
			"ELF Steel: 100-point curve from ELF_MAGIC (CEFC 2020)\n"
			"Custom file: 2-column text file (H [A/m], B [T])")
		self.bh_combo.currentTextChanged.connect(self._on_bh_changed)
		bh_row.addWidget(self.bh_combo)
		self.bh_edit = QLineEdit("")
		self.bh_edit.setReadOnly(True)
		self.bh_edit.setVisible(False)
		bh_row.addWidget(self.bh_edit)
		self.bh_browse = QPushButton("...")
		self.bh_browse.setFixedWidth(30)
		self.bh_browse.setVisible(False)
		self.bh_browse.clicked.connect(self._browse_bh)
		bh_row.addWidget(self.bh_browse)
		mat_layout.addLayout(bh_row, 2, 1)

		# BH curve info display
		self.bh_info = QLabel("ELF Steel: 100 points, Bmax=2.6T (CEFC 2020)")
		self.bh_info.setStyleSheet("color: #666; font-size: 9pt;")
		self.bh_info.setWordWrap(True)
		mat_layout.addWidget(self.bh_info, 3, 0, 1, 2)

		# .hys file (Hysteresis mode)
		self.hys_label = QLabel(".hys file:")
		mat_layout.addWidget(self.hys_label, 4, 0)
		hys_row = QHBoxLayout()
		self.hys_edit = QLineEdit("")
		self.hys_edit.setPlaceholderText("JMAG .hys file")
		hys_row.addWidget(self.hys_edit)
		self.hys_browse = QPushButton("...")
		self.hys_browse.setFixedWidth(30)
		self.hys_browse.clicked.connect(self._browse_hys)
		hys_row.addWidget(self.hys_browse)
		mat_layout.addLayout(hys_row, 4, 1)

		self.hys_info = QLabel("")
		self.hys_info.setStyleSheet("color: #666; font-size: 9pt;")
		self.hys_info.setWordWrap(True)
		mat_layout.addWidget(self.hys_info, 5, 0, 1, 2)

		# Quasi-static steps (Hysteresis mode)
		self.nsteps_label = QLabel("Steps (0->I->0):")
		mat_layout.addWidget(self.nsteps_label, 6, 0)
		self.nsteps_spin = QSpinBox()
		self.nsteps_spin.setRange(1, 200)
		self.nsteps_spin.setValue(10)
		self.nsteps_spin.setToolTip(
			"Number of quasi-static steps.\n"
			"1 = single excitation (no ramp).\n"
			"N>1 = ramp up then down (hysteresis loop).")
		mat_layout.addWidget(self.nsteps_spin, 6, 1)

		self._on_material_changed(self.mat_combo.currentText())
		layout.addWidget(mat_group)

		# --- Picard iteration ---
		iter_group = QGridLayout()
		iter_group.addWidget(QLabel("Max iterations:"), 0, 0)
		self.maxiter_spin = QSpinBox()
		self.maxiter_spin.setRange(1, 200)
		self.maxiter_spin.setValue(30)
		iter_group.addWidget(self.maxiter_spin, 0, 1)

		iter_group.addWidget(QLabel("Relaxation:"), 0, 2)
		self.relax_edit = QLineEdit("0.3")
		self.relax_edit.setToolTip("Under-relaxation (0=full, 0.5=half)")
		iter_group.addWidget(self.relax_edit, 0, 3)

		self.newton_check = QCheckBox("Newton")
		self.newton_check.setToolTip(
			"Newton: anisotropic tangent (quadratic convergence)\n"
			"Picard: isotropic chord (linear convergence, more robust)")
		self.newton_check.setChecked(True)
		iter_group.addWidget(self.newton_check, 0, 4)
		layout.addLayout(iter_group)

		# --- Result ---
		self.result_label = QLabel("")
		self.result_label.setWordWrap(True)
		layout.addWidget(self.result_label)

		# --- Buttons ---
		btn_row = QHBoxLayout()
		btn_row.addStretch()
		self.solve_btn = QPushButton("Solve")
		self.solve_btn.clicked.connect(self._solve)
		btn_row.addWidget(self.solve_btn)
		self.open_btn = QPushButton("Open Result")
		self.open_btn.clicked.connect(self._open_result)
		self.open_btn.setEnabled(False)
		btn_row.addWidget(self.open_btn)
		close_btn = QPushButton("Close")
		close_btn.clicked.connect(self.accept)
		btn_row.addWidget(close_btn)
		layout.addLayout(btn_row)

		# Python info
		py_label = QLabel(f"Python: {self._ext_python or 'Not found'}")
		py_label.setStyleSheet("color: green;" if self._ext_python else "color: red;")
		layout.addWidget(py_label)

	@staticmethod
	def _default_yoke_journal():
		"""Default journal: simple C-dipole yoke + air sphere."""
		lines = [
			"# Accelerator dipole yoke (simplified C-magnet)",
			"# Units: meters",
			"reset",
			"",
			"# === Iron yoke (C-shape) ===",
			"# Outer frame",
			"brick x 0.300 y 0.200 z 0.400",
			"# Aperture (beam pipe hole)",
			"brick x 0.100 y 0.200 z 0.060",
			"move volume 2 x 0 y 0 z 0",
			"subtract volume 2 from volume 1",
			"",
			"# === Air sphere ===",
			"create sphere radius 0.300",
			"subtract volume 1 from volume 3 keep",
			"# Volume 1 = yoke, Volume 3 = air (sphere - yoke)",
			"",
			"# === Imprint and merge ===",
			"imprint volume all",
			"merge volume all",
			"",
			"# === Mesh ===",
			"volume all scheme tetmesh",
			"volume 1 size 0.020",
			"volume 3 size 0.030",
			"mesh volume 1 3",
			"",
			"# === Blocks ===",
			"set duplicate block elements on",
			"block 1 add volume 1",
			'block 1 name "yoke"',
			"block 2 add volume 3",
			'block 2 name "air"',
			"",
			"# Kelvin: use 'Add Kelvin' button in dialog",
			"# outer: auto-created by dialog",
		]
		return "\n".join(lines) + "\n"

	def _on_material_changed(self, text):
		is_linear = "Linear" in text
		is_hys = "Hysteresis" in text
		is_bh = not is_linear and not is_hys
		self.mur_label.setVisible(is_linear)
		self.mur_edit.setVisible(is_linear)
		self.bh_label.setVisible(is_bh)
		self.bh_edit.setVisible(is_bh)
		self.bh_browse.setVisible(is_bh)
		self.bh_info.setVisible(is_bh)
		self.hys_label.setVisible(is_hys)
		self.hys_edit.setVisible(is_hys)
		self.hys_browse.setVisible(is_hys)
		self.hys_info.setVisible(is_hys)
		self.nsteps_label.setVisible(is_hys)
		self.nsteps_spin.setVisible(is_hys)
		# Hysteresis works with both formulations (Hantila decomposition)

	def _browse_coil(self):
		path, _ = QFileDialog.getOpenFileName(
			self, "Select Coil Script", "",
			"Python Files (*.py);;All Files (*)")
		if path:
			self.coil_edit.setText(path)

	def _preview_coil(self):
		"""Run coil script -> STEP -> import into Cubit for visualization."""
		coil_path = self.coil_edit.text().strip()
		if not coil_path or not os.path.exists(coil_path):
			QMessageBox.warning(self, "Error", "Set a valid coil script path first.")
			return

		if not self._ext_python:
			QMessageBox.warning(self, "Error", "External Python not found.")
			return

		# Generate STEP via external Python (needs scipy for CoilBuilder)
		step_file = os.path.join(
			tempfile.gettempdir(), "radia_coil_preview.step").replace("\\", "/")
		code = (
			f"import sys, os; "
			f"sys.path.insert(0, os.path.dirname(r'{coil_path}')); "
			f"sys.path.insert(0, r'{_pkg_root}'); "
			f"import importlib.util; "
			f"spec = importlib.util.spec_from_file_location('cm', r'{coil_path}'); "
			f"mod = importlib.util.module_from_spec(spec); "
			f"spec.loader.exec_module(mod); "
			f"coil = mod.build_coil(); "
			f"coil.write_step(r'{step_file}'); "
			f"print('OK')")

		try:
			r = subprocess.run(
				[self._ext_python, "-c", code],
				capture_output=True, text=True, timeout=30)
			if r.returncode != 0 or "OK" not in r.stdout:
				err = r.stderr.strip().split("\n")[-1] if r.stderr else "Unknown error"
				QMessageBox.warning(self, "Coil Error", f"STEP generation failed:\n{err}")
				return
		except subprocess.TimeoutExpired:
			QMessageBox.warning(self, "Error", "Coil script timed out.")
			return

		if not os.path.exists(step_file):
			QMessageBox.warning(self, "Error", "STEP file not generated.")
			return

		# Import STEP as free body (not added to any block -> not meshed)
		try:
			cubit.cmd(f'import step "{step_file}" heal')
			# Find newly imported volumes (not in any block)
			coil_vids = []
			block_vids = set()
			for bid in cubit.get_block_id_list():
				for v in cubit.get_block_volumes(bid):
					block_vids.add(v)
			for v in cubit.get_entities("volume"):
				if v not in block_vids:
					coil_vids.append(v)
			# Style: lightblue, transparent, no mesh scheme
			for v in coil_vids:
				cubit.cmd(f"color volume {v} lightblue")
				cubit.cmd(f"volume {v} visibility on")
		except Exception as e:
			QMessageBox.warning(self, "Import Error", str(e))

	def _on_bh_changed(self, text):
		"""Update BH curve info and show/hide custom file widgets."""
		is_custom = "Custom" in text
		self.bh_edit.setVisible(is_custom)
		self.bh_browse.setVisible(is_custom)
		if "ELF" in text:
			self.bh_info.setText(
				"ELF Steel: 100 points, Bmax=2.6T (CEFC 2020)")
		elif "Simple" in text:
			self.bh_info.setText(
				"Simple Steel: 12 points, Bmax=2.0T")
		elif is_custom and self.bh_edit.text():
			pass  # keep current info from _browse_bh
		else:
			self.bh_info.setText("")

	def _browse_bh(self):
		path, _ = QFileDialog.getOpenFileName(
			self, "Load BH Curve", "",
			"Text files (*.txt *.csv *.dat);;All Files (*)")
		if path:
			self.bh_edit.setText(path)
			try:
				import numpy as _np
				data = _np.loadtxt(path)
				if data.ndim == 2 and data.shape[1] >= 2:
					self.bh_info.setText(
						f"{os.path.basename(path)}: {data.shape[0]} pts, "
						f"Hmax={data[-1,0]:.0f} A/m, Bmax={data[-1,1]:.2f} T")
			except Exception:
				self.bh_info.setText(f"Loaded: {os.path.basename(path)}")

	def _browse_hys(self):
		path, _ = QFileDialog.getOpenFileName(
			self, "Load Hysteresis Data", "",
			"HYS files (*.hys);;MAT files (*.mat);;All Files (*)")
		if path:
			self.hys_edit.setText(path)
			try:
				# Read .hys header to show info
				with open(path, 'r') as f:
					_model = f.readline().strip()
					dims = f.readline().split()
					nx, ny = int(dims[1]), int(dims[2])
				self.hys_info.setText(
					f"{os.path.basename(path)}: {ny} loops, "
					f"{nx} pts/loop")
			except Exception:
				self.hys_info.setText(f"Loaded: {os.path.basename(path)}")

	def _load_journal(self):
		path, _ = QFileDialog.getOpenFileName(
			self, "Load Cubit Journal", "",
			"Cubit Journal (*.jou);;All Files (*)")
		if path:
			try:
				with open(path, "r", encoding="utf-8") as f:
					self.jou_edit.setPlainText(f.read())
			except Exception as e:
				QMessageBox.warning(self, "Error", f"Failed to load: {e}")

	def _run_journal(self):
		text = self.jou_edit.toPlainText().strip()
		if not text:
			return
		for line in text.splitlines():
			line = line.strip()
			if not line or line.startswith("#"):
				continue
			cubit.cmd(line)
		self._detect_blocks()

	def _on_formulation_changed(self, text):
		"""Enable/disable UI elements based on formulation selection.

		FEM (Omega/A): needs Kelvin, air block, curve/fes order.
		Radia (MMM/MSC): no air, no Kelvin, no FES order.
		"""
		is_radia = "Radia" in text
		self.add_kelvin_btn.setEnabled(not is_radia)
		self.kelvin_radius_edit.setEnabled(not is_radia)
		self.kelvin_sym_combo.setEnabled(not is_radia)
		self.curve_spin.setEnabled(not is_radia)
		self.fes_spin.setEnabled(not is_radia)
		if is_radia:
			self.kelvin_status.setText(
				"Radia solver: no air/Kelvin needed (natural open BC)")
			self.kelvin_status.setStyleSheet("color: gray;")

	def _add_kelvin(self):
		"""Create Kelvin sphere pair with symmetry support.

		Symmetry modes:
		  Full:          2 hemispheres (z-webcut), sphere at model centroid
		  Half (Z):      1 hemisphere, sphere centered on z-symmetry plane
		  Quarter (X,Z): 1 quadrant, sphere centered on x,z-symmetry planes
		  Eighth (X,Y,Z):1 octant, sphere centered on x,y,z-symmetry planes
		"""
		import math

		symmetry = self.kelvin_sym_combo.currentText()

		# Check blocks
		found = {}
		try:
			for bid in cubit.get_block_id_list():
				name = cubit.get_exodus_entity_name("block", bid)
				found[name] = bid
		except Exception:
			pass

		if "kelvin" in found:
			self.kelvin_status.setText("Kelvin block already exists.")
			self.kelvin_status.setStyleSheet("color: orange;")
			return

		air_bid = found.get("air")
		if not air_bid:
			QMessageBox.warning(self, "Error",
				'No "air" block found. Run journal first to create blocks.')
			return

		all_vids = set()
		for bid in cubit.get_block_id_list():
			for v in cubit.get_block_volumes(bid):
				all_vids.add(v)

		if not all_vids:
			QMessageBox.warning(self, "Error", "No meshed volumes found.")
			return

		# Bounding box of all model volumes
		bb_min = [1e30, 1e30, 1e30]
		bb_max = [-1e30, -1e30, -1e30]
		for vid in all_vids:
			bb = cubit.get_bounding_box("volume", vid)
			bb_min[0] = min(bb_min[0], bb[0])
			bb_max[0] = max(bb_max[0], bb[1])
			bb_min[1] = min(bb_min[1], bb[3])
			bb_max[1] = max(bb_max[1], bb[4])
			bb_min[2] = min(bb_min[2], bb[6])
			bb_max[2] = max(bb_max[2], bb[7])

		cx = (bb_min[0] + bb_max[0]) / 2
		cy = (bb_min[1] + bb_max[1]) / 2
		cz = (bb_min[2] + bb_max[2]) / 2

		# Sphere center: at symmetry plane intersection
		# Symmetry plane = bounding box face closest to 0
		def _sym_coord(lo, hi):
			return lo if abs(lo) <= abs(hi) else hi

		if "Quarter" in symmetry or "Eighth" in symmetry:
			sx = _sym_coord(bb_min[0], bb_max[0])
		else:
			sx = cx
		if "Eighth" in symmetry:
			sy = _sym_coord(bb_min[1], bb_max[1])
		else:
			sy = cy
		if symmetry != "Full":
			sz = _sym_coord(bb_min[2], bb_max[2])
		else:
			sz = cz

		# R_kelvin: must enclose model from sphere center
		corners = [(x, y, z) for x in (bb_min[0], bb_max[0])
		                      for y in (bb_min[1], bb_max[1])
		                      for z in (bb_min[2], bb_max[2])]
		max_dist = max(math.sqrt((x - sx)**2 + (y - sy)**2 + (z - sz)**2)
		              for x, y, z in corners)

		r_text = self.kelvin_radius_edit.text().strip()
		if r_text:
			try:
				R_kelvin = float(r_text)
			except ValueError:
				QMessageBox.warning(self, "Error", "Invalid radius value.")
				return
		else:
			R_kelvin = max_dist * 1.3
			self.kelvin_radius_edit.setText(f"{R_kelvin:.6f}")

		if R_kelvin <= max_dist * 0.9:
			QMessageBox.warning(self, "Error",
				f"Kelvin radius ({R_kelvin:.4f}) smaller than "
				f"model extent ({max_dist:.4f}) from sphere center.")
			return

		# Offset for exterior sphere (along x, away from model)
		offset = R_kelvin * 10
		ext_sx = sx + offset

		self.kelvin_status.setText(f"Creating Kelvin sphere ({symmetry})...")
		self.kelvin_status.setStyleSheet("color: blue;")
		QApplication.processEvents()

		try:
			cubit.cmd("set journal off")
			air_vids = list(cubit.get_block_volumes(air_bid))

			# 1. Create interior sphere at symmetry center
			cubit.cmd(f"create sphere radius {R_kelvin}")
			sphere_id = cubit.get_last_id("volume")
			if abs(sx) > 1e-10 or abs(sy) > 1e-10 or abs(sz) > 1e-10:
				cubit.cmd(f"move volume {sphere_id} x {sx} y {sy} z {sz}")

			# 2. Webcut air with sphere (trim air to sphere boundary)
			air_bb = cubit.get_bounding_box("volume", air_vids[0])
			air_diag = math.sqrt(
				air_bb[2]**2 + air_bb[5]**2 + air_bb[8]**2) / 2
			if air_diag > R_kelvin * 1.05:
				cubit.cmd(
					f"webcut volume {' '.join(str(v) for v in air_vids)} "
					f"with sheet body {sphere_id}")
				cubit.cmd(f"delete volume {sphere_id}")
				new_all_vids = set(cubit.get_entities("volume"))
				new_vids = new_all_vids - all_vids
				for nv in new_vids:
					bb = cubit.get_bounding_box("volume", nv)
					nc = ((bb[0]+bb[1])/2, (bb[3]+bb[4])/2, (bb[6]+bb[7])/2)
					dist = math.sqrt(
						(nc[0]-sx)**2 + (nc[1]-sy)**2 + (nc[2]-sz)**2)
					if dist > R_kelvin * 0.8:
						cubit.cmd(f"delete volume {nv}")
			else:
				cubit.cmd(f"delete volume {sphere_id}")

			# 3. Webcut air with z-plane (creates topology curves for mesh copy)
			air_vids_new = list(cubit.get_block_volumes(air_bid))
			if not air_vids_new:
				air_vids_new = air_vids
			cubit.cmd(
				f"webcut volume {' '.join(str(v) for v in air_vids_new)} "
				f"with plane zplane offset {sz}")

			# Track all known model volumes (before exterior sphere)
			all_known = set()
			for bid in cubit.get_block_id_list():
				for v in cubit.get_block_volumes(bid):
					all_known.add(v)

			# 4. Create exterior sphere at offset
			cubit.cmd(f"create sphere radius {R_kelvin}")
			ext_sphere = cubit.get_last_id("volume")
			cubit.cmd(
				f"move volume {ext_sphere} x {ext_sx} y {sy} z {sz}")

			# 5. Webcut exterior sphere: z-plane (always)
			cubit.cmd(
				f"webcut volume {ext_sphere} "
				f"with plane zplane offset {sz}")

			# 6. Additional webcuts for symmetry (exterior only)
			if "Quarter" in symmetry or "Eighth" in symmetry:
				ext_now = [v for v in cubit.get_entities("volume")
				           if v not in all_known]
				if ext_now:
					cubit.cmd(
						f"webcut volume "
						f"{' '.join(str(v) for v in ext_now)} "
						f"with plane xplane offset {ext_sx}")

			if "Eighth" in symmetry:
				ext_now = [v for v in cubit.get_entities("volume")
				           if v not in all_known]
				if ext_now:
					cubit.cmd(
						f"webcut volume "
						f"{' '.join(str(v) for v in ext_now)} "
						f"with plane yplane offset {sy}")

			# 7. Find all exterior volumes (near ext_sx)
			kelvin_vids = []
			for vid in cubit.get_entities("volume"):
				if vid not in all_known:
					bb = cubit.get_bounding_box("volume", vid)
					nc_x = (bb[0] + bb[1]) / 2
					if abs(nc_x - ext_sx) < R_kelvin * 1.1:
						kelvin_vids.append(vid)

			if not kelvin_vids:
				all_now = sorted(cubit.get_entities("volume"))
				kelvin_vids = all_now[-2:]

			# 8. For symmetry modes, keep only matching piece(s)
			if symmetry != "Full":
				dir_x = cx - sx  # model direction from sphere center
				dir_y = cy - sy
				dir_z = cz - sz
				eps = R_kelvin * 0.01

				keep = []
				for kv in kelvin_vids:
					bb = cubit.get_bounding_box("volume", kv)
					nc = ((bb[0]+bb[1])/2, (bb[3]+bb[4])/2,
					      (bb[6]+bb[7])/2)
					ok = True
					# z-check (all symmetry modes)
					if dir_z >= 0 and nc[2] < sz - eps:
						ok = False
					elif dir_z < 0 and nc[2] > sz + eps:
						ok = False
					# x-check (Quarter, Eighth)
					if ok and ("Quarter" in symmetry
					           or "Eighth" in symmetry):
						if dir_x >= 0 and nc[0] < ext_sx - eps:
							ok = False
						elif dir_x < 0 and nc[0] > ext_sx + eps:
							ok = False
					# y-check (Eighth)
					if ok and "Eighth" in symmetry:
						if dir_y >= 0 and nc[1] < sy - eps:
							ok = False
						elif dir_y < 0 and nc[1] > sy + eps:
							ok = False

					if ok:
						keep.append(kv)
					else:
						cubit.cmd(f"delete volume {kv}")
				kelvin_vids = keep

			# 9. Set mesh scheme and size
			kelvin_size = R_kelvin / 3
			for kv in kelvin_vids:
				cubit.cmd(f"volume {kv} scheme tetmesh")
				cubit.cmd(f"volume {kv} size {kelvin_size}")

			# 10. Find interior sphere surfaces (free surfaces on air)
			air_vids_final = []
			for bid in cubit.get_block_id_list():
				name = cubit.get_exodus_entity_name("block", bid)
				if name == "air":
					air_vids_final.extend(cubit.get_block_volumes(bid))

			full_area = 4 * math.pi * R_kelvin**2
			if symmetry == "Full":
				min_area = full_area * 0.15  # hemisphere ~2*pi*R^2
			elif "Half" in symmetry:
				min_area = full_area * 0.15
			elif "Quarter" in symmetry:
				min_area = full_area * 0.05  # quadrant ~pi*R^2
			else:  # Eighth
				min_area = full_area * 0.02  # octant ~pi*R^2/2

			int_hemis = []
			for av in air_vids_final:
				for sid in cubit.get_relatives("volume", av, "surface"):
					adj = cubit.get_relatives("surface", sid, "volume")
					if len(adj) == 1:  # free surface (sphere boundary)
						area = cubit.surface(sid).area()
						if area > min_area:
							int_hemis.append(sid)

			# 11. Find exterior sphere surfaces (largest per kelvin vol)
			ext_hemis = []
			for kv in kelvin_vids:
				best_sid, best_area = None, 0
				for sid in cubit.get_relatives("volume", kv, "surface"):
					area = cubit.surface(sid).area()
					if area > best_area:
						best_area = area
						best_sid = sid
				if best_sid is not None:
					ext_hemis.append(best_sid)

			# Sort by (z, x) for consistent matching
			def _surf_key(sid):
				pts = cubit.get_relatives("surface", sid, "vertex")
				if pts:
					cs = [cubit.vertex(v).coordinates() for v in pts]
					z = sum(c[2] for c in cs) / len(cs)
					x = sum(c[0] for c in cs) / len(cs)
					return (z, x)
				return (0, 0)

			int_hemis.sort(key=_surf_key)
			ext_hemis.sort(key=_surf_key)

			# 12. Copy mesh from interior to exterior
			for src_sid, dst_sid in zip(int_hemis, ext_hemis):
				src_curves = cubit.get_relatives(
					"surface", src_sid, "curve")
				dst_curves = cubit.get_relatives(
					"surface", dst_sid, "curve")
				if not src_curves or not dst_curves:
					continue
				src_c = src_curves[0]
				dst_c = dst_curves[0]
				try:
					src_v = cubit.get_relatives(
						"curve", src_c, "vertex")[0]
					dst_v = cubit.get_relatives(
						"curve", dst_c, "vertex")[0]
					cubit.cmd(
						f"copy mesh surface {src_sid} onto surface "
						f"{dst_sid} source curve {src_c} source vertex "
						f"{src_v} target curve {dst_c} target vertex "
						f"{dst_v}")
				except Exception:
					cubit.cmd(
						f"copy mesh surface {src_sid} onto surface "
						f"{dst_sid} source curve {src_c} "
						f"target curve {dst_c}")

			# 13. Tet mesh kelvin volumes
			for kv in kelvin_vids:
				if len(cubit.get_volume_tets(kv)) == 0:
					cubit.cmd(f"mesh volume {kv}")

			# 14. Create blocks: kelvin, kelvin_int, kelvin_ext
			cubit.cmd("set duplicate block elements on")
			bid_next = max(cubit.get_block_id_list()) + 1
			for kv in kelvin_vids:
				cubit.cmd(f"block {bid_next} add volume {kv}")
			cubit.cmd(f'block {bid_next} name "kelvin"')
			bid_next += 1
			for sid in int_hemis:
				tris = cubit.get_surface_tris(sid)
				if tris:
					cubit.cmd(
						f"block {bid_next} add tri in surface {sid}")
			cubit.cmd(f'block {bid_next} name "kelvin_int"')
			bid_next += 1
			for sid in ext_hemis:
				tris = cubit.get_surface_tris(sid)
				if tris:
					cubit.cmd(
						f"block {bid_next} add tri in surface {sid}")
			cubit.cmd(f'block {bid_next} name "kelvin_ext"')

			n_tets = sum(
				len(cubit.get_volume_tets(kv)) for kv in kelvin_vids)
			self.kelvin_status.setText(
				f"Kelvin added ({symmetry}): R={R_kelvin:.4f}m, "
				f"{len(kelvin_vids)} vol, {n_tets} tets, "
				f"{len(int_hemis)} int + {len(ext_hemis)} ext surf")
			self.kelvin_status.setStyleSheet("color: green;")

		except Exception as e:
			import traceback
			traceback.print_exc()
			self.kelvin_status.setText(f"Error: {e}")
			self.kelvin_status.setStyleSheet("color: red;")
		finally:
			cubit.cmd("set journal on")

		self._detect_blocks()

	def _detect_blocks(self):
		"""Detect yoke/air/kelvin blocks and auto-create outer."""
		found = {}
		try:
			for bid in cubit.get_block_id_list():
				name = cubit.get_exodus_entity_name("block", bid)
				found[name] = bid
		except Exception:
			pass

		# Auto-create outer from air/kelvin free surfaces
		if "outer" not in found:
			self._auto_create_outer(found)
			for bid in cubit.get_block_id_list():
				name = cubit.get_exodus_entity_name("block", bid)
				found[name] = bid

		# Setup kelvin mesh if needed
		if "kelvin" in found and "kelvin_int" not in found:
			self._setup_kelvin(found)
			for bid in cubit.get_block_id_list():
				name = cubit.get_exodus_entity_name("block", bid)
				found[name] = bid

		_optional = {"sym_tangential", "sym_normal", "kelvin", "outer"}
		for key, attr in [("yoke", "_yoke_block"),
		                   ("air", "_air_block"),
		                   ("sym_tangential", "_sym_tan"),
		                   ("sym_normal", "_sym_nor"),
		                   ("kelvin", "_kelvin_block"),
		                   ("outer", "_outer_surface")]:
			label_w = getattr(self, f"label{attr}")
			if key in found:
				label_w.setText(key)
				label_w.setStyleSheet("font-weight: bold; color: green;")
			elif key in _optional:
				label_w.setText("(optional)")
				label_w.setStyleSheet("font-weight: bold; color: orange;")
			else:
				label_w.setText("(not found)")
				label_w.setStyleSheet("font-weight: bold; color: red;")

		has_all = all(k in found for k in ("yoke", "air"))
		has_coil = bool(self.coil_edit.text().strip())
		self.solve_btn.setEnabled(
			has_all and has_coil and self._ext_python is not None)

		if has_all and not has_coil:
			self.result_label.setText("Set coil script path to enable Solve.")

	def _setup_kelvin(self, found):
		"""Delegate to IHFEMDialog._setup_kelvin pattern."""
		try:
			kelvin_bid = found.get("kelvin")
			air_bid = found.get("air")
			if not kelvin_bid or not air_bid:
				return
			kelvin_vids = list(cubit.get_block_volumes(kelvin_bid))
			air_vids = list(cubit.get_block_volumes(air_bid))
			if not kelvin_vids or not air_vids:
				return
			all_meshed = all(
				len(cubit.get_volume_tets(v)) + len(cubit.get_volume_hexes(v)) > 0
				for v in kelvin_vids)
			if all_meshed:
				return

			int_hemis = []
			for av in air_vids:
				best_sid, best_area = None, 0
				for sid in cubit.get_relatives("volume", av, "surface"):
					adj = cubit.get_relatives("surface", sid, "volume")
					if len(adj) == 1:
						area = cubit.surface(sid).area()
						if area > best_area:
							best_area = area
							best_sid = sid
				if best_sid is not None:
					int_hemis.append(best_sid)

			ext_hemis = []
			for kv in kelvin_vids:
				best_sid, best_area = None, 0
				for sid in cubit.get_relatives("volume", kv, "surface"):
					area = cubit.surface(sid).area()
					if area > best_area:
						best_area = area
						best_sid = sid
				if best_sid is not None:
					ext_hemis.append(best_sid)

			if not int_hemis or not ext_hemis:
				return

			cubit.cmd("set duplicate block elements on")
			for src_sid, dst_sid in zip(int_hemis, ext_hemis):
				src_curves = cubit.get_relatives("surface", src_sid, "curve")
				dst_curves = cubit.get_relatives("surface", dst_sid, "curve")
				if not src_curves or not dst_curves:
					continue
				src_c = src_curves[0]
				dst_c = dst_curves[0]
				src_v = int(cubit.curve(src_c).vertices()[0]
				            .entity_name().split()[1])
				dst_v = int(cubit.curve(dst_c).vertices()[0]
				            .entity_name().split()[1])
				cubit.cmd(
					f"copy mesh surface {src_sid} onto surface {dst_sid} "
					f"source curve {src_c} source vertex {src_v} "
					f"target curve {dst_c} target vertex {dst_v}")

			for kv in kelvin_vids:
				if len(cubit.get_volume_tets(kv)) == 0:
					cubit.cmd(f"volume {kv} scheme tetmesh")
					cubit.cmd(f"volume {kv} size 0.015")
					cubit.cmd(f"mesh volume {kv}")

			bid_next = max(cubit.get_block_id_list()) + 1
			for sid in int_hemis:
				cubit.cmd(f"block {bid_next} add tri in surface {sid}")
			cubit.cmd(f'block {bid_next} name "kelvin_int"')
			bid_next += 1
			for sid in ext_hemis:
				cubit.cmd(f"block {bid_next} add tri in surface {sid}")
			cubit.cmd(f'block {bid_next} name "kelvin_ext"')

		except Exception:
			import traceback
			traceback.print_exc()

	def _auto_create_outer(self, found):
		"""Auto-detect outer surface on air/kelvin free surfaces."""
		try:
			# Prefer kelvin outer surfaces, then air
			for block_name in ("kelvin", "air"):
				bid = found.get(block_name)
				if not bid:
					continue
				vids = list(cubit.get_block_volumes(bid))
				if not vids:
					continue
				outer_sids = []
				for vid in vids:
					for sid in cubit.get_relatives("volume", vid, "surface"):
						adj = cubit.get_relatives("surface", sid, "volume")
						if len(adj) == 1:
							outer_sids.append(sid)
				if outer_sids:
					bid_next = max(cubit.get_block_id_list()) + 1
					cubit.cmd("set duplicate block elements on")
					for sid in outer_sids:
						tris = cubit.get_surface_tris(sid)
						if tris:
							cubit.cmd(f"block {bid_next} add tri in surface {sid}")
					cubit.cmd(f'block {bid_next} name "outer"')
					return
		except Exception:
			pass

	def _solve(self):
		"""Run accelerator magnet solver via external Python."""
		if not self._ext_python:
			QMessageBox.warning(self, "Error", "External Python not found.")
			return

		coil_path = self.coil_edit.text().strip()
		if not coil_path:
			QMessageBox.warning(self, "Error", "Coil script path required.")
			return

		self.solve_btn.setEnabled(False)
		self.solve_btn.setText("Solving...")
		self.result_label.setText("Computing...")

		# Save current Cubit model
		tmpdir = tempfile.mkdtemp(prefix="radia_accel_")
		cub5_file = os.path.join(tmpdir, "model.cub5").replace("\\", "/")
		cubit.cmd(f'save cub5 "{cub5_file}" overwrite')

		self._accel_json = os.path.join(tmpdir, "result.json").replace("\\", "/")
		msh_output = os.path.join(tmpdir, "result.msh").replace("\\", "/")
		calc_script = os.path.join(_this_dir, "calc_accel_magnet.py")

		form_text = self.form_combo.currentText()
		if "Omega" in form_text:
			formulation = "omega"
		elif "HCurl" in form_text:
			formulation = "a"
		elif "MMM" in form_text:
			formulation = "mmm"
		else:
			formulation = "msc"
		mat_text = self.mat_combo.currentText()
		is_linear = "Linear" in mat_text
		is_hys = "Hysteresis" in mat_text

		args = [
			calc_script,
			"--coil-script", coil_path,
			"--cub5", cub5_file,
			"--formulation", formulation,
			"--order", str(self.curve_spin.value()),
			"--fes-order", str(self.fes_spin.value()),
			"--max-iter", str(self.maxiter_spin.value()),
			"--relax", self.relax_edit.text().strip() or "0.3",
			] + (["--newton"] if self.newton_check.isChecked() else []) + [
			"--msh-output", msh_output,
			"--output", self._accel_json,
		]

		if is_linear:
			args += ["--material", "linear",
			         "--mu-r", self.mur_edit.text().strip() or "1000"]
		elif is_hys:
			hys = self.hys_edit.text().strip()
			if not hys or not os.path.exists(hys):
				QMessageBox.warning(self, "Error", ".hys file required.")
				self.solve_btn.setEnabled(True)
				self.solve_btn.setText("Solve")
				return
			args += ["--material", "hysteresis",
			         "--hys-file", hys,
			         "--n-steps", str(self.nsteps_spin.value())]
		else:
			bh_sel = self.bh_combo.currentText()
			args += ["--material", "steel"]
			if "Custom" in bh_sel:
				bh = self.bh_edit.text().strip()
				if bh:
					args += ["--bh-file", bh]

		self._process = QProcess(self)
		self._process.readyReadStandardError.connect(self._on_stderr)
		self._process.finished.connect(self._on_finished)
		self._process.start(self._ext_python, args)

	def _on_stderr(self):
		"""Read progress from subprocess stderr."""
		if self._process is None:
			return
		data = self._process.readAllStandardError()
		text = bytes(data).decode("utf-8", errors="replace").strip()
		if text:
			# Show last line of progress
			last_line = text.strip().split("\n")[-1]
			self.result_label.setText(last_line)

	def _on_finished(self, exit_code, exit_status):
		self.solve_btn.setEnabled(True)
		self.solve_btn.setText("Solve")
		self._process = None

		data = None
		json_path = getattr(self, '_accel_json', '')
		if json_path and os.path.exists(json_path):
			try:
				with open(json_path, "r") as f:
					data = json.load(f)
			except Exception:
				pass

		if data is None:
			self.result_label.setText(f"Error (exit code {exit_code})")
			return
		if "error" in data:
			self.result_label.setText(f"Error: {data['error']}")
			return

		B0 = data.get("B_origin_mag", 0)
		L = data.get("L", 0)
		form = data.get("formulation", "?")
		n_iter = data.get("iterations", 0)
		conv = "Yes" if data.get("converged") else "No"
		n_steps = data.get("n_steps", 1)
		hys = "Yes" if data.get("hysteresis") else "No"
		info = (f"|B(0)|={B0:.6f} T, L={L*1e6:.4f} uH\n"
		        f"Formulation={form}, iter={n_iter}, converged={conv}\n"
		        f"DOFs={data.get('ndof', 0)}, t={data.get('t_total', 0):.1f}s")
		if data.get("hysteresis"):
			info += f"\nHysteresis: {n_steps} steps"
			steps = data.get("step_results", [])
			if steps:
				info += " | B(0)=[" + ", ".join(
					f"{s['B_center']:.4f}" for s in steps[-3:]) + "]"
		self.result_label.setText(info)
		self.open_btn.setEnabled(True)

	def _open_result(self):
		"""Open GMSH result file."""
		json_path = getattr(self, '_accel_json', '')
		if not json_path or not os.path.exists(json_path):
			return
		try:
			with open(json_path, "r") as f:
				data = json.load(f)
			msh = data.get("msh_file", "")
			if not msh or not os.path.exists(msh):
				return
			import subprocess as _sp
			ext_py = self._ext_python
			if ext_py:
				ext_pyw = ext_py.replace("python.exe", "pythonw.exe")
				if not os.path.exists(ext_pyw):
					ext_pyw = ext_py
				code = (
					"import sys, gmsh; "
					"gmsh.initialize(sys.argv, run=True); "
					"gmsh.finalize()")
				_sp.Popen([ext_pyw, "-c", code, msh])
			else:
				os.startfile(msh)
		except Exception:
			pass


# ================================================================
# Menu Registration
# ================================================================

def register_menu():
	"""Add 'Radia' submenu to Cubit's Tools menu."""
	main_window = _find_main_window()
	if main_window is None:
		print("WARNING: Could not find Cubit main window. Menu not registered.")
		return

	menu_bar = main_window.menuBar()

	# Find Tools menu
	tools_menu = None
	for action in menu_bar.actions():
		if action.text() == "Tools":
			tools_menu = action.menu()
			break

	if tools_menu is None:
		# Fallback: add to menu bar directly
		tools_menu = menu_bar

	# Remove old menu on re-play (so new class definitions take effect)
	is_reload = False
	for action in list(tools_menu.actions()):
		if action.text() == "Radia-NGSolve":
			sub = action.menu()
			tools_menu.removeAction(action)
			if sub:
				sub.deleteLater()
			# Close stale modeless dialogs
			for attr in ('_radia_ind_dlg', '_radia_fem_dlg', '_radia_accel_dlg'):
				dlg = getattr(main_window, attr, None)
				if dlg is not None:
					dlg.close()
					setattr(main_window, attr, None)
			is_reload = True
			break

	# Create Radia submenu under Tools
	radia_menu = tools_menu.addMenu("Radia-NGSolve")

	# Export Gmsh action
	action_gmsh = QAction("Export Gmsh...", main_window)
	action_gmsh.setStatusTip("Export mesh to Gmsh format (.msh)")
	action_gmsh.triggered.connect(lambda: ExportGmshDialog(main_window).exec())
	radia_menu.addAction(action_gmsh)

	# Volume Calculator action
	action_vol = QAction("Volume Calculator...", main_window)
	action_vol.setStatusTip("Calculate volume of selected volumes (CAD + NGSolve)")
	action_vol.triggered.connect(lambda: VolumeCalculatorDialog(main_window).exec())
	radia_menu.addAction(action_vol)

	# Surface Area Calculator action
	action_surf = QAction("Surface Area...", main_window)
	action_surf.setStatusTip("Calculate surface area (CAD + NGSolve)")
	action_surf.triggered.connect(lambda: SurfaceAreaDialog(main_window).exec())
	radia_menu.addAction(action_surf)

	# IH (BEM): BEM inductance + SIBC heating (surface mesh only)
	action_ih_bem = QAction("IH (BEM)...", main_window)
	action_ih_bem.setStatusTip("Induction heating: BEM inductance + SIBC (surface mesh)")
	def _show_ih_bem():
		if not hasattr(main_window, '_radia_ind_dlg') or main_window._radia_ind_dlg is None:
			main_window._radia_ind_dlg = InductanceDialog(main_window)
		main_window._radia_ind_dlg.show()
		main_window._radia_ind_dlg.raise_()
	action_ih_bem.triggered.connect(_show_ih_bem)
	radia_menu.addAction(action_ih_bem)

	# IH (FEM): FEM-ESIM with Kelvin (volume mesh + SIBC Robin BC)
	action_ih_fem = QAction("IH (FEM)...", main_window)
	action_ih_fem.setStatusTip("Induction heating: FEM + Kelvin + ESIM (volume mesh)")
	def _show_ih_fem():
		if not hasattr(main_window, '_radia_fem_dlg') or main_window._radia_fem_dlg is None:
			main_window._radia_fem_dlg = IHFEMDialog(main_window)
		main_window._radia_fem_dlg.show()
		main_window._radia_fem_dlg.raise_()
	action_ih_fem.triggered.connect(_show_ih_fem)
	radia_menu.addAction(action_ih_fem)

	# Accelerator Magnet: Omega-reduced / A-formulation
	action_accel = QAction("Accelerator Magnet...", main_window)
	action_accel.setStatusTip("Accelerator magnet: Omega-reduced or A-formulation + Kelvin")
	def _show_accel():
		if not hasattr(main_window, '_radia_accel_dlg') or main_window._radia_accel_dlg is None:
			main_window._radia_accel_dlg = AccelMagnetDialog(main_window)
		main_window._radia_accel_dlg.show()
		main_window._radia_accel_dlg.raise_()
	action_accel.triggered.connect(_show_accel)
	radia_menu.addAction(action_accel)

	# Separator + Reload
	radia_menu.addSeparator()
	action_reload = QAction("Reload Panels", main_window)
	action_reload.setStatusTip("Re-read register_toolbar.py from disk (development)")
	def _reload_panels():
		# Use cubit.cmd("play ...") to reload safely outside Qt signal handling.
		# Direct exec() during signal processing crashes Qt (deletes active QMenu).
		startup = os.path.join(os.path.dirname(os.path.abspath(__file__)),
		                       "startup.py").replace("\\", "/")
		cubit.cmd(f'play "{startup}"')
	action_reload.triggered.connect(_reload_panels)
	radia_menu.addAction(action_reload)

	if is_reload:
		print("Radia-NGSolve menu re-registered (code updated).")
	else:
		print("Radia-NGSolve menu registered under Tools.")


# Auto-register when this script is executed
register_menu()
