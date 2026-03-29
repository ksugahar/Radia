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
		QLabel, QLineEdit, QComboBox, QPushButton, QSpinBox,
		QFileDialog, QMainWindow, QMessageBox,
		QTableWidget, QTableWidgetItem, QHeaderView, QGroupBox,
		QPlainTextEdit,
	)
	from PySide6.QtGui import QAction, QFont
	from PySide6.QtCore import Qt, QProcess
except ImportError:
	from PyQt5.QtWidgets import (
		QApplication, QDialog, QVBoxLayout, QHBoxLayout, QGridLayout,
		QLabel, QLineEdit, QComboBox, QPushButton, QSpinBox,
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

		esim_layout.addWidget(QLabel("Material:"), 1, 0)
		self.material_combo = QComboBox()
		self.material_combo.addItems(["Steel", "Copper", "Aluminum"])
		esim_layout.addWidget(self.material_combo, 1, 1)

		esim_layout.addWidget(QLabel("Frequency [Hz]:"), 2, 0)
		self.freq_edit = QLineEdit("50000")
		esim_layout.addWidget(self.freq_edit, 2, 1)

		esim_layout.addWidget(QLabel("Sigma [S/m]:"), 3, 0)
		self.sigma_edit = QLineEdit("2.0e6")
		self.material_combo.currentTextChanged.connect(self._on_material_changed)
		esim_layout.addWidget(self.sigma_edit, 3, 1)

		# mu_r row (SIBC/Dowell only)
		self.mur_label = QLabel("mu_r:")
		esim_layout.addWidget(self.mur_label, 4, 0)
		self.mur_edit = QLineEdit("100")
		esim_layout.addWidget(self.mur_edit, 4, 1)

		# BH curve row (ESIM only)
		self.bh_label = QLabel("BH curve:")
		esim_layout.addWidget(self.bh_label, 5, 0)
		bh_row = QHBoxLayout()
		self.bh_edit = QLineEdit("(built-in Steel)")
		self.bh_edit.setReadOnly(True)
		self.bh_edit.setToolTip("BH curve file (2-column: H[A/m] B[T])")
		bh_row.addWidget(self.bh_edit)
		self.bh_browse = QPushButton("...")
		self.bh_browse.setFixedWidth(30)
		self.bh_browse.clicked.connect(self._browse_bh)
		bh_row.addWidget(self.bh_browse)
		esim_layout.addLayout(bh_row, 5, 1)

		# Curvature row (ESIM only)
		self.curv_label = QLabel("Curvature:")
		esim_layout.addWidget(self.curv_label, 6, 0)
		self.curvature_combo = QComboBox()
		self.curvature_combo.addItems(["Local curvature", "None (flat)"])
		self.curvature_combo.setToolTip(
			"Local curvature: cylindrical cell problem (Bessel I0/I1),\n"
			"None (flat): planar slab cell problem (cosh/sinh).")
		esim_layout.addWidget(self.curvature_combo, 6, 1)

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
		"""Default FEM journal: Periodic Kelvin (2 spheres) + coil + wp hole.

		Interior sphere at origin: air + coil + workpiece hole
		Exterior sphere at offset: Kelvin-mapped domain (same radius)
		Periodic identification via matching surface labels.
		"""
		lines = [
			"# IH (FEM): Periodic Kelvin (2-sphere) + coil + workpiece",
			"# Interior: R=60mm sphere at origin (air + coil + workpiece)",
			"# Exterior: R=60mm sphere at x=150mm (Kelvin domain)",
			"# Kelvin weight: nu = nu0 * (a/r')^2 (conformal symmetry)",
			"# Coil: R=30mm, a=3mm torus (355deg gap)",
			"# Workpiece: R=10mm, H=20mm cylinder (SIBC on interface)",
			"reset",
			"",
			"# === Workpiece cylinder ===",
			"create cylinder height 0.020 radius 0.010",
			"# -> workpiece (volume 1)",
			"",
			"# === Coil torus ===",
			"create torus major radius 0.030 minor radius 0.003",
			"# -> coil (volume 2)",
			"",
			"# === Air sphere (subtract coil and wp, then recreate them) ===",
			"create sphere radius 0.060",
			"create cylinder height 0.020 radius 0.010",
			"create torus major radius 0.030 minor radius 0.003",
			"subtract volume 4 5 from volume 3",
			"# -> air (volume 6)",
			"",
			"# === Exterior domain (Kelvin, offset x=150mm) ===",
			"create sphere radius 0.060",
			"move volume 7 x 0.150 include_merged",
			"# -> kelvin (volume 7)",
			"",
			"# === Imprint and merge (conformal interfaces) ===",
			"imprint volume all",
			"merge volume all",
			"",
			"# === Mesh ===",
			"volume all scheme tetmesh",
			"volume 1 size 0.003",
			"volume 2 size 0.003",
			"volume 6 size 0.008",
			"volume 7 size 0.015",
			"mesh volume all",
			"",
			"# === Blocks ===",
			"set duplicate block elements on",
			"block 1 add volume 1",
			'block 1 name "workpiece"',
			"block 2 add volume 2",
			'block 2 name "coil"',
			"block 3 add volume 6",
			'block 3 name "air"',
			"block 4 add volume 7",
			'block 4 name "kelvin"',
			"",
			"# wp_surface and outer blocks are auto-created by _detect_blocks()",
			"",
			"# === Hide non-essential volumes ===",
			"volume 6 visibility off",
			"volume 7 visibility off",
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

	def _on_material_changed(self, text):
		"""Update sigma, mu_r, and BH curve when material combo changes."""
		sigma_map = {"Steel": "2.0e6", "Copper": "5.8e7", "Aluminum": "3.5e7"}
		mur_map = {"Steel": "100", "Copper": "1", "Aluminum": "1"}
		self.sigma_edit.setText(sigma_map.get(text, "2.0e6"))
		self.mur_edit.setText(mur_map.get(text, "1"))
		if text == "Steel":
			self.bh_edit.setText("(built-in Steel)")
		else:
			self.bh_edit.setText("(linear, mu_r=1)")

	def _browse_bh(self):
		"""Browse for BH curve file (2-column: H[A/m] B[T])."""
		path, _ = QFileDialog.getOpenFileName(
			self, "Load BH Curve", "",
			"Text files (*.txt *.csv *.dat);;All Files (*)")
		if path:
			self.bh_edit.setText(path)
			self.material_combo.setCurrentText("Steel")

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
				"--material", self.material_combo.currentText().lower(),
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
		material = self.material_combo.currentText().lower()

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
		self.setWindowTitle("IH (FEM): Kelvin + ESIM")
		self.setMinimumWidth(550)
		self.setMinimumHeight(500)
		self._ext_python = _find_external_python()
		self._process = None
		self._setup_ui()

	def _setup_ui(self):
		layout = QVBoxLayout(self)

		# --- Journal editor ---
		jou_label = QLabel("Cubit Journal (FEM: coil + air + Kelvin):")
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
			("workpiece:", "_workpiece_block"),
			("coil:", "_coil_block"),
			("air:", "_air_block"),
			("kelvin:", "_kelvin_block"),
			("wp_surface:", "_wp_surface"),
			("outer:", "_outer_surface"),
		]):
			block_group.addWidget(QLabel(lbl), row, 0)
			w = QLabel("(not found)")
			w.setStyleSheet("font-weight: bold; color: red;")
			setattr(self, f"label{attr}", w)
			block_group.addWidget(w, row, 1)
		layout.addLayout(block_group)

		# --- Workpiece settings ---
		wp_group = QGroupBox("Workpiece (SIBC)")
		wp_layout = QGridLayout(wp_group)

		wp_layout.addWidget(QLabel("Impedance:"), 0, 0)
		self.model_combo = QComboBox()
		self.model_combo.addItems(["SIBC", "ESIM"])
		wp_layout.addWidget(self.model_combo, 0, 1)

		wp_layout.addWidget(QLabel("Material:"), 1, 0)
		self.material_combo = QComboBox()
		self.material_combo.addItems(["Steel", "Copper", "Aluminum"])
		wp_layout.addWidget(self.material_combo, 1, 1)

		wp_layout.addWidget(QLabel("Frequency [Hz]:"), 2, 0)
		self.freq_edit = QLineEdit("7000")
		wp_layout.addWidget(self.freq_edit, 2, 1)

		wp_layout.addWidget(QLabel("Sigma [S/m]:"), 3, 0)
		self.sigma_edit = QLineEdit("2.0e6")
		self.material_combo.currentTextChanged.connect(self._on_material_changed)
		wp_layout.addWidget(self.sigma_edit, 3, 1)

		# mu_r (SIBC only, linear)
		self.mur_label = QLabel("mu_r:")
		wp_layout.addWidget(self.mur_label, 4, 0)
		self.mur_edit = QLineEdit("100")
		wp_layout.addWidget(self.mur_edit, 4, 1)

		# BH curve (ESIM only, nonlinear)
		self.bh_label = QLabel("BH curve:")
		wp_layout.addWidget(self.bh_label, 5, 0)
		bh_row = QHBoxLayout()
		self.bh_edit = QLineEdit("(built-in Steel)")
		self.bh_edit.setReadOnly(True)
		bh_row.addWidget(self.bh_edit)
		self.bh_browse = QPushButton("...")
		self.bh_browse.setFixedWidth(30)
		self.bh_browse.clicked.connect(self._browse_bh)
		bh_row.addWidget(self.bh_browse)
		wp_layout.addLayout(bh_row, 5, 1)

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
		self.solve_btn.setToolTip("FEM-ESIM: Kelvin + HCurl + SIBC Karl iteration")
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

	def _on_material_changed(self, text):
		sigma_map = {"Steel": "2.0e6", "Copper": "5.8e7", "Aluminum": "3.5e7"}
		mur_map = {"Steel": "100", "Copper": "1", "Aluminum": "1"}
		self.sigma_edit.setText(sigma_map.get(text, "2.0e6"))
		self.mur_edit.setText(mur_map.get(text, "1"))
		if text == "Steel":
			self.bh_edit.setText("(built-in Steel)")
		else:
			self.bh_edit.setText("(linear, mu_r=1)")

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

		if "outer" not in found and "kelvin" in found:
			self._auto_create_outer(found)
			for bid in cubit.get_block_id_list():
				name = cubit.get_exodus_entity_name("block", bid)
				found[name] = bid

		for key, attr in [("workpiece", "_workpiece_block"),
		                   ("coil", "_coil_block"), ("air", "_air_block"),
		                   ("kelvin", "_kelvin_block"),
		                   ("wp_surface", "_wp_surface"),
		                   ("outer", "_outer_surface")]:
			label_w = getattr(self, f"label{attr}")
			if key in found:
				label_w.setText(key)
				label_w.setStyleSheet("font-weight: bold; color: green;")
			else:
				label_w.setText("(not found)")
				label_w.setStyleSheet("font-weight: bold; color: red;")

		has_all = all(k in found for k in ("coil", "air", "wp_surface"))
		self.solve_btn.setEnabled(has_all and self._ext_python is not None)

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

	def _solve(self):
		"""Run FEM-ESIM via external Python (calc_heating.py)."""
		if not self._ext_python:
			QMessageBox.warning(self, "Error", "External Python not found.")
			return

		self.solve_btn.setEnabled(False)
		self.solve_btn.setText("Solving...")
		self.result_label.setText("Computing FEM-ESIM (Kelvin + SIBC)...")

		# Save current Cubit model
		tmpdir = tempfile.mkdtemp(prefix="radia_fem_")
		cub5_file = os.path.join(tmpdir, "model.cub5").replace("\\", "/")
		cubit.cmd(f'save cub5 "{cub5_file}" overwrite')

		freq = self.freq_edit.text().strip() or "7000"
		sigma = self.sigma_edit.text().strip() or "2e6"
		material = self.material_combo.currentText().lower()

		self._heat_json = os.path.join(tmpdir, "result.json").replace("\\", "/")
		calc_script = os.path.join(_this_dir, "calc_heating.py")

		# Get geometry from Cubit bounding boxes
		try:
			r_coil, a_coil, r_wp, h_wp = 0.03, 0.003, 0.01, 0.02
			for bid in cubit.get_block_id_list():
				name = cubit.get_exodus_entity_name("block", bid)
				if name == "coil":
					vids = list(cubit.get_block_volumes(bid))
					if vids:
						bb = cubit.get_bounding_box("volume", vids[0])
						r_coil = (abs(bb[1]) + abs(bb[0])) / 2
						a_coil = min(abs(bb[1] - bb[0]),
						             abs(bb[5] - bb[4])) / 4
		except Exception:
			pass

		args = [
			calc_script,
			"--r-coil", str(round(r_coil, 6)),
			"--a-coil", str(round(a_coil, 6)),
			"--r-wp", str(round(r_wp, 6)),
			"--h-wp", str(round(h_wp, 6)),
			"--frequency", freq,
			"--sigma", sigma,
			"--material", material,
			"--output", self._heat_json,
		]

		self._process = QProcess(self)
		self._process.finished.connect(self._on_finished)
		self._process.start(self._ext_python, args)

	def _on_finished(self, exit_code, exit_status):
		self.solve_btn.setEnabled(True)
		self.solve_btn.setText("Solve")
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
			self.result_label.setText(f"Error (exit code {exit_code})")
			return
		if "error" in data:
			self.result_label.setText(f"Error: {data['error']}")
			return

		P = data.get("P_total", 0)
		L = data.get("L_coil", 0)
		self.result_label.setText(
			f"P = {P:.4e} W, L = {L*1e9:.2f} nH\n"
			f"H_t = {data.get('H_t_rms', 0):.2f} A/m, "
			f"|Z_s| = {abs(complex(data.get('Z_s', '0'))):.4e}")
		self.open_btn.setEnabled(True)

	def _open_result(self):
		"""Open GMSH result file."""
		json_path = getattr(self, '_heat_json', '')
		if json_path and os.path.exists(json_path):
			try:
				with open(json_path, "r") as f:
					data = json.load(f)
				msh = data.get("msh_file", "")
				if msh and os.path.exists(msh):
					import subprocess
					subprocess.Popen(["gmsh", msh])
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
			for attr in ('_radia_ind_dlg', '_radia_fem_dlg'):
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
