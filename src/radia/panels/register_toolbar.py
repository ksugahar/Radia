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
	)
	from PySide6.QtGui import QAction
	from PySide6.QtCore import Qt, QProcess
except ImportError:
	from PyQt5.QtWidgets import (
		QApplication, QDialog, QVBoxLayout, QHBoxLayout, QGridLayout,
		QLabel, QLineEdit, QComboBox, QPushButton, QSpinBox,
		QFileDialog, QMainWindow, QMessageBox,
		QTableWidget, QTableWidgetItem, QHeaderView, QGroupBox,
		QAction,
	)
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


def _get_selected_volume_ids():
	"""Get currently selected volume IDs. Returns all if none selected."""
	try:
		selected = cubit.parse_cubit_list("volume", "selected")
		if selected:
			return list(selected)
	except Exception:
		pass
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

		ngsolve_layout.addWidget(QLabel("STEP file:"), 1, 0)
		self.step_edit = QLineEdit()
		self.step_edit.setPlaceholderText("(optional, required for order > 1)")
		ngsolve_layout.addWidget(self.step_edit, 1, 1)
		self.step_browse = QPushButton("...")
		self.step_browse.setFixedWidth(30)
		self.step_browse.clicked.connect(self._browse_step)
		ngsolve_layout.addWidget(self.step_browse, 1, 2)

		ngsolve_layout.addWidget(QLabel("Python:"), 2, 0)
		self.python_label = QLabel(self._ext_python or "Not found")
		self.python_label.setStyleSheet(
			"color: green;" if self._ext_python else "color: red;"
		)
		ngsolve_layout.addWidget(self.python_label, 2, 1)

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

	def _browse_step(self):
		path, _ = QFileDialog.getOpenFileName(
			self, "Select STEP File", "", "STEP Files (*.step *.stp);;All Files (*)"
		)
		if path:
			self.step_edit.setText(path)

	def _calculate_cad_volumes(self):
		"""Calculate CAD volumes for selected volumes."""
		vol_ids = _get_selected_volume_ids()
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
		step_file = self.step_edit.text().strip() or None

		self.calc_btn.setEnabled(False)
		self.calc_btn.setText("Calculating...")

		# Save current model to temp cub5
		tmpdir = tempfile.mkdtemp(prefix="radia_vol_")
		cub5_file = os.path.join(tmpdir, "model.cub5").replace("\\", "/")
		cubit.cmd(f'save cub5 "{cub5_file}" overwrite')

		# Build command
		calc_script = os.path.join(_this_dir, "calc_volume.py")
		args = ["--cub5", cub5_file, "--order", str(order)]
		if step_file:
			args += ["--step", step_file]

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
		vol_ids = _get_selected_volume_ids()
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
	"""Dialog for ngsolve.bem inductance extraction."""

	def __init__(self, parent=None):
		super().__init__(parent)
		self.setWindowTitle("Inductance Extractor (ngsolve.bem)")
		self.setMinimumWidth(500)
		self._ext_python = _find_external_python()
		self._result_file = os.path.join(tempfile.gettempdir(), "radia_inductance_result.json")
		self._setup_ui()
		self._populate_blocks()
		self._try_load_existing_result()

	def _setup_ui(self):
		layout = QVBoxLayout(self)

		# Block selection
		grid = QGridLayout()

		grid.addWidget(QLabel("Source block:"), 0, 0)
		self.source_combo = QComboBox()
		self.source_combo.setEditable(True)
		grid.addWidget(self.source_combo, 0, 1)

		grid.addWidget(QLabel("Sink block:"), 1, 0)
		self.sink_combo = QComboBox()
		self.sink_combo.setEditable(True)
		self.sink_combo.addItem("(closed loop)")
		grid.addWidget(self.sink_combo, 1, 1)

		grid.addWidget(QLabel("Conductivity:"), 2, 0)
		self.sigma_edit = QLineEdit("5.8e7")
		grid.addWidget(self.sigma_edit, 2, 1)
		grid.addWidget(QLabel("S/m"), 2, 2)

		grid.addWidget(QLabel("Curve order:"), 3, 0)
		self.order_spin = QSpinBox()
		self.order_spin.setRange(1, 5)
		self.order_spin.setValue(1)
		grid.addWidget(self.order_spin, 3, 1)

		grid.addWidget(QLabel("Frequency:"), 4, 0)
		self.freq_edit = QLineEdit("0")
		grid.addWidget(self.freq_edit, 4, 1)
		grid.addWidget(QLabel("Hz (0 = DC)"), 4, 2)

		layout.addLayout(grid)

		# Result table
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

		# Buttons
		btn_row = QHBoxLayout()
		btn_row.addStretch()
		self.extract_btn = QPushButton("Extract")
		self.extract_btn.clicked.connect(self._extract)
		self.extract_btn.setEnabled(self._ext_python is not None)
		btn_row.addWidget(self.extract_btn)
		self.load_btn = QPushButton("Load Results")
		self.load_btn.clicked.connect(self._load_results)
		self.load_btn.setEnabled(False)
		btn_row.addWidget(self.load_btn)
		close_btn = QPushButton("Close")
		close_btn.clicked.connect(self.accept)
		btn_row.addWidget(close_btn)
		layout.addLayout(btn_row)

	def _set_result(self, param, value):
		"""Set a single-row result in the table."""
		self.result_table.setRowCount(1)
		self.result_table.setItem(0, 0, QTableWidgetItem(param))
		self.result_table.setItem(0, 1, QTableWidgetItem(value))

	def _load_results(self):
		"""Load saved results and set nodal variable for visualization."""
		debug_lines = []
		try:
			# Read result file
			result_file = getattr(self, '_result_file', None)
			if not result_file or not os.path.exists(result_file):
				debug_lines.append("No result file found. Run Extract first.")
				self.debug_text.setText("\n".join(debug_lines))
				return

			with open(result_file, "r") as f:
				data = json.load(f)
			debug_lines.append(f"Loaded: {result_file}")

			node_J = data.get("node_J")
			if not node_J:
				debug_lines.append("ERROR: node_J not in result data")
				self.debug_text.setText("\n".join(debug_lines))
				return

			debug_lines.append(f"node_J: len={len(node_J)}, min={min(node_J):.2e}, max={max(node_J):.2e}")
			debug_lines.append(f"nonzero: {sum(1 for v in node_J if v > 0)}/{len(node_J)}")

			node_ids = list(cubit.get_entities("node"))
			debug_lines.append(f"Cubit nodes: {len(node_ids)}")

			if len(node_J) != len(node_ids):
				debug_lines.append(f"ERROR: mismatch node_J({len(node_J)}) vs cubit({len(node_ids)})")
				self.debug_text.setText("\n".join(debug_lines))
				return

			cubit.set_nodal_variable(node_ids, "J_magnitude", node_J)
			debug_lines.append("set_nodal_variable OK")

			# Export Exodus with nodal variable, then reimport for results display
			exo_file = os.path.join(tempfile.gettempdir(), "inductance_result.exo").replace("\\", "/")
			cubit.cmd(f'export mesh "{exo_file}" overwrite')
			debug_lines.append(f"Exported: {exo_file}")

			# Reimport with nodal variable for contour display
			cubit.cmd("reset")
			cubit.cmd(f'import mesh "{exo_file}" no_geom')
			debug_lines.append("Reimported Exodus with nodal_var")

		except Exception as e:
			debug_lines.append(f"ERROR: {e}")

		self.debug_text.setText("\n".join(debug_lines))

	def _try_load_existing_result(self):
		"""Load existing result JSON if available (from previous Extract)."""
		if os.path.exists(self._result_file):
			try:
				with open(self._result_file, "r") as f:
					data = json.load(f)
				if "inductance_H" in data:
					self._display_result(data)
					self.load_btn.setEnabled(True)
					self.debug_text.setText(f"Previous result loaded: {self._result_file}")
			except Exception:
				pass

	def _populate_blocks(self):
		"""Populate combo boxes with Cubit block names."""
		try:
			block_count = cubit.get_block_count()
			for bid in range(1, block_count + 1):
				try:
					name = cubit.get_exodus_entity_name("block", bid)
					if name:
						self.source_combo.addItem(name)
						self.sink_combo.addItem(name)
				except Exception:
					pass
		except Exception:
			pass

	def _extract(self):
		"""Run inductance extraction via external Python."""
		if not self._ext_python:
			QMessageBox.warning(self, "Error", "External Python not found.")
			return

		source = self.source_combo.currentText().strip()
		sink = self.sink_combo.currentText().strip()
		if sink == "(closed loop)":
			sink = ""

		try:
			sigma = float(self.sigma_edit.text())
		except ValueError:
			QMessageBox.warning(self, "Error", "Invalid conductivity value.")
			return

		try:
			freq = float(self.freq_edit.text())
		except ValueError:
			QMessageBox.warning(self, "Error", "Invalid frequency value.")
			return

		order = self.order_spin.value()

		self.extract_btn.setEnabled(False)
		self.extract_btn.setText("Extracting...")
		self._set_result("Status", "Computing...")

		# Save cub5
		tmpdir = tempfile.mkdtemp(prefix="radia_ind_")
		cub5_file = os.path.join(tmpdir, "model.cub5").replace("\\", "/")
		cubit.cmd(f'save cub5 "{cub5_file}" overwrite')

		# Build command
		calc_script = os.path.join(_this_dir, "calc_inductance.py")
		args = [
			calc_script,
			"--cub5", cub5_file,
			"--source", source,
			"--sink", sink,
			"--sigma", str(sigma),
			"--order", str(order),
			"--freq", str(freq),
		]

		# Run async via QProcess (non-blocking)
		self._process = QProcess(self)
		self._process.finished.connect(self._on_extract_finished)
		self._process.start(self._ext_python, args)

	def _on_extract_finished(self, exit_code, exit_status):
		"""Handle async extraction result."""
		self.extract_btn.setEnabled(True)
		self.extract_btn.setText("Extract")
		data = _parse_json_output(self._process)
		self._process = None
		if data is None:
			self._set_result("Status", "Error")
			return

		if "error" in data:
			QMessageBox.critical(self, "Error", data["error"])
			self._set_result("Status", "Error")
			return

		self._display_result(data)

		# Save result data for Load button
		with open(self._result_file, "w") as f:
			json.dump(data, f)
		self.load_btn.setEnabled(True)
		self.debug_text.setText(f"Result saved: {self._result_file}")

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
			("Inductance", L_str),
			("DOFs", str(data.get("n_free_dofs", ""))),
			("Neg diag", str(data.get("neg_diag", ""))),
			("Surface area", f"{data.get('surface_area', 0):.4e}"),
			("Curve order", str(data.get("order", ""))),
		]

		self.result_table.setRowCount(len(rows))
		for i, (param, val) in enumerate(rows):
			self.result_table.setItem(i, 0, QTableWidgetItem(param))
			item = QTableWidgetItem(val)
			item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
			self.result_table.setItem(i, 1, item)


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

	# Check if already registered (avoid duplicates on re-play)
	for action in tools_menu.actions():
		if action.text() == "Radia-NGSolve":
			print("Radia menu already registered.")
			return

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

	# Inductance Extractor action
	action_ind = QAction("Inductance...", main_window)
	action_ind.setStatusTip("Extract inductance using ngsolve.bem LaplaceSL BEM")
	action_ind.triggered.connect(lambda: InductanceDialog(main_window).exec())
	radia_menu.addAction(action_ind)

	print("Radia-NGSolve menu registered under Tools.")


# Auto-register when this script is executed
register_menu()
