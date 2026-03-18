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
	from PySide6.QtCore import Qt
except ImportError:
	from PyQt5.QtWidgets import (
		QApplication, QDialog, QVBoxLayout, QHBoxLayout, QGridLayout,
		QLabel, QLineEdit, QComboBox, QPushButton, QSpinBox,
		QFileDialog, QMainWindow, QMessageBox,
		QTableWidget, QTableWidgetItem, QHeaderView, QGroupBox,
		QAction,
	)
	from PyQt5.QtCore import Qt


def _find_main_window():
	"""Find Cubit's main QMainWindow."""
	app = QApplication.instance()
	if app is None:
		return None
	for widget in app.topLevelWidgets():
		if isinstance(widget, QMainWindow):
			return widget
	return None


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


def _get_selected_volume_ids():
	"""Get currently selected volume IDs. Returns all if none selected."""
	selected = cubit.parse_cubit_list("volume", "selected")
	if selected:
		return list(selected)
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
		QApplication.processEvents()

		try:
			# Save current model to temp cub5
			tmpdir = tempfile.mkdtemp(prefix="radia_vol_")
			cub5_file = os.path.join(tmpdir, "model.cub5").replace("\\", "/")
			cubit.cmd(f'save cub5 "{cub5_file}" overwrite')

			# Build command
			calc_script = os.path.join(_this_dir, "calc_volume.py")
			cmd = [self._ext_python, calc_script,
			       "--cub5", cub5_file, "--order", str(order)]
			if step_file:
				cmd += ["--step", step_file]

			# Run external Python
			result = subprocess.run(
				cmd, capture_output=True, text=True, timeout=120
			)

			if result.returncode != 0:
				# Filter out NGSolve/Netgen version warnings from stderr
				stderr = result.stderr.strip()
				stderr_lines = [l for l in stderr.split("\n")
				                if l.strip() and "WARNING" not in l
				                and "=====" not in l and "version" not in l]
				error_msg = "\n".join(stderr_lines) or result.stdout.strip()
				if error_msg:
					QMessageBox.critical(self, "NGSolve Error", error_msg[:1000])
					return

			# Parse JSON from last line of stdout (skip cubit/ngsolve banner)
			stdout_lines = result.stdout.strip().split("\n")
			json_line = stdout_lines[-1]
			data = json.loads(json_line)

			if "error" in data:
				QMessageBox.critical(self, "NGSolve Error", data["error"])
				return

			if "warning" in data:
				QMessageBox.warning(self, "Warning", data["warning"])

			# Update table with NGSolve results
			ng_total = data["ngsolve_total"]
			ng_volumes = data.get("volumes", [])

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

		except subprocess.TimeoutExpired:
			QMessageBox.critical(self, "Timeout", "NGSolve calculation timed out (120s).")
		except json.JSONDecodeError as e:
			QMessageBox.critical(self, "Error", f"Failed to parse NGSolve output:\n{e}")
		except Exception as e:
			QMessageBox.critical(self, "Error", str(e))
		finally:
			self.calc_btn.setEnabled(True)
			self.calc_btn.setText("Calculate")


# ================================================================
# Menu Registration
# ================================================================

def register_menu():
	"""Add 'Radia' menu to Cubit's menu bar."""
	main_window = _find_main_window()
	if main_window is None:
		print("WARNING: Could not find Cubit main window. Menu not registered.")
		return

	menu_bar = main_window.menuBar()

	# Check if already registered (avoid duplicates on re-play)
	for action in menu_bar.actions():
		if action.text() == "Radia":
			print("Radia menu already registered.")
			return

	# Create Radia menu
	radia_menu = menu_bar.addMenu("Radia")

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

	print("Radia menu registered.")


# Auto-register when this script is executed
register_menu()
