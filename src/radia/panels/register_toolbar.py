"""
Register custom toolbar buttons in Coreform Cubit.

This script runs inside Cubit on startup (via ~/.cubit) and adds
toolbar buttons using PySide6 for mesh export functions.
"""

import sys
import os

# Add package root to path
_pkg_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _pkg_root not in sys.path:
	sys.path.insert(0, _pkg_root)

import cubit

from PySide6.QtWidgets import (
	QApplication, QDialog, QVBoxLayout, QHBoxLayout,
	QLabel, QLineEdit, QComboBox, QPushButton,
	QFileDialog, QMainWindow, QToolBar, QMessageBox,
)
from PySide6.QtGui import QAction


def _find_main_window():
	"""Find Cubit's main QMainWindow."""
	app = QApplication.instance()
	if app is None:
		return None
	for widget in app.topLevelWidgets():
		if isinstance(widget, QMainWindow):
			return widget
	return None


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


def register_toolbar():
	"""Add mesh export toolbar to Cubit's main window."""
	main_window = _find_main_window()
	if main_window is None:
		print("WARNING: Could not find Cubit main window. Toolbar not registered.")
		return

	# Create toolbar
	toolbar = QToolBar("Mesh Export")
	toolbar.setObjectName("MeshExportToolbar")
	main_window.addToolBar(toolbar)

	# Export Gmsh action
	action_gmsh = QAction("Export Gmsh", main_window)
	action_gmsh.setToolTip("Export mesh to Gmsh format (.msh)")
	action_gmsh.triggered.connect(lambda: ExportGmshDialog(main_window).exec())
	toolbar.addAction(action_gmsh)

	print("Mesh Export toolbar registered.")


# Auto-register when this script is executed
register_toolbar()
