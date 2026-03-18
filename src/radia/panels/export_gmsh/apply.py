"""
Export Gmsh - Command Panel Apply Script

Exports the current Cubit mesh to Gmsh format (.msh) using cubit_mesh_export.
"""

import cubit
from cubit import cubitgui

import sys
import os

# Add repository root to path for cubit_mesh_export
script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(script_dir))
if repo_root not in sys.path:
	sys.path.insert(0, repo_root)

import cubit_mesh_export

# Get values from panel widgets
file_name = cubitgui.UserPanel.line_edit_text("file_name")
version = cubitgui.UserPanel.combo_box_text("version")
dim = cubitgui.UserPanel.combo_box_text("dim")

# Validate file name
if not file_name or file_name.strip() == "":
	print("ERROR: Please specify an output file name.")
else:
	# Add .msh extension if not present
	if not file_name.endswith(".msh"):
		file_name += ".msh"

	print(f"Exporting to Gmsh {version}: {file_name}")

	if version == "4.1":
		cubit_mesh_export.export_Gmesh(cubit, file_name, version=version, DIM=dim)
	else:
		cubit_mesh_export.export_Gmesh(cubit, file_name, version=version)

	print(f"Done: {file_name}")
