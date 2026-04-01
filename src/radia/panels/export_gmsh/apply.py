"""
Export Gmsh - Command Panel Apply Script

Exports the current Cubit mesh to Gmsh format (.msh) using the radia plugin command.
"""

import cubit
from cubit import cubitgui

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

	# Map version string to command parameter
	ver = "4" if version == "4.1" else "2"

	# Build command
	cmd = f'radia export gmsh "{file_name}" version {ver}'
	if version == "4.1" and dim in ("2D", "3D"):
		d = "2" if dim == "2D" else "3"
		cmd += f" dimension {d}"
	cmd += " overwrite"

	cubit.cmd(cmd)
	print(f"Done: {file_name}")
