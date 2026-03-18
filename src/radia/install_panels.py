"""
Panel installer for Coreform Cubit Mesh Export.

Registers the toolbar script in Cubit's startup file (~/.cubit)
so that custom toolbar buttons are loaded automatically on startup.

The toolbar script is used in-place (no copying):
  - Development (pip install -e .): points to the repo's panels/
  - Distribution (pip install):     points to site-packages/panels/

Usage:
    # From command line (after pip install):
    cubit-mesh-export-install-panels

    # From Python:
    from radia.install_panels import install_panels
    install_panels()
"""

import os
import sys


# Marker comments to identify our block in .cubit file
_MARKER_BEGIN = "## BEGIN cubit_mesh_export toolbar"
_MARKER_END = "## END cubit_mesh_export toolbar"


def _get_cubit_startup_file():
	"""Get the path to the user's .cubit startup file."""
	home = os.path.expanduser("~")
	return os.path.join(home, ".cubit")


def _get_toolbar_script():
	"""Get the absolute path to the toolbar registration script."""
	return os.path.join(
		os.path.dirname(os.path.abspath(__file__)),
		"panels", "register_toolbar.py"
	)


def _build_startup_block(script_path):
	"""Build the Python block to insert into .cubit file."""
	# Use forward slashes for Cubit compatibility
	script_path = script_path.replace("\\", "/")
	return (
		f"{_MARKER_BEGIN}\n"
		f"#{{python\n"
		f"exec(open(r\"{script_path}\").read())\n"
		f"#}}python\n"
		f"{_MARKER_END}\n"
	)


def _remove_existing_block(lines):
	"""Remove existing cubit_mesh_export block from lines."""
	result = []
	inside_block = False
	for line in lines:
		if _MARKER_BEGIN in line:
			inside_block = True
			continue
		if _MARKER_END in line:
			inside_block = False
			continue
		if not inside_block:
			result.append(line)
	return result


def install_panels():
	"""Register custom toolbar for Coreform Cubit.

	Adds a startup block to ~/.cubit that loads the toolbar script.
	The toolbar script location is detected automatically:
	  - pip install -e . (editable): repo's panels/ directory
	  - pip install (normal):        site-packages panels/ directory
	"""
	print("=== Coreform Cubit Mesh Export - Toolbar Installer ===\n")

	# Step 1: Locate toolbar script
	script_path = _get_toolbar_script()
	if not os.path.isfile(script_path):
		print(f"ERROR: Toolbar script not found: {script_path}")
		return False

	print(f"Toolbar script: {script_path}")

	# Step 2: Update ~/.cubit
	cubit_file = _get_cubit_startup_file()
	print(f"Startup file:   {cubit_file}")
	print()

	# Read existing content
	lines = []
	if os.path.isfile(cubit_file):
		with open(cubit_file, "r", encoding="utf-8") as f:
			lines = f.readlines()

	# Remove old block if present
	lines = _remove_existing_block(lines)

	# Append new block
	block = _build_startup_block(script_path)
	lines.append("\n" + block)

	# Write back
	with open(cubit_file, "w", encoding="utf-8") as f:
		f.writelines(lines)

	print("=== Installation Complete ===")
	print(f"Updated: {cubit_file}")
	print("Restart Cubit to load the toolbar.")
	return True


def uninstall_panels():
	"""Remove the custom toolbar registration from ~/.cubit."""
	print("=== Coreform Cubit Mesh Export - Toolbar Uninstaller ===\n")

	cubit_file = _get_cubit_startup_file()
	if not os.path.isfile(cubit_file):
		print("Nothing to uninstall (.cubit file not found).")
		return True

	with open(cubit_file, "r", encoding="utf-8") as f:
		lines = f.readlines()

	new_lines = _remove_existing_block(lines)

	if len(new_lines) == len(lines):
		print("Nothing to uninstall (no toolbar block found).")
		return True

	with open(cubit_file, "w", encoding="utf-8") as f:
		f.writelines(new_lines)

	print(f"Removed toolbar from: {cubit_file}")
	print("Restart Cubit to apply changes.")
	return True


def main():
	"""Console script entry point."""
	if len(sys.argv) > 1 and sys.argv[1] == "--uninstall":
		success = uninstall_panels()
	else:
		success = install_panels()
	sys.exit(0 if success else 1)


if __name__ == "__main__":
	main()
