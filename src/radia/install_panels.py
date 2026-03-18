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


def _get_panels_dir():
	"""Get the absolute path to the panels/ directory."""
	return os.path.join(
		os.path.dirname(os.path.abspath(__file__)),
		"panels"
	)


def _generate_startup_script(panels_dir):
	"""Generate startup.py with paths baked in (avoids __file__ issues in Cubit play)."""
	register_path = os.path.join(panels_dir, "register_toolbar.py").replace("\\", "/")
	startup_path = os.path.join(panels_dir, "startup.py")

	# Find Cubit's site-packages for PyQt5
	cubit_path = os.environ.get("CUBIT_PATH", "")
	if not cubit_path:
		import glob
		base = os.environ.get("ProgramFiles", "C:/Program Files")
		candidates = sorted(glob.glob(os.path.join(base, "Coreform Cubit *", "bin")),
		                    reverse=True)
		cubit_path = candidates[0] if candidates else ""
	cubit_site = os.path.join(cubit_path, "python3", "lib", "site-packages").replace("\\", "/")

	# Write startup.py as a single-line script (Cubit play executes line by line)
	content = (
		f'import sys, os; '
		f'_cs = r"{cubit_site}"; '
		f'sys.path.insert(0, _cs) if _cs not in sys.path else None; '
		f'__file__ = r"{register_path}"; '
		f'exec(open(r"{register_path}").read())\n'
	)

	with open(startup_path, "w", encoding="utf-8") as f:
		f.write(content)

	return startup_path


def _build_startup_block(startup_script_path):
	"""Build the block to insert into .cubit file."""
	startup_script_path = startup_script_path.replace("\\", "/")
	return (
		f"{_MARKER_BEGIN}\n"
		f"play \"{startup_script_path}\"\n"
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

	# Step 1: Locate panels directory and generate startup.py
	panels_dir = _get_panels_dir()
	register_script = os.path.join(panels_dir, "register_toolbar.py")
	if not os.path.isfile(register_script):
		print(f"ERROR: Toolbar script not found: {register_script}")
		return False

	startup_script = _generate_startup_script(panels_dir)
	print(f"Toolbar script: {register_script}")
	print(f"Startup script: {startup_script}")

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
	block = _build_startup_block(startup_script)
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
