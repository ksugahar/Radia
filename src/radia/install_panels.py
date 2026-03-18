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

import glob
import os
import sys


# Marker comments to identify our block in .cubit file
_MARKER_BEGIN = "## BEGIN cubit_mesh_export toolbar"
_MARKER_END = "## END cubit_mesh_export toolbar"


def find_cubit_bin():
	"""Find Cubit bin directory (cross-platform).

	Search order:
	  1. CUBIT_PATH environment variable
	  2. Platform-specific common install locations (newest version first)

	Returns:
	  Path to Cubit bin/ directory, or None if not found.
	"""
	# 1. Explicit env var (highest priority, all platforms)
	cubit_path = os.environ.get("CUBIT_PATH")
	if cubit_path and os.path.isdir(cubit_path):
		return cubit_path

	# 2. Platform-specific search
	search_patterns = []
	if sys.platform == "win32":
		for base in [os.environ.get("ProgramFiles", ""),
		             os.environ.get("ProgramFiles(x86)", "")]:
			if base:
				search_patterns.append(os.path.join(base, "Coreform Cubit *", "bin"))
	elif sys.platform == "darwin":
		search_patterns += [
			"/Applications/Coreform-Cubit-*/Coreform Cubit.app/Contents/MacOS",
			"/Applications/Coreform Cubit */bin",
		]
	else:  # Linux
		search_patterns += [
			"/opt/Coreform-Cubit-*/bin",
			"/opt/coreform/cubit-*/bin",
			"/usr/local/Coreform-Cubit-*/bin",
		]

	for pattern in search_patterns:
		candidates = sorted(glob.glob(pattern), reverse=True)  # newest first
		for c in candidates:
			if os.path.isfile(os.path.join(c, "cubit.py")):
				return c

	return None


def find_cubit_site_packages(cubit_bin=None):
	"""Find Cubit's bundled Python site-packages directory.

	Args:
	  cubit_bin: Cubit bin/ path (auto-detected if None)

	Returns:
	  Path to site-packages, or None if not found.
	"""
	if cubit_bin is None:
		cubit_bin = find_cubit_bin()
	if not cubit_bin:
		return None

	# Search for site-packages under Cubit's bundled Python
	# Structure varies: python3/lib/site-packages (Windows/Linux)
	#                   python3/lib/python3.X/site-packages (some Linux)
	candidates = glob.glob(os.path.join(cubit_bin, "python*", "lib", "site-packages"))
	candidates += glob.glob(os.path.join(cubit_bin, "python*", "lib", "python*", "site-packages"))
	if candidates:
		return candidates[0]

	return None


def _get_cubit_startup_file():
	"""Get the path to the user's .cubit startup file."""
	return os.path.join(os.path.expanduser("~"), ".cubit")


def _get_panels_dir():
	"""Get the absolute path to the panels/ directory."""
	return os.path.join(
		os.path.dirname(os.path.abspath(__file__)),
		"panels"
	)


def _generate_startup_script(panels_dir):
	"""Generate startup.py with paths baked in (avoids __file__ issues in Cubit play).

	Paths are determined at install time, not hardcoded in source.
	"""
	register_path = os.path.join(panels_dir, "register_toolbar.py").replace("\\", "/")
	startup_path = os.path.join(panels_dir, "startup.py")

	# Find Cubit's site-packages for PyQt5/PySide6
	cubit_site = find_cubit_site_packages()
	if cubit_site:
		cubit_site = cubit_site.replace("\\", "/")
		site_line = f'_cs = r"{cubit_site}"; sys.path.insert(0, _cs) if _cs not in sys.path else None; '
	else:
		site_line = ""

	# Write startup.py as a single-line script (Cubit play executes line by line)
	content = (
		f'import sys, os; '
		f'{site_line}'
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

	cubit_bin = find_cubit_bin()
	print(f"Cubit bin:      {cubit_bin or 'not found (set CUBIT_PATH)'}")

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
