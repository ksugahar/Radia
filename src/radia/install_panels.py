"""
Panel installer for Coreform Cubit.

Registers the toolbar script in Cubit's startup file (~/.cubit)
so that custom panels are loaded automatically on startup.

The toolbar script is used in-place (no copying):
  - Development (pip install -e .): points to the repo's panels/
  - Distribution (pip install):     points to site-packages/panels/

Usage:
    # From command line (after pip install):
    cubit-install-panels                    # current user + Default profile
    cubit-install-panels --all-users        # all existing users (admin)
    cubit-install-panels --uninstall        # remove from current user
    cubit-install-panels --uninstall --all-users  # remove from all users

    # From Python:
    from radia.install_panels import install_panels
    install_panels()
    install_panels(all_users=True)
"""

import glob
import os
import sys


# Marker comments to identify our block in .cubit file
_MARKER_BEGIN = "## BEGIN radia toolbar"
_MARKER_END = "## END radia toolbar"

# Legacy markers (from older cubit_mesh_export package) — cleaned up on install
_LEGACY_MARKERS = [
	("## BEGIN cubit_mesh_export toolbar", "## END cubit_mesh_export toolbar"),
]


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


def _get_cubit_startup_files(all_users=False):
	"""Get paths to all .cubit startup files to install.

	Args:
	  all_users: If True, include all existing user profiles (Windows).

	Returns list of paths:
	  - Current user: ~/.cubit
	  - Default profile: C:\\Users\\Default\\.cubit (all future users, Windows)
	  - --all-users: C:\\Users\\*\\.cubit (all existing users, Windows)
	"""
	paths = [os.path.join(os.path.expanduser("~"), ".cubit")]
	if sys.platform == "win32":
		users_dir = os.path.join(os.environ.get("SystemDrive", "C:"),
		                         os.sep, "Users")
		# Default profile (future users)
		default = os.path.join(users_dir, "Default", ".cubit")
		if os.path.isdir(os.path.dirname(default)):
			paths.append(default)
		# All existing user profiles
		if all_users:
			skip = {"Default", "Public", "All Users", "Default User"}
			for entry in os.listdir(users_dir):
				if entry in skip:
					continue
				user_dir = os.path.join(users_dir, entry)
				if not os.path.isdir(user_dir):
					continue
				cubit_file = os.path.join(user_dir, ".cubit")
				if cubit_file not in paths:
					paths.append(cubit_file)
	return paths


def _get_panels_dir():
	"""Get the absolute path to the panels/ directory."""
	return os.path.join(
		os.path.dirname(os.path.abspath(__file__)),
		"panels"
	)


def _generate_startup_script(panels_dir):
	"""Generate startup.py with paths baked in (avoids __file__ issues in Cubit play).

	The register_toolbar.py path is baked in at install time.
	Cubit's site-packages path is detected dynamically at runtime
	(version-independent: works across Cubit upgrades without reinstalling).
	"""
	register_path = os.path.join(panels_dir, "register_toolbar.py").replace("\\", "/")
	startup_path = os.path.join(panels_dir, "startup.py")

	# Write startup.py as a single-line script (Cubit play executes line by line)
	# Uses try/except to silently skip if Qt is not available
	# (e.g., when external Python opens a .cub5 and triggers .cubit replay)
	#
	# site-packages detection: find python*/lib/site-packages relative to
	# cubit.py's directory.  This avoids baking in version-specific paths
	# like ".../Cubit 2025.3/bin/python3/lib/site-packages".
	content = (
		f'#!python\n'
		f'import sys, os, glob; '
		f'_cb = os.path.dirname(os.path.abspath(os.path.join(os.path.dirname(sys.executable), "cubit.py"))) if not hasattr(sys, "_cubit_bin") else sys._cubit_bin; '
		f'_sp = glob.glob(os.path.join(_cb, "python*", "lib", "site-packages")) + glob.glob(os.path.join(_cb, "python*", "lib", "python*", "site-packages")); '
		f'sys.path.insert(0, _sp[0]) if _sp and _sp[0] not in sys.path else None; '
		f'__file__ = r"{register_path}"; '
		f"exec(\"try:\\n"
		f" exec(open(r'{register_path}').read())\\n"
		f"except Exception as e:\\n"
		f" import traceback; traceback.print_exc()\")\n"
	)

	with open(startup_path, "w", encoding="utf-8") as f:
		f.write(content)

	return startup_path


def _build_startup_block(startup_script_path, cubit_bin=None):
	"""Build the block to insert into .cubit file."""
	startup_script_path = startup_script_path.replace("\\", "/")
	return (
		f"{_MARKER_BEGIN}\n"
		f"set journal off\n"
		f"play \"{startup_script_path}\"\n"
		f"{_MARKER_END}\n"
	)


def _create_launcher(cubit_bin):
	"""Create a launcher batch file for Cubit with radia plugin support.

	Sets CUBIT_PLUGIN_DIR (for Python API) and passes
	-commandplugindir (for GUI) so export commands work in both modes.
	"""
	if sys.platform != "win32":
		return None

	cubit_exe = os.path.join(cubit_bin, "coreform_cubit.exe")
	if not os.path.isfile(cubit_exe):
		return None

	plugin_dir = os.path.join(cubit_bin, "plugins")
	launcher_path = os.path.join(cubit_bin, "cubit_radia.bat")
	# 'start' treats the first quoted string as window title, making it
	# impossible to quote paths with spaces as the executable. We use
	# a temporary variable to hold the plugin dir (may contain spaces).
	content = (
		f'@set CUBIT_PLUGIN_DIR={plugin_dir}\n'
		f'@set _P={plugin_dir}\n'
		f'@start "" "{cubit_exe}" -commandplugindir "%_P%" %*\n'
	)

	try:
		with open(launcher_path, "w", encoding="utf-8") as f:
			f.write(content)
		return launcher_path
	except PermissionError:
		return None


def _remove_existing_block(lines):
	"""Remove existing toolbar block (current and legacy markers) from lines."""
	# Collect all begin/end marker pairs to remove
	markers = [(_MARKER_BEGIN, _MARKER_END)] + _LEGACY_MARKERS
	result = []
	inside_block = False
	for line in lines:
		if any(begin in line for begin, _ in markers):
			inside_block = True
			continue
		if any(end in line for _, end in markers):
			inside_block = False
			continue
		if not inside_block:
			result.append(line)
	return result


def _get_cubit_ini_paths(all_users=False):
	"""Get paths to all Cubit.ini files (Qt settings).

	Returns list of paths:
	  - Current user: %APPDATA%/Coreform/Cubit.ini
	  - --all-users: all existing user profiles
	"""
	paths = []
	if sys.platform == "win32":
		appdata = os.environ.get("APPDATA", "")
		if appdata:
			paths.append(os.path.join(appdata, "Coreform", "Cubit.ini"))
		if all_users:
			users_dir = os.path.join(os.environ.get("SystemDrive", "C:"),
			                         os.sep, "Users")
			skip = {"Default", "Public", "All Users", "Default User"}
			for entry in os.listdir(users_dir):
				if entry in skip:
					continue
				user_dir = os.path.join(users_dir, entry)
				ini = os.path.join(user_dir, "AppData", "Roaming",
				                   "Coreform", "Cubit.ini")
				if ini not in paths and os.path.isfile(ini):
					paths.append(ini)
	return paths


def _ensure_plugin_path_in_ini(ini_path, plugin_dir):
	"""Ensure plugin_dir is registered in Cubit.ini plugin\\paths.

	Cubit only loads third-party .ccm plugins from directories listed
	in plugin\\paths (Cubit.ini [clarofw] section).  Without this,
	the plugin DLL is found but cubit_plugin_instance() is never called.
	"""
	plugin_dir = plugin_dir.replace("\\", "/")

	if not os.path.isfile(ini_path):
		# Create minimal Cubit.ini with plugin path
		os.makedirs(os.path.dirname(ini_path), exist_ok=True)
		with open(ini_path, "w", encoding="utf-8") as f:
			f.write("[clarofw]\n")
			f.write(f"plugin\\paths={plugin_dir}\n")
		return True

	with open(ini_path, "r", encoding="utf-8") as f:
		content = f.read()
	lines = content.splitlines(keepends=True)

	found = False
	updated = False
	for i, line in enumerate(lines):
		if line.startswith("plugin\\paths="):
			found = True
			current = line.split("=", 1)[1].strip()
			# Check if our path is already there
			existing = [p.strip() for p in current.split(",") if p.strip()]
			normalized = [p.replace("\\", "/") for p in existing]
			if plugin_dir not in normalized:
				if current:
					lines[i] = f"plugin\\paths={current}, {plugin_dir}\n"
				else:
					lines[i] = f"plugin\\paths={plugin_dir}\n"
				updated = True
			break

	if not found:
		# Add under [clarofw] section
		for i, line in enumerate(lines):
			if line.strip() == "[clarofw]":
				lines.insert(i + 1, f"plugin\\paths={plugin_dir}\n")
				updated = True
				break
		else:
			# No [clarofw] section, append one
			lines.append("\n[clarofw]\n")
			lines.append(f"plugin\\paths={plugin_dir}\n")
			updated = True

	if updated:
		with open(ini_path, "w", encoding="utf-8") as f:
			f.writelines(lines)
	return updated


def install_panels(all_users=False):
	"""Register custom toolbar for Coreform Cubit.

	Adds a startup block to ~/.cubit that loads the toolbar script.
	The toolbar script location is detected automatically:
	  - pip install -e . (editable): repo's panels/ directory
	  - pip install (normal):        site-packages panels/ directory

	Also sets CUBIT_PLUGIN_DIR so the radia_cubit.ccm plugin commands
	are available in journal files.

	Args:
	  all_users: If True, install to all existing user profiles (admin).
	"""
	print("=== Coreform Cubit - Panel Installer ===\n")

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

	if all_users:
		print("Mode:           --all-users (all existing profiles)")

	# Step 2: Update .cubit startup files
	cubit_files = _get_cubit_startup_files(all_users=all_users)
	block = _build_startup_block(startup_script)

	for cubit_file in cubit_files:
		lines = []
		if os.path.isfile(cubit_file):
			with open(cubit_file, "r", encoding="utf-8") as f:
				lines = f.readlines()
		lines = _remove_existing_block(lines)
		lines.append("\n" + block)
		try:
			with open(cubit_file, "w", encoding="utf-8") as f:
				f.writelines(lines)
			print(f"Updated: {cubit_file}")
		except PermissionError:
			print(f"SKIP (no permission): {cubit_file}")

	# Step 3: Register plugin directory in Cubit.ini + CUBIT_PLUGIN_DIR
	if cubit_bin:
		plugin_dir = os.path.join(cubit_bin, "plugins")

		# 3a. Cubit.ini (GUI plugin dialog)
		ini_paths = _get_cubit_ini_paths(all_users=all_users)
		for ini_path in ini_paths:
			try:
				if _ensure_plugin_path_in_ini(ini_path, plugin_dir):
					print(f"Plugin path registered: {ini_path}")
				else:
					print(f"Plugin path OK: {ini_path}")
			except PermissionError:
				print(f"SKIP (no permission): {ini_path}")

		# 3b. Launcher batch file (sets CUBIT_PLUGIN_DIR before Cubit starts)
		launcher = _create_launcher(cubit_bin)
		if launcher:
			print(f"Launcher created: {launcher}")
			print(f"  -> Start Cubit via this launcher for plugin support")

	print()
	print("=== Installation Complete ===")
	print("Restart Cubit to load the toolbar.")
	return True


def uninstall_panels(all_users=False):
	"""Remove the custom toolbar registration from .cubit files.

	Args:
	  all_users: If True, uninstall from all existing user profiles (admin).
	"""
	print("=== Coreform Cubit - Panel Uninstaller ===\n")

	for cubit_file in _get_cubit_startup_files(all_users=all_users):
		if not os.path.isfile(cubit_file):
			continue
		with open(cubit_file, "r", encoding="utf-8") as f:
			lines = f.readlines()
		new_lines = _remove_existing_block(lines)
		if len(new_lines) < len(lines):
			try:
				with open(cubit_file, "w", encoding="utf-8") as f:
					f.writelines(new_lines)
				print(f"Removed toolbar from: {cubit_file}")
			except PermissionError:
				print(f"SKIP (no permission): {cubit_file}")

	print("Restart Cubit to apply changes.")
	return True


def main():
	"""Console script entry point."""
	all_users = "--all-users" in sys.argv
	if "--uninstall" in sys.argv:
		success = uninstall_panels(all_users=all_users)
	else:
		success = install_panels(all_users=all_users)
	sys.exit(0 if success else 1)


if __name__ == "__main__":
	main()
