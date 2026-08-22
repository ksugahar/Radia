"""
Panel installer for Coreform Cubit 2025.12+.

Registers the Radia toolbar script in Cubit's startup file so that the
Radia-NGSolve menu is loaded automatically on Cubit startup.

The startup shim is generated outside the Python package:
  - current-user install: %LOCALAPPDATA%/Radia/Cubit/radia_startup.py
  - all-users install:   %ProgramData%/Radia/Cubit/radia_startup.py

This avoids mutating ``site-packages/radia/panels/startup.py`` or the
editable source tree during install, while still baking in the absolute
``register_toolbar.py`` path Cubit needs.
"""

from __future__ import annotations

import glob
import os
import re
import shutil
import sys
import tarfile
import tempfile
from pathlib import Path


_MARKER_BEGIN = "## BEGIN radia toolbar"
_MARKER_END = "## END radia toolbar"

_LEGACY_MARKERS = [
    ("## BEGIN cubit_mesh_export toolbar", "## END cubit_mesh_export toolbar"),
]

_MIN_CUBIT_VERSION = (2025, 12)
_MIN_CUBIT_VERSION_TEXT = "2025.12"


def _parse_cubit_version(path: str | os.PathLike[str]) -> tuple[int, ...]:
    """Return the Cubit version tuple encoded in an install path."""
    p = Path(path)
    parts = [p.name]
    if p.name.lower() == "bin":
        parts.append(p.parent.name)
    text = " ".join(parts)
    match = re.search(r"(\d{4})[.\-_ ](\d+)", text)
    if not match:
        return (0,)
    return tuple(int(g) for g in match.groups())


def _version_key(path: str | os.PathLike[str]) -> tuple[int, ...]:
    version = _parse_cubit_version(path)
    return version if version else (0,)


def _is_supported_cubit_bin(path: str | os.PathLike[str]) -> bool:
    cubit_bin = Path(path)
    return (
        (cubit_bin / "cubit.py").is_file()
        and _parse_cubit_version(cubit_bin) >= _MIN_CUBIT_VERSION
    )


def find_cubit_bin():
    """Find the Coreform Cubit 2025.12+ bin directory.

    Search order:
      1. CUBIT_PATH environment variable
      2. Platform-specific common install locations, newest version first

    Returns:
      Path to Cubit ``bin/`` directory, or None if no supported Cubit is found.
    """
    cubit_path = os.environ.get("CUBIT_PATH")
    if cubit_path:
        candidates = [Path(cubit_path)]
        if (Path(cubit_path) / "bin").is_dir():
            candidates.insert(0, Path(cubit_path) / "bin")
        for candidate in candidates:
            if _is_supported_cubit_bin(candidate):
                return str(candidate)
        return None

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
    else:
        search_patterns += [
            "/opt/Coreform-Cubit-*/bin",
            "/opt/coreform/cubit-*/bin",
            "/usr/local/Coreform-Cubit-*/bin",
        ]

    for pattern in search_patterns:
        candidates = sorted(glob.glob(pattern), key=_version_key, reverse=True)
        for candidate in candidates:
            if _is_supported_cubit_bin(candidate):
                return candidate
    return None


def find_cubit_site_packages(cubit_bin=None):
    """Find Cubit's bundled Python site-packages directory."""
    if cubit_bin is None:
        cubit_bin = find_cubit_bin()
    if not cubit_bin:
        return None

    candidates = glob.glob(os.path.join(cubit_bin, "python*", "lib", "site-packages"))
    candidates += glob.glob(
        os.path.join(cubit_bin, "python*", "lib", "python*", "site-packages")
    )
    return candidates[0] if candidates else None


def _windows_users_dir() -> str:
    return os.path.join(os.environ.get("SystemDrive", "C:") + os.sep, "Users")


def _iter_windows_user_dirs(include_default: bool):
    users_dir = _windows_users_dir()
    if not os.path.isdir(users_dir):
        return

    skip = {"All Users", "Default User", "Public"}
    if include_default:
        default_dir = os.path.join(users_dir, "Default")
        if os.path.isdir(default_dir):
            yield default_dir

    for entry in sorted(os.listdir(users_dir)):
        if entry in skip or (entry == "Default" and include_default):
            continue
        user_dir = os.path.join(users_dir, entry)
        if os.path.isdir(user_dir):
            yield user_dir


def _get_cubit_startup_files(all_users=False):
    """Return .cubit startup files targeted by this install."""
    paths = [os.path.join(os.path.expanduser("~"), ".cubit")]

    if sys.platform == "win32" and all_users:
        for user_dir in _iter_windows_user_dirs(include_default=True):
            cubit_file = os.path.join(user_dir, ".cubit")
            if cubit_file not in paths:
                paths.append(cubit_file)
    return paths


def _get_panels_dir():
    """Get the absolute path to the panels/ directory."""
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "panels")


def _startup_dir(all_users=False):
    if sys.platform == "win32":
        if all_users:
            base = os.environ.get("ProgramData", r"C:\ProgramData")
        else:
            base = os.environ.get("LOCALAPPDATA") or os.path.expanduser("~")
        return os.path.join(base, "Radia", "Cubit")
    return os.path.join(os.path.expanduser("~"), ".radia", "cubit")


def build_official_toolbar_package(output_dir=None):
    """Build Coreform's supported ``WorkflowToolbar`` import package.

    The archive layout and ``.mappings`` contract match Coreform's official
    ``tire_cross_section_tool`` and ``cubit-dagmc-toolbar`` examples.  The
    export implementation is copied into the package so a Cubit toolbar never
    follows a stale editable-install or release-worktree path.
    """
    panels_dir = Path(_get_panels_dir())
    source = panels_dir / "cubit_toolbar"
    template = source / "toolbars" / "radia_export_toolbar.ttb.tmpl"
    menu_source = panels_dir / "radia_export_menu.py"
    if not template.is_file():
        raise FileNotFoundError(f"official toolbar template missing: {template}")
    if not menu_source.is_file():
        raise FileNotFoundError(f"export implementation missing: {menu_source}")

    destination = Path(output_dir or _startup_dir(all_users=False))
    destination.mkdir(parents=True, exist_ok=True)
    package_path = destination / "radia_export_toolbar.tar.gz"

    with tempfile.TemporaryDirectory(
            prefix="radia-export-toolbar-", dir=str(destination)) as temp_dir:
        staging = Path(temp_dir)
        shutil.copytree(
            source / "scripts",
            staging / "scripts",
            ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
        )
        shutil.copytree(source / "icons", staging / "icons")
        (staging / "toolbars").mkdir()
        shutil.copy2(menu_source, staging / "scripts" / "radia_export_menu.py")

        install_dir = staging.as_posix()
        toolbar_text = template.read_text(encoding="utf-8").replace(
            "@TOOLBAR_INSTALL_DIR@", install_dir)
        toolbar_rel = Path("toolbars") / "radia_export_toolbar.ttb"
        (staging / toolbar_rel).write_text(toolbar_text, encoding="utf-8")

        mapped_files = sorted(
            path.relative_to(staging).as_posix()
            for path in staging.rglob("*")
            if path.is_file()
        )
        mappings = "\n".join(
            f"{install_dir}/{relative} => {relative}"
            for relative in mapped_files
        ) + "\n"
        (staging / ".mappings").write_text(mappings, encoding="utf-8")

        with tarfile.open(package_path, "w:gz") as archive:
            for relative in ("scripts", "toolbars", "icons", ".mappings"):
                archive.add(staging / relative, arcname=relative)

    return str(package_path)


def _generate_startup_script(panels_dir, *, all_users=False):
    """Generate the Cubit-played Python shim outside the package tree."""
    register_path = os.path.join(panels_dir, "register_toolbar.py").replace("\\", "/")
    startup_root = _startup_dir(all_users=all_users)
    os.makedirs(startup_root, exist_ok=True)
    startup_path = os.path.join(startup_root, "radia_startup.py")

    content = (
        "#!python\n"
        "import sys, os, glob; "
        "_cb = os.path.dirname(os.path.abspath(os.path.join("
        "os.path.dirname(sys.executable), \"cubit.py\"))) "
        "if not hasattr(sys, \"_cubit_bin\") else sys._cubit_bin; "
        "_sp = glob.glob(os.path.join(_cb, \"python*\", \"lib\", "
        "\"site-packages\")) + glob.glob(os.path.join(_cb, \"python*\", "
        "\"lib\", \"python*\", \"site-packages\")); "
        "sys.path.insert(0, _sp[0]) if _sp and _sp[0] not in sys.path else None; "
        f"__file__ = r\"{register_path}\"; "
        "exec(\"try:\\n"
        f" exec(open(r'{register_path}', encoding='utf-8').read())\\n"
        "except Exception as e:\\n"
        " import traceback; traceback.print_exc()\")\n"
    )

    with open(startup_path, "w", encoding="utf-8") as f:
        f.write(content)
    return startup_path


def _build_startup_block(startup_script_path):
    """Build the block to insert into a .cubit startup file."""
    startup_script_path = startup_script_path.replace("\\", "/")
    return (
        f"{_MARKER_BEGIN}\n"
        "set journal off\n"
        f"play \"{startup_script_path}\"\n"
        f"{_MARKER_END}\n"
    )


def _remove_existing_block(lines):
    """Remove current and legacy toolbar blocks from a .cubit file."""
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
    """Return Cubit.ini paths targeted by this install."""
    paths = []
    if sys.platform != "win32":
        return paths

    appdata = os.environ.get("APPDATA", "")
    if appdata:
        paths.append(os.path.join(appdata, "Coreform", "Cubit.ini"))

    if all_users:
        for user_dir in _iter_windows_user_dirs(include_default=True):
            ini = os.path.join(user_dir, "AppData", "Roaming", "Coreform", "Cubit.ini")
            if ini not in paths:
                paths.append(ini)
    return paths


def _read_text(path):
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        return f.read()


def _ensure_plugin_path_in_ini(ini_path, plugin_dir):
    """Ensure ``plugin_dir`` is registered in Cubit.ini plugin\\paths."""
    plugin_dir = plugin_dir.replace("\\", "/")

    if not os.path.isfile(ini_path):
        os.makedirs(os.path.dirname(ini_path), exist_ok=True)
        with open(ini_path, "w", encoding="utf-8") as f:
            f.write("[clarofw]\n")
            f.write(f"plugin\\paths={plugin_dir}\n")
        return True

    content = _read_text(ini_path)
    lines = content.splitlines(keepends=True)

    found = False
    updated = False
    for i, line in enumerate(lines):
        if line.startswith("plugin\\paths="):
            found = True
            current = line.split("=", 1)[1].strip()
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
        for i, line in enumerate(lines):
            if line.strip() == "[clarofw]":
                lines.insert(i + 1, f"plugin\\paths={plugin_dir}\n")
                updated = True
                break
        else:
            lines.append("\n[clarofw]\n")
            lines.append(f"plugin\\paths={plugin_dir}\n")
            updated = True

    if updated:
        with open(ini_path, "w", encoding="utf-8") as f:
            f.writelines(lines)
    return updated


def _startup_play_path(cubit_file):
    if not os.path.isfile(cubit_file):
        return None
    content = _read_text(cubit_file)
    begin = content.find(_MARKER_BEGIN)
    end = content.find(_MARKER_END)
    if begin < 0 or end < begin:
        return None
    block = content[begin:end]
    match = re.search(r'play\s+"([^"]+)"', block)
    if not match:
        return None
    return match.group(1).replace("/", os.sep)


def _ini_has_plugin_path(ini_path, plugin_dir):
    if not os.path.isfile(ini_path):
        return False
    content = _read_text(ini_path)
    plugin_dir = plugin_dir.replace("\\", "/")
    for line in content.splitlines():
        if not line.startswith("plugin\\paths="):
            continue
        current = line.split("=", 1)[1].strip()
        existing = [p.strip().replace("\\", "/") for p in current.split(",")]
        if plugin_dir in existing:
            return True
    return False


def verify_panel_installation(all_users=False, verbose=True):
    """Verify that Cubit startup and plugin path registration are in place."""
    issues = []
    panels_dir = _get_panels_dir()
    register_script = os.path.join(panels_dir, "register_toolbar.py")
    register_norm = register_script.replace("\\", "/")

    if not os.path.isfile(register_script):
        issues.append(f"toolbar script missing: {register_script}")

    cubit_bin = find_cubit_bin()
    if not cubit_bin:
        issues.append(
            f"Coreform Cubit {_MIN_CUBIT_VERSION_TEXT}+ not found; "
            "set CUBIT_PATH to the 2025.12 bin directory"
        )

    for cubit_file in _get_cubit_startup_files(all_users=all_users):
        play_path = _startup_play_path(cubit_file)
        if not play_path:
            issues.append(f"startup block missing from {cubit_file}")
            continue
        if not os.path.isfile(play_path):
            issues.append(f"startup script missing: {play_path}")
            continue
        startup_text = _read_text(play_path)
        if register_norm not in startup_text.replace("\\", "/"):
            issues.append(
                f"startup script {play_path} does not load {register_norm}"
            )
        cubit_text = _read_text(cubit_file)
        for begin, end in _LEGACY_MARKERS:
            if begin in cubit_text or end in cubit_text:
                issues.append(f"legacy toolbar marker remains in {cubit_file}")

    if cubit_bin:
        plugin_dir = os.path.join(cubit_bin, "plugins")
        for ini_path in _get_cubit_ini_paths(all_users=all_users):
            if not _ini_has_plugin_path(ini_path, plugin_dir):
                issues.append(f"plugin path missing from {ini_path}: {plugin_dir}")

    if verbose:
        print("Panel registration verification:")
        if issues:
            for issue in issues:
                print(f"  [FAIL] {issue}")
        else:
            print("  [OK] .cubit startup and Cubit.ini plugin paths are registered")

    return (not issues), issues


def install_panels(all_users=False):
    """Register the Radia toolbar for Coreform Cubit 2025.12+."""
    print("=== Coreform Cubit - Panel Installer ===\n")

    panels_dir = _get_panels_dir()
    register_script = os.path.join(panels_dir, "register_toolbar.py")
    if not os.path.isfile(register_script):
        print(f"ERROR: Toolbar script not found: {register_script}")
        return False

    try:
        startup_script = _generate_startup_script(panels_dir, all_users=all_users)
        toolbar_package = build_official_toolbar_package(
            output_dir=_startup_dir(all_users=all_users))
    except OSError as exc:
        print(f"ERROR: could not create Cubit startup assets: {exc}")
        return False

    print(f"Toolbar script: {register_script}")
    print(f"Startup script: {startup_script}")
    print(f"Official toolbar package: {toolbar_package}")

    cubit_bin = find_cubit_bin()
    if not cubit_bin:
        print(f"ERROR: Coreform Cubit {_MIN_CUBIT_VERSION_TEXT}+ not found.")
        print("       Set CUBIT_PATH to the Cubit 2025.12 bin directory.")
        return False
    print(f"Cubit bin:      {cubit_bin}")

    if all_users:
        print("Mode:           --all-users (all existing profiles + Default)")

    errors = []
    block = _build_startup_block(startup_script)
    for cubit_file in _get_cubit_startup_files(all_users=all_users):
        try:
            os.makedirs(os.path.dirname(cubit_file), exist_ok=True)
            if os.path.isfile(cubit_file):
                with open(cubit_file, "r", encoding="utf-8", errors="replace") as f:
                    lines = f.readlines()
            else:
                lines = []
            lines = _remove_existing_block(lines)
            lines.append("\n" + block)
            with open(cubit_file, "w", encoding="utf-8") as f:
                f.writelines(lines)
            print(f"Updated: {cubit_file}")
        except OSError as exc:
            msg = f"could not update {cubit_file}: {exc}"
            errors.append(msg)
            print(f"ERROR: {msg}")

    plugin_dir = os.path.join(cubit_bin, "plugins")
    for ini_path in _get_cubit_ini_paths(all_users=all_users):
        try:
            if _ensure_plugin_path_in_ini(ini_path, plugin_dir):
                print(f"Plugin path registered: {ini_path}")
            else:
                print(f"Plugin path OK: {ini_path}")
        except OSError as exc:
            msg = f"could not update {ini_path}: {exc}"
            errors.append(msg)
            print(f"ERROR: {msg}")

    ok, issues = verify_panel_installation(all_users=all_users, verbose=True)
    if errors or not ok:
        print()
        print("=== Installation FAILED ===")
        for err in errors:
            print(f"  - {err}")
        for issue in issues:
            print(f"  - {issue}")
        return False

    print()
    print("=== Installation Complete ===")
    print("Import the official toolbar package once in Cubit via:")
    print("  Tools > Custom Toolbar Editor > Import > Package")
    print(f"  {toolbar_package}")
    print("Then restart Cubit 2025.12 to verify persistent loading.")
    return True


def uninstall_panels(all_users=False):
    """Remove the Radia toolbar registration from .cubit files."""
    print("=== Coreform Cubit - Panel Uninstaller ===\n")

    for cubit_file in _get_cubit_startup_files(all_users=all_users):
        if not os.path.isfile(cubit_file):
            continue
        with open(cubit_file, "r", encoding="utf-8", errors="replace") as f:
            lines = f.readlines()
        new_lines = _remove_existing_block(lines)
        if len(new_lines) < len(lines):
            try:
                with open(cubit_file, "w", encoding="utf-8") as f:
                    f.writelines(new_lines)
                print(f"Removed toolbar from: {cubit_file}")
            except OSError as exc:
                print(f"ERROR: could not update {cubit_file}: {exc}")

    startup_script = os.path.join(_startup_dir(all_users=all_users), "radia_startup.py")
    if os.path.isfile(startup_script):
        try:
            os.remove(startup_script)
            print(f"Removed startup script: {startup_script}")
        except OSError as exc:
            print(f"ERROR: could not remove {startup_script}: {exc}")

    print("Restart Cubit to apply changes.")
    return True


def main():
    """Console script entry point."""
    all_users = "--all-users" in sys.argv
    if "--build-toolbar-only" in sys.argv:
        package = build_official_toolbar_package()
        print(package)
        success = True
    elif "--uninstall" in sys.argv:
        success = uninstall_panels(all_users=all_users)
    elif "--verify-only" in sys.argv:
        success, _ = verify_panel_installation(all_users=all_users, verbose=True)
    else:
        success = install_panels(all_users=all_users)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
