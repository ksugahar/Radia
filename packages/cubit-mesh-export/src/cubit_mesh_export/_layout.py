"""Installed Cubit layout, including the package-based 2026.8 runtime."""
from pathlib import Path
import os


def has_python_binding(bin_dir: Path) -> bool:
    return any((bin_dir / name).is_file() for name in (
        'cubit.py', '_cubit3.pyd', '_cubit3.so', 'cubit/__init__.py'))


def plugin_directory(install_dir: Path) -> Path:
    root = Path(install_dir)
    if root.name.lower() == 'bin':
        root = root.parent
    if (root / 'bin/cubitx.exe').is_file():
        return root / 'plugins'
    return root / 'bin/plugins'


def console_executable(bin_dir: Path) -> Path | None:
    for name in ('cubitx.exe', 'coreform_cubit.com', 'coreform_cubit', 'cubit'):
        candidate = Path(bin_dir) / name
        if candidate.is_file():
            return candidate
    return None


def command_plugin_path(install_dir: Path) -> Path:
    root = Path(install_dir)
    if root.name.lower() == 'bin':
        root = root.parent
    # The 2026.8 SDK loader scans Windows DLLs. Keep .ccm as the package's
    # canonical payload, and use the loader's required name at deployment.
    suffix = '.dll' if (root / 'bin/cubitx.exe').is_file() else '.ccm'
    return plugin_directory(root) / ('cubit_mesh_export' + suffix)


def runtime_environment(bin_dir: Path) -> dict[str, str]:
    env = os.environ.copy()
    acis = Path(bin_dir).parent / 'acis/code/bin'
    if acis.is_dir():
        env['PATH'] = str(acis) + os.pathsep + env.get('PATH', '')
    if (Path(bin_dir) / 'cubitx.exe').is_file():
        env.setdefault('CUBIT_PLUGIN_DIR', str(plugin_directory(Path(bin_dir))))
    return env
