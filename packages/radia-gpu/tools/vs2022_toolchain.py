"""Locate the supported installed VS2022/MSVC14.44 without a fixed host path."""
from pathlib import Path
import os
import re
import subprocess


def discover_vs():
    locator = Path(os.environ.get('ProgramFiles(x86)', 'C:/Program Files (x86)')) / 'Microsoft Visual Studio/Installer/vswhere.exe'
    value = subprocess.check_output([
        str(locator), '-latest', '-version', '[17.0,18.0)', '-products', '*',
        '-requires', 'Microsoft.VisualStudio.Component.VC.Tools.x86.x64',
        '-property', 'installationPath'], text=True).strip()
    if not value:
        raise RuntimeError('VS2022 required; no compiler fallback')
    vs = Path(value)
    compiler_bin(vs)
    if not (vs / 'VC/Auxiliary/Build/vcvars64.bat').is_file():
        raise RuntimeError('VS2022 vcvars64.bat missing')
    return vs


def compiler_bin(vs):
    version = (vs / 'VC/Auxiliary/Build/Microsoft.VCToolsVersion.default.txt').read_text().strip()
    if not re.fullmatch(r'14\.44\.\d+', version):
        raise RuntimeError(f'MSVC14.44 required, got {version}')
    compiler = vs / 'VC/Tools/MSVC' / version / 'bin/Hostx64/x64'
    if not (compiler / 'cl.exe').is_file():
        raise RuntimeError('MSVC14.44 cl.exe missing')
    return compiler
