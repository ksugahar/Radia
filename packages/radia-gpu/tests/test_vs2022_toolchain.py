"""Toolchain discovery rejects unsupported versions and incomplete installs."""
import importlib.util
from pathlib import Path
import pytest

TOOLS = Path(__file__).resolve().parents[1] / 'tools'
spec = importlib.util.spec_from_file_location('vs2022_toolchain', TOOLS / 'vs2022_toolchain.py')
selector = importlib.util.module_from_spec(spec)
spec.loader.exec_module(selector)


def install(tmp_path, version='14.44.35207'):
    vs = tmp_path / 'arbitrary VS2022 location'
    aux = vs / 'VC/Auxiliary/Build'
    aux.mkdir(parents=True)
    (aux / 'Microsoft.VCToolsVersion.default.txt').write_text(version)
    (aux / 'vcvars64.bat').touch()
    compiler = vs / 'VC/Tools/MSVC' / version / 'bin/Hostx64/x64'
    compiler.mkdir(parents=True)
    (compiler / 'cl.exe').touch()
    return vs, compiler


def test_nonstandard_installation_and_supported_version(tmp_path, monkeypatch):
    vs, compiler = install(tmp_path)
    calls = []
    def locate(command, **kwargs):
        calls.append(command)
        return str(vs) + '\n'
    monkeypatch.setattr(selector.subprocess, 'check_output', locate)
    assert selector.discover_vs() == vs
    assert selector.compiler_bin(vs) == compiler
    assert calls[0][calls[0].index('-version') + 1] == '[17.0,18.0)'
    assert '-requires' in calls[0]


def test_missing_vs_rejected(monkeypatch):
    monkeypatch.setattr(selector.subprocess, 'check_output', lambda *a, **k: '')
    with pytest.raises(RuntimeError, match='VS2022 required'):
        selector.discover_vs()


@pytest.mark.parametrize('version', ['14.43.1', '14.50.1', '14.44.1junk'])
def test_unsupported_default_compiler_rejected(tmp_path, monkeypatch, version):
    vs, _ = install(tmp_path, version)
    monkeypatch.setattr(selector.subprocess, 'check_output', lambda *a, **k: str(vs))
    with pytest.raises(RuntimeError, match='MSVC14.44 required'):
        selector.discover_vs()


@pytest.mark.parametrize('component', ['compiler', 'vcvars'])
def test_incomplete_installation_rejected(tmp_path, monkeypatch, component):
    vs, compiler = install(tmp_path)
    missing = compiler / 'cl.exe' if component == 'compiler' else vs / 'VC/Auxiliary/Build/vcvars64.bat'
    missing.unlink()
    monkeypatch.setattr(selector.subprocess, 'check_output', lambda *a, **k: str(vs))
    with pytest.raises(RuntimeError, match='missing'):
        selector.discover_vs()
