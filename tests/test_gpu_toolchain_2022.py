"""Compiler selection fails before native compilation on an incompatible host."""
import importlib.util
import sys
from pathlib import Path
import pytest

P = Path(__file__).parents[1] / 'packages/radia-simulink/tools/matlab_cuda_toolchain.py'
sys.path.insert(0, str(P.parent))
spec = importlib.util.spec_from_file_location('gpu_toolchain_2022', P)
toolchain = importlib.util.module_from_spec(spec)
spec.loader.exec_module(toolchain)


@pytest.fixture(autouse=True)
def sibling_tool_import(monkeypatch):
    monkeypatch.syspath_prepend(str(P.parent))


@pytest.mark.parametrize('key,value', [
    ('CUDA_PATH', 'C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v12.8'),
    ('RADIA_GPU_CUDA_ROOT', 'C:/temp/cuda134'),
    ('RADIA_MEXCUDA_CONFIG', 'private.xml'),
    ('NVCC_APPEND_FLAGS', '--allow-unsupported-compiler'),
    ('NVCC_PREPEND_FLAGS', '-ccbin other'),
])
def test_external_overrides_fail(key, value):
    with pytest.raises(RuntimeError):
        toolchain.reject_overrides({key: value})


def test_bundled_root_allowed():
    toolchain.reject_overrides({'CUDA_PATH': 'C:/Program Files/MATLAB/R2026a/sys/cuda/win64/cuda'})


def test_vs2026_only_is_not_fallback(monkeypatch):
    monkeypatch.setattr(toolchain.os, 'environ', {})
    def locate(command, **kwargs):
        assert command[command.index('-version') + 1] == '[17.0,18.0)'
        return ''
    monkeypatch.setattr(toolchain.subprocess, 'check_output', locate)
    with pytest.raises(RuntimeError, match='no compiler fallback'):
        toolchain.discover_vs()


@pytest.mark.parametrize('version,accepted', [('14.44.35207', True), ('14.51.36231', False)])
def test_default_toolset_boundary(tmp_path, monkeypatch, version, accepted):
    monkeypatch.setattr(toolchain.os, 'environ', {})
    auxiliary = tmp_path / 'VC/Auxiliary/Build'
    auxiliary.mkdir(parents=True)
    (auxiliary / 'Microsoft.VCToolsVersion.default.txt').write_text(version)
    (auxiliary / 'vcvars64.bat').touch()
    compiler = tmp_path / 'VC/Tools/MSVC' / version / 'bin/Hostx64/x64'
    compiler.mkdir(parents=True)
    (compiler / 'cl.exe').touch()
    monkeypatch.setattr(toolchain.subprocess, 'check_output', lambda *a, **k: str(tmp_path))
    if accepted:
        assert toolchain.discover_vs() == tmp_path
    else:
        with pytest.raises(RuntimeError, match='14.44'):
            toolchain.discover_vs()


PROFILE = r'''<config Name="NVIDIA CUDA Compiler" ShortName="nvcc"
HostCompiler="MSVCPP170" Language="CUDA" Version="12.2"
Location="$CUDA_ROOT\bin\nvcc.exe"><vars CMDLINE100="cl file && echo /Fo > NUL"/></config>'''


def test_mathworks_command_template_is_not_generic_xml():
    assert toolchain.validate_stock_profile(PROFILE)['Version'] == '12.2'


@pytest.mark.parametrize('text', [PROFILE.replace('HostCompiler=', 'OtherCompiler='),
                                 PROFILE.replace('MSVCPP170', 'MSVCPP180'),
                                 PROFILE + '<vars FLAGS="--allow-unsupported-compiler"/>'])
def test_profile_identity_and_unsupported_override_fail(text):
    with pytest.raises(RuntimeError):
        toolchain.validate_stock_profile(text)


STOCK = Path('C:/Program Files/MATLAB/R2026a/toolbox/parallel/gpu/extern/src/mex/win64/nvcc_msvcpp2022.xml')


@pytest.mark.skipif(not STOCK.is_file(), reason='Optional installed MATLAB compiler-profile integration check')
def test_installed_stock_profile_read_only():
    before = STOCK.read_bytes()
    assert toolchain.validate_stock_profile(before.decode('utf-8'))['HostCompiler'] == 'MSVCPP170'
    assert STOCK.read_bytes() == before
