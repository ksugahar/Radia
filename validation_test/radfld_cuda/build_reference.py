"""Build an independent CPU pybind oracle from an explicit pre-CUDA Git commit."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys
import sysconfig
import socket


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--commit', default='69301918a')
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    if socket.gethostname().lower() != 't7920':
        raise RuntimeError('Native candidate builds belong to 100')
    import netgen
    root = Path(__file__).resolve().parents[2]
    revision = subprocess.check_output(['git', '-C', str(root), 'rev-parse', args.commit], text=True).strip()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    hashes = {}
    for name in ('rad_hdiv_field_evaluator.cpp', 'rad_hdiv_field_evaluator.h',
                 'rad_hdiv_vim.cpp', 'rad_hdiv_vim.h', 'rad_parallel.h'):
        data = subprocess.check_output(['git', '-C', str(root), 'show', f'{revision}:src/core/{name}'])
        (out/name).write_bytes(data)
        hashes[name] = hashlib.sha256(data).hexdigest()
    net = Path(netgen.__file__).parent
    vs = Path('C:/BuildTools/VS2022-CUDA')
    binary = out/'_radia_field_reference.pyd'
    adapter = Path(__file__).with_name('field_reference.cpp').read_bytes()
    source = out/'field_reference.cpp'
    source.write_bytes(adapter)
    python_lib = Path(sys.base_prefix)/'libs'/f'python{sys.version_info.major}{sys.version_info.minor}.lib'
    command = ['cl', '/nologo', '/LD', '/EHsc', '/MD', '/O2', '/fp:precise', '/std:c++20', '/bigobj', '/utf-8',
               '/DWIN32', '/DWNT', '/DMSVC_EXPRESS', '/DMAX_SYS_DIM=3', '/DLAPACK', '/DNDEBUG', '/DNOMINMAX',
               '/DNGS_PYTHON', '/DNETGEN_PYTHON', '/DNG_PYTHON', '/DHAVE_NETGEN_SOURCES',
               '/DHAVE_STRUCT_TIMESPEC', '/DPYBIND11_SIMPLE_GIL_MANAGEMENT', '/DPy_NO_LINK_LIB',
               f'/I{net / "include"}', f'/I{net / "include/include"}', f'/I{sysconfig.get_path("include")}',
               f'/I{out}', f'/Fo{out}/', str(source), str(out/'rad_hdiv_field_evaluator.cpp'),
               str(out/'rad_hdiv_vim.cpp'), '/link', f'/OUT:{binary}', f'/IMPLIB:{out/"reference.lib"}',
               str(python_lib), str(net/'lib/libngsolve.lib'), str(net/'lib/ngcore.lib'), str(net/'lib/nglib.lib')]
    batch = out/'build.cmd'
    batch.write_text('@echo off\ncall "'+str(vs/'VC/Auxiliary/Build/vcvars64.bat')+
                     '" >nul\nif errorlevel 1 exit /b 1\n'+subprocess.list2cmdline(command)+
                     '\nexit /b %errorlevel%\n', encoding='utf-8')
    with (out/'build.log').open('w') as log:
        subprocess.run(['cmd', '/c', str(batch)], stdout=log, stderr=subprocess.STDOUT, check=True)
    receipt = dict(schema='radia.field-cpu-reference.v1', source_commit=revision, source_sha256=hashes,
                   adapter_sha256=hashlib.sha256(adapter).hexdigest(),
                   binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                   ngsolve=__import__('ngsolve').__version__, python=sys.version)
    (out/'reference.json').write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8')


if __name__ == '__main__':
    main()
