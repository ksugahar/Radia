"""Build and install the patched Netgen library consumed by the HDiv installed-wheel lane.

NGSolve and Netgen stay at their pinned wheels. Pinned Netgen scales every pyramid z by
1 - 1e-10; tools/patches/netgen-pyramid-apex.patch guards the apex only. This tool builds
that source once and replaces nglib.dll, and nothing else, in a run-local environment. A
receipt binds the library to the upstream commit and the patch. Never point it at a shared
runtime; a stock runtime keeps the stock library.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
PATCH = ROOT / 'tools/patches/netgen-pyramid-apex.patch'
UPSTREAM = 'https://github.com/NGSolve/netgen'
TAG = 'v6.2.2608'
UPSTREAM_SHA = '96e5682f6ea43ba77ba3bb2ae4bb4bd1791f506e'
LIBRARY = 'nglib.dll'
RECEIPT = 'netgen-overlay.json'
SCHEMA = 'radia.netgen-overlay.v1'
BUILD_PINS = ('cmake==3.31.10', 'ninja==1.13.2', 'numpy==2.5.4', 'packaging==26.3', 'pathspec==1.1.1',
              'scikit-build-core==1.1.2', 'pybind11-stubgen==3.0.0', 'netgen-occt==7.8.1', 'netgen-occt-devel==7.8.1')
SETTINGS = ('cmake.define.USE_GUI=OFF', 'cmake.define.USE_CCACHE=OFF', 'cmake.define.BUILD_STUB_FILES=OFF')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(command, **kwargs):
    print(' '.join(str(part) for part in command), flush=True)
    subprocess.run([str(part) for part in command], check=True, **kwargs)


def build(work, out):
    """Compile in an environment that already holds MSVC; work is deleted by its owner."""
    work.mkdir(parents=True, exist_ok=False)
    out.mkdir(parents=True, exist_ok=False)
    source = work / 'netgen'
    run(['git', 'clone', '--branch', TAG, '--recurse-submodules', '--shallow-submodules', '--depth', '1', UPSTREAM, source])
    head = subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip()
    if head != UPSTREAM_SHA:
        raise ValueError(f'{TAG} resolves to {head}, expected {UPSTREAM_SHA}')
    run(['git', '-C', source, 'apply', '--check', PATCH])
    run(['git', '-C', source, 'apply', PATCH])
    venv = work / 'venv'
    run([sys.executable, '-m', 'venv', venv])
    python = venv / 'Scripts/python.exe'
    run([python, '-I', '-m', 'pip', 'install', '--only-binary=:all:', *BUILD_PINS])
    env = dict(os.environ, PYTHONNOUSERSITE='1', NETGEN_ARCH='avx2', CMAKE_GENERATOR='Ninja', CC='cl', CXX='cl',
               OpenCascade_DIR=str(venv / 'Lib/cmake/opencascade'), CMAKE_PREFIX_PATH=str(venv / 'Lib'),
               PATH=str(venv / 'Scripts') + os.pathsep + os.environ['PATH'])
    wheels = work / 'wheelhouse'
    run([python, '-I', '-m', 'pip', 'wheel', source, '--no-deps', '--no-build-isolation', '--wheel-dir', wheels,
         *[f'--config-settings={item}' for item in SETTINGS], f'--config-settings=build-dir={work / "build"}',
         f'--config-settings=cmake.define.CMAKE_MAKE_PROGRAM={(venv / "Scripts/ninja.exe").as_posix()}'],
        env=env, cwd=work)
    built = list(wheels.glob('netgen_mesher-*.whl'))
    if len(built) != 1:
        raise ValueError(f'Expected one Netgen wheel, found {[path.name for path in built]}')
    with zipfile.ZipFile(built[0]) as wheel:
        (out / LIBRARY).write_bytes(wheel.read(f'netgen/{LIBRARY}'))
    (out / RECEIPT).write_text(json.dumps(dict(
        schema=SCHEMA, upstream=UPSTREAM, tag=TAG, upstream_sha=UPSTREAM_SHA, patch_sha256=digest(PATCH),
        library=LIBRARY, library_sha256=digest(out / LIBRARY), build_pins=list(BUILD_PINS), settings=list(SETTINGS),
        netgen_arch='avx2'), indent=1) + '\n', encoding='utf-8')


def checked_receipt(overlay):
    receipt = json.loads((Path(overlay) / RECEIPT).read_text(encoding='utf-8'))
    if receipt.get('schema') != SCHEMA or receipt.get('upstream_sha') != UPSTREAM_SHA:
        raise ValueError(f"Overlay is {receipt.get('schema')} at {receipt.get('upstream_sha')}, expected {SCHEMA} at {UPSTREAM_SHA}")
    if receipt.get('patch_sha256') != digest(PATCH):
        raise ValueError(f"Overlay patch {receipt.get('patch_sha256')}, this checkout {digest(PATCH)}")
    actual = digest(Path(overlay) / LIBRARY)
    if actual != receipt.get('library_sha256'):
        raise ValueError(f"Overlay library {actual}, receipt {receipt.get('library_sha256')}")
    return receipt


def locate(python):
    """The netgen package of that interpreter, which must live under its own environment."""
    probe = ('import importlib.metadata as m, importlib.util as u, json, pathlib, sys; '
             's = u.find_spec("netgen"); root = pathlib.Path(s.submodule_search_locations[0]).resolve(); '
             'print(json.dumps(dict(package=str(root), prefix=str(pathlib.Path(sys.prefix).resolve()), '
             'base=str(pathlib.Path(sys.base_prefix).resolve()), version=m.version("netgen-mesher"))))')
    found = json.loads(subprocess.check_output([str(python), '-I', '-c', probe], text=True))
    if found['prefix'] == found['base']:
        raise ValueError(f"{python} is not a virtual environment; the overlay installs only into a run-local one")
    if not Path(found['package']).is_relative_to(found['prefix']):
        raise ValueError(f"netgen at {found['package']} is outside {found['prefix']}")
    if found['version'] != TAG[1:]:
        raise ValueError(f"netgen-mesher {found['version']}, overlay built from {TAG}")
    return Path(found['package'])


def install(overlay, python):
    receipt = checked_receipt(overlay)
    target = locate(python) / LIBRARY
    previous = digest(target)
    shutil.copyfile(Path(overlay) / LIBRARY, target)
    if digest(target) != receipt['library_sha256']:
        raise ValueError(f'{target} does not hold the overlay library after the copy')
    return dict(previous_sha256=previous, library_sha256=receipt['library_sha256'], target=target.name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='mode', required=True)
    compile_ = commands.add_parser('build')
    compile_.add_argument('--work', type=Path, required=True)
    compile_.add_argument('--out', type=Path, required=True)
    place = commands.add_parser('install')
    place.add_argument('--overlay', type=Path, required=True)
    place.add_argument('--python', type=Path, required=True)
    args = parser.parse_args()
    if args.mode == 'build':
        build(args.work.resolve(), args.out.resolve())
    else:
        print(json.dumps(install(args.overlay, args.python)))


if __name__ == '__main__':
    main()
