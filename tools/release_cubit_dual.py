"""Independent Cubit LAB/100 gate, invoked only through release_quad.py.

The published wheel is the byte oracle; only its owning editable package is
installed. Radia/MCP metadata and their selected source remain untouched.
"""
from __future__ import annotations

import base64
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import urllib.request
import zipfile
from email.parser import BytesParser

SCHEMA = 'cubit-mesh-export.release-dual.v2'
TARGETS = ('lab', '100')


def digest(data):
    return hashlib.sha256(data).hexdigest()


def wheel_contract(path, distribution='cubit-mesh-export'):
    prefix = distribution.replace('-', '_') + '/'
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        metadata = BytesParser().parsebytes(archive.read(next(
            n for n in names if n.endswith('.dist-info/METADATA'))))
        if metadata['Name'] != distribution:
            raise ValueError('Not a ' + distribution + ' wheel')
        files = {}
        for name in names:
            if name.startswith(prefix) and not name.endswith('/'):
                content = archive.read(name)
                text = name.endswith(('.py', '.jou', '.svg', '.tmpl', '.json', '.md'))
                if text:
                    content = content.replace(b'\r\n', b'\n').replace(b'\r', b'\n')
                files[name] = {'sha256': digest(content), 'text': text}
    required_files = ('__init__.py', 'common/status.py') if distribution == 'cae-mcp-core' else (
        'cubit_mesh_export.ccm', 'cubit_mesh_curver.pyd', 'toolbar_smoke.py',
        'cubit_gui/toolbar_probe.py', 'mcp/server.py')
    for required in required_files:
        if prefix + required not in files:
            raise ValueError('Wheel lacks ' + required)
    return {'version': metadata['Version'], 'wheel_sha256': digest(Path(path).read_bytes()),
            'files': files, 'requires': metadata.get_all('Requires-Dist', [])}


def verify_published(distribution, contract):
    with urllib.request.urlopen('https://pypi.org/pypi/' + distribution + '/' +
                                contract['version'] + '/json', timeout=30) as response:
        published = json.load(response)
    if contract['wheel_sha256'] not in [u['digests']['sha256'] for u in published['urls']]:
        raise ValueError(distribution + ': wheel bytes are not the published PyPI artifact')


# Same worker executes locally or through the existing Administrator SSH token.
# No GUI process is killed to make installation possible; active Cubit refuses it.
WORKER = r'''
import base64, getpass, hashlib, importlib.metadata as md, json, os
from pathlib import Path
import re, shutil, socket, subprocess, sys, zipfile
cfg = json.loads(base64.b64decode(sys.argv[1]))
root = Path(cfg['source_root']).resolve()
package = root / 'packages/cubit-mesh-export'
core_package = root / 'packages/cae-mcp-core'
out = Path(cfg['output']).resolve()
out.mkdir(parents=True, exist_ok=True)
events = []
result = dict(schema=cfg['schema'], target=cfg['target'], hostname=socket.gethostname(),
              user=getpass.getuser(), interpreter=sys.executable, source_root=str(root),
              source_sha=cfg['source_sha'], version=cfg['version'],
              wheel_sha256=cfg['wheel_sha256'], core=cfg['core'], passed=False)

def command(args, timeout=120):
    p = subprocess.run(args, capture_output=True, text=True, encoding='utf-8',
                       errors='replace', timeout=timeout)
    events.append(dict(command=args, returncode=p.returncode, stdout=p.stdout, stderr=p.stderr))
    if p.returncode:
        raise RuntimeError(str(args) + '\n' + p.stdout[-5000:] + p.stderr[-2000:])
    return p.stdout

def preserved():
    state = {}
    for name in ('radia', 'radia-mcp', 'netgen-mesher', 'ngsolve', 'numpy', 'mcp'):
        try:
            d = md.distribution(name)
            state[name] = dict(version=d.version, direct_url=d.read_text('direct_url.json'))
        except md.PackageNotFoundError:
            state[name] = None
    return state

def verify_files(base, files):
    for relative, expected in files.items():
        data = (base / relative).read_bytes()
        if expected['text']:
            data = data.replace(b'\r\n', b'\n').replace(b'\r', b'\n')
        if hashlib.sha256(data).hexdigest() != expected['sha256']:
            raise RuntimeError('Published wheel/source mismatch: ' + relative)

def recover_scratch(text):
    # Only the exact task-created directory reported by our command is eligible.
    for value in re.findall(r'^\s*Work:\s*(.+)$', text, re.M):
        scratch = Path(value.strip()).resolve()
        if scratch.parent != Path('C:/temp').resolve() or not scratch.name.startswith(
                ('cubit-smoke-', 'cubit-toolbar-smoke-')):
            raise RuntimeError('Unsafe smoke scratch path: ' + str(scratch))
        entries = [scratch, *scratch.rglob('*')]
        if any(p.is_symlink() or p.is_junction() for p in entries):
            raise RuntimeError('Smoke scratch contains a link')
        destination = out / scratch.name
        shutil.copytree(scratch, destination)
        for p in entries:
            if p.is_file() and hashlib.sha256(p.read_bytes()).digest() != hashlib.sha256(
                    (destination / p.relative_to(scratch)).read_bytes()).digest():
                raise RuntimeError('Evidence copy mismatch')
        shutil.rmtree(scratch)
        if scratch.exists():
            raise RuntimeError('Smoke scratch cleanup failed')

before = preserved()
try:
    git = ['git', '-c', 'safe.directory=' + str(root), '-C', str(root)]
    if command(git + ['rev-parse', 'HEAD']).strip() != cfg['source_sha']:
        raise RuntimeError('Wrong release checkout SHA')
    if command(git + ['status', '--porcelain', '--untracked-files=no']).strip():
        raise RuntimeError('Release checkout has tracked modifications')
    verify_files(package / 'src', cfg['files'])
    verify_files(core_package / 'src', cfg['core']['files'])
    from packaging.requirements import Requirement
    for raw in cfg['requires'] + cfg['core']['requires']:
        requirement = Requirement(raw)
        if requirement.marker and not requirement.marker.evaluate():
            continue
        version = cfg['core']['version'] if requirement.name == 'cae-mcp-core' else md.version(requirement.name)
        if version not in requirement.specifier:
            raise RuntimeError('Unsatisfied dependency; update explicitly before deploy: ' + raw)
    for name in ('netgen-mesher', 'ngsolve'):
        if md.version(name) != cfg['dependencies'][name]:
            raise RuntimeError('Dependency mismatch; do not mutate shared solver runtime: ' + name)
    # Load the selected standalone installer directly; avoid old editable imports.
    import importlib.util
    spec = importlib.util.spec_from_file_location('selected_cubit_install', package / 'src/cubit_mesh_export/install.py')
    installer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(installer)
    cubit_dir = installer._find_cubit_dir()
    if cubit_dir is None:
        raise RuntimeError('Supported Cubit installation missing')
    ok, issues = installer.preflight(cubit_dir, verbose=False)
    if not ok:
        raise RuntimeError(str(issues))
    if cfg['action'] == 'deploy':
        command([sys.executable, '-m', 'pip', 'install', '--no-deps', '--no-build-isolation',
                 '-e', str(core_package), '-e', str(package)], 180)
        command([sys.executable, '-m', 'cubit_mesh_export.install'], 180)
    if cfg['action'] != 'preflight':
        probe = "import json,pathlib,importlib.metadata as m,cubit_mesh_export as c;d=m.distribution('cubit-mesh-export');print(json.dumps(dict(version=c.__version__,file=str(pathlib.Path(c.__file__).resolve()),direct_url=json.loads(d.read_text('direct_url.json')))))"
        identity = json.loads(command([sys.executable, '-c', probe]).strip())
        expected_file = package / 'src/cubit_mesh_export/__init__.py'
        if (identity['version'] != cfg['version'] or Path(identity['file']).resolve() != expected_file.resolve()
                or identity['direct_url'].get('dir_info', {}).get('editable') is not True):
            raise RuntimeError('Wrong actual editable import: ' + str(identity))
        result['installed'] = identity
        core_probe = "import json,pathlib,importlib.metadata as m,cae_mcp_core as c;d=m.distribution('cae-mcp-core');print(json.dumps(dict(version=c.__version__,file=str(pathlib.Path(c.__file__).resolve()),direct_url=json.loads(d.read_text('direct_url.json')))))"
        core_identity = json.loads(command([sys.executable, '-c', core_probe]).strip())
        if (core_identity['version'] != cfg['core']['version'] or
                Path(core_identity['file']).resolve() != (core_package / 'src/cae_mcp_core/__init__.py').resolve() or
                core_identity['direct_url'].get('dir_info', {}).get('editable') is not True):
            raise RuntimeError('Wrong shared-core editable import: ' + str(core_identity))
        result['core_installed'] = core_identity
        command([sys.executable, '-m', 'cubit_mesh_export.mcp.server', '--selftest'], 120)
        result['mcp_selftest'] = True
        command([sys.executable, '-m', 'cubit_mesh_export.install', '--verify-only'], 120)
    if cfg['action'] == 'deploy':
        for module, arguments in [('smoke_test', ['--keep']),
                ('toolbar_smoke', ['--restarts', '2', '--keep', '--report-json', str(out / 'gui.json')])]:
            try:
                text = command([sys.executable, '-m', 'cubit_mesh_export.' + module, *arguments], 750)
            finally:
                if events:
                    recover_scratch(events[-1]['stdout'])
            result[module] = True
    result['passed'] = True
except Exception as exc:
    result['error'] = repr(exc)
finally:
    result['preserved_before'] = before
    result['preserved_after'] = preserved()
    result['unrelated_packages_unchanged'] = before == result['preserved_after']
    result['passed'] = result['passed'] and result['unrelated_packages_unchanged']
    result['events'] = events
    (out / (cfg['action'] + '.json')).write_text(json.dumps(result, indent=2), encoding='utf-8')
    print('CUBIT_DUAL_RESULT=' + json.dumps(result))
sys.exit(0 if result['passed'] else 1)
'''


def check_receipt(receipt, contract, target):
    expected = dict(schema=SCHEMA, target=target, source_sha=contract['source_sha'],
                    version=contract['version'], wheel_sha256=contract['wheel_sha256'],
                    core=contract['core'])
    return (all(receipt.get(k) == v for k, v in expected.items())
            and all(receipt.get(k) is True for k in
                    ('passed', 'unrelated_packages_unchanged', 'smoke_test', 'toolbar_smoke', 'mcp_selftest')))


def run(args):
    import tomllib
    contract = wheel_contract(args.wheel)
    contract['core'] = wheel_contract(args.core_wheel, 'cae-mcp-core')
    verify_published('cae-mcp-core', contract['core'])
    contract.update(schema=SCHEMA, source_sha=args.source_sha)
    root = Path(args.source_root_lab)
    metadata = tomllib.loads((root / 'packages/cubit-mesh-export/pyproject.toml').read_text(encoding='utf-8'))
    contract['dependencies'] = {name: next(d.split('==')[1] for d in metadata['project']['dependencies']
                                              if d.startswith(name + '=='))
                                for name in ('netgen-mesher', 'ngsolve')}
    if args.action != 'preflight':
        tag = 'cubit-mesh-export-v' + contract['version']
        p = subprocess.run(['git', '-c', 'safe.directory=' + str(root), '-C', str(root),
                            'rev-parse', tag + '^{}'], capture_output=True, text=True, check=True)
        if p.stdout.strip() != args.source_sha:
            raise ValueError('Tag does not identify selected source')
        verify_published('cubit-mesh-export', contract)
    output = Path(args.evidence_lab)
    output.mkdir(parents=True, exist_ok=True)
    # Both preflights finish before either install starts.
    phases = ['preflight', 'deploy'] if args.action == 'deploy' else [args.action]
    for phase in phases:
        for target in TARGETS:
            if phase == 'done':
                receipt = json.loads((output / target / 'deploy.json').read_text(encoding='utf-8'))
                if not check_receipt(receipt, contract, target):
                    raise ValueError('Missing or stale dual receipt: ' + target)
            cfg = dict(contract, target=target, action='verify' if phase == 'done' else phase,
                       source_root=args.source_root_lab if target == 'lab' else args.source_root_100,
                       output=str(Path(args.evidence_lab if target == 'lab' else args.evidence_100) / target))
            encoded = base64.b64encode(json.dumps(cfg).encode()).decode()
            command = [sys.executable, '-', encoded] if target == 'lab' else [
                'ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=10', '100', 'python', '-', encoded]
            p = subprocess.run(command, input=WORKER, text=True, encoding='utf-8',
                               errors='replace', capture_output=True, timeout=1200)
            lines = [line for line in p.stdout.splitlines() if line.startswith('CUBIT_DUAL_RESULT=')]
            if not lines:
                raise RuntimeError(target + ': no worker receipt\n' + p.stderr[-3000:])
            result = json.loads(lines[-1].partition('=')[2])
            print(f"{phase} {target}: {'PASS' if result['passed'] else 'FAIL'}", flush=True)
            if p.returncode or not result['passed']:
                raise RuntimeError(result.get('error', 'worker failed'))
    if args.action == 'done':
        (output / 'done.json').write_text(json.dumps(dict(contract, passed=True, targets=list(TARGETS)), indent=2))
        print('PASS release-dual: exact published exporter/core wheels, LAB/100 editable and GUI/export/MCP gates; Radia/radia-mcp unchanged')
    return 0
