"""Fail-closed workflow capability routing; stage3 exceptions are explicit debt."""
from __future__ import annotations
import argparse
import ast
import json
import hashlib
import contextlib
import io
from pathlib import Path
import re
import yaml

CAPABILITIES = {'ci-light', 'ci-native', 'gpu-rtx', 'gpu-a100'}
FORBIDDEN = re.compile(r'(?i)(?:matlab(?:engine)?|mexcuda|nvcc|simulink|cuda|cmake|msbuild|cl\.exe|build\.ps1)')
OWNERS = {'radia', 'radia-simulink', 'radia-mcp', 'radia-optuna', 'cubit-mesh-export', 'eqnedit64', 'shared'}
DENIED_IMPORTS = {'matlab', 'cupy', 'torch', 'ngsolve', 'netgen', 'cubit', '_radia_pybind'}


def workflows(root):
    result = {}
    for path in sorted((root / '.github/workflows').glob('*.y*ml')):
        for name, job in yaml.safe_load(path.read_text(encoding='utf-8')).get('jobs', {}).items():
            result[f'{path.name}:{name}'] = job
    return result


def light_sources(root, declaration):
    """Include ALL selector alternatives, nested conftests and static local imports.

    Dynamic file loaders in tests use literal .py paths; resolve those as well.
    This is conservative source inspection, not a Python sandbox.
    """
    tier = json.loads((root / 'tests/test_tier_manifest.json').read_text(encoding='utf-8'))
    paths = set(declaration['entrypoints'])
    name = declaration['test_profile']
    visited = set()
    while name:
        if name in visited:
            raise ValueError('cyclic profile')
        visited.add(name)
        profile = tier['profiles'][name]
        paths.update(profile['paths'])
        name = profile.get('extends')
    import importlib.util
    spec = importlib.util.spec_from_file_location('routing_tier', root / 'tools/run_test_tier.py')
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    runner.ROOT = root
    runner.MANIFEST = root / 'tests/test_tier_manifest.json'
    with contextlib.redirect_stdout(io.StringIO()):
        paths.update(runner.select_impact_tests(list(paths), None, profile_name='fast-contracts'))
    paths.update(['tests/conftest.py', 'tools/ci_runner_routing.py'])
    # Test collection may load intermediate conftests.
    for value in list(paths):
        for parent in (root / value).parents:
            if parent == root: break
            if (parent / 'conftest.py').exists():
                paths.add((parent / 'conftest.py').relative_to(root).as_posix())
    paths.update(declaration.get("audited_dynamic_sources", []))
    queue = list(paths)
    imports = set()
    while queue:
        value = queue.pop()
        file = root / value
        if not file.is_file(): raise ValueError(f'missing light dependency: {value}')
        if file.suffix != '.py': continue
        tree = ast.parse(file.read_text(encoding='utf-8-sig'), filename=value)
        for node in ast.walk(tree):
            names = []
            if isinstance(node, ast.Import): names = [item.name for item in node.names]
            if isinstance(node, ast.ImportFrom) and node.module: names = [node.module]
            for module in names:
                imports.add((value, module))
                relative = module.replace('.', '/')
                candidates = [root / (relative + '.py'), root / 'src' / (relative + '.py'), file.parent / (relative + '.py')]
                if value in declaration.get("deferred_or_mocked_imports", []):
                    continue  # Audited deferred imports / test doubles, bound by closure hash.
                for candidate in candidates:
                    if candidate.is_file():
                        path = candidate.relative_to(root).as_posix()
                        if path not in paths: paths.add(path); queue.append(path)
    return paths, imports


def audit_digest(root, job, paths):
    # Audited closure includes selector and source-only conftest. Any executable
    # change requires a deliberate F2 review, including dynamically loaded code.
    values = {path: (root / path).read_text(encoding='utf-8-sig').replace('\r\n', '\n')
              for path in sorted(paths)}
    values['test_selector_manifest'] = (root / 'tests/test_tier_manifest.json').read_text(encoding='utf-8')
    values['workflow_execution'] = json.dumps({'steps': job['steps'], 'env': job.get('env')}, sort_keys=True)
    return hashlib.sha256(json.dumps(values, sort_keys=True).encode()).hexdigest()


def validate(root, manifest=None, mode=None):
    manifest = manifest or json.loads((root / 'tools/ci_runner_routing.json').read_text(encoding='utf-8'))
    mode = mode or manifest['mode']
    if mode not in ('stage3', 'final'): raise ValueError('unknown routing mode')
    actual = workflows(root)
    if actual.keys() != manifest['jobs'].keys():
        raise ValueError(f'unclassified or removed jobs: {actual.keys() ^ manifest["jobs"].keys()}')
    for key, job in actual.items():
        rule = manifest['jobs'][key]
        if rule.get('owner') not in OWNERS: raise ValueError(f'missing owner: {key}')
        labels = job.get('runs-on')
        local = isinstance(labels, list) and 'self-hosted' in labels
        if rule['target'] == 'github-hosted':
            if local or not isinstance(labels, str) or '${{' in labels: raise ValueError(f'hosted route changed: {key}')
            continue
        if not local: raise ValueError(f'local route changed: {key}')
        caps = CAPABILITIES.intersection(labels)
        if mode == 'stage3' and rule.get('legacy_100'):
            if labels != ['self-hosted', 'Windows', 'X64', '100']: raise ValueError(f'legacy exception mismatch: {key}')
            continue
        expected = {'self-hosted', 'Windows', 'X64', rule['target']}
        if rule['target'] == 'gpu-rtx' or (mode == 'final' and rule['target'] != 'gpu-a100'): expected.add('100')
        if len(caps) != 1 or set(labels) != expected or len(labels) != len(expected):
            raise ValueError(f'capability mismatch: {key}')
        if rule['target'] == 'ci-light':
            declaration = manifest['light'].get(key)
            if not declaration: raise ValueError(f'missing light audit: {key}')
            if any(FORBIDDEN.search(x) for x in declaration['requirements']): raise ValueError(f'heavy requirement: {key}')
            for step in job.get('steps', []):
                run = step.get('run', '')
                if FORBIDDEN.search(run): raise ValueError(f'heavy light command: {key}')
                uses = step.get('uses', '')
                if uses and not uses.startswith(('actions/checkout@', 'actions/upload-artifact@')):
                    raise ValueError(f'unaudited light action: {uses}')
            paths, imports = light_sources(root, declaration)
            denied = [(path, module) for path, module in imports if module.split('.')[0] in DENIED_IMPORTS]
            # tests/conftest owns one explicit source-only guard, tested separately.
            denied = [(path, module) for path, module in denied
                      if module not in declaration.get('deferred_or_mocked_imports', {}).get(path, [])
                      and path != 'tests/conftest.py']
            if job.get('env', {}).get('RADIA_SOURCE_ONLY_CONTRACTS') != '1':
                raise ValueError('source-only import guard disabled')
            if denied: raise ValueError(f'heavy light dependency: {denied}')
            if declaration.get('audit_sha256') != audit_digest(root, job, paths):
                raise ValueError('light execution closure changed: repeat F2 audit')
    return {'jobs': len(actual), 'mode': mode}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['stage3', 'final'])
    args = parser.parse_args()
    try: print(json.dumps(validate(Path(__file__).resolve().parents[1], mode=args.mode)))
    except (ValueError, KeyError) as error: parser.exit(1, str(error) + '\n')
