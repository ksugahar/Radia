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
OWNERS = {'radia', 'radia-simulink', 'radia-mcp', 'radia-optuna', 'cubit-mesh-export', 'shared'}
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
    import importlib.util
    spec = importlib.util.spec_from_file_location('routing_tier', root / 'tools/run_test_tier.py')
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    runner.ROOT = root
    runner.MANIFEST = root / 'tests/test_tier_manifest.json'
    tier = runner.validate_tier_data(root, runner.read_tier_data(runner.MANIFEST))
    # Ownership DATA has its own existing contract; its reader remains sealed.
    spec = importlib.util.spec_from_file_location('routing_ownership', root / 'tools/ci_distribution_contract.py')
    ownership = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ownership)
    owner_data = runner.read_tier_data(root / 'tools/ci_distribution_manifest.json')
    try:
        expected_keys = {'schema', 'distributions', 'infrastructure_owner', 'rule_semantics', 'root_content_diff', 'rules', 'shared_include_closures', 'lanes', 'migration_exceptions'}
        if not isinstance(owner_data, dict) or set(owner_data) != expected_keys:
            raise ValueError('expected supported ownership DATA keys only')
        ownership.validate_schema(owner_data)
        assigned = owner_data['root_content_diff']['attribution']['tests/test_tier_manifest.json']
        if not set(tier['profiles']) <= set(assigned):
            raise ValueError(f"unclassified profile owners: {sorted(set(tier['profiles']) - set(assigned))}")
        for profile in tier['profiles'].values():
            for path in profile['paths']:
                ownership.owner_rule(owner_data, path)
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"ownership DATA: {exc}") from exc
    required = {'tests/test_test_tier_policy.py', 'tests/test_run_test_tier.py'}
    if not required <= set(runner.profile_paths(tier, 'repo-contracts')):
        raise ValueError(f"required tier/coverage contracts missing: {sorted(required - set(runner.profile_paths(tier, 'repo-contracts')))}")
    paths = set(declaration['entrypoints'])
    selected, _ = runner.load_profile(declaration['test_profile'])
    paths.update(selected)
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


# Only these explicit, separately validated DATA documents are excluded.
VALIDATED_DATA = {'tests/test_tier_manifest.json', 'tools/ci_distribution_manifest.json'}


def audit_digest(root, job, paths, declaration):
    sources = {path: (root / path).read_text(encoding='utf-8-sig').replace('\r\n', '\n')
              for path in sorted(paths) if path not in VALIDATED_DATA}
    tier = json.loads((root / 'tests/test_tier_manifest.json').read_text(encoding='utf-8'))
    policy = {key: value for key, value in declaration.items() if key != 'audit_sha256'}
    policy['audited_dynamic_sources'] = [path for path in policy.get('audited_dynamic_sources', []) if path not in VALIDATED_DATA]
    values = {'sources': sources}
    values['audit_format'] = 'executable-closure-v2'
    values['light_declaration'] = policy
    values['active_profile_budget'] = tier['profiles'][declaration['test_profile']].get('max_elapsed_seconds')
    values['workflow_execution'] = job
    return hashlib.sha256(json.dumps(values, sort_keys=True, allow_nan=False).encode()).hexdigest()


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
            if declaration.get('audit_sha256') != audit_digest(root, job, paths, declaration):
                raise ValueError('light execution closure changed: repeat F2 audit')
    return {'jobs': len(actual), 'mode': mode}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['stage3', 'final'])
    args = parser.parse_args()
    try: print(json.dumps(validate(Path(__file__).resolve().parents[1], mode=args.mode)))
    except (ValueError, KeyError) as error: parser.exit(1, str(error) + '\n')
