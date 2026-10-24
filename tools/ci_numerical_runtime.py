"""Verify a numerical plan before it can select executable pytest paths."""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess

from ci_numerical_plan import Tree, git, plan


def verified_plan(root, path):
    root = Path(root).resolve()
    value = json.loads(Path(path).read_text(encoding='utf-8'))
    head = git(root, 'rev-parse', 'HEAD').decode('ascii').strip()
    if value.get('head') != head:
        raise ValueError('numerical plan is not for the checked out HEAD')
    if os.environ.get('GITHUB_EVENT_NAME') and value.get('event') != os.environ['GITHUB_EVENT_NAME']:
        raise ValueError('numerical plan event differs from runner event')
    if os.environ.get('GITHUB_EVENT_PATH'):
        event = json.loads(Path(os.environ['GITHUB_EVENT_PATH']).read_text(encoding='utf-8'))
        if value['event'] == 'pull_request':
            expected_base = git(root, 'merge-base', event['pull_request']['base']['sha'], head).decode('ascii').strip()
            if value['base'] != expected_base: raise ValueError('plan base differs from runner event')
        elif value['event'] == 'push' and event.get('before') and set(event['before']) != {'0'}:
            if value['base'] != event['before']: raise ValueError('plan base differs from push event')
    expected = plan(root, value['base'], head, value['event'])
    if value != expected:
        raise ValueError('numerical plan differs from recomputed immutable selection')
    tree = Tree(root, head)
    paths = {p for row in value['cases'] if row['selected']
             for p in row['dependencies'] if p in tree.files}
    paths.update(p for p in ('tools/ci_numerical_plan.py', 'tools/ci_numerical_runtime.py',
                            'tools/ci_numerical_inputs.json', 'tools/run_test_tier.py',
                            'tools/validation_source_gate.py') if p in tree.files)
    paths = sorted(paths)
    # Hash every current input with Git's clean filters in one process. This
    # accepts a normal CRLF checkout but does not trust stat-cache/skip-worktree.
    for path in paths:
        actual = root / path
        if (tree.files[path][0] not in {'100644', '100755'} or actual.is_symlink()
                or not actual.resolve().is_relative_to(root) or not actual.is_file()):
            raise ValueError('missing/linked numerical input: ' + path)
    request = ''.join(json.dumps(str(root / p), ensure_ascii=False) + '\n' for p in paths).encode('utf-8')
    raw = subprocess.run(['git', '-c', f'safe.directory={root}', 'hash-object', '--stdin-paths'],
                         cwd=root, input=request, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, check=True).stdout
    observed = raw.decode('ascii').splitlines()
    if len(observed) != len(paths): raise ValueError('incomplete numerical source hash batch')
    for path, oid in zip(paths, observed):
        if oid != tree.files[path][2]:
            raise ValueError('numerical input differs from its checked out commit: ' + path)
    return value


def selected_paths(value, profile):
    paths = sorted({row['path'] for row in value['cases']
                    if row['profile'] == profile and row['selected']})
    if not paths:
        raise ValueError('empty numerical execution is not a successful gate: ' + profile)
    return paths


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--profile', required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    try:
        value = verified_plan(args.repo, args.plan)
        result = {'status': 'verified', 'paths': selected_paths(value, args.profile)}
        code = 0
    except Exception as error:
        result = {'status': 'blocked', 'error': str(error)}; code = 2
    args.out.write_text(json.dumps(result, indent=2), encoding='utf-8')
    raise SystemExit(code)
