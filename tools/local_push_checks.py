"""Bounded local push checks. Never imports numerical packages or runs native jobs."""
from __future__ import annotations
import ast
import importlib.util
import json
import re
from pathlib import Path
import subprocess
import sys

# Per-command budgets in seconds.  The shared precheck measured 255.76 s for the
# policy tests alone, so the former fixed 180 s limit failed healthy pushes.
COMMAND_BUDGETS_S = {'policy_lint': 600, 'compact_tests': 1500}

COMPACT_TESTS = ['tests/policy', 'tests/test_ci_execution_policy.py', 'tests/test_sparsesolv_ci_contract.py']

def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def positive_literal_conflicts(text, pattern, edited_lines=None):
    """Only positive membership assertions: negative/synthetic fixtures aren't demands."""
    tree = ast.parse(text)
    result = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Assert) and isinstance(node.test, ast.Compare):
            test = node.test
            if len(test.ops) == 1 and isinstance(test.ops[0], ast.In):
                value = test.left
                if isinstance(value, ast.Constant) and isinstance(value.value, str) and pattern.search(value.value):
                    if edited_lines is None or edited_lines.intersection(range(node.lineno, node.end_lineno + 1)):
                        result.append(node.lineno)
    return result

def minimal_import_findings(text, edited_lines=None):
    """Conservative all-scope direct imports, not a claim of transitive proof.

    When edited_lines is given, statically named imports count only if they touch an
    edited line, like positive_literal_conflicts.  Unresolved dynamic imports stay
    file-wide: an edit elsewhere can change the name they import.
    """
    tree = ast.parse(text)
    findings = []
    for node in ast.walk(tree):
        if edited_lines is not None and isinstance(node, (ast.Import, ast.ImportFrom, ast.Call)) and not edited_lines.intersection(range(node.lineno, node.end_lineno + 1)):
            if isinstance(node, (ast.Import, ast.ImportFrom)):
                continue
            if isinstance(node, ast.Call) and node.args and isinstance(node.args[0], ast.Constant) and isinstance(node.args[0].value, str):
                continue
        modules = []
        if isinstance(node, ast.Import): modules = [a.name for a in node.names]
        if isinstance(node, ast.ImportFrom): modules = [node.module or '']
        if any(m.split('.')[0] == 'radia' for m in modules):
            findings.append({'line': node.lineno, 'reason': 'CPU radia dependency in minimal MCP test'})
        if isinstance(node, ast.Call) and (isinstance(node.func, ast.Name) and node.func.id == '__import__' or isinstance(node.func, ast.Attribute) and node.func.attr == 'import_module'):
            if not node.args or not isinstance(node.args[0], ast.Constant) or not isinstance(node.args[0].value, str):
                findings.append({'line': node.lineno, 'reason': 'unresolved dynamic import requires lane review'})
            elif node.args[0].value.split('.')[0] == 'radia':
                findings.append({'line': node.lineno, 'reason': 'dynamic CPU radia dependency in minimal MCP test'})
    return findings

def git(root, *args):
    return subprocess.check_output(['git', '-c', f'safe.directory={root}', *args], cwd=root, text=True, encoding='utf-8').strip()

def changed_lines(root, base, candidate, path):
    """Do not promote historical assertions elsewhere in an edited file to new debt."""
    diff = git(root, 'diff', '--unified=0', f'{base}...{candidate}', '--', path)
    lines = set()
    for start, count in re.findall(r'^@@ -\d+(?:,\d+)? \+(\d+)(?:,(\d+))? @@', diff, re.M):
        lines.update(range(int(start), int(start) + (int(count) if count else 1)))
    return lines


def plan(root, base, candidate):
    if not re.fullmatch(r'[0-9a-f]{40}', candidate):
        raise ValueError('candidate must be a full immutable Git SHA')
    actual = git(root, 'rev-parse', 'HEAD')
    if actual != candidate or git(root, 'status', '--porcelain', '--untracked-files=no'):
        raise ValueError('local push checks require a clean tracked checkout at exact candidate SHA')
    impact = load(root / 'tools/sparsesolv_ci_impact.py', 'local_impact')
    try:
        git(root, 'cat-file', '-e', f'{base}^{{commit}}')
        raw = git(root, 'diff', '--name-status', '-z', '-M', f'{base}...{candidate}')
    except subprocess.CalledProcessError as exc:
        raise ValueError(f'cannot compare base {base}; fetch the remote base commit/history and retry') from exc
    paths = sorted(impact.changed_paths(raw)) if raw else []
    data = impact.manifest()
    def blobs(path):
        return git(root, 'show', f'{base}:{path}'), git(root, 'show', f'{candidate}:{path}')
    lanes = {lane: dict(zip(('required', 'reason'), impact.lane_required(lane, paths, blobs, data))) for lane in data['lanes']}
    selector = subprocess.run([sys.executable, str(root/'packages/radia-mcp/tools/select_ci_tests.py'), '--changed-files-json', '-', '--base', base, '--head', candidate], input=json.dumps(paths), cwd=root, text=True, encoding='utf-8', capture_output=True, check=True)
    mcp = json.loads(selector.stdout)
    policy = load(root/'tools/policy_lint.py', 'local_policy')
    findings = []
    for path in paths:
        file = root/path
        if not file.is_file() or not path.endswith('.py'): continue
        text = file.read_text(encoding='utf-8-sig')
        added_lines = changed_lines(root, base, candidate, path)
        for line in positive_literal_conflicts(text, policy.PUBLIC_RUNTIME_IDENTITY, added_lines):
            findings.append({'path': path, 'line': line, 'reason': 'positive assertion requires policy-banned runtime literal; review required'})
        if path.startswith('packages/radia-mcp/tests/'):
            findings.extend(dict(path=path, **item) for item in minimal_import_findings(text, added_lines))
    return {'candidate': candidate, 'base': base, 'changed_paths': paths, 'lanes': lanes, 'local_tests': COMPACT_TESTS, 'mcp_selection': mcp, 'static_findings': findings, 'import_audit_scope': 'changed MCP tests, direct/all-scope imports only; no transitive proof', 'heavy_execution': False}

def run(root, base, candidate):
    result = plan(root, base, candidate)
    print(json.dumps(result, ensure_ascii=True), flush=True)
    if result['static_findings']: return 1
    commands = [('policy_lint', [sys.executable, 'tools/policy_lint.py', '--quiet']),
                ('compact_tests', [sys.executable, '-m', 'pytest', *COMPACT_TESTS, '-q', '-o', 'addopts=', '-p', 'no:cacheprovider'])]
    for name, command in commands:
        budget = COMMAND_BUDGETS_S[name]
        try:
            completed = subprocess.run(command, cwd=root, timeout=budget)
        except subprocess.TimeoutExpired:
            print(f'local push check {name} exceeded its {budget} s budget', file=sys.stderr, flush=True)
            return 124
        if completed.returncode: return completed.returncode
    if git(root, 'rev-parse', 'HEAD') != candidate or git(root, 'status', '--porcelain', '--untracked-files=no'):
        raise ValueError('candidate checkout changed during local checks')
    return 0
