"""Executable input closure of the compact AMS PR gate (not the nightly suite).

The native build compiles Radia, SparseSolv, PEEC and axifem. Their source,
headers and toolchain remain indivisible inputs. Python imports are followed
from the actual test profiles; explicit seeds cover file/AST-based loaders.
"""
import ast
import json
from pathlib import Path


def python_inputs(root):
    profiles = json.loads((root / 'tests/test_tier_manifest.json').read_text(encoding='utf-8'))['profiles']
    pending = profiles['sparsesolv']['paths'] + profiles['solver-pr']['paths'] + profiles['solver-heavy']['paths'] + [
        'tests/conftest.py', 'src/radia/__init__.py', 'src/radia/bh_law.py',
        'src/radia/vim/_nonlinear.py',
        'validation_test/induction_heating/aphi_round_wire_validation.py',
    ]
    source_roots = [root / 'src', *(root / 'packages' / package / 'src' for package in
        ('radia-simulink', 'radia-mcp', 'cubit-mesh-export'))]
    seen = set()
    while pending:
        value = pending.pop()
        if value in seen:
            continue
        seen.add(value)
        file = root / value
        if file.suffix != '.py' or not file.is_file():
            continue
        tree = ast.parse(file.read_text(encoding='utf-8-sig'))
        for node in ast.walk(tree):
            # CLI subprocesses and file-based test loaders do not use imports.
            if isinstance(node, ast.Constant) and isinstance(node.value, str) and node.value.endswith('.py'):
                for base in (root, file.parent, root / 'src/radia/panels', root / 'validation_test/panels'):
                    candidate = base / node.value
                    if candidate.is_file() and candidate.is_relative_to(root):
                        pending.append(candidate.relative_to(root).as_posix())
            names = []
            if isinstance(node, ast.Import):
                names = [n.name for n in node.names]
            elif isinstance(node, ast.ImportFrom):
                module = node.module or ''
                if node.level and value.startswith('src/radia/'):
                    parts = value[:-3].split('/')[:-node.level]
                    module = '.'.join(parts[1:] + ([module] if module else []))
                names = [module] + [module + '.' + n.name for n in node.names]
            for name in names:
                for base in (*source_roots, file.parent, root / 'src/radia/panels'):
                    candidate = base / name.replace('.', '/')
                    for path in (candidate.with_suffix('.py'), candidate / '__init__.py'):
                        if path.is_file() and path.is_relative_to(root):
                            pending.append(path.relative_to(root).as_posix())
    return seen


def required(root, changed, contract, data):
    native = data['lanes']['sparsesolv']['pr_heavy_inputs']
    python = python_inputs(root)
    return any(path in python or any(contract.matches(path, pattern) for pattern in native) for path in changed)
