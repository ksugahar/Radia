"""Executable input closure of the full AMS numerical gate.

The native build compiles Radia, SparseSolv, PEEC and axifem. Their source,
headers and toolchain remain indivisible inputs. Python imports are followed
from the actual test profiles; explicit seeds cover file/AST-based loaders.
"""
from functools import lru_cache
import ast
import json
import re
from pathlib import Path

def python_inputs(root):
    return set(_python_inputs(root))


@lru_cache(maxsize=8)
def _python_inputs(root):
    profiles = json.loads((root / 'tests/test_tier_manifest.json').read_text(encoding='utf-8'))['profiles']
    pending = ['tests/conftest.py', 'src/radia/__init__.py', 'src/radia/bh_law.py',
               'src/radia/vim/_nonlinear.py',
               'validation_test/induction_heating/aphi_round_wire_validation.py']
    for profile in ('sparsesolv', 'solver-numerics', 'solver-heavy', 'native-build-contracts'):
        visited = set()
        while profile:
            if profile in visited:
                raise ValueError('cyclic numerical profile inheritance')
            visited.add(profile)
            pending.extend(profiles[profile].get('paths', []))
            profile = profiles[profile].get('extends')
    source_roots = [root / 'src', *(root / 'packages' / package / 'src' for package in
        ('radia-simulink', 'radia-mcp', 'cubit-mesh-export')), root]
    # Conservative floor: dynamically imported solver modules and file-loaded data
    # cannot all be proven by static Python import traversal.
    floor = {p.relative_to(root).as_posix() for p in (root / 'src/radia').rglob('*.py')}
    for folder in ('tests/fixtures', 'validation_test/panels/fixtures'):
        floor.update(p.relative_to(root).as_posix() for p in (root / folder).rglob('*')
                     if p.is_file() and p.suffix not in ('.md', '.pyc') and '__pycache__' not in p.parts)
    seen = set()
    while pending:
        value = pending.pop()
        if value in seen:
            continue
        seen.add(value)
        file = root / value
        if file.suffix != '.py' or not file.is_file():
            continue
        # Importing a module executes every ancestor package initializer too.
        # Traverse those files as inputs, including their own imports.
        for source_root in source_roots:
            if file.is_relative_to(source_root):
                parent = file.parent
                while parent != source_root:
                    initializer = parent / '__init__.py'
                    if initializer.is_file():
                        pending.append(initializer.relative_to(root).as_posix())
                    parent = parent.parent
                break
        tree = ast.parse(file.read_text(encoding='utf-8-sig'))
        for node in ast.walk(tree):
            # CLI subprocesses and file-based test loaders do not use imports.
            if (isinstance(node, ast.Constant) and isinstance(node.value, str)
                    and len(node.value) < 240 and '\n' not in node.value
                    and re.search(r'\.[A-Za-z][A-Za-z0-9_-]{0,15}$', node.value)):
                for base in (root, file.parent, root / 'src/radia/panels', root / 'validation_test/panels'):
                    candidate = base / node.value
                    if (candidate.is_relative_to(root) and candidate.is_file()
                            and '.git' not in candidate.relative_to(root).parts
                            and candidate != root / '.gitignore'):
                        pending.append(candidate.relative_to(root).as_posix())
            names = []
            if isinstance(node, ast.Import):
                names = [n.name for n in node.names]
            elif isinstance(node, ast.ImportFrom):
                module = node.module or ''
                if node.level:
                    for source_root in source_roots:
                        if file.is_relative_to(source_root):
                            parts = list(file.relative_to(source_root).with_suffix('').parts[:-1])
                            if node.level > 1:
                                parts = parts[:-(node.level - 1)]
                            module = '.'.join(parts + ([module] if module else []))
                            break
                names = [module] + [module + '.' + n.name for n in node.names]
            for name in names:
                for base in (*source_roots, file.parent, root / 'src/radia/panels'):
                    candidate = base / name.replace('.', '/')
                    for path in (candidate.with_suffix('.py'), candidate / '__init__.py'):
                        if path.is_file() and path.is_relative_to(root):
                            pending.append(path.relative_to(root).as_posix())
    return frozenset(seen | floor)


def required(root, changed, contract, data):
    native = data['lanes']['sparsesolv']['pr_heavy_inputs']
    if any(contract.numerical_environment_input(path)
           or path.startswith(('src/radia/', 'tests/fixtures/', 'validation_test/panels/fixtures/'))
           or any(contract.matches(path, pattern) for pattern in native) for path in changed):
        return True
    return bool(set(changed) & python_inputs(root))
