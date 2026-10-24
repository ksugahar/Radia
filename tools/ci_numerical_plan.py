"""Plan numerical CI consumers from immutable Git inputs; never execute them.

S1 is diagnostic only. Existing workflow execution and source receipts remain
authoritative until the separately reviewed S2 connects this plan to runners.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
from functools import lru_cache
import json
import os
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
import tomllib

SCHEMA = 'radia.numerical-ci-plan.v1'
TIER = 'tests/test_tier_manifest.json'
DECLARATIONS = 'tools/ci_numerical_inputs.json'
PROFILES = ('sparsesolv', 'solver-numerics', 'solver-heavy', 'native-build-contracts')
GMSH_TESTS = (
    'validation_test/gmsh_post/test_msh42_reference.py',
    'tests/gmsh_post/test_native_field_handoff.py',
    'tests/test_gmsh_post_export_roundtrip.py',
    'tests/test_gmsh_post_export_orientation.py',
    'tests/test_export_filaments_msh.py', 'tests/test_particle_tracks_msh_export.py',
    'tests/test_gmsh_only_output_policy.py', 'tests/test_motor_hdiv_reduced.py',
)
NATIVE_PREFIXES = ('src/core/', 'src/lib/', 'src/ext/', 'cmake/')
GLOBAL_INPUTS = {'CMakeLists.txt', 'Build.ps1', 'MANIFEST.in', TIER,
                 'tools/ci_numerical_plan.py', 'tools/ci_numerical_runtime.py', DECLARATIONS, 'tools/run_test_tier.py',
                 'tools/ci_ams_inputs.py', 'tools/sparsesolv_ci_impact.py',
                 'tools/ci_distribution_manifest.json', 'tools/ci_distribution_contract.py',
                 'tools/native_build_provenance.ps1', 'tools/validation_source_gate.py'}
EXTERNAL = set(sys.stdlib_module_names) | {
    'numpy', 'scipy', 'ngsolve', 'netgen', 'pytest', 'gmsh', 'h5py', 'psutil',
    'matplotlib', 'sympy', 'yaml', 'pefile', 'packaging', 'setuptools', 'pybind11',
    'threadpoolctl', 'pandas', 'pydantic', 'mcp', 'fastmcp', 'sklearn', 'numba',
    'torch', 'kan', 'build123d', 'OCP', 'cubit', 'netCDF4', 'PySide6',
    'skimage', 'trimesh',  # Upstream optional CAD/inspection packages; not local code.
}
NATIVE_MODULES = {'_radia_pybind', 'sparsesolv_ngsolve', 'peec_matrices', 'axifem',
                  'cubit_mesh_curver'}


def git(root, *args):
    return subprocess.check_output(['git', '-c', f'safe.directory={root}', *args],
                                   cwd=root, stderr=subprocess.PIPE)


def unique_pairs(items):
    result = {}
    for key, value in items:
        if key in result:
            raise ValueError('duplicate data key: ' + key)
        result[key] = value
    return result


class Tree:
    def __init__(self, root, rev, declarations=None):
        self.root = root
        self.declarations = declarations or {}
        self.sha = git(root, 'rev-parse', '--verify', '--end-of-options', rev + '^{commit}').decode('utf-8').strip()
        self.files = {}
        for row in git(root, 'ls-tree', '-r', '-z', self.sha).split(b'\0'):
            if row:
                metadata, name = row.split(b'\t', 1)
                mode, kind, oid = metadata.decode('utf-8').split()
                self.files[name.decode('utf-8')] = (mode, kind, oid)
        self.blobs = {}
        self.parsed_nodes = 0
        self.roots = ['src', 'src/radia', *sorted({p.split('/src/', 1)[0] + '/src'
                                     for p in self.files if p.startswith('packages/') and '/src/' in p}), '']

    def load_blobs(self):
        """One immutable batch per side; exclude binary payloads and unrelated data."""
        oids = sorted({oid for path, (mode, kind, oid) in self.files.items()
                       if kind == 'blob' and mode in {'100644', '100755'}
                       and (path.endswith('.py') or path.endswith('pyproject.toml') or path == TIER)})
        request = ''.join(oid + '\n' for oid in oids).encode('ascii')
        raw = subprocess.run(['git', '-c', f'safe.directory={self.root}', 'cat-file', '--batch'],
                             input=request, cwd=self.root, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, check=True).stdout
        offset = 0
        for expected in oids:
            end = raw.index(b'\n', offset)
            oid, kind, size = raw[offset:end].decode('ascii').split()
            if oid != expected or kind != 'blob': raise ValueError('unexpected batch object')
            start = end + 1; finish = start + int(size)
            if raw[finish:finish + 1] != b'\n': raise ValueError('truncated batch object')
            self.blobs[oid] = raw[start:finish]; offset = finish + 1
        if offset != len(raw): raise ValueError('unexpected trailing batch bytes')

    def read(self, path):
        mode, kind, oid = self.files[path]
        if mode not in {'100644', '100755'} or kind != 'blob':
            raise ValueError('non-regular source: ' + path)
        if not self.blobs: self.load_blobs()
        if oid not in self.blobs: raise ValueError('source outside batch contract: ' + path)
        return self.blobs[oid]

    def profiles(self):
        data = json.loads(self.read(TIER), object_pairs_hook=unique_pairs)
        if set(data) != {'schema', 'profiles', 'impact_rules'} or data['schema'] != 'radia.test-tier-manifest.v1':
            raise ValueError('unknown tier schema')
        def resolve(name, seen=()):
            if name in seen: raise ValueError('cyclic profile')
            item = data['profiles'][name]
            paths = resolve(item['extends'], (*seen, name)) if item.get('extends') else []
            for path in item['paths']:
                if not isinstance(path, str) or '\\' in path or '..' in PurePosixPath(path).parts or path not in self.files:
                    raise ValueError('missing/unsafe tier path: ' + str(path))
            return paths + item['paths']
        result = {name: resolve(name) for name in PROFILES}
        if any(not paths or len(paths) != len(set(paths)) for paths in result.values()):
            raise ValueError('empty/duplicate profile paths')
        return result

    def module(self, name, local):
        relative = name.replace('.', '/')
        candidates = []
        for prefix in [*self.roots, str(PurePosixPath(local).parent), 'src/radia/panels']:
            base = (prefix + '/' if prefix and prefix != '.' else '') + relative
            candidates += [base + '.py', base + '/__init__.py']
        return next((p for p in candidates if p in self.files), None)

    def literal_path(self, value, local):
        for base in ['', str(PurePosixPath(local).parent), 'src/radia/panels', 'validation_test/panels']:
            candidate = ((base + '/') if base != '.' and base else '') + value
            if candidate in self.files: return candidate
        return None

    @lru_cache(maxsize=None)
    def edges(self, path):
        """Parse each reached module once into the shared graph for this side."""
        pending = []; unknown = set(); native = set(); components = set()
        if path not in self.files:
            if path != 'tests/conftest.py': unknown.add('missing dependency: ' + path)
            return pending, unknown, native, components
        if not path.endswith('.py'): return pending, unknown, native, components
        self.parsed_nodes += 1
        try: tree = ast.parse(self.read(path).decode('utf-8-sig'))
        except (ValueError, SyntaxError, UnicodeError) as exc:
            unknown.add(path + ': ' + type(exc).__name__); return pending, unknown, native, components
        assignments = {}
        for assignment in ast.walk(tree):
            if isinstance(assignment, ast.Assign):
                for target in assignment.targets:
                    if isinstance(target, ast.Name): assignments.setdefault(target.id, []).append(assignment.value)
        parent = PurePosixPath(path).parent
        for ancestor in [parent, *parent.parents]:
            initializer = str(ancestor / '__init__.py')
            conftest = str(ancestor / 'conftest.py')
            if initializer in self.files: pending.append(initializer)
            if conftest in self.files: pending.append(conftest)
        for node in ast.walk(tree):
            names = []
            if isinstance(node, ast.Import): names = [alias.name for alias in node.names]
            if isinstance(node, ast.ImportFrom):
                module = node.module or ''
                if node.level:
                    base = parent
                    for _ in range(node.level - 1): base = base.parent
                    names = [str(base / module.replace('.', '/')).replace('/', '.')]
                    # Resolve relative imports without assuming a particular package root.
                    stem = str(base / module.replace('.', '/'))
                    found = False
                    for candidate in [stem + '.py', stem + '/__init__.py']:
                        if candidate in self.files: pending.append(candidate); found = True
                    if module.rsplit('.', 1)[-1] in NATIVE_MODULES:
                        native.add(module.rsplit('.', 1)[-1]); found = True
                    for alias in node.names:
                        candidate = stem.rstrip('/') + '/' + alias.name + '.py'
                        if candidate in self.files: pending.append(candidate); found = True
                        if alias.name in NATIVE_MODULES: native.add(alias.name); found = True
                    if not found: unknown.add('unresolved relative import: ' + path + ':' + module)
                    continue
                names = [module] + [module + '.' + a.name for a in node.names if a.name != '*']
            for name in names:
                candidate = self.module(name, path)
                if candidate: pending.append(candidate)
                elif name.rsplit('.', 1)[-1] in NATIVE_MODULES:
                    native.add(name.rsplit('.', 1)[-1])
                elif name.split('.')[0] not in EXTERNAL:
                    # A from-import can name an attribute of a resolved local module.
                    if not any(self.module('.'.join(name.split('.')[:i]), path)
                               for i in range(1, len(name.split('.')))):
                        unknown.add('unresolved import: ' + name + ' in ' + path)
                if name.split('.')[0] == 'gmsh': components.add('gmsh')
                if name.split('.')[0] in {'ngsolve', 'netgen'}: components.add('solver')
            if isinstance(node, ast.Constant) and isinstance(node.value, str):
                value = node.value
                # Only executable call sites are treated as code loaders below.
                # Docstrings, knowledge and generated template fragments are data.
                if len(value) < 240 and '\n' not in value and not value.endswith('.py'):
                    for base in ['', str(parent), 'src/radia/panels', 'validation_test/panels']:
                        candidate = ((base + '/') if base != '.' and base else '') + value
                        if candidate in self.files: pending.append(candidate)
            if isinstance(node, ast.Call):
                name = ast.unparse(node.func)
                if name.endswith(('import_module', 'importorskip', 'run_module', 'find_spec')) or name == '__import__':
                    if node.args and isinstance(node.args[0], ast.Constant) and isinstance(node.args[0].value, str):
                        import_name = node.args[0].value
                        candidate = self.module(import_name, path)
                        if import_name.startswith('.'):
                            relative_path = str(parent / (import_name.lstrip('.').replace('.', '/') + '.py'))
                            candidate = relative_path if relative_path in self.files else None
                        if candidate: pending.append(candidate)
                        elif node.args[0].value.rsplit('.', 1)[-1] in NATIVE_MODULES:
                            native.add(node.args[0].value.rsplit('.', 1)[-1])
                        elif node.args[0].value.split('.')[0] not in EXTERNAL:
                            unknown.add('unresolved dynamic import: ' + path)
                        if node.args[0].value.split('.')[0] in EXTERNAL: components.add('package:' + node.args[0].value.split('.')[0])
                        if node.args[0].value.split('.')[0] == 'gmsh': components.add('gmsh')
                        if node.args[0].value.split('.')[0] in {'ngsolve', 'netgen'}: components.add('solver')
                    else: unknown.add('dynamic import needs declaration: ' + path)
                if name.endswith(('spec_from_file_location', 'run_path', 'SourceFileLoader')) or name in {'eval', 'exec'}:
                    # File execution is not proved complete by ordinary imports.
                    unknown.add('opaque file/code loader: ' + path + ':' + name)
                if name.startswith('subprocess.') and name.rsplit('.', 1)[-1] in {'run', 'Popen', 'call', 'check_call', 'check_output'}:
                    expression = node.args[0] if node.args else None
                    if isinstance(expression, ast.Name):
                        values = assignments.get(expression.id, [])
                        expression = values[0] if len(values) == 1 else None
                    if isinstance(expression, ast.BinOp) and isinstance(expression.op, ast.Add): expression = expression.left
                    # Only a fixed Python script entry is statically bounded.
                    # Dynamic -c programs and arbitrary executables require review.
                    script = None
                    if isinstance(expression, (ast.List, ast.Tuple)) and len(expression.elts) >= 2:
                        expressions = [expression.elts[1]]; visited_names = set()
                        while expressions:
                            for literal in ast.walk(expressions.pop()):
                                if isinstance(literal, ast.Name) and literal.id not in visited_names:
                                    visited_names.add(literal.id)
                                    expressions.extend(assignments.get(literal.id, []))
                                if isinstance(literal, ast.Constant) and isinstance(literal.value, str) and literal.value.endswith('.py'):
                                    candidate_script = self.literal_path(literal.value, path)
                                    if candidate_script:
                                        script = candidate_script; pending.append(script)
                    if script and not any(kw.arg == 'shell' for kw in node.keywords): pending.append(script)
                    else: unknown.add('subprocess entry needs declaration: ' + path)
        binding = self.declarations.get('bindings', {}).get(path)
        if binding:
            if hashlib.sha256(self.read(path)).hexdigest() not in binding['sha256']:
                unknown.add('stale dynamic declaration: ' + path)
            else:
                if not binding.get('entries'):
                    unknown.difference_update(binding['resolves'])
                for dep in binding['dependencies']:
                    if dep not in self.files: unknown.add('declared dependency missing: ' + dep)
                    pending.append(dep)
                for prefix in binding.get('fixture_prefixes', []):
                    pending.extend(p for p in self.files if p.startswith(prefix))
        pending.extend(self.declarations.get('consumer_inputs', {}).get(path, []))
        return pending, unknown, native, components

    @lru_cache(maxsize=None)
    def closure(self, entry):
        pending = [entry, 'tests/conftest.py']; seen = set()
        unknown = set(); native = set(); components = set()
        while pending:
            path = pending.pop()
            if path in seen: continue
            seen.add(path)
            edges, issues, modules, uses = self.edges(path)
            pending.extend(edges); unknown.update(issues)
            native.update(modules); components.update(uses)
        for path in seen:
            binding = self.declarations.get('bindings', {}).get(path, {})
            if binding.get('entries', {}).get(entry) and entry in self.files:
                if (hashlib.sha256(self.read(path)).hexdigest() in binding['sha256']
                        and hashlib.sha256(self.read(entry)).hexdigest() in binding['entries'][entry]):
                    unknown.difference_update(binding['resolves'])
        return seen, unknown, native, components


def changes(root, base, head):
    fields = git(root, 'diff', '--name-status', '-z', '--find-renames', base, head).decode('utf-8').rstrip('\0').split('\0')
    paths = set(); i = 0
    if fields == ['']: return paths
    while i < len(fields):
        status = fields[i]; i += 1
        if not re.fullmatch(r'[ACDMRTUXB][0-9]*', status): raise ValueError('unknown diff status')
        count = 2 if status[0] in 'RC' else 1
        if i + count > len(fields): raise ValueError('truncated diff')
        paths.update(fields[i:i+count]); i += count
    return paths


def dependency_change(old, new):
    """Prove a Gmsh-only metadata change; everything else stays conservative."""
    left, right = tomllib.loads(old.decode('utf-8')), tomllib.loads(new.decode('utf-8'))
    def extract(data):
        project = data.get('project', {})
        project.pop('version', None)
        deps = project.pop('dependencies', [])
        extras = project.pop('optional-dependencies', {})
        entries = {('base', dep) for dep in deps}
        entries |= {(extra, dep) for extra, values in extras.items() for dep in values}
        return data, entries
    a, deps_a = extract(left); b, deps_b = extract(right)
    delta = deps_a ^ deps_b
    if a == b and delta and all(re.match(r'^gmsh\s*(?:[<>=!~;\[]|$)', dep, re.I) for _, dep in delta):
        return 'gmsh'
    return 'solver-or-unknown'


def plan(root, base, head, event):
    policy_sha = git(root, 'rev-parse', 'HEAD').decode('ascii').strip()
    policy_tree = Tree(root, policy_sha)
    declarations = {}
    declaration_blob = None
    if DECLARATIONS in policy_tree.files:
        declaration_blob = policy_tree.files[DECLARATIONS][2]
        declarations = json.loads(git(root, 'cat-file', 'blob', declaration_blob), object_pairs_hook=unique_pairs)
        if set(declarations) != {'schema', 'bindings', 'consumer_inputs'} or declarations['schema'] != 'radia.numerical-inputs.v1':
            raise ValueError('invalid numerical declaration schema')
        def safe_declared_path(value):
            if (not isinstance(value, str) or not value or value.startswith('/')
                    or ':' in value or '\\' in value or '..' in PurePosixPath(value).parts):
                raise ValueError('unsafe declared dependency')
        for path, binding in declarations['bindings'].items():
            if set(binding) - {'sha256', 'resolves', 'dependencies', 'fixture_prefixes', 'reason', 'entries'} or not binding.get('reason'):
                raise ValueError('invalid dynamic declaration: ' + path)
            if not binding.get('sha256') or any(not re.fullmatch('[0-9a-f]{64}', h) for h in binding['sha256']):
                raise ValueError('unsealed dynamic declaration: ' + path)
            for dep in [path, *binding['dependencies'], *binding.get('fixture_prefixes', [])]:
                safe_declared_path(dep)
            for entry, hashes in binding.get('entries', {}).items():
                safe_declared_path(entry)
                if not hashes or any(not re.fullmatch('[0-9a-f]{64}', h) for h in hashes):
                    raise ValueError('unsealed entry scope')
            for issue in binding['resolves']:
                if not (issue.endswith(': ' + path) or ': ' + path + ':' in issue):
                    raise ValueError('dynamic declaration cannot hide another source issue')
        for consumer, paths in declarations['consumer_inputs'].items():
            safe_declared_path(consumer)
            for path in paths: safe_declared_path(path)
    before, after = Tree(root, base, declarations), Tree(root, head, declarations)
    changed = changes(root, before.sha, after.sha)
    old_profiles, new_profiles = before.profiles(), after.profiles()
    reasons = []; global_change = event != 'pull_request' or not changed; gmsh_change = False
    for path in sorted(changed):
        if path.endswith('pyproject.toml'):
            try: component = dependency_change(before.read(path), after.read(path))
            except (KeyError, ValueError, UnicodeError): component = 'solver-or-unknown'
            if component == 'gmsh': gmsh_change = True
            else: global_change = True; reasons.append('dependency/metadata: ' + path)
        elif path in GLOBAL_INPUTS or path.startswith(NATIVE_PREFIXES) or path.startswith('.github/workflows/'):
            global_change = True; reasons.append('native/build/control input: ' + path)
    roots = {profile: sorted(set(old_profiles[profile]) | set(new_profiles[profile])) for profile in PROFILES}
    for snapshot in [before, after]:
        for path in GMSH_TESTS:
            if path in snapshot.files: roots['solver-numerics'] = sorted(set(roots['solver-numerics']) | {path})
    rows = []; known = set(); native_required = set()
    for profile, entries in roots.items():
        for entry in entries:
            dependencies = set(); unknown = set(); native = set(); components = set()
            for snapshot in [before, after]:
                if entry in snapshot.files:
                    seen, issues, modules, uses = snapshot.closure(entry)
                    dependencies |= seen; unknown |= issues; native |= modules; components |= uses
            known |= dependencies
            touched = sorted(changed & dependencies)
            selected = global_change or bool(touched) or bool(unknown) or (gmsh_change and ('gmsh' in components or entry in GMSH_TESTS))
            why = reasons[:] if global_change else []
            if event != 'pull_request': why.append('preserve main/manual/release execution')
            if not changed: why.append('empty diff; conservative execution')
            if touched: why.append('changed dependency: ' + ', '.join(touched))
            if unknown: why.append('unresolved case closure; conservative execution')
            if gmsh_change and ('gmsh' in components or entry in GMSH_TESTS): why.append('Gmsh runtime consumer')
            if selected: native_required |= native
            rows.append({'profile': profile, 'path': entry, 'selected': bool(selected),
                         'reasons': why, 'unknowns': sorted(unknown),
                         'dependencies': sorted(dependencies), 'native_modules': sorted(native), 'components': sorted(components),
                         'dependency_blobs': {p: sorted({s.files[p][2] for s in (before, after) if p in s.files})
                                              for p in sorted(dependencies)},
                         'closure_status': 'unknown' if unknown else 'static graph; not execution proof'})
    unclassified = sorted(p for p in changed if p not in known and not p.endswith(('.md', '.rst'))
                          and p not in GLOBAL_INPUTS and not p.startswith(NATIVE_PREFIXES)
                          and not p.endswith('pyproject.toml') and not p.startswith('.github/workflows/'))
    if unclassified:
        for row in rows:
            row['selected'] = True; row['reasons'].append('unclassified input; full fallback')
            native_required.update(row['native_modules'])
    if native_required:
        for row in rows:
            if row['profile'] == 'native-build-contracts':
                row['selected'] = True
                row['reasons'].append('verify selected native build')
    return {'schema': SCHEMA, 'status': 'planned', 'qualification': False, 'execution_wired': False,
            'base': before.sha, 'head': after.sha, 'event': event,
            'declaration_identity': {'commit': policy_sha, 'blob': declaration_blob}, 'changed_paths': sorted(changed),
            'unclassified_inputs': unclassified, 'cases': rows,
            'native_build_targets': sorted(native_required),
            'ams_paths': sorted({r['path'] for r in rows if r['selected'] and r['profile'] == 'sparsesolv'}),
            'python_test_paths': sorted({r['path'] for r in rows if r['selected'] and r['profile'] in ('solver-numerics', 'native-build-contracts')}),
            'heavy_case_paths': sorted({r['path'] for r in rows if r['selected'] and r['profile'] == 'solver-heavy'}),
            'gmsh_reference_paths': sorted({r['path'] for r in rows if r['selected'] and r['path'] in GMSH_TESTS}),
            'graph_statistics': {label: {'parsed_modules': tree.parsed_nodes, 'blob_count': len(tree.blobs)}
                                 for label, tree in [('base', before), ('head', after)]},
            'note': 'S1 plan only, not test success or authority to skip existing workflow gates'}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--base'); parser.add_argument('--head', default='HEAD')
    parser.add_argument('--event', default='pull_request', choices=['pull_request', 'push', 'workflow_dispatch', 'release', 'schedule'])
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--github-event', action='store_true')
    parser.add_argument('--github-output', type=Path)
    args = parser.parse_args(argv)
    try:
        root = args.repo.resolve()
        if args.github_event:
            event = json.loads(Path(os.environ['GITHUB_EVENT_PATH']).read_text(encoding='utf-8'))
            args.event = os.environ['GITHUB_EVENT_NAME']
            args.head = 'HEAD'
            if args.event == 'pull_request':
                args.base = git(root, 'merge-base', event['pull_request']['base']['sha'], 'HEAD').decode('ascii').strip()
            elif args.event == 'push' and event.get('before') and set(event['before']) != {'0'}:
                args.base = event['before']
            else:
                args.base = 'HEAD'
        if not args.base: raise ValueError('comparison base is required')
        result = plan(root, args.base, args.head, args.event); code = 0
        if args.github_output:
            fields = {'required': any(row['selected'] for row in result['cases']),
                      'native_required': bool(result['native_build_targets']),
                      'ams_required': bool(result['ams_paths']),
                      'python_required': any(row['selected'] and row['profile'] == 'solver-numerics'
                                             for row in result['cases']),
                      'heavy_required': bool(result['heavy_case_paths']),
                      'gmsh_required': bool(result['gmsh_reference_paths']) or any(
                          row['selected'] and 'gmsh' in row['components'] for row in result['cases'])}
            if args.event != 'pull_request': fields['native_required'] = True
            with args.github_output.open('a', encoding='utf-8') as output:
                for key, value in fields.items(): output.write(key + '=' + str(value).lower() + '\n')
    except Exception as exc:
        result = {'schema': SCHEMA, 'status': 'blocked', 'qualification': False, 'execution_wired': False,
                  'error': type(exc).__name__ + ': ' + str(exc), 'required_action': 'do not skip gates; repair inputs or execute full existing gate'}; code = 2
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({key: result[key] for key in ['schema', 'status', 'qualification', 'execution_wired']}))
    return code


if __name__ == '__main__':
    raise SystemExit(main())
