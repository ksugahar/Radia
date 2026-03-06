"""
AST Static Analysis: Extract all gmsh.* API calls from gmsh_builder source
and compare argument counts against actual GMSH Python API signatures.

This catches:
- Wrong number of positional arguments
- Calls to nonexistent GMSH API functions
- Missing required arguments
"""

import ast
import inspect
import os
import sys
import importlib

# Ensure gmsh is available
try:
    import gmsh
except ImportError:
    print("ERROR: gmsh not installed")
    sys.exit(1)


# Functions available in newer GMSH versions that may not exist in the
# installed version.  These are skipped instead of flagged as errors.
VERSION_DEPENDENT_FUNCTIONS = {
    'gmsh.model.getChangedEntities',  # GMSH >= 4.12
}


def resolve_gmsh_function(call_str):
    """Resolve a dotted gmsh call string to the actual function object.
    e.g. 'gmsh.model.occ.addBox' -> gmsh.model.occ.addBox
    """
    parts = call_str.split('.')
    if parts[0] != 'gmsh':
        return None
    obj = gmsh
    for part in parts[1:]:
        try:
            obj = getattr(obj, part)
        except AttributeError:
            return None
    return obj


def get_signature_info(func):
    """Get min/max positional args for a function."""
    try:
        sig = inspect.signature(func)
    except (ValueError, TypeError):
        return None

    min_args = 0
    max_args = 0
    has_var_positional = False
    has_var_keyword = False

    for name, param in sig.parameters.items():
        if param.kind == inspect.Parameter.VAR_POSITIONAL:
            has_var_positional = True
        elif param.kind == inspect.Parameter.VAR_KEYWORD:
            has_var_keyword = True
        elif param.kind in (inspect.Parameter.POSITIONAL_ONLY,
                            inspect.Parameter.POSITIONAL_OR_KEYWORD):
            max_args += 1
            if param.default is inspect.Parameter.empty:
                min_args += 1
        elif param.kind == inspect.Parameter.KEYWORD_ONLY:
            # keyword-only args don't count as positional
            pass

    if has_var_positional:
        max_args = float('inf')

    return {
        'min_args': min_args,
        'max_args': max_args,
        'has_var_positional': has_var_positional,
        'signature': str(sig) if sig else '?',
    }


class GmshCallVisitor(ast.NodeVisitor):
    """AST visitor that extracts all gmsh.* function calls."""

    def __init__(self, filename):
        self.filename = filename
        self.calls = []

    def _get_dotted_name(self, node):
        """Recursively build dotted name from Attribute/Name nodes."""
        if isinstance(node, ast.Name):
            return node.id
        elif isinstance(node, ast.Attribute):
            parent = self._get_dotted_name(node.value)
            if parent:
                return f"{parent}.{node.attr}"
        return None

    def visit_Call(self, node):
        name = self._get_dotted_name(node.func)
        if name and name.startswith('gmsh.'):
            # Count positional args (excluding *args unpacking)
            n_positional = 0
            has_starargs = False
            for arg in node.args:
                if isinstance(arg, ast.Starred):
                    has_starargs = True
                else:
                    n_positional += 1

            # Count keyword args
            keyword_names = []
            has_kwargs = False
            for kw in node.keywords:
                if kw.arg is None:
                    has_kwargs = True
                else:
                    keyword_names.append(kw.arg)

            self.calls.append({
                'function': name,
                'n_positional': n_positional,
                'keyword_names': keyword_names,
                'has_starargs': has_starargs,
                'has_kwargs': has_kwargs,
                'line': node.lineno,
                'file': self.filename,
            })

        self.generic_visit(node)


def analyze_file(filepath):
    """Parse a Python file and extract all gmsh.* calls."""
    with open(filepath, 'r', encoding='utf-8') as f:
        source = f.read()

    try:
        tree = ast.parse(source, filename=filepath)
    except SyntaxError as e:
        print(f"  SYNTAX ERROR in {filepath}: {e}")
        return []

    visitor = GmshCallVisitor(os.path.basename(filepath))
    visitor.visit(tree)
    return visitor.calls


def main():
    gmsh_builder_dir = os.path.join('src', 'radia', 'gmsh_builder')

    # Initialize gmsh to make all APIs available
    gmsh.initialize()

    # Collect all calls
    all_calls = []
    py_files = sorted(f for f in os.listdir(gmsh_builder_dir) if f.endswith('.py'))

    print(f"Scanning {len(py_files)} files in gmsh_builder/\n")

    for filename in py_files:
        filepath = os.path.join(gmsh_builder_dir, filename)
        calls = analyze_file(filepath)
        all_calls.extend(calls)
        if calls:
            # Count unique functions
            unique_fns = set(c['function'] for c in calls)
            print(f"  {filename}: {len(calls)} calls ({len(unique_fns)} unique functions)")

    print(f"\nTotal gmsh.* calls found: {len(all_calls)}")

    # Deduplicate by (function, n_positional, keyword_names) for analysis
    unique_functions = set(c['function'] for c in all_calls)
    print(f"Unique gmsh functions called: {len(unique_functions)}\n")

    # Check each call against actual API
    errors = []
    warnings = []
    not_found = []

    for call in all_calls:
        func_name = call['function']
        func = resolve_gmsh_function(func_name)

        if func is None:
            if func_name in VERSION_DEPENDENT_FUNCTIONS:
                continue  # skip version-dependent API
            not_found.append(call)
            continue

        if not callable(func):
            errors.append({
                **call,
                'error': f"'{func_name}' is not callable (it's a {type(func).__name__})"
            })
            continue

        sig_info = get_signature_info(func)
        if sig_info is None:
            warnings.append({
                **call,
                'warning': f"Could not inspect signature of '{func_name}'"
            })
            continue

        # Check positional arg count
        n_pos = call['n_positional']
        n_kw = len(call['keyword_names'])

        if call['has_starargs'] or call['has_kwargs']:
            # Can't verify dynamic calls
            continue

        # Total args passed = positional + keyword that might be positional
        # But keyword args could be keyword-only, so just check positional
        if n_pos < sig_info['min_args']:
            errors.append({
                **call,
                'error': f"Too few positional args: {n_pos} given, minimum {sig_info['min_args']} required",
                'signature': sig_info['signature'],
            })
        elif n_pos + n_kw > sig_info['max_args'] and sig_info['max_args'] != float('inf'):
            errors.append({
                **call,
                'error': f"Too many args: {n_pos} positional + {n_kw} keyword = {n_pos+n_kw}, maximum {sig_info['max_args']}",
                'signature': sig_info['signature'],
            })

        # Check if keyword names are valid
        if func and not call['has_kwargs']:
            try:
                sig = inspect.signature(func)
                valid_params = set(sig.parameters.keys())
                has_var_kw = any(
                    p.kind == inspect.Parameter.VAR_KEYWORD
                    for p in sig.parameters.values()
                )
                if not has_var_kw:
                    for kw_name in call['keyword_names']:
                        if kw_name not in valid_params:
                            errors.append({
                                **call,
                                'error': f"Unknown keyword argument '{kw_name}'",
                                'signature': sig_info['signature'],
                            })
            except (ValueError, TypeError):
                pass

    # Report
    print("=" * 80)
    print("RESULTS")
    print("=" * 80)

    if not_found:
        # Deduplicate
        nf_unique = {}
        for nf in not_found:
            key = nf['function']
            if key not in nf_unique:
                nf_unique[key] = []
            nf_unique[key].append(f"{nf['file']}:{nf['line']}")

        print(f"\n--- NONEXISTENT FUNCTIONS ({len(nf_unique)} unique) ---")
        for func_name, locations in sorted(nf_unique.items()):
            print(f"  {func_name}")
            for loc in locations[:3]:
                print(f"    at {loc}")
            if len(locations) > 3:
                print(f"    ... and {len(locations)-3} more")

    if errors:
        # Deduplicate errors by function+error
        err_unique = {}
        for err in errors:
            key = (err['function'], err['error'])
            if key not in err_unique:
                err_unique[key] = []
            err_unique[key].append(f"{err['file']}:{err['line']}")

        print(f"\n--- ARGUMENT ERRORS ({len(err_unique)} unique) ---")
        for (func_name, error_msg), locations in sorted(err_unique.items()):
            print(f"\n  ERROR: {func_name}")
            print(f"    {error_msg}")
            if 'signature' in errors[0]:
                # Find the original error to get signature
                for e in errors:
                    if e['function'] == func_name and e['error'] == error_msg and 'signature' in e:
                        print(f"    Actual signature: {e['signature']}")
                        break
            for loc in locations[:5]:
                print(f"    at {loc}")
            if len(locations) > 5:
                print(f"    ... and {len(locations)-5} more")

    if warnings:
        print(f"\n--- WARNINGS ({len(warnings)}) ---")
        for w in warnings[:10]:
            print(f"  {w['function']} at {w['file']}:{w['line']}: {w['warning']}")

    if not errors and not not_found:
        print("\nAll gmsh.* API calls appear to have correct argument counts!")

    # Summary
    print(f"\n{'=' * 80}")
    print(f"SUMMARY: {len(all_calls)} total calls, {len(unique_functions)} unique functions")
    print(f"  Nonexistent: {len(set(nf['function'] for nf in not_found)) if not_found else 0}")
    print(f"  Arg errors:  {len(err_unique) if errors else 0}")
    print(f"  Warnings:    {len(warnings)}")
    print(f"{'=' * 80}")

    gmsh.finalize()

    return len(errors) + len(not_found)


if __name__ == '__main__':
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    exit_code = main()
    sys.exit(1 if exit_code > 0 else 0)
