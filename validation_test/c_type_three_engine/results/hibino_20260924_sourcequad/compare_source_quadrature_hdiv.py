"""Experimental source-load quadrature only; never a three-engine pass."""
import hashlib
import importlib.util
import json
import os
import platform
import sys
from pathlib import Path

# Patch only a disposable venv's source BEFORE importing Radia. The numerical
# change is one exact expression; Gram/mass quadrature remains untouched.
root = Path(__file__).resolve().parent
assert Path(sys.prefix).resolve() == (root / 'venv').resolve(), sys.prefix
source = Path(sys.prefix) / 'Lib/site-packages/radia/vim/_solve.py'
original = source.read_bytes()
needle = b'H_ext * fes.TestFunction() * ng.dx'
assert original.count(needle) == 1
replacement = b'H_ext * fes.TestFunction() * ng.dx(bonus_intorder=int(__import__("os").environ["HDIV_SOURCE_BONUS"]))'
modified = original.replace(needle, replacement)
(root / 'original_solve.py').write_bytes(original)
(root / 'experimental_solve.py').write_bytes(modified)
source.write_bytes(modified)
try:
    import numpy as np
    import run_three_engine as d
    original_solve = d.vim.Solve
    original_adapter = d.solve_hdiv
    captured = {}
    rows = []
    class DiagnosticDone(Exception):
        pass
    def solve(*args, **kwargs):
        kwargs['tol'] = 1e-10
        result = original_solve(*args, **kwargs)
        captured.clear()
        captured.update({k: v for k, v in result.items()
                         if isinstance(v, (str, int, float, bool, type(None)))})
        return result
    def adapter(*args, **kwargs):
        field, diag = original_adapter(*args, **kwargs)
        rows.append(dict(bonus_intorder=bonus, field_T=field.tolist(),
                         diagnostics=diag, solve_scalars=dict(captured)))
        raise DiagnosticDone()
    d.vim.Solve = solve
    d.solve_hdiv = adapter
    for bonus in (0, 4, 8):
        os.environ['HDIV_SOURCE_BONUS'] = str(bonus)
        sys.argv = ['run_three_engine.py', '--mesh-dir', str(root/'meshes'),
                    '--mode', 'linear', '--hdiv-order', '1', '--fem-order', '2',
                    '--threads', '8', '--output', str(root/'not_a_three_engine_result.json')]
        try:
            d.main()
        except DiagnosticDone:
            pass
    baseline = np.asarray(json.loads((root/'baseline.json').read_text())['fields_T']['hdiv_mmm'])
    for row in rows:
        row['relative_vector_difference_from_previous_fine'] = float(
            np.linalg.norm(np.asarray(row['field_T'])-baseline)/np.linalg.norm(baseline))
    result = dict(scope=__doc__, host=platform.node(), python=sys.version, rows=rows,
                  original_source_sha256=hashlib.sha256(original).hexdigest(),
                  experimental_source_sha256=hashlib.sha256(modified).hexdigest(),
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (root/'source_quadrature_comparison.json').write_text(json.dumps(result, indent=2))
    print(json.dumps(result), flush=True)
finally:
    source.write_bytes(original)
