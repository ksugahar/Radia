"""Runtime-independent gates for the repaired ESRF6 acceptance driver."""
import hashlib
import math
from pathlib import Path, PureWindowsPath


def require_accepted_hdiv(result):
    """Require literal success flags and a finite, nonnegative true residual."""
    if (result.get('completed') is not True or result.get('accepted') is not True
            or type(result.get('case')) is not int or result['case'] != 6):
        raise ValueError('An accepted, completed ESRF6 HDiv result is required')
    try:
        stats = result['nonlinear_stats']
        residual = stats['nonlinear_final_relative_residual']
        tolerance = result['options']['nl_tol']
        valid = (type(residual) in (int, float) and math.isfinite(residual)
                 and residual >= 0 and type(tolerance) in (int, float)
                 and math.isfinite(tolerance) and tolerance > 0
                 and residual <= tolerance
                 and stats['nonlinear_converged_final_stage'] is True)
    except (KeyError, TypeError, OverflowError) as exc:
        raise ValueError('Missing or invalid HDiv residual contract') from exc
    if not valid:
        raise ValueError('HDiv true residual did not pass')


def _sha256(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def verify_relocated_identity(implementation, package, iron_mesh, wheel):
    """Verify content, not historical absolute paths; return current provenance."""
    package = Path(package).resolve(strict=True)
    paths = {'native': package / '_radia_pybind.pyd',
             'iron': Path(iron_mesh), 'wheel': Path(wheel)}
    expected = implementation['identities']
    if set(expected) != set(paths):
        raise ValueError('Expected exactly native, iron and wheel identities')
    identities = {}
    for name, path in paths.items():
        path = path.resolve(strict=True)
        digest = _sha256(path)
        if digest != expected[name]['sha256']:
            raise RuntimeError(f'{name} differs from the HDiv input/runtime')
        identities[name] = {'path': str(path), 'sha256': digest}
    actual = {p.relative_to(package).as_posix(): _sha256(p)
              for p in sorted(package.rglob('*.py'))}
    sources = implementation['python_sources']
    normalized = {PureWindowsPath(name).as_posix(): digest
                  for name, digest in sources.items()}
    if not sources or len(normalized) != len(sources) or actual != normalized:
        raise RuntimeError('Installed Python sources differ from the HDiv wheel')
    return {'package_path': str(package), 'identities': identities,
            'python_sources': actual}
