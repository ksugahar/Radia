"""Explicit GPU-only resident AMG adapter; import is numerical-runtime free."""
from importlib import import_module
import hashlib
import json
import math
from numbers import Integral, Real
from pathlib import Path
import re
import sysconfig


def _qualified_native():
    # Verify installed bytes before opening any CUDA search directory.
    from . import verify
    identity = verify()
    root = Path(__file__).resolve().parent
    name = '_resident_amg' + sysconfig.get_config_var('EXT_SUFFIX')
    if not (root / 'resident-amg-build.json').is_file() or not (root / name).is_file():
        raise ImportError('Qualified installed resident AMG receipt and binary required')
    receipt = json.loads((root / 'resident-amg-build.json').read_bytes())
    if (receipt.get('schema') != 'radia.resident-amg-python-build.v1'
            or receipt.get('source_commit') != identity['source_commit']
            or receipt.get('source_dirty') is not False
            or receipt.get('binary') != name
            or receipt.get('extension_suffix') != sysconfig.get_config_var('EXT_SUFFIX')
            or receipt.get('binary_sha256') != hashlib.sha256((root / name).read_bytes()).hexdigest()
            or any(receipt.get(key) is not False for key in (
                'fault_injection', 'device_fault_injection', 'lifecycle_diagnostic',
                'matlab_engine_user_dependency', 'matlab_dlls_loaded'))
            or any(receipt.get(key) is not True for key in (
                'runtime_import_checked', 'redistribution_accepted', 'matlab_dll_search_excluded'))
            or any(not isinstance(receipt.get(key), str) or not re.fullmatch(r'[0-9a-f]{64}', receipt[key]) for key in (
                'runtime_validation_sha256', 'loaded_runtime_sha256',
                'redistribution_reference_sha256'))):
        raise ImportError('Qualified installed resident AMG receipt required')
    from ._cuda_runtime import prepare_cuda_runtime
    prepare_cuda_runtime()
    return import_module('._resident_amg', __package__)


class ResidentAMG:
    """Own a native hierarchy. CPU workflows remain independent in radia.

    ``levels`` and ``coarse_inverse`` follow the checked native row-major CSR
    contract. No hierarchy construction or physical formulation is performed.
    """

    def __init__(self, levels, coarse_inverse, *, backend):
        if backend != 'gpu':
            raise ValueError('ResidentAMG requires explicit backend="gpu"')
        self._owner = _qualified_native().Resident(levels, coarse_inverse)

    def solve(self, rhs, *, tolerance=1e-12, max_iterations=10000):
        if (isinstance(tolerance, bool) or not isinstance(tolerance, Real)
                or not math.isfinite(tolerance) or not 0 < tolerance <= 1e-6):
            raise ValueError('AMG tolerance must be finite and in (0, 1e-6]')
        if (isinstance(max_iterations, bool) or not isinstance(max_iterations, Integral)
                or not 1 <= max_iterations <= 100000):
            raise ValueError('AMG iteration limit must be an integer in [1, 100000]')
        if self._owner is None:
            raise RuntimeError('Resident AMG owner is closed')
        return self._owner.solve(rhs, float(tolerance), int(max_iterations))

    def close(self):
        owner = getattr(self, '_owner', None)
        if owner is not None:
            owner.close()
            self._owner = None

    def __enter__(self):
        if self._owner is None:
            raise RuntimeError('Resident AMG owner is closed')
        return self

    def __exit__(self, *_):
        self.close()
