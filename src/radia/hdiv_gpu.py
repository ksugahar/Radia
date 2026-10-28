"""Headless HDiv-MMM workflow over the shared CUDA material/Newton core.

NGSolve owns FE assembly; Radia owns charge-Gram preparation. The caller owns
the TaskManager region. Construction prepares operators, never a CPU solution.
"""
from __future__ import annotations

__all__ = ["PreparedHDivNonlinear", "prepare_hdiv_host_setup"]

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import threading
import tempfile
import math
import numbers
from contextlib import contextmanager


class PreparedHDivNonlinear:
    """Persistent GPU solve for a checked body mesh and explicit inverse-BH law.

    ``law`` is the shared inverse-BH polynomial payload, not a new table fit.
    Only the explicit GPU route is implemented here. CPU workflows remain in
    radia.vim; unsupported geometry or a missing native binding raises.
    """

    @classmethod
    def from_bh_table(cls, mesh_path, bh_table, *, backend, label_contract, order=1):
        """Use the existing Radia CPU polynomial construction, without a CPU solve.

        Rows are [H in A/m, B in tesla]; one uniform soft-iron law is supported.
        """
        if backend != "gpu":
            raise ValueError("This workflow requires explicit backend='gpu'")
        import numpy as np
        from radia.vim._nonlinear import _validate_bh_table, _bh_inverse_payload
        table = np.asarray(bh_table, dtype=np.float64)
        _validate_bh_table(table, context="shared HDiv GPU")
        if table.shape[0] > 4096:
            raise ValueError("BH table exceeds the shared material budget")
        law = _bh_inverse_payload(table[:, 0], table[:, 1])
        return cls(mesh_path, law, backend=backend, label_contract=label_contract, order=order)

    def __init__(self, mesh_path, law, *, backend, label_contract, order=1):
        if backend != "gpu":
            raise ValueError("This workflow requires explicit backend='gpu'; CPU uses radia.vim")
        if isinstance(order, bool) or order not in (1, 2):
            raise ValueError("HDiv order must be 1 or 2")
        with _checked_mesh_snapshot(mesh_path, label_contract) as (snapshot, mesh_hash, contract_hash):
            self.mesh_sha256 = mesh_hash
            self.label_contract_sha256 = contract_hash
            self._initialize(snapshot, law, order)

    def _initialize(self, path, law, order):
        import numpy as np
        from .gpu import _native
        _hdiv_cuda = _native("_hdiv_cuda")
        self._lock = threading.RLock()
        self._operator = self._material = None
        self._closed = False
        self.order = order
        try:
            setup = _prepare_host_components(path, order)
            self.dofs = setup["dofs"]
            B, snapshot, E = setup["charge_map"], setup["snapshot"], setup["evaluation"]
            self._mass, self._rhs = setup["mass"], setup["rhs_operator"]
            self._operator = _hdiv_cuda._ResidentHDiv(snapshot["leaves"],
                snapshot["permutation"], snapshot["scale"], self.dofs,
                B.indptr.astype(np.uint64), B.indices.astype(np.uint64),
                np.ascontiguousarray(B.data), backend="gpu")
            native_law = dict(law)
            for name in ("knots", "coefficients", "integral_base"):
                native_law[name] = np.ascontiguousarray(native_law[name], dtype=np.float64).ravel()
            self._material = _hdiv_cuda._ResidentBHMaterial(E.indptr.astype(np.uint64),
                E.indices.astype(np.uint64), np.ascontiguousarray(E.data),
                np.ascontiguousarray(setup["weights"]), self.dofs, native_law, backend="gpu")
        except BaseException:
            self.close()
            raise

    def solve(self, applied_field, *, initial=None, tolerance=1e-11,
              linear_tolerance=5e-14, max_iterations=40, max_linear_iterations=10000,
              preconditioner="reference-diagonal", degree=32, lower_fraction=.01):
        """Solve a uniform applied H (A/m); return owned coefficients and history.

        Native Newton recertifies the original nonlinear equation. This is not
        an independent CPU certificate or a performance claim.
        """
        if preconditioner not in ("reference-diagonal", "reference-chebyshev", "reference-chebyshev-graph"):
            raise ValueError("Unsupported explicit shared Newton preconditioner")
        for value in (tolerance, linear_tolerance):
            if isinstance(value, bool) or not isinstance(value, numbers.Real) or not math.isfinite(value) or not 0 < value <= 1e-6:
                raise ValueError("Solve tolerances must be finite in (0,1e-6]")
        for value in (degree, max_iterations, max_linear_iterations):
            if isinstance(value, bool) or not isinstance(value, numbers.Integral) or value < 1:
                raise ValueError("Iteration limits and degree must be positive integers")
        if degree > 64 or isinstance(lower_fraction, bool) or not isinstance(lower_fraction, numbers.Real) or not math.isfinite(lower_fraction) or not 0 < lower_fraction < 1:
            raise ValueError("Invalid shared Newton polynomial settings")
        import numpy as np
        with self._lock:
            if self._closed:
                raise RuntimeError("Prepared HDiv workflow has been closed")
            field = np.asarray(applied_field, dtype=np.float64)
            if field.shape != (3,) or not np.isfinite(field).all():
                raise ValueError("applied_field must contain three finite H components in A/m")
            x = np.zeros(self.dofs) if initial is None else np.asarray(initial)
            if x.dtype != np.float64 or x.shape != (self.dofs,) or not np.isfinite(x).all():
                raise ValueError("initial must be finite FP64 in the prepared FE numbering")
            m = self._mass
            result = self._operator.solve_nonlinear_material(self._material,
                m.indptr.astype(np.uint64), m.indices.astype(np.uint64),
                np.ascontiguousarray(m.data), np.ascontiguousarray(self._rhs @ field),
                np.ascontiguousarray(x), tolerance, linear_tolerance, max_iterations,
                max_linear_iterations, preconditioner, degree, lower_fraction)
            return dict(result, backend="gpu", mesh_sha256=self.mesh_sha256,
                        order=self.order, dofs=self.dofs, field_units="A/m")

    def stats(self):
        with self._lock:
            return dict(closed=self._closed, backend="gpu", dofs=getattr(self, "dofs", None),
                operator=self._operator.stats() if self._operator is not None else None,
                material=self._material.stats() if self._material is not None else None)

    def close(self):
        with self._lock:
            try:
                if self._material is not None:
                    self._material.close()
            finally:
                if self._operator is not None:
                    self._operator.close()
                self._closed = True

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


@contextmanager
def _checked_mesh_snapshot(mesh_path, label_contract):
    path = Path(mesh_path).resolve(strict=True)
    contract = Path(label_contract).resolve(strict=True)
    mesh_bytes, contract_bytes = path.read_bytes(), contract.read_bytes()
    with tempfile.TemporaryDirectory(prefix="radia-hdiv-mesh-") as directory:
        suffix = ".vol.gz" if path.name.lower().endswith(".vol.gz") else ".vol"
        if not path.name.lower().endswith((".vol", ".vol.gz")):
            raise ValueError("HDiv mesh input must be .vol or .vol.gz")
        snapshot = Path(directory) / ("mesh" + suffix)
        checked_contract = Path(directory) / "labels.json"
        snapshot.write_bytes(mesh_bytes)
        checked_contract.write_bytes(contract_bytes)
        checked = subprocess.run([sys.executable, "-m", "cubit_mesh_export.check",
            str(snapshot), "--strict-labels", "--contract", str(checked_contract), "--format", "json"],
            capture_output=True, text=True)
        if checked.returncode:
            raise RuntimeError("cubit-mesh-export required for check-vol; mesh validation failed: "
                               + checked.stderr.strip() + "\n" + checked.stdout.strip())
        if not json.loads(checked.stdout).get("passed"):
            raise ValueError("Mesh check did not pass: " + checked.stdout)
        yield snapshot, hashlib.sha256(mesh_bytes).hexdigest(), hashlib.sha256(contract_bytes).hexdigest()


def _prepare_host_components(path, order, *, maximum_dofs=16384):
    """One NG-owned FE/charge preparation shared by Python and MATLAB publishers.

    No resident GPU operator and no CPU solution is created here. The caller
    owns TaskManager; the private native material staging shares the same FE.
    """
    import numpy as np
    import ngsolve as ng
    from scipy.sparse import coo_matrix, csc_matrix
    from radia.vim._vim import build_charge_gram
    from .gpu import _native
    _hdiv_cuda = _native("_hdiv_cuda")
    mesh = ng.Mesh(str(path))
    fes = ng.HDiv(mesh, order=order)
    nd = fes.ndof
    if not 0 < nd <= maximum_dofs:
        raise ValueError("HDiv workflow dimension exceeds bounded scope")
    B, gram, _ = build_charge_gram(fes, eps=1e-10, leafsize=64, eta=2.)
    snapshot = gram._resident_leaf_snapshot(16)
    if snapshot["dimension"] != B.shape[0]:
        raise ValueError("Charge snapshot numbering mismatch")
    u, v = fes.TnT()
    mass = ng.BilinearForm(fes, symmetric=True)
    mass += ng.InnerProduct(u, v)*ng.dx
    mass.Assemble()
    rows, columns, values = mass.mat.COO()
    M = coo_matrix((values, (rows, columns)), shape=(nd, nd)).tocsr()
    M.sort_indices()
    rhs = np.empty((nd, 3), dtype=np.float64)
    for component in range(3):
        direction = tuple(float(i == component) for i in range(3))
        source = ng.LinearForm(fes)
        source += ng.InnerProduct(ng.CF(direction), v)*ng.dx
        source.Assemble()
        rhs[:, component] = source.vec.FV().NumPy()
    staging = _hdiv_cuda._MaterialStagingContext(str(path), order)
    try:
        payload = staging.quadrature()
    finally:
        staging.close()
    if payload["dofs"] != nd:
        raise ValueError("Material and charge FE numbering mismatch")
    E = csc_matrix((payload["evaluation_values"], payload["evaluation_row_indices"],
        payload["evaluation_column_offsets"]), shape=(payload["evaluation_rows"], nd)).tocsr()
    E.sort_indices()
    # Same declared NG default volume rule as the existing CPU mean field.
    volume = float(ng.Integrate(1, mesh))
    if not math.isfinite(volume) or volume <= 0:
        raise ValueError("Invalid mesh volume")
    return dict(dofs=nd, mass=M, rhs_operator=rhs, mean_operator=rhs.T/volume,
                charge_map=B, snapshot=snapshot, evaluation=E,
                weights=np.asarray(payload["weights"]), quadrature_order=payload["integration_order"])


def prepare_hdiv_host_setup(mesh_path, *, backend, label_contract, order=1, maximum_dofs=16384):
    """Return owned zero-based FE/charge/material arrays; no resident GPU solve.

    NGSolve assembly stays under the caller's TaskManager. Consumers publish
    their own ABI layout without reimplementing FE preparation or domain math.
    """
    if backend != "gpu":
        raise ValueError("Host setup requires explicit backend='gpu'")
    if isinstance(order, bool) or order not in (1, 2):
        raise ValueError("HDiv order must be 1 or 2")
    if isinstance(maximum_dofs, bool) or not isinstance(maximum_dofs, int) or not 0 < maximum_dofs <= 16384:
        raise ValueError("Unsupported bounded host setup dimension")
    with _checked_mesh_snapshot(mesh_path, label_contract) as (path, mesh_hash, contract_hash):
        data = _prepare_host_components(path, order, maximum_dofs=maximum_dofs)
    return dict(data, mesh_sha256=mesh_hash, label_contract_sha256=contract_hash, order=order)
