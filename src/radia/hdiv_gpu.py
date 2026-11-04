"""Headless HDiv-MMM workflow over the shared CUDA material/Newton core.

NGSolve owns FE assembly; Radia owns charge-Gram preparation. The caller owns
the TaskManager region. Construction prepares operators, never a CPU solution.
The charge Gram, the charge map B/B^T, the material evaluation map and the
Newton reference mass are generated on the device from affine or P2-curved TET
or Q2-lattice HEX geometry and NGSolve reference tables; no host-built sparse
operator is uploaded.
"""
from __future__ import annotations

__all__ = ["PreparedHDivLinear", "PreparedHDivNonlinear", "prepare_hdiv_device_setup"]

import base64
from importlib import metadata
from pathlib import PureWindowsPath
import re
import hashlib
import json
from pathlib import Path
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

    ``gram_aca`` holds the device ACA parameters of the charge Gram
    (eps, eta, leaf_size, max_rank, rank_cap). Leaf size and admissibility
    are problem dependent and are never chosen automatically; ``None``
    selects the documented CPU charge-Gram preset ``GRAM_ACA_PRESET``.
    ``stats()`` reports the parameters actually used.

    Supported domain (every bound fails loudly before device allocation):

    - affine tetrahedra (``curve_order=1``) or NGSolve ``Curve(2)``
      isoparametric P2 tetrahedra (``curve_order=2``), or pure hexahedra
      (linear or ``Curve(2)``, family ``hex``, detected from the mesh, with the
      CPU graded near family on thin and non-conforming BDM1 pairs), HDiv
      order 1 or 2, one soft-iron law with at most 4096 table rows;
    - material rule: NGSolve ``IntegrationRule(TET or HEX, 2*max(6, 2*order+2))``
      (343 points per TET element), shared with the CPU route and never reduced;
    - material: no fixed element, DOF or point cap. Numbering is int32, and
      the mass triples ``elements*local**2`` are at most 2**31-1 (CUB item
      count). The
      device bytes in ``stats()["material_budget"]`` (resident element
      tables, DOF incidence and its sort, mass triples, mass CSR and frozen
      tangent values, material point state) plus CUB scratch must fit the live
      free memory above the reserve max(2 GiB, total/5); the error names the
      limiting term. E and E^T act matrix-free; each linearization builds the
      Newton tangent on the device on the mass pattern; the 343-point material
      state remains;
    - charge Gram and charge map: no fixed charge, cell or DOF cap. Charge,
      cell, face and DOF numbering is int32, and the charge-map nonzeros are
      at most 2**31-1 (CUB item count). The device bytes in
      ``stats()["charge_budget"]`` (prepared charges, B and B^T, the map
      build scratch, resident O(n) vectors, and the dense n**2 Gram or the
      ACA batch scratch) must fit the live free memory above the same
      reserve; the error names the limiting term. ACA factor storage is
      known only after compression and is checked live where it is
      allocated. ACA leaf_size <= 4096 and max_rank <= 1024;
    - binding host snapshot of the tables: at most 512 MiB.
    """

    @classmethod
    def from_bh_table(cls, mesh_path, bh_table, *, backend, label_contract, mesh_check_receipt,
                      mesh_check_receipt_sha256, mesh_check_provenance, order=1, gram_aca=None, curve_order=1,
                      image=None, image_cyclic=None, image_cyclic_alternating=False, cyclic_periodic_boundaries=(),
                      q2_deformation=None):
        """Use the existing Radia CPU polynomial construction, without a CPU solve.

        Rows are [H in A/m, B in tesla]; one uniform soft-iron law is supported.
        """
        if backend != "gpu":
            raise ValueError("This workflow requires explicit backend='gpu'")
        if isinstance(order, bool) or order not in (1, 2):
            raise ValueError("HDiv order must be 1 or 2")
        curve_order = _checked_curve_order(curve_order, order)
        aca = _checked_gram_aca(gram_aca)
        with _checked_mesh_snapshot(mesh_path, label_contract, mesh_check_receipt,
                                    mesh_check_receipt_sha256, mesh_check_provenance) as (snapshot, mesh_hash, contract_hash):
            import numpy as np
            from radia.vim._nonlinear import _validate_bh_table, _bh_inverse_payload
            table = np.asarray(bh_table, dtype=np.float64)
            _validate_bh_table(table, context="shared HDiv GPU")
            if table.shape[0] > 4096:
                raise ValueError("BH table exceeds the shared material budget")
            law = _bh_inverse_payload(table[:, 0], table[:, 1])
            prepared = cls.__new__(cls)
            prepared.mesh_sha256 = mesh_hash
            prepared.label_contract_sha256 = contract_hash
            prepared._initialize(snapshot, law, order, aca, curve_order,
                                 _DeviceImages(image, image_cyclic, image_cyclic_alternating, cyclic_periodic_boundaries),
                                 q2_deformation)
            return prepared

    def __init__(self, mesh_path, law, *, backend, label_contract, mesh_check_receipt,
                 mesh_check_receipt_sha256, mesh_check_provenance, order=1, gram_aca=None, curve_order=1,
                 image=None, image_cyclic=None, image_cyclic_alternating=False, cyclic_periodic_boundaries=(),
                 q2_deformation=None):
        if backend != "gpu":
            raise ValueError("This workflow requires explicit backend='gpu'; CPU uses radia.vim")
        if isinstance(order, bool) or order not in (1, 2):
            raise ValueError("HDiv order must be 1 or 2")
        curve_order = _checked_curve_order(curve_order, order)
        aca = _checked_gram_aca(gram_aca)
        with _checked_mesh_snapshot(mesh_path, label_contract, mesh_check_receipt,
                                    mesh_check_receipt_sha256, mesh_check_provenance) as (snapshot, mesh_hash, contract_hash):
            self.mesh_sha256 = mesh_hash
            self.label_contract_sha256 = contract_hash
            self._initialize(snapshot, law, order, aca, curve_order,
                             _DeviceImages(image, image_cyclic, image_cyclic_alternating, cyclic_periodic_boundaries),
                             q2_deformation)

    def _initialize(self, path, law, order, aca, curve_order, images=None, q2_deformation=None):
        import numpy as np
        from .gpu import _native
        _hdiv_cuda = _native("_hdiv_cuda")
        self._lock = threading.RLock()
        self._operator = self._material = None
        self._closed = False
        self.order = order
        self.family = None
        self.gram_aca = dict(aca)
        try:
            setup = _prepare_device_components(path, order, _hdiv_cuda, aca=aca, curve_order=curve_order, images=images,
                                               q2_deformation=q2_deformation)
            self._gram_arguments = setup["charge_gram"]
            self.family = setup["family"]
            self.dofs = setup["dofs"]
            self.fespace, self._fe_owner = setup["fespace"], setup["fe_owner"]
            self._rhs = setup["rhs_operator"]
            self._mean = setup["mean_operator"]
            self.material_budget = setup["material_budget"]
            native_law = dict(law)
            for name in ("knots", "coefficients", "integral_base"):
                native_law[name] = np.ascontiguousarray(native_law[name], dtype=np.float64).ravel()
            # The material is built first: its live device budget check runs before the Gram.
            self._material = _hdiv_cuda._ResidentBHMaterial.from_flat_tables(
                **setup["evaluation_tables"], dofs=self.dofs, law=native_law, backend="gpu")
            self.material_assembly = _checked_material_assembly(self._material.assembly(),
                                                                setup["evaluation_tables"])
            _require_budget_cover(self.material_assembly, self._material.stats())
            self._operator = _hdiv_cuda._ResidentHDiv.from_flat_charges(
                **setup["charge_gram"], dofs=self.dofs, backend="gpu")
            self.assembly = _checked_device_assembly(self._operator.assembly(), setup["charge_gram"])
        except BaseException:
            self.close()
            raise

    def solve(self, applied_field, *, initial=None, tolerance=1e-11,
              linear_tolerance=1e-6, max_iterations=40, max_linear_iterations=30000,
              preconditioner="reference-diagonal", degree=32, lower_fraction=.01):
        """Solve an applied H (A/m), uniform or a CoefficientFunction; return coefficients, history and the mean magnetization.

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
            source = _applied_source(self.fespace, self._rhs, applied_field)
            x = np.zeros(self.dofs) if initial is None else np.asarray(initial)
            if x.dtype != np.float64 or x.shape != (self.dofs,) or not np.isfinite(x).all():
                raise ValueError("initial must be finite FP64 in the prepared FE numbering")
            # The reference mass is the material context's device-built E^T W E.
            result = self._operator.solve_nonlinear_material(self._material, source,
                np.ascontiguousarray(x), tolerance, linear_tolerance, max_iterations,
                max_linear_iterations, preconditioner, degree, lower_fraction)
            return dict(result, mean_magnetization=self._mean @ result["x"], backend="gpu",
                        mesh_sha256=self.mesh_sha256, order=self.order, dofs=self.dofs, field_units="A/m")

    def field(self, x, points):
        """H in A/m of a solved state at (points, 3) probes, images folded; a mixed mesh refuses.

        The device route of ``vim.FieldFromSolution(..., algorithm="direct")``: the charges
        are the resident charge map's product and only the state, geometry and probes cross.
        """
        import numpy as np
        with self._lock:
            if self._closed:
                raise RuntimeError("Prepared HDiv workflow has been closed")
            state = np.ascontiguousarray(x, dtype=np.float64)
            probes = np.ascontiguousarray(points, dtype=np.float64)
            if state.shape != (self.dofs,) or probes.ndim != 2 or probes.shape[1] != 3:
                raise ValueError("field requires the prepared state and (points, 3) probes: %r, %r"
                                 % (state.shape, probes.shape))
            return _checked_device_field(self._operator.field(state, probes), self._gram_arguments, state, probes)

    def stats(self):
        with self._lock:
            return dict(closed=self._closed, backend="gpu", dofs=getattr(self, "dofs", None),
                family=getattr(self, "family", None), gram_route=_gram_route(getattr(self, "family", None)),
                gram_compression=_GRAM_COMPRESSION,
                gram_aca=dict(self.gram_aca) if hasattr(self, "gram_aca") else None,
                assembly=getattr(self, "assembly", None),
                material_assembly=getattr(self, "material_assembly", None),
                material_budget=getattr(self, "material_budget", None),
                charge_budget=self.assembly["charge_budget"] if getattr(self, "assembly", None) else None,
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


class PreparedHDivLinear:
    """Persistent GPU linear isotropic material solve on a checked body mesh.

    Solves (B^T G B + inv_chi M) x = f(H0) with the device-built charge Gram,
    charge map and reference mass M (the material quadrature mass of
    ``PreparedHDivNonlinear``); only vectors cross. The supported domain and
    the ``gram_aca`` contract are those of ``PreparedHDivNonlinear`` without
    the inverse-BH law, plus conforming RT0 (``order=0``) and an explicit
    exact ``gram_compression="dense"``. RT0 is affine TET only. ``inv_chi`` is
    per solve, so one prepared operator serves several susceptibilities.
    """

    def __init__(self, mesh_path, *, backend, label_contract, mesh_check_receipt,
                 mesh_check_receipt_sha256, mesh_check_provenance, order=1, gram_aca=None,
                 gram_compression="aca", curve_order=1, image=None, image_cyclic=None,
                 image_cyclic_alternating=False, cyclic_periodic_boundaries=(), q2_deformation=None):
        if backend != "gpu":
            raise ValueError("This workflow requires explicit backend='gpu'; CPU uses radia.vim")
        if isinstance(order, bool) or order not in (0, 1, 2):
            raise ValueError("HDiv linear order must be 0 (RT0), 1 or 2")
        curve_order = _checked_curve_order(curve_order, order)
        aca = _checked_gram_compression(gram_compression, gram_aca)
        import numpy as np
        from .gpu import _native
        _hdiv_cuda = _native("_hdiv_cuda")
        self._lock = threading.RLock()
        self._operator = self._mass = None
        self._closed = False
        self.order = order
        self.family = None
        self.gram_compression = gram_compression
        self.gram_aca = None if aca is None else dict(aca)
        with _checked_mesh_snapshot(mesh_path, label_contract, mesh_check_receipt,
                                    mesh_check_receipt_sha256, mesh_check_provenance) as (path, mesh_hash, contract_hash):
            self.mesh_sha256 = mesh_hash
            self.label_contract_sha256 = contract_hash
            try:
                setup = _prepare_device_components(path, order, _hdiv_cuda, aca=aca, compression=gram_compression,
                                                   curve_order=curve_order, q2_deformation=q2_deformation,
                                                   images=_DeviceImages(image, image_cyclic, image_cyclic_alternating,
                                                                        cyclic_periodic_boundaries))
                self._gram_arguments = setup["charge_gram"]
                self.family = setup["family"]
                self.dofs = setup["dofs"]
                self.fespace, self._fe_owner = setup["fespace"], setup["fe_owner"]
                self._rhs = setup["rhs_operator"]
                self._mean = setup["mean_operator"]
                # The mass is built first: its live device budget check runs before the Gram.
                self._mass = _hdiv_cuda._ResidentMass.from_flat_tables(
                    **setup["evaluation_tables"], dofs=self.dofs, backend="gpu")
                self.mass_assembly = _checked_material_assembly(self._mass.assembly(), setup["evaluation_tables"],
                                                                mass=True)
                self._operator = _hdiv_cuda._ResidentHDiv.from_flat_charges(
                    **setup["charge_gram"], dofs=self.dofs, backend="gpu")
                self.assembly = _checked_device_assembly(self._operator.assembly(), setup["charge_gram"])
            except BaseException:
                self.close()
                raise

    def solve(self, applied_field, inv_chi, *, tolerance=1e-10, max_iterations=1000,
              preconditioner="mass-diagonal", degree=None, lower_fraction=None):
        """Solve an applied H (A/m), uniform or a CoefficientFunction; return owned coefficients and the mean magnetization.

        The native PCG recertifies the true residual of the original linear
        equation; Chebyshev settings are accepted only with a Chebyshev mode.
        """
        if preconditioner not in _LINEAR_PRECONDITIONERS:
            raise ValueError("Unsupported explicit linear material preconditioner")
        for value in (inv_chi, tolerance):
            if isinstance(value, bool) or not isinstance(value, numbers.Real) or not math.isfinite(value) or value <= 0:
                raise ValueError("inv_chi and tolerance must be finite and positive")
        if tolerance > 1e-6:
            raise ValueError("True residual tolerance must be at most 1e-6")
        if isinstance(max_iterations, bool) or not isinstance(max_iterations, numbers.Integral) or max_iterations < 1:
            raise ValueError("max_iterations must be a positive integer")
        import numpy as np
        with self._lock:
            if self._closed:
                raise RuntimeError("Prepared HDiv workflow has been closed")
            result = self._operator.solve_linear(self._mass, _applied_source(self.fespace, self._rhs, applied_field), float(inv_chi),
                                                 float(tolerance), int(max_iterations), preconditioner,
                                                 degree, lower_fraction)
            return dict(result, mean_magnetization=self._mean @ result["x"], backend="gpu",
                        mesh_sha256=self.mesh_sha256, order=self.order, dofs=self.dofs, field_units="A/m")

    def field(self, x, points):
        """H in A/m of a solved state at (points, 3) probes, images folded; a mixed mesh refuses.

        The device route of ``vim.FieldFromSolution(..., algorithm="direct")``: the charges
        are the resident charge map's product and only the state, geometry and probes cross.
        """
        import numpy as np
        with self._lock:
            if self._closed:
                raise RuntimeError("Prepared HDiv workflow has been closed")
            state = np.ascontiguousarray(x, dtype=np.float64)
            probes = np.ascontiguousarray(points, dtype=np.float64)
            if state.shape != (self.dofs,) or probes.ndim != 2 or probes.shape[1] != 3:
                raise ValueError("field requires the prepared state and (points, 3) probes: %r, %r"
                                 % (state.shape, probes.shape))
            return _checked_device_field(self._operator.field(state, probes), self._gram_arguments, state, probes)

    def stats(self):
        with self._lock:
            return dict(closed=self._closed, backend="gpu", dofs=getattr(self, "dofs", None),
                family=self.family, gram_route=_gram_route(self.family), gram_compression=self.gram_compression,
                gram_aca=None if self.gram_aca is None else dict(self.gram_aca),
                assembly=getattr(self, "assembly", None), mass_assembly=getattr(self, "mass_assembly", None),
                charge_budget=self.assembly["charge_budget"] if getattr(self, "assembly", None) else None,
                operator=self._operator.stats() if self._operator is not None else None,
                mass=self._mass.stats() if self._mass is not None else None)

    def close(self):
        with self._lock:
            try:
                if self._mass is not None:
                    self._mass.close()
            finally:
                if self._operator is not None:
                    self._operator.close()
                self._closed = True

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


_LINEAR_PRECONDITIONERS = ("mass-diagonal", "mass-cholesky", "mass-polynomial", "mass-polynomial-graph",
                           "mass-chebyshev", "mass-chebyshev-graph")
# TET geometry family per NGSolve curve order: affine TET, or isoparametric P2 TET
# (Mesh.Curve(2), the CPU curved route's geometry). A pure-HEX mesh is the
# family "hex" at either curve order: its Q2 lattice carries linear and Curve(2) maps.
_FAMILIES = {1: "affine-tet", 2: "curved-tet"}
_HEX = "hex"
_MIXED = "mixed"
_GRAM_ROUTES = {"affine-tet": "device-flat-tet", "curved-tet": "device-curved-tet", _HEX: "device-hex",
                _MIXED: "device-mixed"}


def _checked_curve_order(curve_order, order):
    """Explicit curve order; curved geometry carries HDiv orders 1 and 2 only, as the CPU curved routes."""
    if isinstance(curve_order, bool) or curve_order not in _FAMILIES:
        raise ValueError("curve_order must be 1 (linear geometry) or 2 (NGSolve Curve(2)): %r" % (curve_order,))
    if curve_order == 2 and order not in (1, 2):
        raise ValueError("Curved geometry requires HDiv order 1 or 2")
    return int(curve_order)


def _mesh_family(mesh, curve_order, order):
    """Device geometry family of a loaded TET, HEX or straight tet, pyramid and hex mesh; anything else fails loudly."""
    from radia.vim._vim import _volume_vertex_counts
    counts = _volume_vertex_counts(mesh)
    if counts == {4}:
        return _FAMILIES[curve_order]
    if counts == {8}:
        if order not in (1, 2):
            raise NotImplementedError("Device HEX route requires HDiv order 1 or 2; RT0 is affine TET only")
        return _HEX
    if counts <= {4, 5, 8}:
        if order not in (1, 2) or curve_order != 1:
            raise NotImplementedError("Device tet, pyramid and hex route requires HDiv order 1 or 2 on straight geometry")
        return _MIXED
    raise NotImplementedError("Device HDiv route requires TET, pyramid and HEX elements: element vertex counts %r"
                              % (sorted(counts),))


def _gram_route(family):
    return None if family is None else _GRAM_ROUTES[family]


def _strict_json(data):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError("Duplicate checked-mesh JSON key: " + key)
            result[key] = value
        return result
    def constant(value):
        raise ValueError("Nonfinite checked-mesh JSON: " + value)
    return json.loads(data, object_pairs_hook=pairs, parse_constant=constant)


def _verify_mesh_receipt(mesh_bytes, contract_bytes, receipt_bytes, expected_sha256, provenance):
    """Check pinned producer evidence, not repeat geometry analysis or authenticate consent.

    The caller obtains both pins from an independently approved artifact manifest.
    Package bootstrap may have loaded CPU libraries already. This gate precedes
    workflow-owned FE/native/GPU entry, not the import of the parent package.
    """
    if (not isinstance(expected_sha256, str) or not re.fullmatch(r"[0-9a-f]{64}", expected_sha256)
            or hashlib.sha256(receipt_bytes).hexdigest() != expected_sha256):
        raise ValueError("Checked-mesh external receipt SHA256 mismatch")
    value = _strict_json(receipt_bytes)
    if value.get("schema") != "radia.checked-vol-receipt.v1":
        raise ValueError("Unsupported checked-mesh receipt schema")
    if not isinstance(provenance, dict) or set(provenance) != {"checker", "argv", "ngsolve_version"}:
        raise ValueError("Externally approved checker provenance is required")
    for key in provenance:
        if value.get(key) != provenance[key]:
            raise ValueError("Checked-mesh approved provenance mismatch: " + key)
    checker = value["checker"]
    if (not isinstance(checker, dict) or set(checker) != {"distribution", "version", "record_sha256"}
            or checker["distribution"] != "cubit-mesh-export"
            or not isinstance(checker["version"], str) or checker["version"].lower() in ("", "unknown")
            or not re.fullmatch(r"[0-9a-f]{64}", checker["record_sha256"])):
        raise ValueError("Checked-mesh distribution identity is incomplete")
    argv = value["argv"]
    if (not isinstance(argv, list) or len(argv) != 9 or not all(isinstance(x, str) for x in argv)
            or argv[1:3] != ["-m", "cubit_mesh_export.check"]
            or argv[4:6] != ["--strict-labels", "--contract"] or argv[7:] != ["--format", "json"]
            or PureWindowsPath(argv[3]).name not in ("mesh.vol", "mesh.vol.gz")
            or PureWindowsPath(argv[6]).name != "labels.json"
            or PureWindowsPath(argv[3]).parent != PureWindowsPath(argv[6]).parent):
        raise ValueError("Checked-mesh canonical argv/default options mismatch")
    if value.get("returncode") != 0 or type(value.get("returncode")) is not int:
        raise ValueError("Checked-mesh checker execution failed")
    for key, data in (("mesh_sha256", mesh_bytes), ("contract_sha256", contract_bytes)):
        if value.get(key) != hashlib.sha256(data).hexdigest():
            raise ValueError("Checked-mesh content mismatch: " + key)
    try:
        raw_report = base64.b64decode(value["report_base64"], validate=True)
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("Checked-mesh raw report is invalid") from error
    if value.get("report_sha256") != hashlib.sha256(raw_report).hexdigest():
        raise ValueError("Checked-mesh raw report SHA256 mismatch")
    report = _strict_json(raw_report)
    if (report.get("schema") != "cubit-mesh-export.vol-check.v1" or report.get("passed") is not True
            or report.get("warnings") != [] or report.get("vol_file") != argv[3]):
        raise ValueError("Checked-mesh authoritative report did not pass")
    cad = report.get("cad_reference", {})
    if cad.get("available") is not False or cad.get("auto_discovered") is not False or report.get("json_file") is not None:
        raise ValueError("Checked-mesh receipt must not use unbound sibling CAD")
    for name in ("labels", "quality"):
        part = report.get(name)
        if not isinstance(part, dict) or part.get("passed") is not True or part.get("warnings") != []:
            raise ValueError("Checked-mesh authoritative component failed: " + name)
    labels = report["labels"]
    if labels.get("schema") != "radia.vol-label-contract.v1" or labels.get("strict_labels") is not True:
        raise ValueError("Checked-mesh strict label contract was not executed")
    contract = _strict_json(contract_bytes)
    if contract.get("schema") != "radia.vol-label-contract.v1":
        raise ValueError("Unsupported mesh label contract schema")
    # This compares the declared label sets; only the producer analyzes geometry.
    for category in ("required", "allowed"):
        for kind, names in contract.get(category, {}).items():
            expected = names.split(",") if isinstance(names, str) else names
            if not isinstance(expected, list) or {str(x).strip() for x in expected} != set(labels.get(category, {}).get(kind, [])):
                raise ValueError("Checked-mesh semantic contract mismatch: " + category + "/" + kind)
    ownership = report.get("boundary_domain_ownership")
    if not isinstance(ownership, dict) or ownership.get("issues") != []:
        raise ValueError("Checked-mesh boundary/domain ownership did not pass")
    runtime_ngsolve = metadata.version("ngsolve")  # Metadata only, no FE/native import.
    if (value["ngsolve_version"] != runtime_ngsolve
            or report.get("versions", {}).get("ngsolve") != runtime_ngsolve
            or report.get("versions", {}).get("cubit_mesh_export") != checker["version"]):
        raise ValueError("Checked-mesh checker/runtime ABI mismatch")
    return value


@contextmanager
def _checked_mesh_snapshot(mesh_path, label_contract, mesh_check_receipt,
                           mesh_check_receipt_sha256, mesh_check_provenance):
    path = Path(mesh_path).resolve(strict=True)
    contract = Path(label_contract).resolve(strict=True)
    mesh_bytes, contract_bytes = path.read_bytes(), contract.read_bytes()
    suffix = ".vol.gz" if path.name.lower().endswith(".vol.gz") else ".vol"
    if not path.name.lower().endswith((".vol", ".vol.gz")):
        raise ValueError("HDiv mesh input must be .vol or .vol.gz")
    receipt_bytes = Path(mesh_check_receipt).read_bytes()
    _verify_mesh_receipt(mesh_bytes, contract_bytes, receipt_bytes,
                         mesh_check_receipt_sha256, mesh_check_provenance)
    with tempfile.TemporaryDirectory(prefix="radia-hdiv-mesh-") as directory:
        snapshot = Path(directory) / ("mesh" + suffix)
        snapshot.write_bytes(mesh_bytes)  # Never re-read mutable caller files after verification.
        yield snapshot, hashlib.sha256(mesh_bytes).hexdigest(), hashlib.sha256(contract_bytes).hexdigest()


# Device-generated Gram and charge-map ledger: primitives and NGSolve tables
# in, three validity flags out.
_DEVICE_GRAM_UPLOADS = ("cell_vertices_m", "face_vertices_m", "charge_host", "charge_kind",
                        "charge_exponents", "tet_rule", "tri_rule", "tet_low_rule",
                        "tri_low_rule", "invalid_flag", "charge_map_cell_vertices_m",
                        "charge_map_face_vertices_m", "charge_cell_tables", "charge_face_tables",
                        "charge_cell_class", "charge_cell_dofs", "charge_face_class",
                        "charge_face_dofs", "charge_map_invalid_flag")
_CHARGE_TABLE_KEYS = ("cell_class", "cell_dofs", "cell_tables", "cell_exponents",
                      "face_class", "face_dofs", "face_tables", "face_exponents")
# Curved P2 additions: corner numbers for touching, and the mass rules of the
# per-host device Gram solve of the charge map.
_CURVED_CHARGE_KEYS = ("cell_numbers", "face_numbers", "cell_mass_rule", "face_mass_rule")
# Mixed mesh: charge map tables per class and DOFs per host, in host order.
_MIXED_CHARGE_KEYS = ("host_class", "class_info", "host_dofs", "class_tables")
_CURVED_GRAM_UPLOADS = ("curve_gauss_points", "curve_gauss_weights", "cell_vertex_numbers", "face_vertex_numbers",
                        "charge_cell_mass_rule", "charge_face_mass_rule", "charge_cell_exponents",
                        "charge_face_exponents")
_RULE_ROW_BYTES = 32  # FlatRulePoint: xi[3], weight as float64.
# Device ACA of the charge Gram. The preset carries the CPU charge-Gram H-matrix
# parameters (radia.vim build_charge_gram eps, leaf size, eta); it is a starting
# point, not a tuned choice. eps is the ACA stopping parameter, not an operator
# error bound. A leaf reaching max_rank unconverged is stored as an exact dense
# leaf of the same operator and counted (rank_cap="dense") or raises ("fail").
_GRAM_COMPRESSION = "aca"
GRAM_ACA_PRESET = dict(eps=1e-10, eta=2.0, leaf_size=64, max_rank=200, rank_cap="dense")
# Metadata rows that exist only when the partition has admissible blocks.
_ACA_OPTIONAL_ROWS = (("aca_leaf_spans", "H2D"), ("aca_leaf_status", "D2H"), ("aca_active_count", "D2H"))
_ACA_COUNTERS = ("elements", "clusters", "dense_leaves", "lowrank_leaves", "admissible_blocks",
                 "rank_cap_promotions", "storage_promotions", "zero_rank_leaves", "rank_sum", "max_rank", "aca_steps",
                 "aca_batches", "factor_coefficients", "dense_coefficients")


def _checked_gram_compression(compression, aca):
    """Return the checked ACA parameters, or ``None`` for the explicit exact dense Gram."""
    if compression == "dense":
        if aca is not None:
            raise ValueError("gram_aca is rejected for gram_compression='dense'")
        return None
    if compression != _GRAM_COMPRESSION:
        raise ValueError("gram_compression must be explicit 'aca' or 'dense': %r" % (compression,))
    return _checked_gram_aca(aca)


def _checked_gram_aca(aca):
    """Return an owned copy of explicit device ACA parameters; ``None`` is the preset.

    Ranges match the native partition and compression checks, so a bad value
    fails before any mesh staging or device work.
    """
    if aca is None:
        return dict(GRAM_ACA_PRESET)
    if not isinstance(aca, dict) or set(aca) != set(GRAM_ACA_PRESET):
        raise ValueError("gram_aca requires exactly eps, eta, leaf_size, max_rank, rank_cap: %r" % (aca,))

    def real(name, valid):
        value = aca[name]
        if isinstance(value, bool) or not isinstance(value, numbers.Real) or not math.isfinite(value) or not valid(value):
            raise ValueError("gram_aca %s out of range: %r" % (name, value))
        return float(value)

    def integer(name, upper):
        value = aca[name]
        if isinstance(value, bool) or not isinstance(value, numbers.Integral) or not 1 <= value <= upper:
            raise ValueError("gram_aca %s must be an integer in [1,%d]: %r" % (name, upper, value))
        return int(value)
    if aca["rank_cap"] not in ("dense", "fail"):
        raise ValueError("gram_aca rank_cap must be 'dense' or 'fail': %r" % (aca["rank_cap"],))
    return dict(eps=real("eps", lambda v: 0 < v < 1), eta=real("eta", lambda v: v > 0),
                leaf_size=integer("leaf_size", 4096), max_rank=integer("max_rank", 1024),
                rank_cap=str(aca["rank_cap"]))


def _expected_device_ledger(arguments):
    """Exact (payload, direction) -> bytes of the Gram-construction ledger.

    Covers the Gram and the device charge map only; mass, evaluation and
    preconditioner uploads of the full workflow are outside this ledger.
    Empty low rules (far split disabled) are absent rows, never zero-byte rows.
    The charge map crosses as NGSolve tables, classes and DOF numbers; its
    exponents are host-side layout checks and never uploaded. The resident
    plan size is owned by the native side and required positive. The ACA
    route adds its charge permutation, dense-leaf spans and integer leaf
    rank/status and step counts; no Gram, sigma or factor value crosses.
    The curved P2 family uploads P2 nodes in place of vertices, plus its
    Duffy rule, corner numbers, charge mass rules and exponents. The HEX
    family uploads Q2 lattices, host spans, host-grouped charges and
    exponents (its 1D rules are kernel arguments, its charge map takes no
    geometry) and reads back its pair-table sizes and block flags.
    """
    if arguments.get("family") == _MIXED:
        return _expected_mixed_ledger(arguments)
    if arguments.get("family") == _HEX:
        return _expected_hex_ledger(arguments)
    curved = arguments.get("family") == "curved-tet"
    sizes = {"cell_vertices_m": arguments["cell_vertices"].nbytes,
             "face_vertices_m": arguments["face_vertices"].nbytes,
             "charge_host": arguments["charge_host"].nbytes,
             "charge_kind": arguments["charge_kind"].nbytes,
             "charge_exponents": arguments["charge_exponents"].nbytes,
             "tet_rule": len(arguments["tet_rule"])*_RULE_ROW_BYTES,
             "tri_rule": len(arguments["tri_rule"])*_RULE_ROW_BYTES,
             "tet_low_rule": len(arguments["tet_low_rule"])*_RULE_ROW_BYTES,
             "tri_low_rule": len(arguments["tri_low_rule"])*_RULE_ROW_BYTES,
             "invalid_flag": 4,
             "charge_map_cell_vertices_m": arguments["cell_vertices"].nbytes,
             "charge_map_face_vertices_m": arguments["face_vertices"].nbytes,
             "charge_map_invalid_flag": 4}
    sizes.update(("charge_" + name, arguments[name].nbytes) for name in _CHARGE_TABLE_KEYS
                 if curved or not name.endswith("exponents"))
    uploads = _DEVICE_GRAM_UPLOADS
    if curved:
        sizes.update(cell_p2_nodes_m=sizes.pop("cell_vertices_m"), face_p2_nodes_m=sizes.pop("face_vertices_m"),
                     curve_gauss_points=arguments["curve_gl"].nbytes, curve_gauss_weights=arguments["curve_gw"].nbytes,
                     cell_vertex_numbers=arguments["cell_numbers"].nbytes,
                     face_vertex_numbers=arguments["face_numbers"].nbytes,
                     charge_cell_mass_rule=arguments["cell_mass_rule"].nbytes,
                     charge_face_mass_rule=arguments["face_mass_rule"].nbytes)
        uploads = ("cell_p2_nodes_m", "face_p2_nodes_m") + uploads[2:] + _CURVED_GRAM_UPLOADS
    expected = {(name, "H2D"): sizes[name] for name in uploads if sizes[name]}
    expected.update({("prepare_flag", "D2H"): 4, ("energy_flag", "D2H"): 4,
                     ("charge_map_flag", "D2H"): 4, ("resident_plan", "H2D"): None})
    if arguments.get("image_masks") is not None:
        expected[("image_transforms", "H2D")] = _IMAGE_ROW_BYTES*arguments["image_masks"].size
    if arguments["compression"] == "aca":
        expected.update({("aca_permutation", "H2D"): 8*arguments["charge_host"].size,
                         ("aca_dense_spans", "H2D"): None, ("dense_flag", "D2H"): 4})
    return expected


# One native HexImage row: flips, forward and inverse rotation, sign and flag, all float64.
_IMAGE_ROW_BYTES = 72


def _expected_hex_ledger(arguments):
    hosts = arguments["cell_class"].size + arguments["face_class"].size
    expected = {(name, "H2D"): size for name, size in (
        ("cell_q2_nodes_m", arguments["cell_vertices"].nbytes), ("face_q2_nodes_m", arguments["face_vertices"].nbytes),
        ("hex_host_first", 4*hosts), ("hex_host_count", 4*hosts), ("charge_host", arguments["charge_host"].nbytes),
        ("charge_exponents", arguments["charge_exponents"].nbytes), ("hex_invalid_flag", 4), ("invalid_flag", 4),
        ("charge_map_invalid_flag", 4))}
    expected.update((("charge_" + name, "H2D"), arguments[name].nbytes) for name in _CHARGE_TABLE_KEYS
                    if not name.endswith("exponents"))
    expected.update({("prepare_flag", "D2H"): 4, ("hex_pair_count", "D2H"): 16, ("hex_block_values", "D2H"): 16,
                     ("hex_block_flag", "D2H"): 4, ("hex_symmetrize_flag", "D2H"): 4, ("energy_flag", "D2H"): 4,
                     ("charge_map_flag", "D2H"): 4, ("resident_plan", "H2D"): None})
    if arguments.get("image_masks") is not None:
        expected[("image_transforms", "H2D")] = _IMAGE_ROW_BYTES*arguments["image_masks"].size
    if arguments["compression"] == "aca":
        expected.update({("aca_permutation", "H2D"): 8*arguments["charge_host"].size,
                         ("aca_dense_spans", "H2D"): None, ("dense_flag", "D2H"): 4})
    return expected


def _expected_mixed_ledger(arguments):
    """The HEX route's lattices, spans and flags, plus host types and the per-class charge map tables."""
    tables = arguments["mixed_charges"]
    hosts = arguments["host_types"].size
    expected = {(name, "H2D"): size for name, size in (
        ("cell_q2_nodes_m", arguments["cell_vertices"].nbytes), ("face_q2_nodes_m", arguments["face_vertices"].nbytes),
        ("hex_host_first", 4*hosts), ("hex_host_count", 4*hosts), ("charge_host", arguments["charge_host"].nbytes),
        ("charge_exponents", arguments["charge_exponents"].nbytes), ("host_types", arguments["host_types"].nbytes),
        ("hex_invalid_flag", 4), ("invalid_flag", 4), ("charge_map_invalid_flag", 4),
        ("charge_class_tables", tables["class_tables"].nbytes), ("charge_class_info", tables["class_info"].nbytes),
        ("charge_host_class", tables["host_class"].nbytes), ("charge_host_dofs", tables["host_dofs"].nbytes),
        ("charge_row_host", 4*arguments["charge_host"].size), ("charge_host_first", 4*hosts),
        ("charge_host_dof_start", 4*hosts), ("charge_host_entry", 8*hosts))}
    expected.update({("prepare_flag", "D2H"): 4, ("hex_pair_count", "D2H"): 16, ("hex_block_values", "D2H"): 16,
                     ("hex_block_flag", "D2H"): 4, ("hex_symmetrize_flag", "D2H"): 4, ("energy_flag", "D2H"): 4,
                     ("charge_map_flag", "D2H"): 4, ("resident_plan", "H2D"): None})
    if arguments.get("image_masks") is not None:
        expected[("image_transforms", "H2D")] = _IMAGE_ROW_BYTES*arguments["image_masks"].size
    if arguments["compression"] == "aca":
        expected.update({("aca_permutation", "H2D"): 8*arguments["charge_host"].size,
                         ("aca_dense_spans", "H2D"): None, ("dense_flag", "D2H"): 4})
    return expected


def _checked_device_assembly(assembly, arguments):
    """Reject any Gram assembly ledger other than the exact primitive-only device route.

    A host-built charge-map row is rejected like any other unexpected payload.
    """
    charges = int(arguments["charge_host"].size)
    if assembly.get("route") != _gram_route(arguments["family"]) \
            or assembly.get("charges") != charges:
        raise RuntimeError("Unexpected HDiv Gram assembly route: %r" % (assembly.get("route"),))
    if assembly.get("compression") != arguments["compression"]:
        raise RuntimeError("HDiv Gram compression %r, requested %r"
                           % (assembly.get("compression"), arguments["compression"]))
    counters = assembly.get("hmatrix")
    if arguments["compression"] == "aca" and (set(counters) != set(_ACA_COUNTERS) or not all(
            type(v) is int and v >= 0 for v in counters.values())):
        raise RuntimeError("HDiv ACA Gram counters rejected: %r" % (counters,))
    assembly["charge_budget"] = _checked_charge_budget(assembly.get("charge_budget"), arguments["compression"])
    expected = _expected_device_ledger(arguments)
    optional = set(_ACA_OPTIONAL_ROWS) if arguments["compression"] == "aca" else set()
    seen = {}
    for row in assembly["transfer_ledger"]:
        key, size = (row["payload"], row["direction"]), row["bytes"]
        if key in optional and key not in seen and type(size) is int and size > 0:
            seen[key] = size
            continue
        if key in seen or key not in expected or type(size) is not int or size < 0:
            raise RuntimeError("HDiv device Gram transfer ledger row rejected: %r" % (row,))
        want = expected[key]
        if (want is None and size <= 0) or (want is not None and size != want):
            raise RuntimeError("HDiv device Gram ledger bytes for %r: actual %r, expected %r"
                               % (key, size, "positive" if want is None else want))
        seen[key] = size
    if set(seen) - optional != set(expected):
        raise RuntimeError("HDiv device Gram transfer ledger missing rows: %r"
                           % (sorted(set(expected) - set(seen)),))
    return assembly


# The Gram-construction rows a TET field evaluation repeats: geometry, charge layout and rules, never a matrix.
_FIELD_GRAM_ROWS = ("cell_vertices_m", "face_vertices_m", "cell_p2_nodes_m", "face_p2_nodes_m", "charge_host", "charge_kind",
                    "charge_exponents", "tet_rule", "tri_rule", "tet_low_rule", "tri_low_rule", "invalid_flag",
                    "curve_gauss_points", "curve_gauss_weights", "cell_vertex_numbers", "face_vertex_numbers", "image_transforms")


def _expected_tet_field_ledger(result, arguments):
    curved = arguments["family"] == "curved-tet"
    if result["route"] != ("device-curved-tet-direct" if curved else "device-flat-tet-polynomial"):
        raise RuntimeError("Unexpected device field route: %r" % (result["route"],))
    expected = {key: size for key, size in _expected_device_ledger(arguments).items() if key[1] == "H2D" and key[0] in _FIELD_GRAM_ROWS}
    hosts = (arguments["cell_vertices"].size//(30 if curved else 12)) + (arguments["face_vertices"].size//(18 if curved else 9))
    expected.update({("field_host_first", "H2D"): 8*(hosts + 1), ("field_host_order", "H2D"): 4*arguments["charge_host"].size,
                     ("field_polynomial_flag", "D2H"): 4})
    return expected


def _expected_hex_field_ledger(result, arguments):
    hosts = arguments["cell_class"].size + arguments["face_class"].size
    expected = {(name, "H2D"): size for name, size in (
        ("cell_q2_nodes_m", arguments["cell_vertices"].nbytes), ("face_q2_nodes_m", arguments["face_vertices"].nbytes),
        ("hex_host_first", 4*hosts), ("hex_host_count", 4*hosts), ("charge_host", arguments["charge_host"].nbytes),
        ("charge_exponents", arguments["charge_exponents"].nbytes), ("hex_invalid_flag", 4))}
    if arguments.get("image_masks") is not None:
        expected[("image_transforms", "H2D")] = _IMAGE_ROW_BYTES*arguments["image_masks"].size
    expected[("field_host_affine", "D2H")] = hosts
    if result["route"] == "device-hex-polynomial":
        expected[("field_polynomial_flag", "D2H")] = 4
    elif result["route"] != "device-hex-cloud":
        raise RuntimeError("Unexpected device field route: %r" % (result["route"],))
    return expected


def _checked_device_field(result, arguments, state, probes):
    """Reject a field evaluation whose ledger holds anything but geometry, the state, probes and values."""
    expected = (_expected_hex_field_ledger if arguments["family"] == _HEX else _expected_tet_field_ledger)(result, arguments)
    expected.update({("field_far_gauss", "H2D"): None, ("field_state", "H2D"): state.nbytes, ("field_probes", "H2D"): probes.nbytes,
                     ("prepare_flag", "D2H"): 4, ("field_flag", "D2H"): 4,
                     ("field_values", "D2H"): probes.nbytes})
    seen = {}
    for row in result["transfer_ledger"]:
        key, size = (row["payload"], row["direction"]), row["bytes"]
        want = expected.get(key, -1)
        if key in seen or want == -1 or type(size) is not int or (want is None and size <= 0) \
                or (want is not None and size != want):
            raise RuntimeError("Device field ledger row rejected: %r" % (row,))
        seen[key] = size
    if set(expected) - set(seen):
        raise RuntimeError("Device field ledger missing rows: %r" % (sorted(set(expected) - set(seen)),))
    return result


def _checked_charge_exponents(kind, exponents, family):
    """Bound every exponent before any degree sum (face z = 0).

    TET: total degree, cell <= 1 and face <= 2. HEX and mixed: per-axis exponents <= 2,
    and on mixed meshes the rational pyramid monomial coded (1, 1, -1).
    """
    import numpy as np
    if exponents.shape != (3*kind.size,) or not np.isin(kind, (0, 1)).all():
        raise ValueError("Charge kinds/exponents disagree: kinds %r, exponents %r"
                         % (kind.shape, exponents.shape))
    expo = exponents.reshape(-1, 3).astype(np.int64)
    axis = family in (_HEX, _MIXED)
    limit = np.full((kind.size, 1), 2) if axis else np.where(kind == 0, 1, 2)[:, None]
    bad = (expo < 0) | (expo > limit)
    if family == _MIXED:
        bad[(kind == 0) & (expo == (1, 1, -1)).all(axis=1)] = False
    bad[:, 2] |= (kind == 1) & (expo[:, 2] != 0)
    bad = bad.any(axis=1)
    if not axis:
        bad |= expo.sum(axis=1) > limit[:, 0]
    if bad.any():
        first = int(np.flatnonzero(bad)[0])
        raise ValueError("Device charge Gram supports TET cell degree <= 1, TET face degree <= 2 and HEX axis "
                         "degree <= 2: charge %d kind %d exponents %r" % (first, int(kind[first]), expo[first].tolist()))


def _rule_rows(points, weights, width):
    """(points, 4) rows of xi and weight, the native FlatRulePoint layout."""
    import numpy as np
    points = np.asarray(points, dtype=np.float64).reshape(-1, width)
    weights = np.asarray(weights, dtype=np.float64).reshape(-1)
    if len(points) != len(weights):
        raise ValueError("Reference rule points and weights disagree")
    rows = np.zeros((len(points), 4))
    rows[:, :width], rows[:, 3] = points, weights
    return rows


# Two evaluations of the same NGSolve P2 map (vectorized CPU-route trafo and the
# native staging's mapped points) agree to roundoff, never to a modelling error.
_P2_GEOMETRY_ROUNDOFF = 1e-12


class _DeviceImages:
    """The caller's image request, checked against the mesh with the CPU contract before any staging."""

    def __init__(self, image=None, image_cyclic=None, image_cyclic_alternating=False, cyclic_periodic_boundaries=()):
        self.image, self.cyclic = image, image_cyclic
        self.terms = _device_image_terms(image, image_cyclic, image_cyclic_alternating)
        if isinstance(cyclic_periodic_boundaries, str):
            raise ValueError("cyclic_periodic_boundaries must be a pair of labels: %r" % (cyclic_periodic_boundaries,))
        self.periodic = tuple(str(name) for name in (cyclic_periodic_boundaries or ()))
        if self.periodic:
            if image_cyclic is None:
                raise ValueError("cyclic_periodic_boundaries requires image_cyclic=N")
            if len(self.periodic) != 2 or self.periodic[0] == self.periodic[1]:
                raise ValueError("cyclic_periodic_boundaries must name two distinct rotation-related sector faces")
            if image_cyclic_alternating:
                raise NotImplementedError("Device connected antiperiodic HDiv sectors are not wired")

    def reduced_space(self, mesh, order):
        """The CPU route's ``Compress(Periodic(HDiv))`` of the identified sector faces, after its mesh checks."""
        from radia.vim._solve import _hdiv_space_with_image_constraints
        return _hdiv_space_with_image_constraints(mesh, order, [], self.periodic)[0]

    def checked_terms(self, mesh, order, family):
        """Native Gram terms, or ``None``, after the ``vim.Solve`` mesh checks of a space without a compressed trace."""
        if self.terms is None:
            return None
        if self.periodic and family != _HEX:
            raise NotImplementedError("A connected periodic sector needs the pure HEX family, not %r" % (family,))
        from radia.vim import _image
        from radia.vim._solve import _hdiv_space_with_image_constraints
        if self.cyclic is not None and not self.periodic:
            labels = [name for name in mesh.GetBoundaries() if "periodic" in name.lower() or "cyclic" in name.lower()]
            if labels:
                raise NotImplementedError("Device cyclic route has an unidentified azimuthal seam: the mesh has "
                                          "boundary labels %r; pass cyclic_periodic_boundaries" % (labels,))
        planes = _image.parse_image_string(self.image) if self.image is not None else []
        _, constrained, _ = _hdiv_space_with_image_constraints(mesh, order, list(planes), self.periodic)
        if constrained:
            raise NotImplementedError("Device solve with a '+' mirror plane needs the compressed normal trace of %d "
                                      "DOFs, which the device tables do not carry; the '-' planes and rotations do not"
                                      % (len(constrained),))
        return self.terms


def _device_image_terms(image=None, image_cyclic=None, image_cyclic_alternating=False):
    """Mirror and cyclic terms of the device Gram, in the ``vim.Solve`` order and semantics.

    Returns ``None`` without images, else the int32 masks (bits 1 x, 2 y, 4 z), float64 signs
    and, for a cyclic list, the float64 rotation angles about +z. A mirror and a cyclic list
    together are refused as on the CPU route; the native binding repeats the closure checks.
    """
    import numpy as np
    from radia.vim import _image
    masks, signs, angles = [], [], None
    if image is not None:
        for axes, sign in _image.image_group(_image.parse_image_string(image)):
            masks.append(int(sum(1 << axis for axis in axes)))
            signs.append(float(sign))
    if image_cyclic is not None:
        if isinstance(image_cyclic, bool) or not isinstance(image_cyclic, numbers.Integral) or image_cyclic < 2:
            raise ValueError("image_cyclic must be an N-fold count >= 2: %r" % (image_cyclic,))
        if image_cyclic_alternating and image_cyclic % 2:
            raise ValueError("image_cyclic_alternating needs an even pole count: %r" % (image_cyclic,))
        if masks:
            raise NotImplementedError("Mirror image planes and image_cyclic cannot be combined, as on the CPU route")
        fold = int(image_cyclic)
        masks = [0]*(fold - 1)
        signs = [-1.0 if image_cyclic_alternating and k % 2 else 1.0 for k in range(1, fold)]
        angles = [2.0*math.pi*k/fold for k in range(1, fold)]
    elif image_cyclic_alternating:
        raise ValueError("image_cyclic_alternating requires image_cyclic")
    if not masks:
        return None
    terms = dict(image_masks=np.ascontiguousarray(masks, dtype=np.int32),
                 image_signs=np.ascontiguousarray(signs, dtype=np.float64))
    if angles is not None:
        terms["image_rotations"] = np.ascontiguousarray(angles, dtype=np.float64)
    return terms


def _device_gram_arguments(tables, inputs, geometry, aca, compression=_GRAM_COMPRESSION, images=None):
    """Map the staged TET charge inputs and NGSolve tables to the device binding.

    The Gram is compressed by device ACA with the caller's checked parameters
    (``_checked_gram_aca``); only the geometric partition is built on the
    host. The inner subtraction rule is unused by the exact analytic BDM1/BDM2
    potentials, as on the CPU; higher degrees fail loudly on the device.
    The table staging must reproduce the Gram geometry and charge order
    (``geometry`` is its cell and exterior-face coordinates): exactly for
    affine vertices; for curved P2 nodes to roundoff, with identical corner
    numbers. The curved family carries the CPU curved route's rules
    (``_curved_gram_rules``), far factor and Duffy rule, and the device
    consumes the staged P2 nodes that also build its charge map. The HEX
    family carries the CPU hex route's 1D rules and bands
    (``_prepare_gpu_hex_charge_inputs``) and consumes the staged Q2
    lattices, which must equal the CPU lattices to roundoff. The mixed
    family (``_prepare_gpu_mixed_charge_inputs``) consumes the staged host
    lattices the same way, with their host types and per-class charge tables.
    ``images`` (``_device_image_terms``) folds mirror and cyclic terms on every family.
    """
    import numpy as np
    family = inputs.get("family", "affine-tet")
    if len(inputs.get("image_masks", ())):
        raise NotImplementedError("Device charge Gram takes its image terms from the request, not the CPU inputs: family %r" % (family,))
    curved = family == "curved-tet"
    common = dict(charge_host=np.ascontiguousarray(inputs["charge_host"], dtype=np.int32),
                  charge_kind=np.ascontiguousarray(inputs["charge_kind"], dtype=np.int32),
                  charge_exponents=np.ascontiguousarray(inputs["charge_expo"], dtype=np.int32),
                  compression=compression, aca=_checked_gram_compression(compression, aca))
    cells_m, faces_m = geometry
    if curved or family in (_HEX, _MIXED):
        cells_py = np.ascontiguousarray(inputs["cell_nodes"], dtype=np.float64).reshape(-1)
        faces_py = np.ascontiguousarray(inputs["face_nodes"], dtype=np.float64).reshape(-1)
        if cells_py.size != (30 if curved else 81)*int(inputs["n_el"]) or cells_m.shape != cells_py.shape \
                or faces_m.shape != faces_py.shape:
            raise ValueError("Staged and CPU %s charge geometry disagree in size" % family)
        scale = max(float(np.abs(cells_py).max()), float(np.abs(faces_py).max()), np.finfo(float).tiny)
        deviation = max(float(np.abs(cells_m - cells_py).max()), float(np.abs(faces_m - faces_py).max()))
        if deviation > _P2_GEOMETRY_ROUNDOFF*scale:
            raise ValueError("Staged %s geometry and the CPU route disagree: %r m at scale %r m"
                             % (family, deviation, scale))
    if family in (_HEX, _MIXED):
        rules = inputs["rules"]
        rows = np.concatenate([np.c_[np.full(len(rules[name][0]), index, dtype=np.float64),
                                     np.asarray(rules[name][0], dtype=np.float64),
                                     np.asarray(rules[name][1], dtype=np.float64)]
                               for index, name in enumerate(("out", "inner", "pair", "pair_affine", "near", "inner_self"))
                               if name in rules])
        empty = np.zeros((0, 4))
        arguments = dict(common, cell_vertices=cells_m, face_vertices=faces_m, tet_rule=empty, tri_rule=empty,
                         tet_low_rule=empty, tri_low_rule=empty, far_factor=1e30, family=family,
                         hex_rules=np.ascontiguousarray(rows),
                         hex_bands=np.array([inputs["near_grade"], inputs["far_inner"]], dtype=np.float64))
    elif curved:
        numbers = (np.ravel(tables["cell_numbers"]), np.ravel(tables["face_numbers"]))
        if not (np.array_equal(numbers[0], inputs["cell_vertices"]) and np.array_equal(numbers[1], inputs["face_vertices"])):
            raise ValueError("Staged and CPU curved corner numbers disagree")
        rules = inputs["rules"]
        arguments = dict(common, cell_vertices=cells_m, face_vertices=faces_m,
                         tet_rule=_rule_rows(*rules["tet"], 3), tri_rule=_rule_rows(*rules["tri"], 2),
                         tet_low_rule=_rule_rows(*rules["tet_lo"], 3), tri_low_rule=_rule_rows(*rules["tri_lo"], 2),
                         far_factor=float(inputs["ho_far_factor"]), family="curved-tet",
                         cell_numbers=np.ascontiguousarray(numbers[0], dtype=np.int32),
                         face_numbers=np.ascontiguousarray(numbers[1], dtype=np.int32),
                         curve_gl=np.ascontiguousarray(rules["gl"], dtype=np.float64),
                         curve_gw=np.ascontiguousarray(rules["gw"], dtype=np.float64),
                         cell_mass_rule=tables["cell_mass_rule"], face_mass_rule=tables["face_mass_rule"])
    else:
        cells = np.ascontiguousarray(inputs["cell_verts"], dtype=np.float64).reshape(-1)
        faces = np.ascontiguousarray(inputs["face_verts"], dtype=np.float64).reshape(-1)
        if cells.size != 12*int(inputs["n_el"]):
            raise ValueError("Charge cell geometry does not match n_el")

        def rule(kind, suffix, width):
            return _rule_rows(inputs["ref_%s_pts%s" % (kind, suffix)], inputs["ref_%s_w%s" % (kind, suffix)], width)

        far = "ho_far_factor" in inputs and bool(np.isfinite(inputs["ho_far_factor"]))
        arguments = dict(common, cell_vertices=cells, face_vertices=faces, family=family,
                         tet_rule=rule("tet", "", 3), tri_rule=rule("tri", "", 2),
                         tet_low_rule=rule("tet", "_lo", 3) if far else np.zeros((0, 4)),
                         tri_low_rule=rule("tri", "_lo", 2) if far else np.zeros((0, 4)),
                         far_factor=float(inputs["ho_far_factor"]) if far else 1e30)
        if not (np.array_equal(cells_m, cells) and np.array_equal(faces_m, faces)):
            raise ValueError("Charge tables and Gram geometry disagree")
    _checked_charge_exponents(arguments["charge_kind"], arguments["charge_exponents"], family)
    if images is not None:
        arguments.update(images)
    if family == _MIXED:
        staged = (tables["charge_host"], tables["charge_kind"], tables["charge_exponents"], tables["host_types"])
        cpu = (arguments["charge_host"], arguments["charge_kind"], arguments["charge_exponents"], inputs["host_types"])
        if not all(np.array_equal(a, b) for a, b in zip(staged, cpu)):
            raise ValueError("Charge tables and Gram charge order disagree")
        arguments.update(_EMPTY_CHARGE_TABLES(), host_types=tables["host_types"],
                         mixed_charges={name: tables[name] for name in _MIXED_CHARGE_KEYS})
        return arguments
    ne, nf = tables["cell_class"].size, tables["face_class"].size
    cm, fm = len(tables["cell_exponents"]), len(tables["face_exponents"])
    face_expo = np.c_[tables["face_exponents"], np.zeros(fm, dtype=np.int32)]
    layout = (np.r_[np.repeat(np.arange(ne), cm), np.repeat(np.arange(nf), fm)],
              np.r_[np.zeros(ne*cm, dtype=np.int32), np.ones(nf*fm, dtype=np.int32)],
              np.r_[np.tile(np.ravel(tables["cell_exponents"]), ne), np.tile(face_expo.ravel(), nf)])
    staged = (arguments["charge_host"], arguments["charge_kind"], arguments["charge_exponents"])
    if not all(np.array_equal(a, b) for a, b in zip(layout, staged)):
        raise ValueError("Charge tables and Gram charge order disagree")
    arguments.update((name, tables[name]) for name in _CHARGE_TABLE_KEYS)
    return arguments


def _EMPTY_CHARGE_TABLES():
    """Uniform-host tables, unused on a mixed mesh: its hosts take the per-class tables."""
    import numpy as np
    return dict(cell_class=np.zeros(0, dtype=np.int32), cell_dofs=np.zeros((0, 0), dtype=np.int32),
                cell_tables=np.zeros((0, 0, 0)), cell_exponents=np.zeros((0, 3), dtype=np.int32),
                face_class=np.zeros(0, dtype=np.int32), face_dofs=np.zeros((0, 0), dtype=np.int32),
                face_tables=np.zeros((0, 0, 0)), face_exponents=np.zeros((0, 2), dtype=np.int32))


def _prepare_device_components(path, order, native, *, aca, compression=_GRAM_COMPRESSION, curve_order=1, images=None,
                               q2_deformation=None):
    """FE operators plus the primitive device-Gram payload; no CPU Gram or H-matrix.

    ``curve_order=2`` curves the FE mesh and the native staging mesh with
    NGSolve ``Curve(2)``; sources, mean field, tables and geometry then all
    refer to the curved mesh. The geometry family (``family``) follows from
    the mesh: pure TET by curve order, pure HEX as ``hex``, and a straight mesh
    with pyramids, or with both tets and hexes, as ``mixed``.

    ``native`` is the caller's verified ``_hdiv_cuda`` module: one preparation
    verifies the native and CUDA library closure once, not per staging helper.
    The material tables pass the native host-only scope check before any
    charge staging or Gram construction, so an unsupported material rule
    fails before the costly device work. The check returns the material
    device byte budget (``material_budget``); the live free-memory check
    runs when the material context is constructed.

    ``images`` carries the mirror and cyclic terms and, for a connected sector,
    the labels of its identified faces: the tables then take the periodic
    wrapper's element orientation and the ``Compress(Periodic(HDiv))`` numbering,
    and those faces leave the charge skin. ``q2_deformation`` maps a HEX mesh
    with the caller's Curve(2) + SetDeformation; its lattices replace the
    staged ones and sources and mean field refer to the deformed mesh.
    """
    import ngsolve as ng
    from radia.vim._vim import (_curve_mesh, _prepare_gpu_charge_inputs, _prepare_gpu_curved_charge_inputs,
                                _prepare_gpu_hex_charge_inputs, _prepare_gpu_mixed_charge_inputs)
    aca = _checked_gram_compression(compression, aca)
    curve_order = _checked_curve_order(curve_order, order)
    mesh = ng.Mesh(str(path))
    family = _mesh_family(mesh, curve_order, order)
    image_terms = None if images is None else images.checked_terms(mesh, order, family)
    periodic = () if images is None else images.periodic
    if q2_deformation is not None and (family != _HEX or curve_order != 2):
        raise NotImplementedError("A Q2 deformation requires the pure HEX family and curve_order=2")
    # A deformed mesh stages on its straight lattices: only the geometry is replaced below.
    staged_curve = 1 if q2_deformation is not None else curve_order
    if staged_curve == 2:
        _curve_mesh(mesh, 2)
    fes = ng.HDiv(mesh, order=order)
    evaluation = _device_evaluation_tables(path, fes, native, staged_curve, family, periodic)
    numbering = _reduced_numbering(fes, images.reduced_space(mesh, order), mesh) if periodic else None
    if periodic:
        for group in evaluation["tables"]:
            group.update(element_dofs=numbering[group["element_dofs"]], periodic_trace=True)
    budget = _checked_material_budget(native._ResidentBHMaterial.validate_flat_tables(
        **evaluation, dofs=fes.ndof if numbering is None else int(numbering.max()) + 1))
    if family == _MIXED:
        inputs = _prepare_gpu_mixed_charge_inputs(fes, backend="gpu")
    elif family == _HEX:
        inputs = _prepare_gpu_hex_charge_inputs(fes, backend="gpu", excluded_boundaries=periodic)
    elif curve_order == 2:
        inputs = _prepare_gpu_curved_charge_inputs(fes, backend="gpu")
    else:
        _, inputs = _prepare_gpu_charge_inputs(fes, backend="gpu")
    tables, geometry = _device_charge_tables(path, fes, native, staged_curve, family, periodic)
    arguments = _device_gram_arguments(tables, inputs, geometry, aca, compression, image_terms)
    if periodic:
        arguments.update(cell_dofs=numbering[arguments["cell_dofs"]], face_dofs=numbering[arguments["face_dofs"]],
                         periodic_trace=True)
    deformation = None if q2_deformation is None else _deformed_lattices(mesh, order, periodic, q2_deformation, arguments, evaluation)
    space = images.reduced_space(mesh, order) if periodic else (fes if deformation is None else ng.HDiv(mesh, order=order))
    return dict(_fe_sources(mesh, space), family=family, evaluation_tables=evaluation, material_budget=budget,
                charge_gram=arguments, fespace=space, fe_owner=(mesh, deformation))


def _require_budget_cover(assembly, stats):
    """Exact persistent plus setup-transient allocations must not exceed the native preflight bound."""
    used = assembly["persistent_device_bytes"] + assembly["setup_transient_device_bytes"]
    bound = stats.get("device_payload_bound_bytes")
    if type(bound) is not int or used > bound:
        raise RuntimeError("HDiv material allocations %r bytes exceed the preflight budget %r bytes" % (used, bound))


_CHARGE_BUDGET_TERMS = ("charge_prepare", "charge_map", "charge_map_scratch", "resident_vectors", "gram_dense",
                        "aca_batch_scratch")


def _checked_charge_budget(budget, compression):
    """Owned copy of the native pre-allocation charge budget; exactly one Gram storage term is nonzero."""
    terms = budget.get("terms") if isinstance(budget, dict) else None
    if (budget.get("schema") if isinstance(budget, dict) else None) != "radia.hdiv-flat-charge-budget.v1" \
            or not isinstance(terms, dict) or tuple(terms) != _CHARGE_BUDGET_TERMS:
        raise RuntimeError("Unexpected HDiv charge budget: %r" % (budget,))
    values = {name: terms[name] for name in _CHARGE_BUDGET_TERMS}
    storage = "gram_dense" if compression == "dense" else "aca_batch_scratch"
    if any(isinstance(v, bool) or not isinstance(v, int) or v < 0 for v in values.values()) \
            or any(values[name] <= 0 for name in _CHARGE_BUDGET_TERMS[:4] + (storage,)) \
            or values["aca_batch_scratch" if storage == "gram_dense" else "gram_dense"] != 0 \
            or budget.get("total_bytes") != sum(values.values()) \
            or budget.get("limiting_term") != max(values, key=values.get):
        raise RuntimeError("Inconsistent HDiv charge budget: %r" % (budget,))
    return dict(terms=values, total_bytes=budget["total_bytes"], limiting_term=budget["limiting_term"],
                excludes=budget.get("excludes"))


_MATERIAL_BUDGET_TERMS = ("reference_tables", "incidence", "incidence_sort", "mass_triples", "mass_csr",
                          "material_state")


def _checked_material_budget(budget):
    """Owned copy of the native material byte budget; the limiting term is the largest."""
    terms = budget.get("terms") if isinstance(budget, dict) else None
    if (budget.get("schema") if isinstance(budget, dict) else None) != "radia.hdiv-flat-material-budget.v1" \
            or not isinstance(terms, dict) or tuple(terms) != _MATERIAL_BUDGET_TERMS:
        raise RuntimeError("Unexpected HDiv material budget: %r" % (budget,))
    values = {name: terms[name] for name in _MATERIAL_BUDGET_TERMS}
    if any(isinstance(v, bool) or not isinstance(v, int) or v <= 0 for v in values.values()) \
            or budget.get("total_bytes") != sum(values.values()) \
            or budget.get("limiting_term") != max(values, key=values.get):
        raise RuntimeError("Inconsistent HDiv material budget: %r" % (budget,))
    return dict(terms=values, total_bytes=budget["total_bytes"], limiting_term=budget["limiting_term"],
                excludes=budget.get("excludes"))


def _fe_sources(mesh, fes):
    """NG-owned uniform-source and mean-field operators (vectors, not matrices)."""
    import numpy as np
    import ngsolve as ng
    nd = fes.ndof
    v = fes.TestFunction()
    rhs = np.empty((nd, 3), dtype=np.float64)
    for component in range(3):
        direction = tuple(float(i == component) for i in range(3))
        source = ng.LinearForm(fes)
        source += ng.InnerProduct(ng.CF(direction), v)*ng.dx
        source.Assemble()
        rhs[:, component] = source.vec.FV().NumPy()
    # Same declared NG default volume rule as the existing CPU mean field.
    volume = float(ng.Integrate(1, mesh))
    if not math.isfinite(volume) or volume <= 0:
        raise ValueError("Invalid mesh volume")
    return dict(dofs=nd, rhs_operator=rhs, mean_operator=rhs.T/volume)


def _applied_source(fespace, uniform, applied_field):
    """The weak source: ``uniform @ H`` of a constant H in A/m, or the ``vim.Solve`` form of a CoefficientFunction."""
    import numpy as np
    import ngsolve as ng
    if isinstance(applied_field, ng.CoefficientFunction):
        if applied_field.dim != 3:
            raise ValueError("applied_field must be a 3-component CoefficientFunction: dim %r" % (applied_field.dim,))
        source = ng.LinearForm(fespace)
        source += applied_field*fespace.TestFunction()*ng.dx
        source.Assemble()
        vector = np.array(source.vec.FV().NumPy(), dtype=np.float64)
    else:
        field = np.asarray(applied_field, dtype=np.float64)
        if field.shape != (3,) or not np.isfinite(field).all():
            raise ValueError("applied_field must contain three finite H components in A/m or be a CoefficientFunction")
        vector = uniform @ field
    if not np.isfinite(vector).all():
        raise ValueError("applied_field must be finite")
    return np.ascontiguousarray(vector)


def _reduced_numbering(plain, reduced, mesh):
    """Plain HDiv DOF -> ``Compress(Periodic(HDiv))`` DOF, read element by element from NGSolve; fails on any conflict."""
    import numpy as np
    import ngsolve as ng
    numbering = np.full(plain.ndof, -1, dtype=np.int32)
    for region in (ng.VOL, ng.BND):
        for element in mesh.Elements(region):
            source, target = plain.GetDofNrs(element), reduced.GetDofNrs(element)
            if len(source) != len(target) or min(target) < 0:
                raise NotImplementedError("Device tables carry no removed DOF: element %r maps %r to %r" % (element, source, target))
            for a, b in zip(source, target):
                if numbering[a] not in (-1, b):
                    raise RuntimeError("NGSolve reduced numbering is not a map: DOF %d to %d and %d" % (a, numbering[a], b))
                numbering[a] = b
    if numbering.min() < 0 or not np.array_equal(np.unique(numbering), np.arange(reduced.ndof)):
        raise RuntimeError("Reduced numbering does not cover the %d reduced DOFs" % (reduced.ndof,))
    return numbering


def _hex_lattice_signs(nodes):
    """Sign of det J at each Q2 cell centre (node ix + 3 iy + 9 iz)."""
    import numpy as np
    ex, ey, ez = nodes[:, 14] - nodes[:, 12], nodes[:, 16] - nodes[:, 10], nodes[:, 22] - nodes[:, 4]
    return np.sign(np.einsum("ij,ij->i", ex, np.cross(ey, ez)))


def _deformed_lattices(mesh, order, periodic, q2_deformation, arguments, evaluation):
    """Replace the staged straight lattices by those of the caller's Curve(2) + SetDeformation map.

    Geometry is a non-matrix input: B rows do not depend on it beyond the sign of det J, which must hold, and the
    device builds E, the mass and the Gram from these nodes. Returns the deformation for its owner to keep alive.
    """
    import numpy as np
    import ngsolve as ng
    from radia.vim._vim import _prepare_gpu_hex_charge_inputs
    if not callable(q2_deformation):
        raise ValueError("q2_deformation must be a callable applied to the loaded NGSolve mesh")
    deformation = q2_deformation(mesh)
    if mesh.deformation is None or mesh.GetCurveOrder() != 2:
        raise ValueError("q2_deformation must leave a Curve(2) mesh with a live deformation")
    moved = _prepare_gpu_hex_charge_inputs(ng.HDiv(mesh, order=order), backend="gpu", excluded_boundaries=periodic,
                                           live_deformation=True)
    cells = np.ascontiguousarray(moved["cell_nodes"], dtype=np.float64).reshape(-1)
    faces = np.ascontiguousarray(moved["face_nodes"], dtype=np.float64).reshape(-1)
    if cells.shape != arguments["cell_vertices"].shape or faces.shape != arguments["face_vertices"].shape:
        raise ValueError("Deformed and staged HEX lattices disagree in size")
    if not np.array_equal(_hex_lattice_signs(cells.reshape(-1, 27, 3)), _hex_lattice_signs(arguments["cell_vertices"].reshape(-1, 27, 3))):
        raise ValueError("q2_deformation changes the sign of a cell Jacobian")
    arguments.update(cell_vertices=cells, face_vertices=faces)
    evaluation["tables"][0]["vertices"] = cells.reshape(-1, 81)
    return deformation


def _device_evaluation_tables(path, fes, native, curve_order, family, periodic=()):
    """NGSolve reference-shape tables for device E construction; no evaluation matrix.

    ``tables`` holds one table set per element type, a single one on a pure mesh.
    """
    import numpy as np
    staging = native._MaterialStagingContext(str(path), fes.globalorder, curve_order, list(periodic))
    try:
        staged = list(staging.mixed_tables()) if family == _MIXED else [staging.flat_tables()]
    finally:
        staging.close()
    groups = []
    for tables in staged:
        if tables.get("schema") != "radia.hdiv-flat-evaluation-tables.v1" or tables["dofs"] != fes.ndof:
            raise ValueError("Material and charge FE numbering mismatch")
        if tables.get("family") != family:
            raise ValueError("Evaluation tables family %r, requested %r" % (tables.get("family"), family))
        groups.append(dict(rule=np.ascontiguousarray(tables["rule"], dtype=np.float64),
                           reference_shapes=np.ascontiguousarray(tables["reference_shapes"], dtype=np.float64),
                           vertices=np.ascontiguousarray(tables["vertices"], dtype=np.float64),
                           element_class=np.ascontiguousarray(tables["element_class"], dtype=np.int32),
                           element_dofs=np.ascontiguousarray(tables["element_dofs"], dtype=np.int32),
                           element_points=np.ascontiguousarray(tables["element_points"], dtype=np.int32)))
    if sum(group["element_points"].size for group in groups) != staged[0]["total_points"]:
        raise ValueError("Device evaluation tables do not cover every material point")
    return dict(tables=groups)


def _device_charge_tables(path, fes, native, curve_order, family, periodic=()):
    """NGSolve charge tables for device B construction; no charge map.

    Returns the binding tables and the staged cell and exterior-face
    geometry: affine vertices (host-side consistency check only), curved P2
    nodes or HEX Q2 lattices, which the curved and HEX device Grams consume.
    Curved tables add the corner numbers and the charge mass rules.
    """
    import numpy as np
    staging = native._MaterialStagingContext(str(path), fes.globalorder, curve_order, list(periodic))
    try:
        staged = staging.charge_tables()
    finally:
        staging.close()
    if staged.get("schema") != "radia.hdiv-flat-charge-tables.v1" or staged["dofs"] != fes.ndof:
        raise ValueError("Charge tables and charge FE numbering mismatch")
    if staged.get("family") != family:
        raise ValueError("Charge tables family %r, requested %r" % (staged.get("family"), family))
    keys = _CHARGE_TABLE_KEYS + (_CURVED_CHARGE_KEYS if family == "curved-tet" else ())
    if family == _MIXED:
        keys = _MIXED_CHARGE_KEYS + ("host_types", "charge_host", "charge_kind", "charge_exponents")
    tables = {name: np.ascontiguousarray(staged[name], dtype=np.float64
                                         if name.endswith(("tables", "rule")) else np.int32)
              for name in keys}
    geometry = tuple(np.ascontiguousarray(staged[name], dtype=np.float64).reshape(-1)
                     for name in ("cell_vertices", "face_vertices"))
    return tables, geometry


# Device-generated evaluation and mass ledger: reference tables in; validity flags,
# the mass sparsity count and its Jacobi bound out.

def _checked_material_assembly(assembly, tables, *, mass=False):
    """Reject any material ledger other than the exact device-built evaluation route.

    Covers material setup only; per-solve state uploads are outside this ledger.
    Persistent and setup-transient device bytes must be reported as positive integers.
    The mass-only context (``mass=True``) runs the same build without the
    inverse-BH law. Curved P2 tables (10 nodes per element) upload P2 nodes
    and HEX tables (27-node Q2 lattices) their lattices in place of vertices,
    and report their own routes. The element-type groups of a mixed mesh add
    up per payload, with pyramid vertices as their own payload.
    """
    groups = tables["tables"]
    names = {12: "element_vertices_m", 15: "element_pyramid_vertices_m", 30: "element_p2_nodes_m", 81: "element_q2_nodes_m"}
    family = {30: "curved-tet", 81: "hex", 12: "flat-tet"}.get(groups[0]["vertices"].shape[-1])
    if len(groups) > 1 or groups[0]["vertices"].shape[-1] == 15:
        family = _MIXED
    route = "device-%s-%s" % (family, "mass" if mass else "matrix-free")
    if assembly.get("route") != route:
        raise RuntimeError("Unexpected HDiv material assembly route: %r" % (assembly.get("route"),))
    expected = {("evaluation_invalid_flag", "H2D"): 4}
    for group in groups:
        for name, size in ((names[group["vertices"].shape[-1]], group["vertices"].nbytes),
                           ("reference_shapes", group["reference_shapes"].nbytes), ("material_rule", group["rule"].nbytes),
                           ("element_class", group["element_class"].nbytes), ("element_dofs", group["element_dofs"].nbytes),
                           ("element_points", group["element_points"].nbytes)):
            expected[(name, "H2D")] = expected.get((name, "H2D"), 0) + size
    expected.update({("evaluation_flag", "D2H"): 4, ("mass_nonzeros", "D2H"): 4, ("mass_flag", "D2H"): 4,
                     ("mass_spectral_bound", "D2H"): 8})
    if not mass:
        expected[("inverse_bh_law", "H2D")] = None
    for name in ("persistent_device_bytes", "setup_transient_device_bytes"):
        value = assembly.get(name)
        if type(value) is not int or value <= 0:
            raise RuntimeError("HDiv material %s rejected: %r" % (name, value))
    seen = {}
    for row in assembly["transfer_ledger"]:
        key, size = (row["payload"], row["direction"]), row["bytes"]
        if key in seen or key not in expected or type(size) is not int or size < 0:
            raise RuntimeError("HDiv material transfer ledger row rejected: %r" % (row,))
        want = expected[key]
        if (want is None and size <= 0) or (want is not None and size != want):
            raise RuntimeError("HDiv material ledger bytes for %r: actual %r, expected %r"
                               % (key, size, "positive" if want is None else want))
        seen[key] = size
    if set(seen) != set(expected):
        raise RuntimeError("HDiv material transfer ledger missing rows: %r" % (sorted(set(expected) - set(seen)),))
    return assembly


def prepare_hdiv_device_setup(mesh_path, *, backend, label_contract, mesh_check_receipt,
                              mesh_check_receipt_sha256, mesh_check_provenance, order=1, gram_aca=None,
                              gram_compression="aca", curve_order=1):
    """Return the owned primitive device-construction payload of a checked mesh.

    Schema ``radia.hdiv-device-setup.v1``: NGSolve reference tables and charge
    geometry for the device Gram, charge map, evaluation map and reference
    mass, plus the NG-assembled source and mean-field vectors. No matrix and
    no resident GPU object is created; the binding consumers (pybind11 and
    the plain MEX) build every operator on the device from these arrays.
    Order 0 is conforming RT0 (affine TET, linear consumers only).
    ``curve_order=2`` selects the isoparametric P2 family ``curved-tet``
    (orders 1 and 2); its payload adds the curved Gram and charge-map tables.
    A pure-HEX mesh is the family ``hex`` (orders 1 and 2, linear or
    ``Curve(2)``); its payload adds the 1D Gram rules and bands. The caller
    owns TaskManager.
    """
    if backend != "gpu":
        raise ValueError("Device setup requires explicit backend='gpu'")
    if isinstance(order, bool) or order not in (0, 1, 2):
        raise ValueError("HDiv device setup order must be 0 (RT0), 1 or 2")
    curve_order = _checked_curve_order(curve_order, order)
    aca = _checked_gram_compression(gram_compression, gram_aca)
    with _checked_mesh_snapshot(mesh_path, label_contract, mesh_check_receipt,
                                mesh_check_receipt_sha256, mesh_check_provenance) as (path, mesh_hash, contract_hash):
        from .gpu import _native
        native = _native("_hdiv_cuda")
        setup = _prepare_device_components(path, order, native, aca=aca, compression=gram_compression,
                                           curve_order=curve_order)
    return dict(schema="radia.hdiv-device-setup.v1", backend="gpu", family=setup["family"], order=order,
                dofs=setup["dofs"], mesh_sha256=mesh_hash, label_contract_sha256=contract_hash,
                evaluation_tables=setup["evaluation_tables"], charge_gram=setup["charge_gram"],
                material_budget=setup["material_budget"], rhs_operator=setup["rhs_operator"],
                mean_operator=setup["mean_operator"])
