"""Headless HDiv-MMM workflow over the shared CUDA material/Newton core.

NGSolve owns FE assembly; Radia owns charge-Gram preparation. The caller owns
the TaskManager region. Construction prepares operators, never a CPU solution.
The charge Gram, the charge map B/B^T, the material evaluation map and the
Newton reference mass are generated on the device from flat-TET geometry and
NGSolve reference tables; no host-built sparse operator is uploaded.
"""
from __future__ import annotations

__all__ = ["PreparedHDivNonlinear", "prepare_hdiv_host_setup"]

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
    """

    @classmethod
    def from_bh_table(cls, mesh_path, bh_table, *, backend, label_contract, mesh_check_receipt,
                      mesh_check_receipt_sha256, mesh_check_provenance, order=1):
        """Use the existing Radia CPU polynomial construction, without a CPU solve.

        Rows are [H in A/m, B in tesla]; one uniform soft-iron law is supported.
        """
        if backend != "gpu":
            raise ValueError("This workflow requires explicit backend='gpu'")
        if isinstance(order, bool) or order not in (1, 2):
            raise ValueError("HDiv order must be 1 or 2")
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
            prepared._initialize(snapshot, law, order)
            return prepared

    def __init__(self, mesh_path, law, *, backend, label_contract, mesh_check_receipt,
                 mesh_check_receipt_sha256, mesh_check_provenance, order=1):
        if backend != "gpu":
            raise ValueError("This workflow requires explicit backend='gpu'; CPU uses radia.vim")
        if isinstance(order, bool) or order not in (1, 2):
            raise ValueError("HDiv order must be 1 or 2")
        with _checked_mesh_snapshot(mesh_path, label_contract, mesh_check_receipt,
                                    mesh_check_receipt_sha256, mesh_check_provenance) as (snapshot, mesh_hash, contract_hash):
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
            setup = _prepare_device_components(path, order, _hdiv_cuda)
            self.dofs = setup["dofs"]
            self._rhs = setup["rhs_operator"]
            self._operator = _hdiv_cuda._ResidentHDiv.from_flat_charges(
                **setup["charge_gram"], dofs=self.dofs, backend="gpu")
            self.assembly = _checked_device_assembly(self._operator.assembly(), setup["charge_gram"])
            native_law = dict(law)
            for name in ("knots", "coefficients", "integral_base"):
                native_law[name] = np.ascontiguousarray(native_law[name], dtype=np.float64).ravel()
            self._material = _hdiv_cuda._ResidentBHMaterial.from_flat_tables(
                **setup["evaluation_tables"], dofs=self.dofs, law=native_law, backend="gpu")
            self.material_assembly = _checked_material_assembly(self._material.assembly(),
                                                                setup["evaluation_tables"])
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
            # The reference mass is the material context's device-built E^T W E.
            result = self._operator.solve_nonlinear_material(self._material,
                np.ascontiguousarray(self._rhs @ field),
                np.ascontiguousarray(x), tolerance, linear_tolerance, max_iterations,
                max_linear_iterations, preconditioner, degree, lower_fraction)
            return dict(result, backend="gpu", mesh_sha256=self.mesh_sha256,
                        order=self.order, dofs=self.dofs, field_units="A/m")

    def stats(self):
        with self._lock:
            return dict(closed=self._closed, backend="gpu", dofs=getattr(self, "dofs", None),
                gram_route="device-flat-tet", gram_compression=_GRAM_COMPRESSION,
                assembly=getattr(self, "assembly", None),
                material_assembly=getattr(self, "material_assembly", None),
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
_RULE_ROW_BYTES = 32  # FlatRulePoint: xi[3], weight as float64.
# Device ACA of the charge Gram with the CPU charge-Gram H-matrix parameters
# (radia.vim build_charge_gram eps, leaf size, eta). eps is the ACA stopping
# parameter, not an operator error bound. A leaf reaching max_rank unconverged
# is stored as an exact dense leaf of the same operator and counted.
_GRAM_COMPRESSION = "aca"
_GRAM_ACA = dict(eps=1e-10, eta=2.0, leaf_size=64, max_rank=200, rank_cap="dense")
# Metadata rows that exist only when the partition has admissible blocks.
_ACA_OPTIONAL_ROWS = (("aca_leaf_spans", "H2D"), ("aca_leaf_status", "D2H"), ("aca_active_count", "D2H"))
_ACA_COUNTERS = ("elements", "clusters", "dense_leaves", "lowrank_leaves", "admissible_blocks",
                 "rank_cap_promotions", "storage_promotions", "rank_sum", "max_rank", "aca_steps",
                 "aca_batches", "factor_coefficients", "dense_coefficients")


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
    """
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
                 if not name.endswith("exponents"))
    expected = {(name, "H2D"): sizes[name] for name in _DEVICE_GRAM_UPLOADS if sizes[name]}
    expected.update({("prepare_flag", "D2H"): 4, ("energy_flag", "D2H"): 4,
                     ("charge_map_flag", "D2H"): 4, ("resident_plan", "H2D"): None})
    if arguments["compression"] == "aca":
        expected.update({("aca_permutation", "H2D"): 8*arguments["charge_host"].size,
                         ("aca_dense_spans", "H2D"): None, ("dense_flag", "D2H"): 4})
    return expected


def _checked_device_assembly(assembly, arguments):
    """Reject any Gram assembly ledger other than the exact primitive-only device route.

    A host-built charge-map row is rejected like any other unexpected payload.
    """
    charges = int(arguments["charge_host"].size)
    if assembly.get("route") != "device-flat-tet" or assembly.get("charges") != charges:
        raise RuntimeError("Unexpected HDiv Gram assembly route: %r" % (assembly.get("route"),))
    if assembly.get("compression") != arguments["compression"]:
        raise RuntimeError("HDiv Gram compression %r, requested %r"
                           % (assembly.get("compression"), arguments["compression"]))
    counters = assembly.get("hmatrix")
    if arguments["compression"] == "aca" and (set(counters) != set(_ACA_COUNTERS) or not all(
            type(v) is int and v >= 0 for v in counters.values())):
        raise RuntimeError("HDiv ACA Gram counters rejected: %r" % (counters,))
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


def _checked_charge_exponents(kind, exponents):
    """Bound every exponent before any degree sum (cell <= 1, face <= 2, face z = 0)."""
    import numpy as np
    if exponents.shape != (3*kind.size,) or not np.isin(kind, (0, 1)).all():
        raise ValueError("Charge kinds/exponents disagree: kinds %r, exponents %r"
                         % (kind.shape, exponents.shape))
    expo = exponents.reshape(-1, 3).astype(np.int64)
    limit = np.where(kind == 0, 1, 2)[:, None]
    bad = (expo < 0) | (expo > limit)
    bad[:, 2] |= (kind == 1) & (expo[:, 2] != 0)
    bad = bad.any(axis=1) | (expo.sum(axis=1) > limit[:, 0])
    if bad.any():
        first = int(np.flatnonzero(bad)[0])
        raise ValueError("Device charge Gram supports cell degree <= 1 and face degree <= 2: "
                         "charge %d kind %d exponents %r" % (first, int(kind[first]), expo[first].tolist()))


def _device_gram_arguments(tables, inputs, geometry):
    """Map the staged flat-TET charge inputs and NGSolve tables to the device binding.

    The Gram is compressed by device ACA with the CPU charge-Gram H-matrix
    parameters (``_GRAM_ACA``); only the geometric partition is built on the
    host. The inner subtraction rule is unused by the exact analytic BDM1/BDM2
    potentials, as on the CPU; higher degrees fail loudly on the device.
    The table staging must reproduce the Gram geometry and charge order
    exactly (``geometry`` is its cell and exterior-face vertex coordinates).
    """
    import numpy as np
    if len(inputs.get("image_masks", ())):
        raise NotImplementedError("Device charge Gram has no image route")
    cells = np.ascontiguousarray(inputs["cell_verts"], dtype=np.float64).reshape(-1)
    faces = np.ascontiguousarray(inputs["face_verts"], dtype=np.float64).reshape(-1)
    if cells.size != 12*int(inputs["n_el"]):
        raise ValueError("Charge cell geometry does not match n_el")

    def rule(kind, suffix, width):
        points = np.asarray(inputs["ref_%s_pts%s" % (kind, suffix)], dtype=np.float64).reshape(-1, width)
        weights = np.asarray(inputs["ref_%s_w%s" % (kind, suffix)], dtype=np.float64).reshape(-1)
        if len(points) != len(weights):
            raise ValueError("Reference rule points and weights disagree")
        rows = np.zeros((len(points), 4))
        rows[:, :width], rows[:, 3] = points, weights
        return rows

    far = "ho_far_factor" in inputs and bool(np.isfinite(inputs["ho_far_factor"]))
    arguments = dict(
        cell_vertices=cells, face_vertices=faces,
        charge_host=np.ascontiguousarray(inputs["charge_host"], dtype=np.int32),
        charge_kind=np.ascontiguousarray(inputs["charge_kind"], dtype=np.int32),
        charge_exponents=np.ascontiguousarray(inputs["charge_expo"], dtype=np.int32),
        tet_rule=rule("tet", "", 3), tri_rule=rule("tri", "", 2),
        tet_low_rule=rule("tet", "_lo", 3) if far else np.zeros((0, 4)),
        tri_low_rule=rule("tri", "_lo", 2) if far else np.zeros((0, 4)),
        far_factor=float(inputs["ho_far_factor"]) if far else 1e30,
        compression=_GRAM_COMPRESSION, aca=dict(_GRAM_ACA))
    _checked_charge_exponents(arguments["charge_kind"], arguments["charge_exponents"])
    cells_m, faces_m = geometry
    if not (np.array_equal(cells_m, cells) and np.array_equal(faces_m, faces)):
        raise ValueError("Charge tables and Gram geometry disagree")
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


def _prepare_device_components(path, order, native, *, maximum_dofs=16384):
    """FE operators plus the primitive device-Gram payload; no CPU Gram or H-matrix.

    ``native`` is the caller's verified ``_hdiv_cuda`` module: one preparation
    verifies the native and CUDA library closure once, not per staging helper.
    """
    import ngsolve as ng
    from radia.vim._vim import _prepare_gpu_charge_inputs
    mesh = ng.Mesh(str(path))
    fes = ng.HDiv(mesh, order=order)
    if not 0 < fes.ndof <= maximum_dofs:
        raise ValueError("HDiv workflow dimension exceeds bounded scope")
    _, inputs = _prepare_gpu_charge_inputs(fes, backend="gpu")
    tables, geometry = _device_charge_tables(path, fes, native)
    return dict(_fe_sources(mesh, fes), evaluation_tables=_device_evaluation_tables(path, fes, native),
                charge_gram=_device_gram_arguments(tables, inputs, geometry))


def _prepare_host_components(path, order, *, maximum_dofs=16384):
    """NG-owned FE/charge preparation for the MATLAB publisher.

    No resident GPU operator and no CPU solution is created here. The caller
    owns TaskManager; the private native material staging shares the same FE.
    The publisher still ships host H-matrix leaves; its device route is pending.
    """
    import ngsolve as ng
    from radia.vim._vim import build_charge_gram
    mesh = ng.Mesh(str(path))
    fes = ng.HDiv(mesh, order=order)
    if not 0 < fes.ndof <= maximum_dofs:
        raise ValueError("HDiv workflow dimension exceeds bounded scope")
    B, gram, _ = build_charge_gram(fes, eps=1e-10, leafsize=64, eta=2.)
    snapshot = gram._resident_leaf_snapshot(16)
    if snapshot["dimension"] != B.shape[0]:
        raise ValueError("Charge snapshot numbering mismatch")
    return dict(_fe_components(path, mesh, fes), charge_map=B, snapshot=snapshot)


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


def _fe_forms(mesh, fes):
    """Sources plus the NG-assembled host mass (MATLAB publisher route only)."""
    import ngsolve as ng
    from scipy.sparse import coo_matrix
    nd = fes.ndof
    u, v = fes.TnT()
    mass = ng.BilinearForm(fes, symmetric=True)
    mass += ng.InnerProduct(u, v)*ng.dx
    mass.Assemble()
    rows, columns, values = mass.mat.COO()
    M = coo_matrix((values, (rows, columns)), shape=(nd, nd)).tocsr()
    M.sort_indices()
    return dict(_fe_sources(mesh, fes), mass=M)


def _fe_components(path, mesh, fes):
    """FE forms plus the NG-built evaluation CSR (MATLAB publisher route)."""
    import numpy as np
    from scipy.sparse import csc_matrix
    from .gpu import _native
    _hdiv_cuda = _native("_hdiv_cuda")
    forms = _fe_forms(mesh, fes)
    nd = forms["dofs"]
    staging = _hdiv_cuda._MaterialStagingContext(str(path), fes.globalorder)
    try:
        payload = staging.quadrature()
    finally:
        staging.close()
    if payload["dofs"] != nd:
        raise ValueError("Material and charge FE numbering mismatch")
    E = csc_matrix((payload["evaluation_values"], payload["evaluation_row_indices"],
        payload["evaluation_column_offsets"]), shape=(payload["evaluation_rows"], nd)).tocsr()
    E.sort_indices()
    return dict(forms, evaluation=E, weights=np.asarray(payload["weights"]),
                quadrature_order=payload["integration_order"])


def _device_evaluation_tables(path, fes, native):
    """NGSolve reference-shape tables for device E construction; no evaluation matrix."""
    import numpy as np
    staging = native._MaterialStagingContext(str(path), fes.globalorder)
    try:
        tables = staging.flat_tables()
    finally:
        staging.close()
    if tables.get("schema") != "radia.hdiv-flat-evaluation-tables.v1" or tables["dofs"] != fes.ndof:
        raise ValueError("Material and charge FE numbering mismatch")
    arguments = dict(rule=np.ascontiguousarray(tables["rule"], dtype=np.float64),
                     reference_shapes=np.ascontiguousarray(tables["reference_shapes"], dtype=np.float64),
                     vertices=np.ascontiguousarray(tables["vertices"], dtype=np.float64),
                     element_class=np.ascontiguousarray(tables["element_class"], dtype=np.int32),
                     element_dofs=np.ascontiguousarray(tables["element_dofs"], dtype=np.int32),
                     element_points=np.ascontiguousarray(tables["element_points"], dtype=np.int32))
    if arguments["element_points"].size != tables["total_points"]:
        raise ValueError("Device evaluation tables do not cover every material point")
    return arguments


def _device_charge_tables(path, fes, native):
    """NGSolve charge-projection tables for device B construction; no charge map.

    Returns the binding tables and, for the host-side consistency check only,
    the staged cell and exterior-face vertex coordinates.
    """
    import numpy as np
    staging = native._MaterialStagingContext(str(path), fes.globalorder)
    try:
        staged = staging.charge_tables()
    finally:
        staging.close()
    if staged.get("schema") != "radia.hdiv-flat-charge-tables.v1" or staged["dofs"] != fes.ndof:
        raise ValueError("Charge tables and charge FE numbering mismatch")
    tables = {name: np.ascontiguousarray(staged[name], dtype=np.float64 if name.endswith("tables") else np.int32)
              for name in _CHARGE_TABLE_KEYS}
    geometry = tuple(np.ascontiguousarray(staged[name], dtype=np.float64).reshape(-1)
                     for name in ("cell_vertices", "face_vertices"))
    return tables, geometry


# Device-generated evaluation and mass ledger: reference tables in; validity flags,
# the mass sparsity count and its Jacobi bound out.
_DEVICE_EVALUATION_UPLOADS = ("element_vertices_m", "reference_shapes", "material_rule", "element_class",
                              "element_dofs", "element_points", "evaluation_invalid_flag")


def _checked_material_assembly(assembly, tables):
    """Reject any material ledger other than the exact device-built evaluation route.

    Covers material setup only; per-solve state uploads are outside this ledger.
    Persistent and setup-transient device bytes must be reported as positive integers.
    """
    if assembly.get("route") != "device-flat-tet-evaluation":
        raise RuntimeError("Unexpected HDiv material assembly route: %r" % (assembly.get("route"),))
    sizes = {"element_vertices_m": tables["vertices"].nbytes, "reference_shapes": tables["reference_shapes"].nbytes,
             "material_rule": tables["rule"].nbytes, "element_class": tables["element_class"].nbytes,
             "element_dofs": tables["element_dofs"].nbytes, "element_points": tables["element_points"].nbytes,
             "evaluation_invalid_flag": 4}
    expected = {(name, "H2D"): sizes[name] for name in _DEVICE_EVALUATION_UPLOADS}
    expected.update({("evaluation_flag", "D2H"): 4, ("inverse_bh_law", "H2D"): None,
                     ("mass_nonzeros", "D2H"): 4, ("mass_flag", "D2H"): 4, ("mass_spectral_bound", "D2H"): 8})
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


def prepare_hdiv_host_setup(mesh_path, *, backend, label_contract, mesh_check_receipt,
                            mesh_check_receipt_sha256, mesh_check_provenance, order=1, maximum_dofs=16384):
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
    with _checked_mesh_snapshot(mesh_path, label_contract, mesh_check_receipt,
                                mesh_check_receipt_sha256, mesh_check_provenance) as (path, mesh_hash, contract_hash):
        data = _prepare_host_components(path, order, maximum_dofs=maximum_dofs)
    return dict(data, mesh_sha256=mesh_hash, label_contract_sha256=contract_hash, order=order)
