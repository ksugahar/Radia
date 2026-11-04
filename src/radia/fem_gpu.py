"""FEM solves with every operator built on the selected backend.

``Backend="gpu"`` builds stiffness and mass matrices, the sparsity pattern, the
AMG, AMS and BDDC preconditioners and their coarse operators on the device from
mesh arrays and reference tabulations. ``Backend="cpu"`` runs the NGSolve and
compact-solver route on the same inputs. There is no fallback between them.

What stays on the host for a device solve, and crosses to the device as data:
DOF numbers, vertex ordering classes, shape tabulations at quadrature points,
integration rules, vertex coordinates, element coefficients, curved-element
Jacobians, the right-hand side and the lift. NGSolve supplies these; it builds no
operator of the device route. ``ledger`` in each result lists every crossing, and a
crossing outside the permitted payload families raises.
"""
from __future__ import annotations

from dataclasses import dataclass
import time

import numpy as np

from ._residual_gate import RELATIVE_LIMIT

PERMITTED_PAYLOADS = ("mesh.", "table.", "coef.", "rhs.", "lift.", "result.", "scalar.", "probe.")
AMS_MIN_COARSE = 500
AMG_MIN_COARSE = 50
AMG_THETA = 0.25


@dataclass(frozen=True)
class Coefficients:
    """Per-material and per-boundary constants of one sesquilinear-free form.

    Each field is a float or a mapping from a region name to a float. ``stiff`` and
    ``mass`` act on the volume, ``trace`` on the boundary. The ``_im`` fields are
    the imaginary parts of the system; the ``_aux`` fields define the real operator
    from which the AMG or AMS hierarchy is built.
    """

    stiff: object = 0.0
    stiff_im: object = 0.0
    stiff_aux: object = 0.0
    mass: object = 0.0
    mass_im: object = 0.0
    mass_aux: object = 0.0
    trace: object = 0.0
    trace_im: object = 0.0
    trace_aux: object = 0.0


def _native():
    from .gpu import _native as load
    return load("_fem_cuda")


def _check_backend(Backend):
    if Backend not in ("gpu", "cpu"):
        raise ValueError(f"Backend must be 'gpu' or 'cpu', got {Backend!r}")


def _region_value(value, region):
    if isinstance(value, dict):
        return float(value.get(region, 0.0))
    return float(value)


def _check_regions(value, regions, label):
    if isinstance(value, dict):
        unknown = sorted(set(value) - set(regions))
        if unknown:
            raise ValueError(f"{label} names regions the mesh does not have: {unknown}")


def _space_kind(fes):
    name = type(fes).__name__
    if name not in ("H1", "HCurl"):
        raise ValueError(f"the FEM device route supports H1 and HCurl spaces, got {name}")
    return name


def _element_rows(fes, mesh, coefficients):
    """Group the elements by (VOL or BND, element type); no matrix is touched."""
    import ngsolve as ng

    materials, boundaries = mesh.GetMaterials(), mesh.GetBoundaries()
    for field, regions in (("stiff", materials), ("stiff_im", materials), ("stiff_aux", materials), ("mass", materials),
                           ("mass_im", materials), ("mass_aux", materials), ("trace", boundaries),
                           ("trace_im", boundaries), ("trace_aux", boundaries)):
        _check_regions(getattr(coefficients, field), regions, field)
    volume = {m: (_region_value(coefficients.stiff, m), _region_value(coefficients.stiff_im, m),
                  _region_value(coefficients.stiff_aux, m), _region_value(coefficients.mass, m),
                  _region_value(coefficients.mass_im, m), _region_value(coefficients.mass_aux, m)) for m in set(materials)}
    surface = {b: (0.0, 0.0, 0.0, _region_value(coefficients.trace, b), _region_value(coefficients.trace_im, b),
                   _region_value(coefficients.trace_aux, b)) for b in set(boundaries)}
    groups = {}
    for vb, table in ((ng.VOL, volume), (ng.BND, surface)):
        for element in fes.Elements(vb):
            columns = table[element.mat]
            if vb == ng.BND and not any(columns):
                continue
            row = groups.setdefault((vb, element.type), {"numbers": [], "dofs": [], "vertices": [], "coef": []})
            row["numbers"].append(element.nr)
            row["dofs"].append(element.dofs)
            row["vertices"].append([v.nr for v in element.vertices])
            row["coef"].append(columns)
    return groups


def _derivative_rule_order(fe_order, et, bonus):
    """Rule order of a derivative form: 2 * element order, less two on simplices, plus the bonus."""
    import ngsolve as ng

    simplex = et in (ng.ET.TET, ng.ET.TRIG)
    return max(2 * fe_order - (2 if simplex else 0) + bonus, 0)


def _tabulate(fes, mesh, vb, representatives, points, nl, quantity):
    import ngsolve as ng

    table = np.zeros((len(representatives), len(points), nl, 3))
    for c, number in enumerate(representatives):
        ei = ng.ElementId(vb, int(number))
        fe, trafo = fes.GetFE(ei), mesh.GetTrafo(ei)
        for q, point in enumerate(points):
            x, y, z = (tuple(point) + (0.0, 0.0, 0.0))[:3]
            if quantity == "grad":
                values = np.array(fe.CalcDShape(x, y, z))
                table[c, q, :, :values.shape[1]] = values
            elif quantity == "scalar":
                table[c, q, :, 0] = np.array(fe.CalcShape(x, y, z)).ravel()
            elif quantity == "covariant":
                values = np.array(fe.CalcShape(x, y, z))
                table[c, q, :, :values.shape[1]] = values
            elif quantity == "curl":
                # The binding returns the physical curl J curl_ref / det J.
                mip = trafo(x, y, z)
                J = np.array(mip.jacobi).reshape(3, 3)
                curl = np.array(fe.CalcCurlShape(mip))
                table[c, q] = (np.linalg.det(J) * np.linalg.inv(J) @ curl.T).T
            else:
                raise ValueError(f"unknown tabulation {quantity}")
    return table


def _geometry_derivatives(mesh, vb, representative, points):
    """First-order geometry derivatives at the rule points, from the mesh order-1 space."""
    import ngsolve as ng

    fe = ng.H1(mesh, order=1).GetFE(ng.ElementId(vb, int(representative)))
    tables = []
    for point in points:
        x, y, z = (tuple(point) + (0.0, 0.0, 0.0))[:3]
        tables.append(np.array(fe.CalcDShape(x, y, z), dtype=np.float64))
    return np.array(tables, dtype=np.float64)


def _curved_jacobians(mesh, vb, rule_by_type, numbers):
    """Element Jacobians of a curved mesh, taken from the NGSolve mapping."""
    import ngsolve as ng

    all_types = [el.type for el in mesh.Elements(vb)]
    rules = {et: rule_by_type.get(et) or ng.IntegrationRule(et, 0) for et in set(all_types)}
    points = mesh.MapToAllElements(rules, vb)
    dim = 3 if vb == ng.VOL else 2
    flat = np.array(ng.specialcf.JacobianMatrix(3, dim)(points), dtype=np.float64)
    counts = np.array([len(rules[et].points) for et in all_types], dtype=np.int64)
    starts = np.concatenate(([0], np.cumsum(counts)))
    nq = int(counts[numbers[0]])
    picks = np.concatenate([np.arange(starts[n], starts[n] + nq) for n in numbers])
    return flat[picks].reshape(len(numbers) * nq * 3 * dim)


def _form_input(fes, mesh, vb, et, rule, representatives, numbers, nl, quantity, law, curved, tensor=False):
    import ngsolve as ng

    points = [tuple(p) for p in rule.points]
    weights = np.array(list(rule.weights), dtype=np.float64)
    affine = not curved and et in (ng.ET.TET, ng.ET.TRIG) and not tensor
    table = _tabulate(fes, mesh, vb, representatives, points, nl, quantity)
    form = {"law": law, "nq": len(points), "affine": affine, "weights": weights,
            "table": np.ascontiguousarray(table, dtype=np.float64), "dgeo": None, "jacobian": None, "tensor": None}
    if curved:
        form["jacobian"] = _curved_jacobians(mesh, vb, {et: rule}, numbers)
    else:
        form["dgeo"] = np.ascontiguousarray(_geometry_derivatives(mesh, vb, representatives[0], points[:1] if affine else points))
    return form


def _group_points(mesh, rule, et, numbers):
    """Mapped rule points of one element group, in element then rule order."""
    import ngsolve as ng

    all_types = [el.type for el in mesh.Elements(ng.VOL)]
    rules = {kind: rule if kind == et else ng.IntegrationRule(kind, 0) for kind in set(all_types)}
    points = mesh.MapToAllElements(rules, ng.VOL)
    counts = np.array([len(rules[kind].points) for kind in all_types], dtype=np.int64)
    starts = np.concatenate(([0], np.cumsum(counts)))
    nq = len(rule.points)
    return points[np.concatenate([np.arange(starts[n], starts[n] + nq) for n in numbers])]


def _point_tensor(tangent, points):
    """The caller's symmetric tensor at mapped points: six values per point."""
    values = np.asarray(tangent(points), dtype=np.float64)
    expected = (len(points), 6)
    if values.shape != expected or not np.isfinite(values).all():
        raise ValueError(f"the point tensor must be finite with shape {expected}, got {values.shape}")
    return np.ascontiguousarray(values.reshape(-1))


def _group_inputs(fes, mesh, coefficients, kind, bonus, native, tangent=None):
    import ngsolve as ng

    if mesh.dim != 3:
        raise ValueError(f"the FEM device route needs a three-dimensional mesh, got dimension {mesh.dim}")
    curved = mesh.GetCurveOrder() > 1
    vertices = np.array([v.point for v in mesh.vertices], dtype=np.float64)
    inputs = []
    rows = _element_rows(fes, mesh, coefficients)
    for (vb, et), row in rows.items():
        dofs = np.array(row["dofs"], dtype=np.int32)
        if dofs.min() < 0:
            raise ValueError("an element holds an inactive DOF; the device route does not support it")
        numbers = np.array(row["numbers"], dtype=np.int64)
        topology = np.array(row["vertices"], dtype=np.int32)
        ordering = np.argsort(topology, axis=1)
        unique, first, inverse = np.unique(ordering, axis=0, return_index=True, return_inverse=True)
        representatives = numbers[first]
        nl = dofs.shape[1]
        dim = 3 if vb == ng.VOL else 2
        # The rule follows the element's own order: HCurl on a hexahedron reports p + 1.
        fe_order = fes.GetFE(ng.ElementId(vb, int(representatives[0]))).order
        group = {"dim": dim, "ne": len(numbers), "nl": nl, "ncls": len(unique), "ngeo": topology.shape[1],
                 "dofs": np.ascontiguousarray(dofs, dtype=np.int32),
                 "cls": np.ascontiguousarray(inverse.ravel(), dtype=np.int32),
                 "nodes": None if curved else np.ascontiguousarray(vertices[topology], dtype=np.float64),
                 "coef": np.ascontiguousarray(np.array(row["coef"], dtype=np.float64)), "stiff": None, "mass": None}
        if vb == ng.VOL:
            quantity, law = ("grad", native.LAW_GRAD) if kind == "H1" else ("curl", native.LAW_CURL)
            rule = ng.IntegrationRule(et, _derivative_rule_order(fe_order, et, bonus[0]))
            group["stiff"] = _form_input(fes, mesh, vb, et, rule, representatives, numbers, nl, quantity, law, curved, tensor=tangent is not None)
            if tangent is not None:
                if kind != "H1":
                    raise ValueError("a point tensor applies to a gradient form on H1")
                group["_points"] = _group_points(mesh, rule, et, numbers)
                group["stiff"]["tensor"] = _point_tensor(tangent, group["_points"])
        trace_bonus = bonus[2] if vb == ng.BND else bonus[1]
        mass_rule = ng.IntegrationRule(et, 2 * fe_order + trace_bonus)
        quantity, law = ("scalar", native.LAW_VALUE) if kind == "H1" else ("covariant", native.LAW_COVARIANT)
        group["mass"] = _form_input(fes, mesh, vb, et, mass_rule, representatives, numbers, nl, quantity, law, curved)
        inputs.append(group)
    return inputs, vertices


def _edge_arrays(fes, mesh):
    """Mesh edges with their low-order DOF numbering checked against the space."""
    import ngsolve as ng

    edges = np.array([[v.nr for v in e.vertices] for e in mesh.edges], dtype=np.int32)
    for number in range(mesh.nedge):
        dofs = fes.GetDofNrs(ng.NodeId(ng.EDGE, number))
        if not dofs or dofs[0] != number:
            raise ValueError(f"the lowest-order DOF of edge {number} is {dofs}; AMS needs DOF number equal to edge number")
    return np.ascontiguousarray(edges)


def _bddc_split(fes, mesh, kind, groups):
    """Wirebasket numbering and the uniform element split of every group."""
    import ngsolve as ng

    if kind == "HCurl":
        is_wirebasket = np.arange(fes.ndof) < mesh.nedge
    else:
        is_wirebasket = np.array([fes.CouplingType(d) == ng.COUPLING_TYPE.WIREBASKET_DOF for d in range(fes.ndof)])
    wb_of_dof = np.full(fes.ndof, -1, dtype=np.int32)
    wb_of_dof[is_wirebasket] = np.arange(int(is_wirebasket.sum()), dtype=np.int32)
    split = []
    for group in groups:
        local = is_wirebasket[group["dofs"].reshape(group["ne"], group["nl"])]
        if not (local == local[0]).all():
            raise ValueError("the wirebasket split differs between elements of one group")
        split.append({"nw": int(local[0].sum()), "ni": int((~local[0]).sum()),
                      "wb_local": np.flatnonzero(local[0]).astype(np.int32),
                      "int_local": np.flatnonzero(~local[0]).astype(np.int32)})
    return split, wb_of_dof, int(is_wirebasket.sum())


def _complex_parts(value, ndof, name, complex_values):
    if value is None:
        return None, None
    array = np.asarray(value)
    if array.shape != (ndof,) or not np.isfinite(array).all():
        raise ValueError(f"{name} must be a finite vector of length {ndof}, got shape {array.shape}")
    if not complex_values:
        if np.iscomplexobj(array) and np.any(array.imag != 0):
            raise ValueError(f"{name} is complex but the system is real")
        return np.ascontiguousarray(array.real, dtype=np.float64), None
    return np.ascontiguousarray(array.real, dtype=np.float64), np.ascontiguousarray(array.imag, dtype=np.float64)


class DeviceSystem:
    """Operators of one space on the device. Mesh arrays cross; matrices never do."""

    def __init__(self, fes, coefficients, *, bonus=(0, 0, 0), complex_values=None, tangent=None):
        """``tangent(points)`` returns the symmetric tensor xx xy xz yy yz zz at mapped integration
        points, evaluated on the host; the stiffness column then scales that tensor."""
        native = _native()
        self.tangent = tangent
        self.fes, self.mesh = fes, fes.mesh
        self.kind = _space_kind(fes)
        self.order = int(fes.globalorder)
        self.complex_values = bool(fes.is_complex) if complex_values is None else bool(complex_values)
        self._times = {}
        start = time.perf_counter()
        self.groups, self.vertices = _group_inputs(fes, self.mesh, coefficients, self.kind, bonus, native, tangent)
        self.free_dofs = np.array(list(fes.FreeDofs()), dtype=np.int32)
        self._times["host_setup_seconds"] = time.perf_counter() - start
        self.native = native.System(int(fes.ndof), self.complex_values, self.groups, self.free_dofs)
        self.hierarchy = None
        self.preconditioner = None

    def update_tangent(self):
        """Evaluate the point tensors again and rebuild the operators on the device.

        The preconditioner is dropped; call ``setup`` before the next solve.
        """
        if self.tangent is None:
            raise ValueError("this system was built without a point tensor")
        tensors = [_point_tensor(self.tangent, group["_points"]) if group["stiff"] is not None else None for group in self.groups]
        self.native.update_tensors(tensors)
        self.preconditioner = None
        self.hierarchy = None

    def setup(self, preconditioner, *, cycles=None, min_coarse=None, theta=AMG_THETA):
        if preconditioner == "amg":
            self.hierarchy = [self.native.setup_amg(AMG_MIN_COARSE if min_coarse is None else min_coarse, theta)]
        elif preconditioner == "ams":
            if self.kind != "HCurl" or self.fes.ndof != self.mesh.nedge:
                raise ValueError("AMS needs a lowest-order HCurl space with one DOF per mesh edge, "
                                 f"got {self.fes.ndof} DOFs and {self.mesh.nedge} edges")
            edges = _edge_arrays(self.fes, self.mesh)
            self.hierarchy = self.native.setup_ams(int(self.mesh.nedge), int(self.mesh.nv), edges,
                                                   np.ascontiguousarray(self.vertices), AMS_MIN_COARSE if min_coarse is None else min_coarse, theta)
        elif preconditioner == "bddc":
            split, wb_of_dof, nw = _bddc_split(self.fes, self.mesh, self.kind, self.groups)
            if nw >= self.fes.ndof:
                raise ValueError("BDDC applies to p>1 or to a space with DOFs outside the wirebasket; use AMG or AMS on a pure "
                                 f"lowest-order space, got {self.fes.ndof} DOFs and {nw} wirebasket DOFs")
            if self.kind == "HCurl":
                edges = _edge_arrays(self.fes, self.mesh)
                self.hierarchy = self.native.setup_bddc(split, wb_of_dof, nw, int(self.mesh.nedge), int(self.mesh.nv), edges,
                                                        np.ascontiguousarray(self.vertices), 3 if cycles is None else int(cycles),
                                                        AMS_MIN_COARSE if min_coarse is None else min_coarse, theta)
            else:
                self.hierarchy = self.native.setup_bddc(split, wb_of_dof, nw, 0, 0, None, None, 1 if cycles is None else int(cycles),
                                                        AMG_MIN_COARSE if min_coarse is None else min_coarse, theta)
        else:
            raise ValueError(f"unknown device preconditioner {preconditioner!r}")
        self.preconditioner = preconditioner
        return self.hierarchy

    def action(self, which, x):
        index = {"re": 0, "im": 1, "aux": 2}[which]
        return self.native.probe_action(index, np.ascontiguousarray(x, dtype=np.float64))

    def apply_preconditioner(self, x):
        parts = _complex_parts(x, self.fes.ndof, "probe", self.complex_values)
        y_re, y_im = self.native.probe_preconditioner(parts[0], parts[1])
        return y_re + 1j * y_im if self.complex_values else y_re

    def ledger(self):
        rows = self.native.ledger()
        forbidden = [r for r in rows if not r["payload"].startswith(PERMITTED_PAYLOADS)]
        if forbidden:
            raise RuntimeError(f"a payload outside the permitted families crossed: {forbidden}")
        return rows

    def solve(self, rhs, *, lift=None, tolerance=1e-10, max_iterations=500):
        ndof = int(self.fes.ndof)
        rhs_re, rhs_im = _complex_parts(rhs, ndof, "rhs", self.complex_values)
        lift_re, lift_im = _complex_parts(lift, ndof, "lift", self.complex_values)
        if rhs_re is None:
            raise ValueError("a right-hand side is required")
        result = self.native.solve(rhs_re, rhs_im, lift_re, lift_im, float(tolerance), int(max_iterations), RELATIVE_LIMIT)
        x_re, x_im = result.pop("x_re"), result.pop("x_im")
        result["solution"] = x_re + 1j * x_im if self.complex_values else x_re
        result["backend"] = "gpu"
        result["device"] = _native().device_name()
        result["preconditioner"] = self.preconditioner
        result["hierarchy"] = self.hierarchy
        result["ledger"] = self.ledger()
        result["nonzeros"] = self.native.nonzeros()
        result["device_peak_bytes"] = self.native.peak_device_bytes()
        result["timing_diagnostic"] = {**self._times, **self.native.timings()}
        return result


def _material_cf(mesh, value):
    import ngsolve as ng

    if isinstance(value, dict):
        return mesh.MaterialCF(value, default=0.0)
    return ng.CF(float(value))


def _boundary_cf(mesh, value):
    import ngsolve as ng

    if isinstance(value, dict):
        return mesh.BoundaryCF(value, default=0.0)
    return ng.CF(float(value))


def _is_zero(value):
    return not value if not isinstance(value, dict) else not any(value.values())


def cpu_form(fes, coefficients, *, bonus=(0, 0, 0), part="system"):
    """The same form on the CPU; ``part`` is ``"system"`` or ``"aux"``. Not assembled."""
    import ngsolve as ng

    mesh = fes.mesh
    kind = _space_kind(fes)
    derivative = ng.grad if kind == "H1" else ng.curl
    target_fes = fes
    u, v = target_fes.TnT()
    form = ng.BilinearForm(target_fes)
    if part == "system":
        terms = ((coefficients.stiff, coefficients.stiff_im), (coefficients.mass, coefficients.mass_im),
                 (coefficients.trace, coefficients.trace_im))
    elif part == "aux":
        terms = ((coefficients.stiff_aux, 0.0), (coefficients.mass_aux, 0.0), (coefficients.trace_aux, 0.0))
    else:
        raise ValueError(f"unknown form part {part!r}")
    stiff, mass, trace = terms

    def scalar(pair, region_cf):
        re = region_cf(mesh, pair[0])
        return re if _is_zero(pair[1]) else re + 1j * region_cf(mesh, pair[1])

    if not (_is_zero(stiff[0]) and _is_zero(stiff[1])):
        form += scalar(stiff, _material_cf) * ng.InnerProduct(derivative(u), derivative(v)) * ng.dx(bonus_intorder=bonus[0])
    if not (_is_zero(mass[0]) and _is_zero(mass[1])):
        form += scalar(mass, _material_cf) * ng.InnerProduct(u, v) * ng.dx(bonus_intorder=bonus[1])
    if not (_is_zero(trace[0]) and _is_zero(trace[1])):
        form += scalar(trace, _boundary_cf) * ng.InnerProduct(u.Trace(), v.Trace()) * ng.ds(bonus_intorder=bonus[2])
    return form


def _cpu_cg(matrix, pre, b, free, tolerance, max_iterations):
    x = b.CreateVector()
    r, z, p, q = b.CreateVector(), b.CreateVector(), b.CreateVector(), b.CreateVector()
    x[:] = 0
    r.data = b
    b_norm = float(np.linalg.norm(b.FV().NumPy()[free]))
    if not b_norm > 0:
        raise RuntimeError("CPU CG: the right-hand side is zero on the free DOFs")
    z.data = pre * r
    p.data = z
    rho = float(np.dot(r.FV().NumPy()[free], z.FV().NumPy()[free]))
    history = []
    for iteration in range(1, max_iterations + 1):
        q.data = matrix * p
        curvature = float(np.dot(p.FV().NumPy()[free], q.FV().NumPy()[free]))
        if not curvature > 0 or not rho > 0:
            raise RuntimeError(f"CPU CG breakdown at iteration {iteration}: rho {rho} curvature {curvature}")
        alpha = rho / curvature
        x.data += alpha * p
        r.data -= alpha * q
        norm = float(np.linalg.norm(r.FV().NumPy()[free]))
        history.append(norm / b_norm)
        if norm <= tolerance * b_norm:
            return x, iteration, history
        z.data = pre * r
        rho_next = float(np.dot(r.FV().NumPy()[free], z.FV().NumPy()[free]))
        p.data = z + (rho_next / rho) * p
        rho = rho_next
    raise RuntimeError(f"CPU CG did not converge in {max_iterations} iterations: relative residual {history[-1]}")


def cpu_solve(fes, coefficients, rhs, *, lift=None, preconditioner="direct", bonus=(0, 0, 0), tolerance=1e-10, max_iterations=500,
              cycles=None, min_coarse=None, theta=AMG_THETA):
    """Solve the same system with the CPU route. Direct solves use sparsecholesky only."""
    import ngsolve as ng

    start = time.perf_counter()
    kind = _space_kind(fes)
    mesh = fes.mesh
    complex_values = bool(fes.is_complex)
    ndof = int(fes.ndof)
    free = np.array(list(fes.FreeDofs()), dtype=bool)
    form = cpu_form(fes, coefficients, bonus=bonus)
    if preconditioner == "bddc":
        if kind == "HCurl":
            for dof in range(mesh.nedge, fes.ndof):
                if fes.CouplingType(dof) == ng.COUPLING_TYPE.WIREBASKET_DOF:
                    fes.SetCouplingType(dof, ng.COUPLING_TYPE.INTERFACE_DOF)
            import radia.sparsesolv_ngsolve  # noqa: F401  (registers sparsesolv_ams)
            pre = ng.Preconditioner(form, "bddc", coarsetype="sparsesolv_ams", coarseflags={"cycles": 3 if cycles is None else int(cycles)})
        else:
            pre = ng.Preconditioner(form, "bddc", coarsetype="h1amg")
    form.Assemble()
    matrix = form.mat
    f = matrix.CreateRowVector()
    g = matrix.CreateRowVector()
    f.FV().NumPy()[:] = np.asarray(rhs)
    g.FV().NumPy()[:] = 0 if lift is None else np.asarray(lift)
    residual = f.CreateVector()
    residual.data = f - matrix * g
    residual.FV().NumPy()[~free] = 0
    if preconditioner == "direct":
        inverse = matrix.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
        x = residual.CreateVector()
        x.data = inverse * residual
        iterations = 0
    else:
        import radia.sparsesolv_ngsolve as ssn

        if preconditioner == "amg":
            pre = ssn.CompactAMGPreconditioner(mat=_aux_matrix(fes, coefficients, bonus), freedofs=fes.FreeDofs(), theta=theta, max_levels=25,
                                               min_coarse=AMG_MIN_COARSE if min_coarse is None else min_coarse, num_smooth=1, print_level=0)
        elif preconditioner == "ams":
            pre = _cpu_ams(fes, coefficients, bonus, ssn, ndof)
        elif preconditioner != "bddc":
            raise ValueError(f"unknown CPU preconditioner {preconditioner!r}")
        if complex_values:
            solver = ssn.COCRSolver(matrix, pre, freedofs=fes.FreeDofs(), maxiter=max_iterations, tol=tolerance, printrates=False)
            x = residual.CreateVector()
            x.data = solver * residual
            iterations = int(solver.iterations)
        else:
            x, iterations, _history = _cpu_cg(matrix, pre, residual, free, tolerance, max_iterations)
    check = residual.CreateVector()
    check.data = matrix * x
    difference = check.FV().NumPy()[free] - residual.FV().NumPy()[free]
    true_residual = float(np.linalg.norm(difference) / np.linalg.norm(residual.FV().NumPy()[free]))
    if not true_residual <= RELATIVE_LIMIT:
        raise RuntimeError(f"CPU route: true relative residual {true_residual} exceeds the limit {RELATIVE_LIMIT}")
    if preconditioner != "direct" and not iterations:
        raise RuntimeError("CPU route: the iterative solver reported no iteration")
    solution = np.array(x.FV().NumPy(), copy=True)
    solution[~free] = 0
    solution = solution + np.array(g.FV().NumPy(), copy=True)
    return {"solution": solution, "backend": "cpu", "preconditioner": preconditioner, "iterations": iterations, "true_residual": true_residual,
            "timing_diagnostic": {"total_seconds": time.perf_counter() - start}}


def _aux_matrix(fes, coefficients, bonus):
    form = cpu_form(fes, coefficients, bonus=bonus, part="aux")
    form.Assemble()
    return form.mat


def _cpu_ams(fes, coefficients, bonus, ssn, ndof):
    import ngsolve as ng

    mesh = fes.mesh
    real_space = ng.HCurl(mesh, order=1, nograds=True, complex=False)
    aux = _aux_matrix(real_space, coefficients, bonus)
    gradient, _ = real_space.CreateGradient()
    coordinates = np.array([v.point for v in mesh.vertices], dtype=np.float64)
    arguments = dict(grad_mat=gradient, freedofs=fes.FreeDofs(), coord_x=coordinates[:, 0].tolist(), coord_y=coordinates[:, 1].tolist(),
                     coord_z=coordinates[:, 2].tolist(), cycle_type=1, print_level=0)
    if fes.is_complex:
        return ssn.ComplexCompactAMSPreconditioner(a_real_mat=aux, ndof_complex=ndof, **arguments)
    return ssn.CompactAMSPreconditioner(mat=aux, **arguments)


def _solve(space, coefficients, rhs, *, Backend, lift, preconditioner, bonus, tolerance, max_iterations, cycles):
    _check_backend(Backend)
    if Backend == "cpu":
        return cpu_solve(space, coefficients, rhs, lift=lift, preconditioner=preconditioner, bonus=bonus, tolerance=tolerance,
                         max_iterations=max_iterations, cycles=cycles)
    if preconditioner == "direct":
        raise ValueError("the device route has no direct solver; select amg, ams or bddc")
    system = DeviceSystem(space, coefficients, bonus=bonus)
    system.setup(preconditioner, cycles=cycles)
    return system.solve(rhs, lift=lift, tolerance=tolerance, max_iterations=max_iterations)


def solve_h1(space, coefficients, rhs, *, Backend, lift=None, preconditioner=None, bonus=(0, 0, 0), tolerance=1e-10,
             max_iterations=500, cycles=None):
    """Solve an H1 form: ``amg`` with one DOF per vertex, ``bddc`` otherwise; no implicit backend change."""
    if _space_kind(space) != "H1":
        raise ValueError("solve_h1 needs an H1 space")
    selected = preconditioner or ("amg" if space.ndof == space.mesh.nv else "bddc")
    return _solve(space, coefficients, rhs, Backend=Backend, lift=lift, preconditioner=selected, bonus=bonus, tolerance=tolerance,
                  max_iterations=max_iterations, cycles=cycles)


def solve_hcurl(space, coefficients, rhs, *, Backend, lift=None, preconditioner=None, bonus=(0, 0, 0), tolerance=1e-8,
                max_iterations=500, cycles=None):
    """Solve a curl-curl form: ``ams`` with one DOF per edge, ``bddc`` over AMS otherwise (p>1, or HCurl order 1 on hexes)."""
    if _space_kind(space) != "HCurl":
        raise ValueError("solve_hcurl needs an HCurl space")
    selected = preconditioner or ("ams" if space.ndof == space.mesh.nedge else "bddc")
    return _solve(space, coefficients, rhs, Backend=Backend, lift=lift, preconditioner=selected, bonus=bonus, tolerance=tolerance,
                  max_iterations=max_iterations, cycles=cycles)
