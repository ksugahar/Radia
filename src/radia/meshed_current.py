"""Impressed DC current in a connected closed conductor.

Two discretisations of the same resistive (minimum-dissipation) current:

* :func:`solve_closed_coil_current` -- mixed RT0/P0, exactly divergence-free.
* :func:`solve_closed_coil_current_phi` -- A-phi primal form: an electric scalar
  potential in H1 with one thick cut across the loop. Its piecewise-constant J
  is weakly divergence-free against every P1 function, i.e. orthogonal to the
  gradients inside lowest-order Nedelec spaces, which keeps a curl-curl
  right-hand side compatible without projection.
"""
from __future__ import annotations

import math
import re
import time

import numpy as np
import ngsolve as ng


def solve_closed_coil_current(mesh, drive, *, current_A, drive_length_m,
                              materials=("coil",), inverse="pardiso"):
    """Solve the RT0/P0 minimum-energy current with insulating conductor walls.

    ``drive`` must describe a unit tangential drive along a complete straight
    leg of length ``drive_length_m``. Its volume pairing with a conservative
    current is then length times section current. Geometry/source identity is
    the caller's responsibility. Only one connected, straight-tet conductor
    is supported. Use a caller-owned TaskManager for parallel assembly.

    The field minimizes integral |J|^2 with zero discrete divergence, in
    contrast to a continuous scalar-potential gradient whose normal component
    can jump between cells. No reference-solver field enters this solve.
    """
    current_A, length = float(current_A), float(drive_length_m)
    if not math.isfinite(current_A) or not math.isfinite(length) or length <= 0:
        raise ValueError("finite current and positive finite drive length required")
    names = tuple(materials)
    if not names or not set(names) <= set(mesh.GetMaterials()):
        raise ValueError("conductor material labels are missing")
    if mesh.GetCurveOrder() != 1:
        raise ValueError("closed-coil current currently requires straight tet4 geometry")
    if drive.dim != 3:
        raise ValueError("drive must be a three-component coefficient function")
    region = mesh.Materials("|".join(re.escape(name) for name in names))
    cells = [el for el in mesh.Elements() if el.mat in names]
    if not cells or any(el.type != ng.ET.TET for el in cells):
        raise ValueError("nonempty tetrahedral conductor required")
    started = time.perf_counter()
    owners = {}
    adjacent = {el.nr: [] for el in cells}
    for el in cells:
        for face in el.faces:
            owners.setdefault(face.nr, []).append(el.nr)
    for sharing in owners.values():
        if len(sharing) > 2:
            raise ValueError("nonmanifold conductor face")
        if len(sharing) == 2:
            a, b = sharing
            adjacent[a].append(b)
            adjacent[b].append(a)
    visited, pending = set(), [cells[0].nr]
    while pending:
        cell = pending.pop()
        if cell not in visited:
            visited.add(cell)
            pending.extend(adjacent[cell])
    if len(visited) != len(cells):
        raise ValueError("conductor must be face-connected; split independent coils")
    V = ng.Compress(ng.HDiv(mesh, order=0, RT=True, definedon=region))
    Q = ng.Compress(ng.L2(mesh, order=0, definedon=region))
    X = V * Q
    (u, p), (v, q) = X.TnT()
    a = ng.BilinearForm(X, symmetric=True)
    a += (ng.InnerProduct(u, v) + ng.div(u)*q + ng.div(v)*p)*ng.dx(definedon=region)
    f = ng.LinearForm(X)
    f += ng.InnerProduct(drive, v)*ng.dx(definedon=region)
    free = X.FreeDofs()
    for face, sharing in owners.items():
        if len(sharing) == 1:
            dofs = V.GetDofNrs(ng.NodeId(ng.FACE, face))
            if len(dofs) != 1 or dofs[0] < 0:
                raise RuntimeError("unexpected RT0 wall DOF layout")
            free[dofs[0]] = False
    free[V.ndof] = False  # One pressure gauge for the connected conductor.
    solution = ng.GridFunction(X)
    a.Assemble()
    f.Assemble()
    solution.vec.data = a.mat.Inverse(free, inverse=inverse)*f.vec
    flux = solution.components[0]
    pairing = float(ng.Integrate(ng.InnerProduct(flux, drive), mesh, definedon=region))
    drive_energy = float(ng.Integrate(ng.InnerProduct(drive, drive), mesh, definedon=region))
    if (not math.isfinite(pairing) or not math.isfinite(drive_energy)
            or pairing <= 1e-12*max(drive_energy, 1e-300)):
        raise ValueError("drive does not excite a resolved closed current path")
    current = ng.GridFunction(V, name="conservative_current_density")
    current.vec.data = (current_A*length/pairing)*flux.vec
    div2 = float(ng.Integrate(ng.div(current)**2, mesh, definedon=region))
    energy = float(ng.Integrate(ng.InnerProduct(current, current), mesh, definedon=region))
    relative_divergence = length*math.sqrt(max(div2, 0)/max(energy, 1e-300))
    if not math.isfinite(relative_divergence) or relative_divergence > 1e-8:
        raise RuntimeError("current divergence gate failed")
    return {"current": current, "space": V,
            "stats": {"current_A": current_A, "drive_length_m": length,
                      "conductor_cells": len(cells), "ndof": X.ndof,
                      "relative_divergence": relative_divergence,
                      "section_current_A": float(ng.Integrate(ng.InnerProduct(current, drive),
                                                  mesh, definedon=region))/length,
                      "runtime_s": time.perf_counter()-started}}


class _CoilCurrentResult(dict):
    """Result dict whose ``"phi"`` GridFunction is built on first access."""

    def __init__(self, build_phi, values):
        super().__init__(values)
        self._build_phi = build_phi

    def __missing__(self, key):
        if key != "phi":
            raise KeyError(key)
        self["phi"] = value = self._build_phi()
        return value


def solve_closed_coil_current_phi(mesh, *, current_A, cut_origin, cut_normal,
                                  cut_radius_m, materials=("coil",),
                                  inverse="iccg"):
    """A-phi DC current of a closed coil: J = -sigma (grad phi + h), phi in H1.

    ``h`` is the gradient of a function that jumps by one across a thick cut:
    the conductor faces crossing the plane through ``cut_origin`` with normal
    ``cut_normal``, restricted to cells whose centroid lies within
    ``cut_radius_m`` of the origin in that plane (so a plane crossing the loop
    twice cuts it once, as a current condition with a direction and an origin
    does). The cut must span the conductor cross-section completely; its
    border must lie on the conductor surface, otherwise this fails loudly.
    Conductor walls are insulating (natural J.n = 0). ``current_A`` is the net
    current through the cut in the direction of ``cut_normal``. Uniform
    conductivity; straight tets; one face-connected conductor.

    The computation is one native pass over the conductor only
    (``radia.sparsesolv_ngsolve.ClosedCoilCurrentPhi``): its tets, faces and
    connectivity, the cut and its completeness, the jump function, the P1
    Laplacian on the conductor vertices (one vertex gauged; the default
    ``inverse="iccg"`` is IC(0)-CG to a true relative residual of 1e-12 and
    fails if it does not converge, any NGSolve inverse type such as
    ``"sparsecholesky"`` solves directly), the current, the weak-divergence
    gate and the cut-face flux. Parallel under a caller-owned TaskManager.

    Returns ``current`` (piecewise-constant GridFunction components in A/m^2
    on the conductor, as a vector CoefficientFunction), ``components``,
    ``phi`` and ``stats`` with the weak-divergence residual, the face-averaged
    flux through the cut and the timing. ``phi`` (a GridFunction on an H1
    space over the mesh, restricted to the conductor) is built on first access.
    """
    import radia.sparsesolv_ngsolve as ssn

    current_A = float(current_A)
    origin = np.asarray(cut_origin, dtype=float).reshape(3)
    normal = np.asarray(cut_normal, dtype=float).reshape(3)
    radius = float(cut_radius_m)
    if not (math.isfinite(current_A) and np.all(np.isfinite(origin))
            and np.all(np.isfinite(normal)) and math.isfinite(radius) and radius > 0):
        raise ValueError("finite current, cut origin, cut normal and positive cut radius required")
    if np.linalg.norm(normal) < 1e-12:
        raise ValueError("cut normal must be nonzero")
    normal = normal / np.linalg.norm(normal)
    names = tuple(materials)
    if not names or not set(names) <= set(mesh.GetMaterials()):
        raise ValueError("conductor material labels are missing")
    if mesh.GetCurveOrder() != 1:
        raise ValueError("closed-coil current currently requires straight tet4 geometry")
    started = time.perf_counter()
    indices = [index for index, name in enumerate(mesh.GetMaterials()) if name in set(names)]
    native = ssn.ClosedCoilCurrentPhi(mesh, indices, current_A, origin.tolist(), normal.tolist(),
                                      radius, inverse)
    element_numbers = np.asarray(native["elements"])
    density = np.asarray(native["density"])

    # Order-0 L2 carries one dof per element, numbered by the element.
    scalar = ng.L2(mesh, order=0)
    if (scalar.ndof != mesh.ne or scalar.GetDofNrs(ng.ElementId(ng.VOL, mesh.ne - 1))[0] != mesh.ne - 1):
        raise RuntimeError("unexpected order-0 L2 dof layout")
    components = []
    for k in range(3):
        gf = ng.GridFunction(scalar)
        data = gf.vec.FV().NumPy()
        data[:] = 0.0
        data[element_numbers] = density[:, k]
        components.append(gf)
    current = ng.CF(tuple(components))
    region = mesh.Materials("|".join(re.escape(name) for name in names))
    vertices = np.asarray(native["vertices"], dtype=np.int64)
    potential = np.asarray(native["phi"])

    def build_phi():
        potential_space = ng.H1(mesh, order=1, definedon=region)     # dof = vertex number
        phi = ng.GridFunction(potential_space, name="electric_scalar_potential")
        phi.vec.FV().NumPy()[:] = 0.0
        phi.vec.FV().NumPy()[vertices] = potential
        return phi

    face_flux = float(native["cut_face_flux_A"])
    return _CoilCurrentResult(build_phi, {
            "current": current, "components": components,
            "stats": {"method": "A-phi H1 potential with one thick cut",
                      "current_A": current_A, "conductor_cells": int(len(element_numbers)),
                      "cut_faces": int(native["cut_faces"]),
                      "jump_support_cells": int(native["support_cells"]),
                      "phi_ndof": int(native["phi_ndof"]),
                      "relative_weak_divergence": float(native["relative_weak_divergence"]),
                      "cut_face_flux_A": face_flux,
                      "cut_face_flux_relative_error": (abs(face_flux - current_A) / abs(current_A)
                                                       if current_A else abs(face_flux)),
                      "inverse": inverse, "native_timing_s": dict(native["timing"]),
                      "runtime_s": time.perf_counter() - started}})
