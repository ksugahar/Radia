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


def _conductor_tets(mesh, names):
    """Vertex ids (0-based), element numbers and coordinates of conductor tets."""
    elements = mesh.ngmesh.Elements3D().NumPy()
    materials = mesh.GetMaterials()
    labels = np.asarray([materials[i - 1] for i in elements["index"]])
    mask = np.isin(labels, list(names))
    if not mask.any():
        raise ValueError("nonempty tetrahedral conductor required")
    if np.any(elements["np"][mask] != 4):
        raise ValueError("closed-coil current requires straight tet4 conductor cells")
    tets = elements["nodes"][mask][:, :4].astype(np.int64) - 1
    points = np.array([p.p for p in mesh.ngmesh.Points()], dtype=float)
    return tets, np.flatnonzero(mask), points


def _tet_faces(tets):
    """Unique faces of a tet list: (sorted vertex triples, owner pairs, local ids)."""
    local = np.array([[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]])
    triples = np.sort(tets[:, local].reshape(-1, 3), axis=1)
    owners = np.repeat(np.arange(len(tets)), 4)
    unique, inverse, counts = np.unique(triples, axis=0, return_inverse=True,
                                        return_counts=True)
    if np.any(counts > 2):
        raise ValueError("nonmanifold conductor face")
    order = np.argsort(inverse, kind="stable")
    sorted_inverse = inverse.ravel()[order]
    starts = np.searchsorted(sorted_inverse, np.arange(len(unique)))
    first = owners[order][starts]
    second = np.full(len(unique), -1, dtype=np.int64)
    shared = counts == 2
    second[shared] = owners[order][starts[shared] + 1]
    return unique, first, second


def _connected(n_cells, first, second):
    """True when the face-adjacency graph of the conductor is connected."""
    shared = second >= 0
    parent = np.arange(n_cells)

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for a, b in zip(first[shared], second[shared]):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    return len({find(i) for i in range(n_cells)}) == 1


def solve_closed_coil_current_phi(mesh, *, current_A, cut_origin, cut_normal,
                                  cut_radius_m, materials=("coil",),
                                  inverse="sparsecholesky"):
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
    conductivity; straight tets; one face-connected conductor. The caller owns
    TaskManager for assembly.

    Returns ``current`` (piecewise-constant VectorL2 GridFunction in A/m^2 on
    the conductor), ``phi`` and ``stats`` with the weak-divergence residual,
    the face-averaged flux through the cut and the timing.
    """
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
    tets, element_numbers, points = _conductor_tets(mesh, names)
    faces, first, second = _tet_faces(tets)
    if not _connected(len(tets), first, second):
        raise ValueError("conductor must be face-connected; split independent coils")

    # Cells near the cut: signed distance of the centroid and in-plane radius.
    centroids = points[tets].mean(axis=1)
    offset = centroids - origin
    signed = offset @ normal
    in_plane = np.linalg.norm(offset - np.outer(signed, normal), axis=1)
    near = in_plane < radius
    shared = second >= 0
    a, b = first[shared], second[shared]
    crossing = near[a] & near[b] & ((signed[a] < 0) != (signed[b] < 0))
    cut_faces = faces[shared][crossing]
    if len(cut_faces) == 0:
        raise ValueError("cut plane does not cross the conductor within cut_radius_m")
    # Completeness: every border edge of the cut sheet must lie on the wall.
    edge_pairs = np.sort(np.concatenate([cut_faces[:, [0, 1]], cut_faces[:, [0, 2]],
                                         cut_faces[:, [1, 2]]]), axis=1)
    edges, edge_count = np.unique(edge_pairs, axis=0, return_counts=True)
    border = edges[edge_count == 1]
    wall = faces[~shared]
    wall_edges = np.unique(np.sort(np.concatenate([wall[:, [0, 1]], wall[:, [0, 2]],
                                                   wall[:, [1, 2]]]), axis=1), axis=0)
    wall_set = {tuple(e) for e in wall_edges}
    loose = [tuple(e) for e in border if tuple(e) not in wall_set]
    if loose:
        raise ValueError(f"cut does not span the conductor section ({len(loose)} interior border edges)")


    # Jump function tau: one at cut vertices seen from the positive side, zero
    # elsewhere; h = grad tau is piecewise constant on its support cells.
    cut_vertices = np.unique(cut_faces)
    on_cut = np.isin(tets, cut_vertices)
    support = near & (signed >= 0) & on_cut.any(axis=1)
    edge_matrix = points[tets[:, 1:]] - points[tets[:, :1]]       # rows x_i - x_0
    volume = np.abs(np.linalg.det(edge_matrix)) / 6.0
    if np.any(volume <= 0.0):
        raise ValueError("degenerate conductor tetrahedron")
    lambda_grad = np.transpose(np.linalg.inv(edge_matrix), (0, 2, 1))   # grad lambda_1..3
    grads = np.concatenate([-lambda_grad.sum(axis=1, keepdims=True), lambda_grad], axis=1)
    jump_gradient = np.einsum("nk,nkd->nd", (on_cut & support[:, None]).astype(float), grads)

    # Order-0 L2 carries one dof per element, numbered by the element.
    scalar = ng.L2(mesh, order=0)
    if (scalar.ndof != mesh.ne or scalar.GetDofNrs(ng.ElementId(ng.VOL, mesh.ne - 1))[0] != mesh.ne - 1):
        raise RuntimeError("unexpected order-0 L2 dof layout")

    def piecewise(values):
        components = []
        for k in range(3):
            gf = ng.GridFunction(scalar)
            data = np.zeros(mesh.ne)
            data[element_numbers] = values[:, k]
            gf.vec.FV().NumPy()[:] = data
            components.append(gf)
        return ng.CF(tuple(components)), components

    h_cf, _ = piecewise(jump_gradient)
    region = mesh.Materials("|".join(re.escape(name) for name in names))
    potential_space = ng.H1(mesh, order=1, definedon=region)     # dof = vertex number
    u, v = potential_space.TnT()
    a_form = ng.BilinearForm(potential_space, symmetric=True)
    a_form += ng.grad(u) * ng.grad(v) * ng.dx(definedon=region)
    f_form = ng.LinearForm(potential_space)
    f_form += -ng.InnerProduct(h_cf, ng.grad(v)) * ng.dx(definedon=region)
    a_form.Assemble()
    f_form.Assemble()
    free = potential_space.FreeDofs()
    free[int(tets[0, 0])] = False  # One potential gauge for the connected conductor.
    phi = ng.GridFunction(potential_space, name="electric_scalar_potential")
    phi.vec.data = a_form.mat.Inverse(free, inverse=inverse) * f_form.vec
    phi_values = phi.vec.FV().NumPy()
    field = np.einsum("nk,nkd->nd", phi_values[tets], grads) + jump_gradient
    energy = float(np.sum(volume * np.einsum("nd,nd->n", field, field)))
    if not math.isfinite(energy) or energy <= 0.0:
        raise ValueError("cut does not drive a resolved closed current path")
    # Flux of E through the cut along +n is -energy (tau jumps up across it).
    scale = -current_A / energy
    density = scale * field
    current, current_components = piecewise(density)

    # Weak divergence against every conductor P1 function (the Nedelec gradients).
    divergence = ng.LinearForm(potential_space)
    divergence += ng.InnerProduct(current, ng.grad(v)) * ng.dx(definedon=region)
    drive = ng.LinearForm(potential_space)
    drive += scale * ng.InnerProduct(h_cf, ng.grad(v)) * ng.dx(definedon=region)
    divergence.Assemble()
    drive.Assemble()
    drive_norm = float(np.linalg.norm(drive.vec.FV().NumPy()))
    relative_divergence = (float(np.linalg.norm(divergence.vec.FV().NumPy())) / drive_norm
                           if drive_norm > 0 else 0.0)
    if current_A != 0.0 and not relative_divergence <= 1e-8:
        raise RuntimeError(f"weak current divergence gate failed ({relative_divergence:.2e})")
    # Independent check: face-averaged flux through the cut sheet along +n.
    triangle = points[cut_faces]
    area_normal = 0.5 * np.cross(triangle[:, 1] - triangle[:, 0], triangle[:, 2] - triangle[:, 0])
    area_normal *= np.sign(area_normal @ normal)[:, None]
    mean_density = 0.5 * (density[a[crossing]] + density[b[crossing]])
    face_flux = float(np.sum(np.einsum("nd,nd->n", mean_density, area_normal)))
    return {"current": current, "components": current_components, "phi": phi,
            "stats": {"method": "A-phi H1 potential with one thick cut",
                      "current_A": current_A, "conductor_cells": int(len(tets)),
                      "cut_faces": int(len(cut_faces)), "jump_support_cells": int(support.sum()),
                      "phi_ndof": int(sum(free)), "relative_weak_divergence": relative_divergence,
                      "cut_face_flux_A": face_flux,
                      "cut_face_flux_relative_error": (abs(face_flux - current_A) / abs(current_A)
                                                       if current_A else abs(face_flux)),
                      "runtime_s": time.perf_counter() - started}}
