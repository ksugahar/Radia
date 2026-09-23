"""Locally conservative impressed current in a connected closed conductor."""
from __future__ import annotations

import math
import re
import time

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
