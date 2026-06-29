"""kelvin_validate.py

Layer 4 of the Kelvin helper API (api_plan.md): end-to-end validation
harness. Runs a full-A FEM self-inductance computation on a Sugahara
two-sphere Kelvin mesh and compares it against a Radia analytical
reference (via rad.Fld 'a' energy integral) and, optionally, against a
user-supplied closed-form L_analytical.

Initial M3 cut (2026-04-15): single-coil self-inductance only. The
caller supplies an OCC inner shape list, a matching ``rad`` container
representing the same physical current, and a J_src CF for the FEM
leg. The harness builds the Kelvin geometry, meshes, solves, and
reports the three L values.
"""

from __future__ import annotations

import math
import time

import numpy as np

from ngsolve import (CoefficientFunction as CF, GridFunction,
                      Integrate, InnerProduct, TaskManager, VectorH1,
                      Mesh)
from netgen.occ import OCCGeometry

from radia.kelvin_geometry import add_kelvin_exterior_domain
from radia.kelvin_solver import (solve_full_A_kelvin, inductance_from_energy,
                                 NU_0)


def _per_vertex_A_from_radia(mesh, radia_handle):
    """Evaluate Radia A (Cartesian) at every vertex of ``mesh``.

    Uses the physical (un-pulled-back) A; for the J . A energy integral
    only the coil region contributes and the coil lives in the inner
    physical domain (no Kelvin map needed for those vertices).
    """
    import radia as rad

    nv = mesh.nv
    arr = np.zeros((nv, 3), dtype=float)
    for vi in range(nv):
        p = mesh.vertices[vi].point
        arr[vi] = rad.Fld(radia_handle, 'a', [p[0], p[1], p[2]])
    return arr


def L_radia_reference_energy(mesh, radia_handle, J_source_cf,
                               coil_material, I_total,
                               vh1_order=1, int_order=10):
    """Radia reference inductance via ``W = 0.5 * int J . A dV``.

    A is evaluated at mesh vertices by Radia and loaded into a
    VectorH1 GridFunction; the energy integral is restricted to
    ``coil_material``. Caveat: only vertex dof are filled, so the
    effective A approximation is O(h) regardless of vh1_order > 1
    unless edge/face dof are also populated. Use this leg as a
    cross-check of magnitude and sign, not as a sub-percent truth.
    For a tighter reference leg, supply ``L_analytical`` to
    ``compare_against_radia_self_inductance``.
    """
    A_per_vert = _per_vertex_A_from_radia(mesh, radia_handle)
    fes_vec = VectorH1(mesh, order=vh1_order)
    gfu_A = GridFunction(fes_vec)
    nv = mesh.nv
    data = gfu_A.vec.FV().NumPy()
    data[:nv] = A_per_vert[:, 0]
    data[nv:2 * nv] = A_per_vert[:, 1]
    data[2 * nv:3 * nv] = A_per_vert[:, 2]

    W = Integrate(0.5 * InnerProduct(J_source_cf, gfu_A),
                   mesh, definedon=mesh.Materials(coil_material),
                   order=int_order).real
    return 2.0 * W / (I_total ** 2), W


def compare_against_radia_self_inductance(
        inner_shapes, radia_handle, J_source_cf,
        kelvin_offset=(0.15, 0.0, 0.0), R_K=0.06,
        maxh=12e-3, source_material="coil",
        order=1, I_total=1.0,
        L_analytical=None,
        skip_radia=False,
        verbose=True):
    """Build Kelvin mesh, run FEM, compute Radia reference, report diffs.

    Args:
        inner_shapes: OCC shape or list of sub-shapes for the inner
            domain (must include a face named ``"kelvin_int"`` on the
            outward sphere). The source region must carry material
            name ``source_material``.
        radia_handle: Radia object whose ``rad.Fld(., 'a', pt)`` yields
            the physical A of the coil (same physical current as
            ``J_source_cf``).
        J_source_cf: NGSolve vector CF for J_src (A/m^2), supported on
            ``source_material``.
        kelvin_offset, R_K, maxh: Kelvin geometry parameters.
        source_material: mesh material name for the current volume.
        order: HCurl polynomial order.
        I_total: total coil current (A).
        L_analytical: optional closed-form L for a third leg of the
            comparison.
        verbose: print progress lines.

    Returns:
        dict with keys ``L_FEM``, ``L_Radia``, ``L_analytical`` (or None),
        ``mesh``, ``gfu``, ``nu_cf``, and timing entries
        ``t_mesh``, ``t_fem``, ``t_radia``.
    """
    t_mesh_0 = time.perf_counter()
    geometry, info = add_kelvin_exterior_domain(
        inner_shapes, offset=kelvin_offset, R_K=R_K, inner_maxh=maxh)
    # NOTE: Caller MUST be inside `with TaskManager():` per CLAUDE.md
    # "Caller Wraps, Helper Does NOT" (2026-05-27).
    ngmesh = OCCGeometry(geometry).GenerateMesh(maxh=maxh, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(order + 1)
    t_mesh = time.perf_counter() - t_mesh_0
    if verbose:
        print(f"  mesh: {mesh.ne} tets, {mesh.nv} verts, "
              f"mats={mesh.GetMaterials()} ({t_mesh:.1f}s)")

    t_fem_0 = time.perf_counter()
    res = solve_full_A_kelvin(mesh, J_source_cf,
                                R_K=R_K, offset=kelvin_offset,
                                source_material=source_material,
                                order=order)
    L_FEM, _ = inductance_from_energy(res["gfu"], res["nu_cf"], mesh,
                                        I_total)
    t_fem = time.perf_counter() - t_fem_0
    if verbose:
        print(f"  FEM   L = {L_FEM*1e9:.4f} nH ({t_fem:.1f}s, "
              f"ndof={res['fes'].ndof})")

    if skip_radia or radia_handle is None:
        L_Radia = None
        t_radia = 0.0
    else:
        t_radia_0 = time.perf_counter()
        L_Radia, _ = L_radia_reference_energy(
            mesh, radia_handle, J_source_cf,
            coil_material=source_material, I_total=I_total)
        t_radia = time.perf_counter() - t_radia_0
        if verbose:
            print(f"  Radia L = {L_Radia*1e9:.4f} nH ({t_radia:.1f}s, "
                  f"vertex-projected; O(h) approximation)")
            print(f"    FEM vs Radia: {(L_FEM/L_Radia - 1)*100:+.2f}%")
    if verbose and L_analytical is not None:
        print(f"    Analytical  = {L_analytical*1e9:.4f} nH")
        print(f"    FEM vs ana  : {(L_FEM/L_analytical - 1)*100:+.2f}%")
        if L_Radia is not None:
            print(f"    Rad vs ana  : {(L_Radia/L_analytical - 1)*100:+.2f}%")

    return {
        "L_FEM": L_FEM,
        "L_Radia": L_Radia,
        "L_analytical": L_analytical,
        "mesh": mesh,
        "gfu": res["gfu"],
        "nu_cf": res["nu_cf"],
        "t_mesh": t_mesh,
        "t_fem": t_fem,
        "t_radia": t_radia,
    }
