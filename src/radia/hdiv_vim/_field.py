"""Correct field reconstruction for the HDiv-VIM.

Given a SOLVED magnetization gfM (an HDiv GridFunction, any order), reconstruct the magnetic field (B or H)
at arbitrary points via the EXACT analytic field of the per-element magnetization (radia.Fld -- the
documented "NGSolve magnetization -> Radia open-boundary field" pipeline).

DO NOT reconstruct the demag field as M_mass^{-1} N m (the energy operator applied + mass-inverted): the
demag operator N = B^T G B has a SOLENOIDAL NULLSPACE at order>=2 (divergence-free RT bubbles carry zero
charge, so B maps them to zero), so M_mass^{-1} N m does NOT equal the field -- it gives GARBAGE point
fields at order>=2 (it happens to work only at order 0/1).  The demag FACTOR (the energy Rayleigh quotient
m^T N m / m^T M_mass m) stays correct because it is the energy, not the field.  See
tests/feec/test_hdiv_vim_field_reconstruction.py (uniform-M sphere: B_inside = MU0*(2/3)*M at ALL orders).
"""
import numpy as np
import ngsolve as ng


def reconstruct_field(mesh, gfM, points, quantity="b", units="m"):
    """Magnetic field at `points` from a solved HDiv-VIM magnetization gfM, via the exact analytic field of
    the per-element M (radia.Fld).  CALLER wraps in TaskManager when combining with other NGSolve work.

    mesh     : the NGSolve tet mesh gfM lives on (the magnetized body, e.g. steel-only).
    gfM      : HDiv GridFunction (any order) holding the solved magnetization M(x) [A/m].
    points   : (N,3) array-like of query points (same length units as the mesh).
    quantity : 'b' (Tesla, = MU0*(H_demag + M) inside the body) or 'h' (A/m, the demag/stray H field).
    returns  : (N,3) ndarray of the field FROM THE MAGNETIZATION.  Add any SOURCE/coil field separately --
               radia.Fld here is the magnetization's contribution only.

    Per-element M is the centroid value (== the average for <= linear M; ~ for quadratic).  This is the
    field of the piecewise-constant-M body: EXACT for uniform M at ALL orders, and it carries the
    inter-element M variation the high-order solve resolved.  (Reconstructing a per-element POLYNOMIAL M
    field would need element subdivision or a polynomial-charge field kernel -- a future refinement.)
    """
    import radia as rad
    from radia.netgen_mesh_import import netgen_mesh_to_radia

    Mel = []
    for el in mesh.Elements(ng.VOL):
        c = np.array([mesh[v].point for v in el.vertices]).mean(0)   # element centroid
        mp = mesh(c[0], c[1], c[2])
        Mel.append([float(gfM[i](mp)) for i in range(3)])
    cont = netgen_mesh_to_radia(mesh, material=lambda i: {"magnetization": Mel[i]},
                                units=units, verbose=False)
    pts = np.asarray(points, float).reshape(-1, 3)
    return np.array([rad.Fld(cont, quantity, [float(p[0]), float(p[1]), float(p[2])]) for p in pts], float)
