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


# ---------------------------------------------------------------------------------------------------
# Polynomial-charge field kernel (the order>=2 field: the centroid reconstruct_field above is only
# correct for piecewise-constant M; for a polynomial M the VOLUME charge rho = -div M is non-zero and
# its omission is a 90-230% error at div M != 0 -- measured, see tests/feec/test_hdiv_vim_poly_field.py).
# ---------------------------------------------------------------------------------------------------
def _scalar(val):
    return float(val[0] if isinstance(val, tuple) else val)


def reconstruct_field_polynomial(mesh, gfM, points, quad=4, quantity="h", include_volume=True):
    """EXTERNAL-point magnetic field of a polynomial magnetization gfM (HDiv order p), via its FULL
    polynomial charge -- volume rho = -div M (L2 order p-1) AND surface sigma = M.n (SurfaceL2 order p):

        H(r) = (1/4pi)[ INT_V (-div M)(r') (r-r')/|r-r'|^3 dV' + INT_S (M.n)(r') (r-r')/|r-r'|^3 dS' ]

    This is the order>=2 generalisation of `reconstruct_field` (which uses the per-element CENTROID M ==
    surface charge only; it DROPS the volume charge, a 90-230% error wherever div M != 0).  H = -grad phi_M
    with phi_M the magnetic scalar potential of the charges; verified on the uniform sphere (center -M/3,
    external dipole) and a linear-M body (volume-charge term essential; coarse->fine self-convergent).

    ELEMENT-TYPE AGNOSTIC: the quadrature points / weights / normals come from NGSolve's own
    ElementTransformation + IntegrationRule(el.type) + specialcf.normal, so the SAME code handles
    tetrahedral AND hexahedral (and prism) meshes, flat or curved (mip.weight carries the curved
    Jacobian, specialcf.normal the curved outward normal).

    mesh     : NGSolve mesh (tet or hex) gfM lives on.
    gfM      : HDiv GridFunction (order p) of the solved magnetization M(x) [A/m].
    points   : (N,3) query points -- EXTERNAL to the body (the integrands are 1/r^2-singular; an INTERNAL/
               near point needs the singular-aware kernel, a future step).
    quad     : controls the volume/surface integration order (intorder = 2*quad); raise for near-surface
               observation.
    quantity : 'h' (A/m, the demag/stray H) or 'b' (Tesla); 'b' = MU0*H at an EXTERNAL point (no M there).
    include_volume : drop the rho term when False (== the centroid/surface-only field) -- for diagnostics.
    returns  : (N,3) ndarray of the field FROM THE MAGNETIZATION (add any source/coil field separately).

    CALLER wraps in TaskManager when combining with other NGSolve work.  This is a reference (Python)
    implementation; a C++/H-matrix-accelerated version is the next productionisation step.
    """
    divM = ng.div(gfM)
    nsurf = ng.specialcf.normal(mesh.dim)
    sigma_cf = ng.InnerProduct(gfM.Trace(), nsurf)      # M.n on the boundary (correct outward normal)
    obs = np.asarray(points, float).reshape(-1, 3)
    H = np.zeros((len(obs), 3))
    inv4pi = 1.0 / (4.0 * np.pi)
    iord = 2 * quad

    if include_volume:                                  # volume charge rho = -div M
        for i in range(mesh.GetNE(ng.VOL)):
            ei = ng.ElementId(ng.VOL, i)
            trafo = mesh.GetTrafo(ei)
            for ip in ng.IntegrationRule(mesh[ei].type, iord):
                mip = trafo(ip)
                xq = np.array([mip.point[k] for k in range(3)])
                rho = -_scalar(divM(mip))
                d = obs - xq
                rn = np.linalg.norm(d, axis=1)
                H += inv4pi * (ip.weight * mip.measure) * rho * d / rn[:, None] ** 3

    for i in range(mesh.GetNE(ng.BND)):                 # surface charge sigma = M.n
        ei = ng.ElementId(ng.BND, i)
        trafo = mesh.GetTrafo(ei)
        for ip in ng.IntegrationRule(mesh[ei].type, iord):
            mip = trafo(ip)
            xq = np.array([mip.point[k] for k in range(3)])
            sig = _scalar(sigma_cf(mip))
            d = obs - xq
            rn = np.linalg.norm(d, axis=1)
            H += inv4pi * (ip.weight * mip.measure) * sig * d / rn[:, None] ** 3

    if quantity.lower() == "b":
        return (4e-7 * np.pi) * H                        # external point: B = MU0 * H (no M present)
    return H
