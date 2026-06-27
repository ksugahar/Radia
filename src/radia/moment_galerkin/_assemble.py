"""radia.moment_galerkin._assemble -- moment-basis assembly for the SYMMETRIC moment-Galerkin MMMM demag.

The moment-Galerkin demag operator is N = B^T G B, the SAME symmetric structure HDiv-VIM uses, but on the
MMMM moment basis instead of HDiv's RT0 flux.  G is the EXACT analytic charge-Gram H-matrix
(`radia._radia_pybind._ChargeGramHMatrix`, the same C++ kernel HDiv-VIM ships); B maps the per-element moment
amplitudes to the face surface charges sigma; M_mass is the magnetization mass.  The physical SPD demag system

    ( (1/chi) M_mass + B^T G B ) m = M_mass H_ext

is solved by the existing C++ Krylov kernel.  The per-element moment-basis assembly (B, M_mass, the boundary
triangle soup) is the C++ `RadMomentGalerkinAssemble` parallel kernel (ngcore::ParallelFor, exposed as
`_rp._moment_galerkin_assemble`); this file is the thin Python packaging layer (COO -> csr).

MIXED ELEMENTS (tet / pyramid / wedge / hex -- the old yano-type element zoo): each element is a vertex array
of 4 (tetrahedron), 5 (pyramid), 6 (wedge/prism) or 8 (hexahedron) points, in the radia / netgen corner order.
The per-element DOF count is 3 + max(0, nFace - 4):  tet 3 (dipole only), wedge/pyramid 4 (+1 quad mode),
hex 5 (+2 quad modes).  A single mixed list assembles into one B / M_mass / N -- the DOF layout is the
per-element `dof_off` / `ndof_elem` (variable stride), not a uniform stride.

Two moment orders (`quad=`):
  quad=False (default) -- 3 DOF/element: a CONSTANT magnetization M_0 (the dipole / standard demag).
      B col k = sigma_f = M_0 . n_f = n_f[k];  M_mass = diag(V_e).
  quad=True            -- 3 + (nFace-4) DOF/element: dipole + the geometry-adaptive QUAD residual-eigenmode
      amplitudes (phi_q = the zero-monopole, zero-dipole-moment face-charge modes, the same Gram-Schmidt basis
      the C++ collocation MMMM `momentResidualEigenmodes` uses).  The quad mode's magnetization is the LINEAR
      field M_q(x) = Dm_q.(x-c), so the quad MASS is the magnetization Gram W_quad[q,r] = Dm_q : M2 : Dm_r,
      M2 = INT (x-c)(x-c) dV (a unified star-triangulation sub-tet quadrature, exact on every element type).
      The quad DOF buy per-element accuracy under SKEW / gradient loads; they are ~null for axial / symmetric
      loads (correct physics, never worse), and a tetrahedron has 0 quad modes (nFace = 4).

Symmetry is by construction: G is symmetric (Coulomb reciprocity) -> N = B^T G B symmetric to machine
precision (validated 1e-16 on tet/wedge/pyramid/hex + mixed) -> loop modes field-null -> mu_r-independent /
loop-free convergence with no loop-star.

ACA defaults (leaf=40, eta=0.5): tuned so the H-matvec reproduces the exact analytic Gram to ~1e-4 on
near-heavy geometry (eta is the dominant knob; the charge group must not be split, so leaf~40 keeps the
triangles of a few elements intact).

FAST BUILD defaults (near_factor=2, far_quad=4): the SAME precision-preserving NEAR/FAR split HDiv-VIM ships on
this shared charge-Gram kernel.  NEAR pairs (r <= near_factor*(size_a+size_b)) stay analytic-exact; FAR pairs use
a low-order degree-2 double-quadrature of 1/r (far_quad>0, O((size/r)^4)) that reproduces the all-analytic Gram
at ~monopole cost -- 3-6x faster build at N>=256, growing with N.  Compact/small geometries are ALL-NEAR at
near_factor=2, so the fast build is BIT-IDENTICAL to all-analytic there.  The fast FAR rule leaves the RAW Gram
~1e-5 asymmetric, but (a) the production solve applies G via the EXACTLY-symmetric H-matvec (mass_riesz
symmetric=true default) and (b) demag_factor is a quadratic form c.Gc, insensitive to the antisymmetric part.
Pass near_factor=1e30 to force all-analytic.
"""
import numpy as np
import scipy.sparse as sp
import radia._radia_pybind as _rp

# Hexahedron face -> vertex indices (the bottom 0-3 / top 4-7 corner order the C++ kernel also uses).  The C++
# kernel carries its own tet/pyramid/wedge/hex face tables; this constant is retained for reference / external use.
HEX_FACES = [[0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 5, 4], [3, 2, 6, 7], [0, 3, 7, 4], [1, 2, 6, 5]]

# vertex count -> element name (the supported polyhedra).
_ELEM_NAME = {4: "tetrahedron", 5: "pyramid", 6: "wedge", 8: "hexahedron"}


def assemble_moment_system(elements, *, quad=False, eps=1e-9, leaf=40, eta=0.5, near_factor=2.0, far_quad=4,
                           build=True):
    """Build the moment-Galerkin demag pieces for a MIXED list of polyhedra.

    elements : list of (nv, 3) vertex arrays; nv = 4 (tet) / 5 (pyramid) / 6 (wedge) / 8 (hex), radia corner order.
    quad : False -> 3 DOF/element (dipole); True -> 3 + (nFace-4) DOF/element (+ quad residual eigenmodes).
    eps, leaf, eta : ACA H-matrix parameters for the C++ charge-Gram (validated defaults).
    near_factor : analytic NEAR/FAR entry split (default 2 = fast; 1e30 = all-analytic).
    far_quad : FAR rule (default 4 = precision-preserving degree-2 double-quad; 0 = centroid-monopole).
    build : build the H-matrix now.

    Returns dict(G, B (csr n_charge x ndof_total), M_mass (csr ndof_total x ndof_total), vols, n_elem, n_charge,
                 ndof_total, ndof_elem (per-element DOF), dof_off (per-element global DOF start),
                 ndof_per (uniform DOF/element if all the same type, else None),
                 all_tris ((n_charge,3,3) triangle vertices in B-row order, for field reconstruction)).
    """
    elements = [np.asarray(V, float) for V in elements]
    n = len(elements)
    if n == 0:
        raise ValueError("moment_galerkin: empty element list")
    vcounts = []
    for V in elements:
        if V.ndim != 2 or V.shape[1] != 3 or V.shape[0] not in _ELEM_NAME:
            raise ValueError("moment_galerkin: each element must be an (nv,3) array with nv in {4,5,6,8} "
                             "(tet/pyramid/wedge/hex); got shape %s" % (V.shape,))
        vcounts.append(int(V.shape[0]))
    verts = np.concatenate([V.reshape(-1) for V in elements])
    # Per-element moment-basis assembly in C++ (ngcore::ParallelFor): face geometry, residual eigenmodes, second
    # moment -> B COO triplets, M_mass COO, boundary triangle soup.  This thin layer just packages the arrays.
    out = _rp._moment_galerkin_assemble(verts.tolist(), [int(v) for v in vcounts], bool(quad))
    n_charge = int(out["n_charge"]); ndof_total = int(out["ndof_total"])
    ndof_elem = np.asarray(out["ndof_elem"], int)
    dof_off = np.asarray(out["dof_off"], int)
    vols = np.asarray(out["vols"], float)
    if not np.all(np.isfinite(vols)) or np.any(vols <= 0.0):
        raise ValueError("moment_galerkin: degenerate / non-positive-volume element (check the vertex order)")
    # B and M_mass from the COO triplets (the C++ emits dense per-element blocks incl. structural zeros; filter).
    Br = np.asarray(out["B_rows"]); Bc = np.asarray(out["B_cols"]); Bd = np.asarray(out["B_data"], float)
    nz = Bd != 0.0
    B = sp.csr_matrix((Bd[nz], (Br[nz], Bc[nz])), shape=(n_charge, ndof_total))
    Mr = np.asarray(out["M_rows"]); Mc = np.asarray(out["M_cols"]); Md = np.asarray(out["M_data"], float)
    mz = Md != 0.0
    M_mass = sp.csr_matrix((Md[mz], (Mr[mz], Mc[mz])), shape=(ndof_total, ndof_total))
    face_tris = list(out["face_tris"])
    all_tris = np.asarray(face_tris, float).reshape(n_charge, 3, 3)
    ndof_per = int(ndof_elem[0]) if np.all(ndof_elem == ndof_elem[0]) else None
    G = _rp._ChargeGramHMatrix(cell_verts=[], face_verts=face_tris, n_el=0,
                               eps=eps, leaf=int(leaf), eta=float(eta), near_factor=float(near_factor),
                               far_quad=int(far_quad), build=bool(build))
    return {"G": G, "B": B, "M_mass": M_mass, "vols": vols, "n_elem": n, "n_hex": n, "n_charge": n_charge,
            "ndof_total": ndof_total, "ndof_elem": ndof_elem, "dof_off": dof_off, "ndof_per": ndof_per,
            "all_tris": all_tris}
