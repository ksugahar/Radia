"""radia.moment_galerkin._assemble -- moment-basis assembly for the SYMMETRIC moment-Galerkin MMMM demag.

The moment-Galerkin demag operator is N = B^T G B, the SAME symmetric structure HDiv-VIM uses, but on the
MMMM moment basis instead of HDiv's RT0 flux.  G is the EXACT analytic charge-Gram H-matrix
(`radia._radia_pybind._ChargeGramHMatrix`, the same C++ kernel HDiv-VIM ships); B maps the per-element moment
amplitudes to the face surface charges sigma; M_mass is the magnetization mass.  The physical SPD demag system

    ( (1/chi) M_mass + B^T G B ) m = M_mass H_ext

is solved by the existing C++ Krylov kernel.  The per-hex moment-basis assembly (B, M_mass, the boundary
triangle soup) is the C++ `RadMomentGalerkinAssembleHex` parallel kernel (ngcore::ParallelFor, exposed as
`_rp._moment_galerkin_assemble_hex`); this file is the thin Python packaging layer (COO -> csr, block-diag).

Two moment orders (`quad=`):
  quad=False (default) -- 3 DOF/hex: a CONSTANT magnetization M_0 per hex (the dipole / standard demag).
      B col k = sigma_f = M_0 . n_f = n_f[k];  M_mass = diag(V_e).
  quad=True            -- 5 DOF/hex: 3 dipole + 2 QUAD residual-eigenmode amplitudes (the geometry-adaptive
      quadrupole patterns phi_q = the zero-monopole, zero-dipole-moment face-charge modes, the same
      Gram-Schmidt basis the C++ collocation MMMM `momentResidualEigenmodes` uses).  The quad mode's
      magnetization is the LINEAR field M_q(x) = Dm_q.(x-c), Dm_q = sum_f Ae_f n_f (x) d_f * (phi_q/Ae_f),
      so the quad MASS is the magnetization Gram W_quad[q,r] = Dm_q : M2 : Dm_r, M2 = INT (x-c)(x-c) dV.
      The 2 quad DOF buy per-element accuracy under SKEW / gradient loads (the iron-yoke case); they are
      ~null for axial / symmetric loads (correct physics, never worse).  Validated: a 3-hex bar 5-DOF field
      is strictly closer to a fine-mesh rad.Solve(yano) truth than dipole-only at every near probe under an
      oblique field, ties for axial; the cube quad amplitudes are ~0 by symmetry (the dipole 1/3 demag is
      preserved).  NOTE: on a cube / near-cube the quad mass is near-singular (the residual modes carry ~zero
      second moment) -- benign because the cube also nulls the quad amplitudes, and the diagonal-Jacobi C++
      solve (auto_prec) stays well-defined (its diagonal includes the demag self-energy N_diag).

Symmetry is by construction: G is symmetric (Coulomb reciprocity) -> N = B^T G B symmetric to machine
precision (validated 1e-16 on the loop-heavy C-yoke vs collocation MMMM's 1.56) -> loop modes field-null ->
mu_r-independent / loop-free convergence with no loop-star.

ACA defaults (leaf=40, eta=0.5): tuned so the H-matvec reproduces the exact analytic Gram to ~1e-4 on
near-heavy geometry (eta is the dominant knob; the MMMM 12-face-triangle charge group must not be split, so
leaf~40 keeps ~3 hexes intact).

FAST BUILD defaults (near_factor=2, far_quad=4): the SAME precision-preserving NEAR/FAR split HDiv-VIM ships on
this shared charge-Gram kernel.  NEAR pairs (r <= near_factor*(size_a+size_b)) stay analytic-exact; FAR pairs use
a low-order degree-2 double-quadrature of 1/r (far_quad>0, O((size/r)^4)) that reproduces the all-analytic Gram
at ~monopole cost -- 3-6x faster build at N>=256, growing with N.  Compact/small geometries (cube, few-hex bars)
are ALL-NEAR at near_factor=2, so the fast build is BIT-IDENTICAL to all-analytic there (the goldens are
unchanged).  The fast FAR rule leaves the RAW Gram entries ~1e-5 asymmetric, but (a) the production solve applies
G via the EXACTLY-symmetric H-matvec (mass_riesz symmetric=true default) and (b) demag_factor is a quadratic form
c.Gc, insensitive to the antisymmetric part -- so the moment-Galerkin symmetry/accuracy guarantees are preserved.
Pass near_factor=1e30 to force all-analytic (e.g. to assert the formulation's exact symmetry on a large mesh).
"""
import numpy as np
import scipy.sparse as sp
import radia._radia_pybind as _rp

# Hexahedron face -> vertex indices (the bottom 0-3 / top 4-7 corner order the C++ kernel also uses).
HEX_FACES = [[0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 5, 4], [3, 2, 6, 7], [0, 3, 7, 4], [1, 2, 6, 5]]


def assemble_moment_system(hexes, *, quad=False, eps=1e-9, leaf=40, eta=0.5, near_factor=2.0, far_quad=4,
                           build=True):
    """Build the moment-Galerkin demag pieces for a list of hexahedra.

    quad : False -> 3 DOF/hex (dipole); True -> 5 DOF/hex (3 dipole + 2 quad residual eigenmodes).
    eps, leaf, eta : ACA H-matrix parameters for the C++ charge-Gram (validated defaults).
    near_factor : analytic NEAR/FAR entry split (default 2 = fast; 1e30 = all-analytic).
    far_quad : FAR rule (default 4 = precision-preserving degree-2 double-quad; 0 = centroid-monopole).
    build : build the H-matrix now.

    Returns dict(G, B (csr n_charge x ndof_per*n_hex), M_mass (csr), vols, n_hex, n_charge, ndof_per,
                 all_tris ((n_charge,3,3) triangle vertex array in B-row order, for field reconstruction)).
    """
    hexes = [np.asarray(V, float) for V in hexes]
    n = len(hexes)
    if n == 0:
        raise ValueError("moment_galerkin: empty hex list")
    hexverts = np.ascontiguousarray(np.array(hexes, float)).reshape(-1)
    if hexverts.size != 24 * n:
        raise ValueError("moment_galerkin: each hex needs 8 vertices x 3 coordinates")
    # Per-hex moment-basis assembly in C++ (ngcore::ParallelFor): face geometry, residual eigenmodes, second
    # moment -> B COO triplets, M_mass blocks, boundary triangle soup.  This thin layer just packages the arrays.
    out = _rp._moment_galerkin_assemble_hex(hexverts.tolist(), int(n), bool(quad))
    ndof_per = int(out["ndof_per"]); n_charge = int(out["n_charge"])
    vols = np.asarray(out["vols"], float)
    if not np.all(np.isfinite(vols)) or np.any(vols <= 0.0):
        raise ValueError("moment_galerkin: degenerate / non-positive-volume hex (check the 8-vertex corner order)")
    # B from the COO triplets (the C++ emits all 12*ndof_per slots/hex incl. structural zeros; filter to recover
    # the exact sparsity the dense path expects).
    Br = np.asarray(out["B_rows"]); Bc = np.asarray(out["B_cols"]); Bd = np.asarray(out["B_data"], float)
    nz = Bd != 0.0
    B = sp.csr_matrix((Bd[nz], (Br[nz], Bc[nz])), shape=(n_charge, ndof_per * n))
    # M_mass: dipole -> diag(V_e); quad -> block-diag of the 5x5 magnetization-mass blocks (vectorized).
    if quad:
        Wb = np.asarray(out["Wblocks"], float).reshape(n, 5, 5)
        ii, jj = np.meshgrid(np.arange(5), np.arange(5), indexing="ij")
        base = np.arange(n)[:, None, None] * 5
        Mr = (base + ii[None]).ravel(); Mc = (base + jj[None]).ravel()
        M_mass = sp.csr_matrix((Wb.ravel(), (Mr, Mc)), shape=(5 * n, 5 * n))
    else:
        M_mass = sp.diags(np.repeat(vols, 3)).tocsr()
    face_tris = list(out["face_tris"])
    all_tris = np.asarray(face_tris, float).reshape(n_charge, 3, 3)
    G = _rp._ChargeGramHMatrix(cell_verts=[], face_verts=face_tris, n_el=0,
                               eps=eps, leaf=int(leaf), eta=float(eta), near_factor=float(near_factor),
                               far_quad=int(far_quad), build=bool(build))
    return {"G": G, "B": B, "M_mass": M_mass, "vols": vols, "n_hex": n, "n_charge": n_charge,
            "ndof_per": ndof_per, "all_tris": all_tris}
