"""Acoustic fluid-structure interaction (FSI), DtN exterior.

A solid ELASTIC scatterer (vector NGSolve VectorH1(order=p) elasticity FEM in the
interior) radiating into an unbounded fluid, closed by the EXACT spherical
Helmholtz Dirichlet-to-Neumann operator ("high-order Zs" -- a radiating impedance,
NOT a Kelvin boundary and NOT a PML).  Plane wave exp(ikz); fluid c=1, rho=1;
speeds/densities are ratios to the fluid.  Ported from the matlab-acoustic-fembem
fsiCoupledSolve / sphericalDtnOperator; the P1 interior elasticity is delegated to
NGSolve (Complement NGSolve) so p=1,2,... come from one implementation.

DtN reduced block system (u interior displacement, c = (N+1)^2 harmonic coeffs):
  [ Kdyn                    G' Phi          ] [u]   [ -G' p_incGamma ]
  [ -rhoF w^2 (Phi' G)      Gram .* lam'    ] [c] = [ -Phi' Minc     ]
Kdyn = K - w^2 M (elasticity), G = int_Gamma q (u.n) ds (mixed coupling), Phi =
real spherical harmonics at the boundary vertices, lam_n = k h_n'(kR)/h_n(kR).
Exterior scattered field straight from c: p_s(r,dir) = sum c_col h_n(kr)/h_n(kR)
Y_n^m(dir) -- no dense layer-potential assembly.

Validated (validation_test/acoustics): converges to the analytic Faran elastic
sphere (radia.acoustics.elastic_sphere_scattering) under refinement; P2 converges
~O(h^2) (much faster than P1), e.g. kR=2 vs-Faran 15%->8.6%->3.7% at maxh
0.4/0.28/0.20 (P1 is ~25% at maxh 0.20). The stiff limit reproduces the rigid
sphere (formulation gate).
"""
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from scipy.special import spherical_jn, spherical_yn, lpmv, gammaln


# ===== spherical Hankel + real spherical harmonics + DtN (numpy) ============== #
def _sph_h1(l, x):
    return spherical_jn(l, x) + 1j * spherical_yn(l, x)


def _sph_h1_d(l, x):
    return spherical_jn(l, x, derivative=True) + 1j * spherical_yn(l, x, derivative=True)


def real_spherical_harmonics(dirs, N):
    """Orthonormal real spherical harmonics Y_n^m at unit directions (nP x (N+1)^2).

    Columns: n = 0..N, within n the order m = 0 then (+m cos, -m sin); degree_of
    tags each column with its n so it can carry the DtN eigenvalue Lambda_n.
    """
    npn = len(dirs)
    ct = dirs[:, 2]
    ph = np.arctan2(dirs[:, 1], dirs[:, 0])
    nModes = (N + 1) ** 2
    Phi = np.zeros((npn, nModes))
    degree_of = np.zeros(nModes, int)
    col = 0
    for n in range(N + 1):
        for m in range(n + 1):
            cnm = np.sqrt((2*n + 1) / (4*np.pi)) * np.exp(0.5 * (gammaln(n - m + 1) - gammaln(n + m + 1)))
            Pm = lpmv(m, n, ct)                     # P_n^m(ct), Condon-Shortley phase
            if m == 0:
                Phi[:, col] = cnm * Pm; degree_of[col] = n; col += 1
            else:
                s2 = np.sqrt(2) * cnm
                Phi[:, col] = s2 * Pm * np.cos(m * ph); degree_of[col] = n
                Phi[:, col + 1] = s2 * Pm * np.sin(m * ph); degree_of[col + 1] = n
                col += 2
    return Phi, degree_of


def spherical_dtn(bnd_vtx, k, Mb, degree=-1, tol=3e-2):
    """Exact spherical Helmholtz DtN on the boundary vertices (fail-loud off-sphere)."""
    nB = len(bnd_vtx)
    Afit = np.hstack([2 * bnd_vtx, np.ones((nB, 1))])
    pfit, *_ = np.linalg.lstsq(Afit, (bnd_vtx**2).sum(1), rcond=None)
    center = pfit[:3]
    R = np.sqrt(pfit[3] + center @ center)
    rel = bnd_vtx - center
    radii = np.linalg.norm(rel, axis=1)
    deviation = np.max(np.abs(radii - R)) / R
    if deviation > tol:
        raise ValueError(f"DtN needs a spherical truncation; radius deviation "
                         f"{deviation:.3g} exceeds tol {tol}")
    N = degree
    if N < 0:
        N = int(np.ceil(k * R)) + 8
    if (N + 1)**2 > nB:
        N = min(N, max(1, int(np.floor(np.sqrt(nB))) - 1))
    dirs = rel / radii[:, None]
    Phi, degree_of = real_spherical_harmonics(dirs, N)
    x0 = k * R
    Lambda = np.array([k * _sph_h1_d(n, x0) / _sph_h1(n, x0) for n in range(N + 1)])
    lam_col = Lambda[degree_of]
    MbPhi = Mb @ Phi
    Gram = Phi.T @ MbPhi
    return {"Phi": Phi, "lam_col": lam_col, "Gram": Gram, "degree_of": degree_of,
            "center": center, "radius": R, "degree": N, "num_modes": (N + 1)**2,
            "Lambda": Lambda, "deviation": deviation}


def _exterior_from_coeffs(c, dtn, k, obs):
    center, R = dtn["center"], dtn["radius"]
    rel = np.asarray(obs, float) - center
    r = np.linalg.norm(rel, axis=1)
    dirs = rel / r[:, None]
    N = dtn["degree"]
    Phi_obs, degree_of = real_spherical_harmonics(dirs, N)
    hkR = {n: _sph_h1(n, k * R) for n in range(N + 1)}
    ratio = np.empty((len(obs), dtn["num_modes"]), complex)
    for col in range(dtn["num_modes"]):
        ratio[:, col] = _sph_h1(degree_of[col], k * r) / hkR[degree_of[col]]
    return (Phi_obs * ratio) @ c


# ===== NGSolve interior elasticity + coupling ================================= #
def _ng_to_scipy(mat):
    rows, cols, vals = mat.COO()
    return sp.csr_matrix((np.array(vals, dtype=complex), (np.array(rows), np.array(cols))),
                         shape=(mat.height, mat.width))


def sphere_mesh(radius=1.0, maxh=0.4):
    """A sphere volume tet mesh (boundary face named 'gamma') for the FSI solve."""
    from netgen.occ import Sphere, OCCGeometry
    from ngsolve import Mesh
    solid = Sphere((0, 0, 0), radius)
    solid.faces.name = "gamma"
    return Mesh(OCCGeometry(solid).GenerateMesh(maxh=maxh))


def fsi_dtn_solve(mesh, k, cL=2.0, cT=1.0, rho_s=1.5, rho_f=1.0,
                  order=1, boundary="gamma", obs=None):
    """FSI coupled solve with the spherical-DtN exterior.

    mesh: an NGSolve Mesh whose truncation boundary (``boundary``) is a sphere.
    Returns a dict with the DtN harmonic coefficients ``c``, solve ``residual``,
    the ``dtn`` operator info, and (if ``obs`` given) ``scattered`` / ``incident``
    / ``total`` complex pressures at the exterior observation points.
    """
    from ngsolve import (VectorH1, H1, BilinearForm, LinearForm, TaskManager, BND,
                         InnerProduct, Sym, Grad, div, ds, dx, specialcf, exp, z as Zc)

    omega = k
    mu = rho_s * cT**2
    lam = rho_s * (cL**2 - 2 * cT**2)

    with TaskManager():
        fes_u = VectorH1(mesh, order=order, complex=True)
        u, v = fes_u.TnT()
        a = BilinearForm(fes_u, symmetric=True)
        a += (2 * mu * InnerProduct(Sym(Grad(u)), Sym(Grad(v))) + lam * div(u) * div(v)
              - omega**2 * rho_s * InnerProduct(u, v)) * dx
        a.Assemble()

        fes_p = H1(mesh, order=1, complex=True)
        p, q = fes_p.TnT()
        nrm = specialcf.normal(mesh.dim)
        mb = BilinearForm(fes_p, symmetric=True)
        mb += p * q * ds(boundary)
        mb.Assemble()
        gform = BilinearForm(trialspace=fes_u, testspace=fes_p, check_unused=False)
        gform += InnerProduct(u, nrm) * q * ds(boundary)
        gform.Assemble()
        minc = LinearForm(fes_p)
        minc += (1j * k * exp(1j * k * Zc)) * nrm[2] * q * ds(boundary)
        minc.Assemble()

    Kdyn = _ng_to_scipy(a.mat)
    Mb_full = _ng_to_scipy(mb.mat)
    G_full = _ng_to_scipy(gform.mat)
    Minc_full = np.array(minc.vec, dtype=complex)

    bnd = np.array(sorted({vx.nr for el in mesh.Elements(BND) for vx in el.vertices}), int)
    coords = np.array([list(mesh.vertices[i].point) for i in bnd])
    Mb = Mb_full[np.ix_(bnd, bnd)]
    G = G_full[bnd, :]                                   # nB x ndof_u
    Minc = Minc_full[bnd]
    pincB = np.exp(1j * k * coords[:, 2])

    dtn = spherical_dtn(coords, k, Mb)
    Phi, lam_col = dtn["Phi"], dtn["lam_col"]
    Gram = dtn["Gram"]
    ndof_u = Kdyn.shape[0]

    A12 = sp.csr_matrix(G.T @ Phi)
    A21 = sp.csr_matrix(-rho_f * omega**2 * (Phi.T @ G))
    A22 = sp.csr_matrix(Gram * lam_col[None, :])
    Amat = sp.bmat([[Kdyn, A12], [A21, A22]]).tocsc()
    rhs = np.concatenate([-(G.T @ pincB), -(Phi.T @ Minc)])
    x = spsolve(Amat, rhs)
    c = x[ndof_u:]
    residual = float(np.linalg.norm(Amat @ x - rhs))

    out = {"c": c, "residual": residual, "dtn": dtn, "order": int(order),
           "ndof_u": int(ndof_u), "num_boundary_nodes": int(len(bnd)), "wavenumber": k}
    if obs is not None:
        obs = np.asarray(obs, float)
        out["scattered"] = _exterior_from_coeffs(c, dtn, k, obs)
        out["incident"] = np.exp(1j * k * obs[:, 2])
        out["total"] = out["incident"] + out["scattered"]
    return out
