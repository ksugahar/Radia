"""Shared compute helpers for hex_vs_mmmm_crossvalidation.ipynb -- HDiv-VIM hex (RT1 charge-Gram) vs
collocation-MMMM hex, on the SAME structured hex mesh.

Two INDEPENDENT hex soft-iron demag backends are cross-validated:
  * collocation MMMM (six-face surface charge) via the public rad.Solve(demag_backend='collocation_mmmm')
    on the ObjHexahedron iron built by radia.vim.MeshSoftIron.
  * HDiv-VIM RT1: public radia.vim.Solve and rad.Solve(auto) now accept pure hex.  This historical
    helper still drives the same wired hex Gram with the shipped production mass-Riesz CG directly so the
    notebook remains byte-stable with its 2026-07-04 result sidecar.

Per CLAUDE.md "TaskManager Wrap Policy: Caller Wraps, Helper Does NOT" -- these helpers open NO TaskManager;
the notebook wraps each call in `with ngsolve.TaskManager():`.

The GetTrafo first-touch flake (memory ngsolve-gettrafo-first-touch-garbage) is a KNOWN bursty NGSolve issue
in the hex-charge basis extraction; the determinism contract fail-LOUD raises rather than returning garbage.
hdiv_hex_solve retries the build a few times (each retry re-draws the extraction) so a showcase run is robust.
"""
import numpy as np
import scipy.sparse as sp
import ngsolve as ng
from netgen.meshing import Mesh as _NGMesh, MeshPoint, Element3D, Element2D, FaceDescriptor
from netgen.csg import Pnt

MU0 = 4.0e-7 * np.pi

# standard hex node ordering: bottom quad (z-) CCW then top quad (z+) CCW; 6 faces (outward) for bnd quads
_HEX = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
_FACES = {(-1, 0, 0): (0, 4, 7, 3), (1, 0, 0): (1, 2, 6, 5), (0, -1, 0): (0, 1, 5, 4),
          (0, 1, 0): (3, 7, 6, 2), (0, 0, -1): (0, 3, 2, 1), (0, 0, 1): (4, 5, 6, 7)}


def build_voxel_hex_mesh(iron_pred, x0, y0, z0, h, nx, ny, nz):
    """Genuine structured-hex NGSolve mesh of a voxel predicate on cell CENTERS.  Verified: Integrate(1)
    == n_iron*h^3, positive Jacobians, HDiv(order=1) + MeshSoftIron + ChargeGram all accept it."""
    keep = np.zeros((nx, ny, nz), bool)
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                keep[i, j, k] = iron_pred(x0 + (i + 0.5) * h, y0 + (j + 0.5) * h, z0 + (k + 0.5) * h)
    m = _NGMesh(dim=3)
    m.SetMaterial(1, "iron")
    fd = m.Add(FaceDescriptor(bc=1, domin=1, surfnr=1))
    pids = {}

    def pid(i, j, k):
        key = (i, j, k)
        if key not in pids:
            pids[key] = m.Add(MeshPoint(Pnt(x0 + i * h, y0 + j * h, z0 + k * h)))
        return pids[key]

    n_hex = 0
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if not keep[i, j, k]:
                    continue
                corner = [pid(i + dx, j + dy, k + dz) for (dx, dy, dz) in _HEX]
                m.Add(Element3D(1, corner))
                n_hex += 1
                for (nrm, loc) in _FACES.items():
                    ni, nj, nk = i + nrm[0], j + nrm[1], k + nrm[2]
                    if not (0 <= ni < nx and 0 <= nj < ny and 0 <= nk < nz and keep[ni, nj, nk]):
                        m.Add(Element2D(fd, [corner[t] for t in loc]))
    return ng.Mesh(m), n_hex


def cube_mesh(n, L=0.02):
    """n x n x n structured hex cube, edge L, centred at the origin."""
    from ngsolve.meshes import MakeStructured3DMesh
    return MakeStructured3DMesh(hexes=True, nx=n, ny=n, nz=n,
                                mapping=lambda x, y, z: (L * (x - 0.5), L * (y - 0.5), L * (z - 0.5)))


def cyoke_pred(cx, cy, cz):
    """The OCC cyoke() as a voxel predicate: outer box minus inner cavity minus the +x gap opening (a
    genuine non-convex, reentrant-corner C-yoke; the same shape as the retired
    C-yoke nonlinear prototype inventoried in vim_examples_retirement_results.json)."""
    if not (abs(cx) <= 0.06 and abs(cy) <= 0.06 and abs(cz) <= 0.02):
        return False
    if abs(cx) <= 0.035 and abs(cy) <= 0.035:
        return False
    return cx < 0.018


def cyoke_mesh(h):
    """Voxelized hex C-yoke at spacing h over the [-0.06,0.06]^2 x [-0.02,0.02] bounding box."""
    x0, y0, z0 = -0.06, -0.06, -0.02
    nx = int(round(0.12 / h)); ny = int(round(0.12 / h)); nz = int(round(0.04 / h))
    return build_voxel_hex_mesh(cyoke_pred, x0, y0, z0, h, nx, ny, nz)


def hdiv_hex_solve(mesh, mu_r, H_vec, retries=4):
    """HDiv-VIM hex RT1: wired hex charge Gram (ChargeGram) + the SHIPPED production symmetric
    mass-Riesz CG (_solve_linear_mass_riesz_cpp).  Retries the build on the known GetTrafo first-touch
    flake (fail-loud in the determinism contract).  Returns per-element M, volume-avg M, demag factor."""
    from radia.vim import ChargeGram
    from radia.vim._solve import _solve_linear_mass_riesz_cpp
    n_el = mesh.GetNE(ng.VOL)
    last = None
    for _ in range(max(1, retries)):
        try:
            fes = ng.HDiv(mesh, order=1)
            B, G, M_mass = ChargeGram(fes)          # hex auto-branch (Q1 vol charge + Q2 geometry)
            break
        except RuntimeError as e:
            if "GetTrafo lattice evaluation unstable" in str(e):
                last = e
                continue
            raise
    else:
        raise last
    B = sp.csr_matrix(B); Mm = sp.csr_matrix(M_mass)
    n_face = fes.ndof
    gfH = ng.GridFunction(fes); gfH.Set(ng.CoefficientFunction(tuple(H_vec)))
    h_ext = gfH.vec.FV().NumPy().copy()
    gfMu = ng.GridFunction(fes)
    unit = tuple(1.0 if abs(v) == max(abs(np.array(H_vec))) and v != 0 else 0.0 for v in H_vec)
    gfMu.Set(ng.CoefficientFunction(unit))
    mu = gfMu.vec.FV().NumPy().copy()
    chi = float(mu_r) - 1.0
    m, iters = _solve_linear_mass_riesz_cpp(G, B, Mm, n_face, h_ext, chi, 1e-10, 6000)
    gfM = ng.GridFunction(fes); gfM.vec.FV().NumPy()[:] = m
    gfMc = ng.GridFunction(ng.VectorL2(mesh, order=0)); gfMc.Set(gfM)
    M_el = gfMc.vec.FV().NumPy().reshape(3, n_el).T.copy()
    vol = ng.Integrate(ng.CoefficientFunction(1.0), mesh)
    M_avg = np.array([ng.Integrate(gfM[i], mesh) for i in range(3)]) / vol
    Nmu = B.T @ np.asarray(G.matvec_sym((B @ mu).tolist()), float)
    demag = float((mu @ Nmu) / (mu @ np.asarray(Mm @ mu).ravel()))
    return dict(M=M_el, M_avg=M_avg, demag=demag, iters=int(iters), n_el=n_el,
                n_charge=int(B.shape[0]), ndof=n_face)


def mmmm_hex_solve(mesh, mu_r, H_vec, probe):
    """collocation MMMM hex: MeshSoftIron -> ObjHexahedron, rad.Solve(collocation_mmmm), then the
    iron reaction B at probe points via rad.Fld.  Applied field via ObjBckg B = mu0*H_vec."""
    import radia as rad
    from radia.vim import MeshSoftIron
    rad.UtiDelAll()
    core = MeshSoftIron(mesh, mu_r=float(mu_r))
    bkg = rad.ObjBckg(lambda p: [MU0 * H_vec[0], MU0 * H_vec[1], MU0 * H_vec[2]])
    rad.Solve(rad.ObjCnt([core, bkg]), 1e-6, 3000, 0, demag_backend="collocation_mmmm")   # LU
    M_el = np.array([m for (_c, m) in rad.ObjM(core)], float)
    B_iron = np.array(rad.Fld(core, 'b', list(probe))).reshape(-1, 3)
    return dict(M_avg=M_el.mean(0), B_iron=B_iron, n_el=len(M_el))


def iron_external_field(mesh, M_el, probe):
    """Rebuild a Radia hex iron from a per-element M and evaluate its external B (SAME field kernel)."""
    import radia as rad
    from radia import netgen_mesh_import as nmi
    rad.UtiDelAll()
    iron = nmi.netgen_mesh_to_radia(
        mesh, material=lambda i: {'magnetization': M_el[i].tolist()}, units='m',
        allow_hex=True, verbose=False)
    return np.array(rad.Fld(iron, 'b', list(probe))).reshape(-1, 3)


def agreement(hd_Mavg, mm_Mavg, hd_B, mm_B, axis):
    """Order-independent agreement metrics: rel diff of the axis volume-avg M, and the per-probe vector B diff."""
    dM = float(abs(hd_Mavg[axis] - mm_Mavg[axis]) / (abs(mm_Mavg[axis]) + 1e-30))
    num = np.linalg.norm(np.asarray(hd_B) - np.asarray(mm_B), axis=1)
    den = np.linalg.norm(np.asarray(mm_B), axis=1) + 1e-30
    dB = (num / den)
    return dM, dB.tolist(), float(dB.max())
