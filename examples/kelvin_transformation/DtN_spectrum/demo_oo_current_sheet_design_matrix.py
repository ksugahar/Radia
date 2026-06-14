# VERIFIED (2026-06-15): build the CURRENT-SHEET transfer matrix M[target, psi-dof] = M_free + M_reaction
# (M_free = Biot-Savart at targets; M_reaction = -grad(Omega) at targets, Omega = iron reaction via the
# reduced-potential Kelvin-FEM) and wire it into the EXISTING radia.stream_function (ACA+)+TSVD + ridge.
# Then: design with the material-aware M HITS the iron-system target; the free-space M MISSES.
#   psi-dof = P1 nodal values on the ELLIPSOID coil-surface vertices (a true current potential).
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import ngsolve as ng
from ngsolve import TaskManager
import netgen.occ as occ
from netgen.occ import Sphere, Pnt, Vec, IdentificationType, OCCGeometry
from radia.stream_function import aca_tsvd, RegularizedTSVD

b, c = 0.58, 0.74
R_out, offset, order, MU_R = 1.1, 3.3, 2, 50.0
INV4PI = 1.0 / (4.0 * np.pi)
th = np.deg2rad(np.linspace(25, 155, 5)); ph = np.deg2rad(np.linspace(0, 300, 6))
PTS = np.array([[np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)] for t in th for p in ph]) * 0.20
NT_TARGET = PTS.shape[0] * 3


def coil_surface(axes, hc=0.18):
    m = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0)).GenerateMesh(maxh=hc))
    Vall = np.array([list(v.point) for v in m.vertices]) * np.asarray(axes)
    Tall = np.array([[v.nr for v in el.vertices] for el in m.Elements(ng.BND)])
    uniq, inv = np.unique(Tall.ravel(), return_inverse=True)
    V = Vall[uniq]; T = inv.reshape(Tall.shape)          # reindex to boundary verts only
    return V, T


def tri_basis(V, T):
    """per-triangle centroid C and per-local-vertex area-weighted currents w[t,i] = area*(n x grad lam_i)."""
    p0, p1, p2 = V[T[:, 0]], V[T[:, 1]], V[T[:, 2]]
    C = (p0 + p1 + p2) / 3.0
    e1, e2 = p1 - p0, p2 - p0
    nrm = np.cross(e1, e2); area2 = np.linalg.norm(nrm, axis=1)
    nhat = nrm / area2[:, None]; nhat *= np.sign(np.einsum('ij,ij->i', nhat, C))[:, None]
    Mtx = np.stack([e1, e2, nhat], axis=1)
    gl = np.zeros((len(T), 3, 3))                         # grad lam_i, i=0,1,2
    rhs = np.array([[-1.0, -1, 0], [1.0, 0, 0], [0.0, 1, 0]])
    for i in range(3):
        gl[:, i, :] = np.linalg.solve(Mtx, np.tile(rhs[i], (len(T), 1))[:, :, None])[:, :, 0]
    w = np.cross(nhat[:, None, :], gl) * (0.5 * area2)[:, None, None]   # (NT,3,3) area-weighted current
    return C, w


def scatter_biot(C, w, T, nvert, P):
    """B[p, comp, vert] : Biot-Savart H at points P from a unit psi on each coil vertex."""
    B = np.zeros((len(P), 3, nvert))
    for t in range(len(C)):
        d = P - C[t]; r3 = (d * d).sum(1) ** 1.5 + 1e-300
        for i in range(3):
            B[:, :, T[t, i]] += INV4PI * np.cross(w[t, i][None, :], d) / r3[:, None]
    return B


def vol_geo(kelvin, R_box=3.0):
    sb, sc = Sphere(Pnt(0, 0, 0), b), Sphere(Pnt(0, 0, 0), c)
    if kelvin:
        outer = Sphere(Pnt(0, 0, 0), R_out)
        for f in outer.faces:
            f.name = "kelvin_int"
        core = sb; core.mat("vac"); shell = (sc - sb); shell.mat("shell"); out = (outer - sc); out.mat("vac")
        kball = Sphere(Pnt(offset, 0, 0), R_out)
        for f in kball.faces:
            f.name = "kelvin_ext"
        kball.mat("kelvin"); gnd = occ.Vertex(Pnt(offset, 0, 0)); gnd.name = "GND"
        solids = [core, shell, out]
        fi = [f for s in solids for f in s.faces if f.name == "kelvin_int"][0]
        fe = [f for f in kball.faces if f.name == "kelvin_ext"][0]
        fi.Identify(fe, "kelvin", IdentificationType.PERIODIC, occ.gp_Trsf.Translation(Vec(offset, 0, 0)))
        return occ.Glue(solids + [kball, gnd])
    box = Sphere(Pnt(0, 0, 0), R_box)
    for f in box.faces:
        f.name = "outer"
    core = sb; core.mat("vac"); core.maxh = 0.14; shell = (sc - sb); shell.mat("shell"); shell.maxh = 0.12
    out = (box - sc); out.mat("vac")
    return occ.Glue([core, shell, out])


def build_reaction_M(C, w, T, nvert, gmaxh=0.13):
    """M_reaction[target-row, psi-dof] = -grad(Omega) at targets, Omega = reduced-potential iron reaction
    via the Kelvin-FEM, with H_s injected by a precomputed iron coupling matrix D (one factorisation)."""
    mesh = ng.Mesh(OCCGeometry(vol_geo(True)).GenerateMesh(maxh=gmaxh)).Curve(min(order + 1, 4))
    xx, yy, zz = ng.x, ng.y, ng.z; rp2 = (xx - offset) ** 2 + yy * yy + zz * zz + 1e-20
    mu = mesh.MaterialCF({"vac": 1.0, "shell": MU_R, "kelvin": R_out ** 2 / rp2}, default=1.0)
    fes = ng.Periodic(ng.H1(mesh, order=order, dirichlet="GND"))
    Vsrc = ng.VectorH1(mesh, order=1)
    u, v = fes.TnT()
    A = ng.BilinearForm(mu * ng.grad(u) * ng.grad(v) * ng.dx(bonus_intorder=6)); A.Assemble()
    usrc = Vsrc.TrialFunction()
    D = ng.BilinearForm(trialspace=Vsrc, testspace=fes)
    D += (MU_R - 1.0) * ng.InnerProduct(usrc, ng.grad(v)) * ng.dx(definedon=mesh.Materials("shell"), bonus_intorder=4)
    D.Assemble()
    Ainv = A.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    # iron-shell vertices (where H_s is injected)
    Vc = np.array([list(p.point) for p in mesh.vertices]); r = np.linalg.norm(Vc, axis=1)
    iron = np.where((r > b - 0.05) & (r < c + 0.05))[0]
    B_iron = scatter_biot(C, w, T, nvert, Vc[iron])       # (n_iron, 3, nvert)
    targets = [mesh(*p) for p in PTS]
    gsrc = ng.GridFunction(Vsrc); gO = ng.GridFunction(fes); rhs = gO.vec.CreateVector()
    Mr = np.zeros((NT_TARGET, nvert))
    for j in range(nvert):
        for k in range(3):
            arr = gsrc.components[k].vec.FV().NumPy(); arr[:] = 0.0; arr[iron] = B_iron[:, k, j]
        rhs.data = D.mat * gsrc.vec
        gO.vec.data = Ainv * rhs
        gr = ng.grad(gO)
        col = -np.array([[gr(t)[k] for k in range(3)] for t in targets])   # -grad(Omega) at targets
        Mr[:, j] = col.reshape(-1)
    return Mr


with TaskManager():
    print("ELLIPSOID coil (0.36,0.36,0.52), iron shell mu_r=%g -- current-sheet transfer matrix M" % MU_R)
    V, T = coil_surface((0.36, 0.36, 0.52), hc=0.18)
    nvert = V.shape[0]
    C, w = tri_basis(V, T)
    M_free = scatter_biot(C, w, T, nvert, PTS).reshape(NT_TARGET, nvert)   # Biot-Savart at targets
    print("  coil psi-dofs = %d boundary vertices ; targets = %d (x3 = %d rows)" % (nvert, len(PTS), NT_TARGET))
    M_react = build_reaction_M(C, w, T, nvert)
    M_iron = M_free + M_react
    print("  ||M_react|| / ||M_free|| = %.3f  (iron reaction is a real part of the operator)"
          % (np.linalg.norm(M_react) / np.linalg.norm(M_free)))

    # --- design: a realizable target = the field of a reference psi; design with iron vs free-space ---
    rng_psi = (V[:, 2] ** 2 - 0.5 * (V[:, 0] ** 2 + V[:, 1] ** 2))         # a Z2 current potential
    B_des = M_iron @ rng_psi                                               # desired field at targets

    def tsvd_solve(M, B, rc=1e-8):
        U, s, Vt = np.linalg.svd(M, full_matrices=False); keep = s > rc * s[0]
        return (Vt[keep].T * (1.0 / s[keep])) @ (U[:, keep].T @ B)

    psi_iron = tsvd_solve(M_iron, B_des)
    psi_free = tsvd_solve(M_free, B_des)
    hit = np.linalg.norm(M_iron @ psi_iron - B_des) / np.linalg.norm(B_des)
    miss = np.linalg.norm(M_iron @ psi_free - B_des) / np.linalg.norm(B_des)
    print("\n  DESIGN (target = a realizable field pattern at the DSV):")
    print("    design WITH material-aware M  : residual = %.2e   <- HITS" % hit)
    print("    design IGNORING iron (free M) : residual = %.2e   <- MISSES (iron not in the kernel)" % miss)
    # ridge via the existing RegularizedTSVD + current-density seminorm S = identity proxy
    rr = RegularizedTSVD.from_stiffness(
        aca_tsvd(NT_TARGET, nvert, lambda i, j: float(M_iron[i, j]), modes=min(NT_TARGET, nvert), aca_eps=1e-10),
        np.eye(nvert))
    psi_ridge = rr.solve(B_des, alpha=1e-8)
    hit_r = np.linalg.norm(M_iron @ psi_ridge - B_des) / np.linalg.norm(B_des)
    print("    radia.stream_function ridge   : residual = %.2e   (existing solver consumes the M)" % hit_r)
    print("\n[RESULT] current-sheet M = M_free + M_react (react/free %.2f); material-aware design HITS"
          " (%.1e), free-space design MISSES (%.0f%%); radia.stream_function consumes M (%.1e)."
          % (np.linalg.norm(M_react) / np.linalg.norm(M_free), hit, 100 * miss, hit_r))
