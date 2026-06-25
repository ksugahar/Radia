"""Step-1 verification: the C++ BuildMomentSystem == the Python moment reference, entry-for-entry.

Builds the same C-yoke hexes, assembles the moment system A,rhs (a) in Python (the prototype build logic
using GetFaceGeom + GetCentroidFieldGrad) and (b) in C++ via rad.BuildMomentSystem, and compares.  Pass =
machine precision.

The reference mirrors the current C++ residual-eigenmode quadrupole rows: for each element, the remaining
face-charge modes are the orthonormal complement of monopole + three dipole moment functionals.  This keeps
the verification aligned with the hex/wedge/pyramid-capable BuildMomentSystemCore implementation.
"""
import numpy as np
import radia as rad
from multipole_moment_iter_scaling import build_cyoke_hexes


def _norm(row, rhs):
    n = np.linalg.norm(row)
    return (row / n, rhs / n) if n > 1e-300 else (row, rhs)


def _residual_eigenmodes(Ae, d):
    nF = len(Ae)
    basis = []

    def ortho_add(v):
        v = np.asarray(v, float).copy()
        for b in basis:
            v -= np.dot(b, v) * b
        n = np.linalg.norm(v)
        if n > 1e-9:
            q = v / n
            basis.append(q)
            return q
        return None

    ortho_add(Ae)
    for k in range(3):
        ortho_add(Ae * d[:, k])

    modes = []
    for e in range(nF):
        v = np.zeros(nF)
        v[e] = 1.0
        q = ortho_add(v)
        if q is not None:
            modes.append(q)
        if len(basis) >= nF:
            break
    return modes


def main():
    nxy, nz, chi = 8, 2, 999.0
    Happ = np.array([0.0, 1e3, 0.0])
    hexes = build_cyoke_hexes(nxy, nz)
    rad.UtiDelAll(); rad.set_demag_backend("yano")
    objs = [rad.ObjHexahedron([list(v) for v in V], [0, 0, 0]) for V in hexes]
    for h in objs:
        rad.MatApl(h, rad.MatLin(chi + 1.0))
    hd = rad.BuildMatrix(rad.ObjCnt(objs))

    # (a) Python prototype assembly
    G = np.asarray(rad.GetFaceGeom(hd), float); C = np.asarray(rad.GetCentroidFieldGrad(hd), float)
    dof = G.shape[0]; elem = G[:, 0].astype(int); area = G[:, 1]; fc = G[:, 2:5]; nrm = G[:, 5:8]; ecen = G[:, 8:11]
    n_el = int(elem.max()) + 1; dofs_of = [np.where(elem == e)[0] for e in range(n_el)]
    A_py = np.zeros((dof, dof)); b_py = np.zeros(dof); r = 0
    for e in range(n_el):
        fs = dofs_of[e]; Ae = area[fs]; ne = nrm[fs]; ce = ecen[fs[0]]; d = fc[fs] - ce
        Ve = (1.0 / 3.0) * np.sum(Ae * np.sum(fc[fs] * ne, axis=1)); F0 = C[e, 0:3, :]; Ginv = C[e, 3:9, :]
        dip = np.zeros((3, dof))
        for k in range(3):
            dip[k, fs] = Ae * d[:, k]
        for k in range(3):
            row, rhs = _norm(dip[k, :] / Ve - chi * F0[k, :], chi * Happ[k]); A_py[r, :] = row; b_py[r] = rhs; r += 1
        mono = np.zeros(dof); mono[fs] = Ae
        row, rhs = _norm(mono, 0.0); A_py[r, :] = row; b_py[r] = rhs; r += 1
        for phi in _residual_eigenmodes(Ae, d):
            Bm = np.divide(phi, Ae, out=np.zeros_like(phi), where=Ae > 1e-300)
            row = np.zeros(dof); row[fs] = Ae * Bm
            cm = np.array([np.sum(Ae * ne[:, k] * Bm) for k in range(3)]); row -= (cm @ dip) / Ve
            Dm = np.array([[np.sum(Ae * d[:, jj] * ne[:, ii] * Bm) for jj in range(3)] for ii in range(3)])
            Dvec = np.array([Dm[0, 0], Dm[1, 1], Dm[2, 2], Dm[0, 1]+Dm[1, 0], Dm[0, 2]+Dm[2, 0], Dm[1, 2]+Dm[2, 1]])
            row -= chi * (Dvec @ Ginv); row, rhs = _norm(row, 0.0); A_py[r, :] = row; b_py[r] = rhs; r += 1

    # (b) C++ assembly
    A_cpp, rhs_cpp, dof2 = rad.BuildMomentSystem(hd, chi, float(Happ[0]), float(Happ[1]), float(Happ[2]))
    A_cpp = np.asarray(A_cpp, float); rhs_cpp = np.asarray(rhs_cpp, float)
    rad.UtiDelAll()

    dA = float(np.max(np.abs(A_cpp - A_py))); db = float(np.max(np.abs(rhs_cpp - b_py)))
    print(f"dof: py={dof} cpp={dof2}")
    print(f"max|A_cpp - A_py|   = {dA:.3e}")
    print(f"max|rhs_cpp - b_py| = {db:.3e}")
    ok = (dof == dof2) and dA < 1e-9 and db < 1e-9
    print("RESULT:", "PASS -- C++ moment system == Python prototype" if ok else "FAIL")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
