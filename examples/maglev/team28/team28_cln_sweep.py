"""team28_cln_sweep.py -- CLN-reduced TEAM 28 levitation force vs height.

Sweeps the disk height and, at each height, computes the levitation force
both from the full split-(K+sN) FEM solve and from the N-stage CLN/Cauer
reduced model, and compares to the lab full-FEM ground truth Fz1(dZ).

Confirms the 5-6 stage CLN reproduces the WHOLE force-vs-height curve
(incl. the levitation equilibrium ~dZ=+4mm where lift == disk weight),
not just the dZ=0 point.

Run:  python team28_cln_sweep.py
"""
import numpy as np
import scipy.sparse.linalg as spla
from numpy import pi

from team28_cln_force import (
    setup, to_csr, force_from_vec, aluminium_z, FREQ,
)

# Lab full-FEM ground truth (W:\...\2024_08_TEAM28\50Hz_可動\axisymmetric_mixed.mat)
# dZ in mm -> Fz1 in N (NEGATIVE = upward lift).
REF = {
    -6: -5.6713, -3: -3.5780, 0: -2.1928, 4: -1.0655,
    8: -0.4483, 12: -0.1279, 15: -0.0033,
}
DISK_WEIGHT = 1.055   # N (R=65mm, t=3mm Al disk, 107.5 g)
NSTAGE = 6


def cln_force_at(al_z, nstage=NSTAGE):
    mesh, fes, fesPhi, fesB, sig, K, Nmat, F = setup(al_z)
    ndof = fes.ndof
    s = 2 * pi * FREQ * 1j
    free = np.array([i for i in range(ndof) if fes.FreeDofs()[i]])
    Kf = to_csr(K.mat, ndof)[free][:, free].tocsc()
    Nf = to_csr(Nmat.mat, ndof)[free][:, free].tocsc()
    Ff = np.array(F.vec.FV().NumPy(), dtype=complex)[free]

    Klu = spla.splu(Kf)
    # full
    xf = np.zeros(ndof, dtype=complex)
    xf[free] = spla.spsolve((Kf + s * Nf).tocsc(), Ff)
    fz_full = force_from_vec(xf, fes, fesPhi, fesB, sig)
    # CLN/Cauer Krylov basis from the coil source
    V = []
    v0 = Klu.solve(Ff); v0 /= np.sqrt(abs(v0 @ v0.conj())); V.append(v0)
    for _ in range(nstage - 1):
        w = Klu.solve(Nf @ V[-1])
        for vv in V:
            w = w - (vv.conj() @ w) * vv
        nrm = np.sqrt(abs(w @ w.conj()))
        if nrm < 1e-12:
            break
        V.append(w / nrm)
    Vm = np.array(V).T
    y = np.linalg.solve(Vm.conj().T @ (Kf @ Vm) + s * (Vm.conj().T @ (Nf @ Vm)),
                        Vm.conj().T @ Ff)
    xr = np.zeros(ndof, dtype=complex); xr[free] = Vm @ y
    fz_cln = force_from_vec(xr, fes, fesPhi, fesB, sig)
    return fz_full, fz_cln


def main():
    print(f" CLN stages = {NSTAGE};  disk weight = {DISK_WEIGHT:.3f} N "
          f"(equilibrium where lift == -weight)\n")
    print(" dZ[mm]  full-FEM   CLN(6)    lab-ref    CLN-vs-lab")
    rows = []
    for dz_mm in sorted(REF):
        al_z = aluminium_z + dz_mm * 1e-3
        fz_full, fz_cln = cln_force_at(al_z)
        ref = REF[dz_mm]
        err = abs(fz_cln - ref) / abs(ref) * 100 if ref != 0 else float('nan')
        print(f"  {dz_mm:4d}   {fz_full:+.4f}  {fz_cln:+.4f}  {ref:+.4f}    {err:6.2f} %")
        rows.append((dz_mm, fz_full, fz_cln, ref))

    # Physical levitation equilibrium.  The lab/repo F_z column is the verbatim
    # TEAM 28 integral Re[B_r J_t] = 2x the physical time-averaged Lorentz force
    # (verified ratio 1.9998); the disk floats where the PHYSICAL lift == weight,
    # so balance 0.5*F_z against -DISK_WEIGHT (NOT F_z, which gave a spurious
    # +4mm / 14.9mm height).  Absolute disk-bottom z = 10.8mm + dZ.
    PHYS = 0.5
    dz = np.array([r[0] for r in rows], float)
    fc = PHYS * np.array([r[2] for r in rows], float)
    tgt = -DISK_WEIGHT
    eq = None
    for i in range(len(dz) - 1):
        if (fc[i] - tgt) * (fc[i + 1] - tgt) <= 0:
            t = (tgt - fc[i]) / (fc[i + 1] - fc[i])
            eq = dz[i] + t * (dz[i + 1] - dz[i]); break
    if eq is not None:
        print(f"\n CLN levitation equilibrium (PHYSICAL lift==weight): dZ = "
              f"{eq:.2f} mm  -> absolute z = {10.8 + eq:.2f} mm")
        print(" (published measured steady-state z = 11.5 mm)")
    else:
        print("\n (equilibrium not bracketed)")


if __name__ == "__main__":
    main()
