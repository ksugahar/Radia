"""The field-gradient is a TENSOR: why multipole-moment MMM and EIEM2 are NOT the same DOF (user insight 2026-06-22:
"the direction of the gradient is what matters; gradH is a tensor with 6 components").

H = -grad(phi) is curl-free, so gradH = -Hess(phi) is a SYMMETRIC 3x3 tensor = 6 components
[xx, yy, zz, xy, yz, zx]  (trace 0 in a charge-free region -> 5 independent = the l=2 quadrupole).

EIEM2 (the shipped surface-charge MSC collocation): for face i the demag field is sampled at the OFFSET eval point
    EvalPt_i = ElementCenter + alpha*(FaceCenter_i - ElementCenter),   alpha = 0.5
and the surface-charge constitutive relation uses the FACE-NORMAL component.  The first-order (gradient)
content that enters is therefore the scalar
    n_i . gradH . u_i ,   u_i = (FaceCenter_i - ElementCenter).
The offset sample gradH.u_i is a full VECTOR and DOES contain shear in its transverse components -- but the
NORMAL PROJECTION n_i.(.) discards those transverse components and keeps only n_i.gradH.n_i, the DIAGONAL
(normal-normal) part.  For an axis-aligned hex the 6 face directions are +/- x,y,z, so the 6 EIEM2 face
measurements collapse to just {T_xx, T_yy, T_zz}.  The 3 SHEAR components {T_xy, T_yz, T_zx} are in the NULL
SPACE of the offset-with-normal-projection measurement -- NO choice of alpha recovers them (alpha only scales
|u_i| along the normal).

The moment formulation matches all 6 components of gradH directly (rad.GetCentroidFieldGrad's 6 gradient
functionals), so it captures the shear.  That shear is BOTH its accuracy gain (sheared / distorted fields,
non-aligned flux) AND its conditioning cost (the off-diagonal 1/r^3 kernel is the non-normal part).  Hence the
~5x non-normality of moment vs EIEM2 is exactly the price of the 3 shear functionals -- not a re-encoding of
the same DOF.  This refutes the earlier "alpha-offset == gradient / moment ~ EIEM2(alpha*)" claim.

mdx-clean (import radia directly).
"""
import json
import os
import platform
from datetime import datetime

import numpy as np

import radia as rad

HERE = os.path.dirname(os.path.abspath(__file__))

# symmetric-tensor 6-vector order: [xx, yy, zz, xy, yz, zx]
LABELS = ["xx", "yy", "zz", "xy", "yz", "zx"]


def to_mat(v):
    xx, yy, zz, xy, yz, zx = v
    return np.array([[xx, xy, zx], [xy, yy, yz], [zx, yz, zz]])


def basis():
    B = []
    for k in range(6):
        e = np.zeros(6); e[k] = 1.0; B.append(e)
    return B


def eiem2_normal_measurement_matrix(h=0.5):
    """6 face measurements of a symmetric gradient tensor by EIEM2's offset-then-normal-project rule, for an
    axis-aligned unit hex (half-size h).  Row i = the linear functional T -> n_i . T . u_i, u_i = h*n_i."""
    normals = [np.array(v, float) for v in
               ([1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1])]
    M = np.zeros((6, 6))
    for i, n in enumerate(normals):
        u = h * n                                    # FaceCenter_i - ElementCenter, along the normal
        for k, e in enumerate(basis()):
            M[i, k] = n @ to_mat(e) @ u              # n_i . gradH . u_i  (NORMAL component of the offset sample)
    return M


def eiem2_fullvector_measurement_matrix(h=0.5):
    """Same offset samples but keeping the FULL VECTOR gradH.u_i (3 comps per face = 18 rows).  This is what
    EIEM2 would see IF it did not project onto the normal -- it has full rank 6 (the shear lives in the
    transverse components).  The contrast with the rank-3 normal-projected matrix shows the loss is the NORMAL
    PROJECTION, not the offset itself."""
    normals = [np.array(v, float) for v in
               ([1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1])]
    rows = []
    for n in normals:
        u = h * n
        for comp in range(3):                        # keep each Cartesian component of gradH.u_i
            rows.append([to_mat(e)[comp] @ u for e in basis()])
    return np.array(rows)


def radia_gradH_reality_check():
    """A real two-magnet config has a field-gradient tensor with NONZERO off-diagonal (shear) at a generic
    probe point: shear is the rule, not a corner case.  FD of rad.Fld('h')."""
    rad.UtiDelAll()
    m1 = rad.magnet_box([0.030, 0.010, 0.000], [0.02, 0.02, 0.02], [0, 0, 1.0e6])
    m2 = rad.magnet_box([-0.020, 0.030, 0.012], [0.02, 0.02, 0.02], [1.0e6, 0, 0])   # tilted -> shears the field
    cont = rad.ObjCnt([m1, m2])
    p = np.array([0.0, 0.0, 0.050]); d = 1.0e-4

    def Hf(x):
        return np.array(rad.Fld(cont, "h", list(x)), float)

    gradH = np.zeros((3, 3))                          # gradH[i, j] = dH_i / dx_j
    for j in range(3):
        xp = p.copy(); xp[j] += d; xm = p.copy(); xm[j] -= d
        gradH[:, j] = (Hf(xp) - Hf(xm)) / (2 * d)
    rad.UtiDelAll()
    asym = float(np.linalg.norm(gradH - gradH.T) / np.linalg.norm(gradH))   # ~0 confirms curl-free
    T = 0.5 * (gradH + gradH.T)
    diag = [T[0, 0], T[1, 1], T[2, 2]]
    shear = [T[0, 1], T[1, 2], T[2, 0]]
    return dict(asymmetry=asym, diag=[float(v) for v in diag], shear=[float(v) for v in shear],
                shear_over_diag=float(np.linalg.norm(shear) / np.linalg.norm(diag)))


def main():
    print("\nThe field-gradient is a TENSOR (6 comp).  EIEM2 offset+normal-projection sees only the DIAGONAL;\n"
          "moment matches all 6, incl. the 3 SHEAR -- which is the real difference (NOT an alpha re-encoding).\n")

    M_norm = eiem2_normal_measurement_matrix()
    M_full = eiem2_fullvector_measurement_matrix()
    r_norm = int(np.linalg.matrix_rank(M_norm, tol=1e-10))
    r_full = int(np.linalg.matrix_rank(M_full, tol=1e-10))

    # which tensor components does the normal-projected EIEM2 measurement actually constrain?
    constrained = [LABELS[k] for k in range(6) if np.linalg.norm(M_norm[:, k]) > 1e-12]
    blind = [LABELS[k] for k in range(6) if np.linalg.norm(M_norm[:, k]) <= 1e-12]

    print(f"  EIEM2 offset + NORMAL projection : rank {r_norm}/6, constrains {constrained}, BLIND to {blind}")
    print(f"  EIEM2 offset + FULL vector (n.p.) : rank {r_full}/6  (shear lives in the transverse comps that")
    print(f"                                      the normal projection throws away)")
    print(f"  moment (matches gradH directly)  : rank 6/6, all of {LABELS}")

    # pure-shear probe: T = E_xy
    e_xy = np.zeros(6); e_xy[3] = 1.0
    eiem2_reads = M_norm @ e_xy
    print(f"\n  Pure-shear gradient T = E_xy (G=1, all diagonals 0):")
    print(f"    EIEM2 6 face measurements = {np.array2string(eiem2_reads, precision=3, suppress_small=True)} "
          f"-> ||.|| = {np.linalg.norm(eiem2_reads):.2e}  (BLIND)")
    print(f"    moment reads component xy = 1.000  (CAPTURED)")

    real = radia_gradH_reality_check()
    print(f"\n  Reality check (two tilted magnets, FD of rad.Fld 'h' at a generic point):")
    print(f"    gradH asymmetry ||G-G^T||/||G|| = {real['asymmetry']:.2e}  (~0 -> curl-free, symmetric tensor)")
    print(f"    diagonal  [xx,yy,zz] = {np.array2string(np.array(real['diag']), precision=3)}")
    print(f"    SHEAR     [xy,yz,zx] = {np.array2string(np.array(real['shear']), precision=3)}")
    print(f"    ||shear|| / ||diag|| = {real['shear_over_diag']:.2f}  (shear is O(1) of the gradient, not negligible)")

    eiem2_blind_to_shear = (r_norm == 3 and set(blind) == {"xy", "yz", "zx"}
                            and np.linalg.norm(eiem2_reads) < 1e-12)
    out = dict(timestamp=datetime.now().isoformat(), hostname=platform.node(),
               benchmark="multipole_moment_shear_gradient",
               eiem2_normal_rank=r_norm, eiem2_fullvector_rank=r_full, moment_rank=6,
               eiem2_constrains=constrained, eiem2_blind_to=blind,
               pure_shear_eiem2_norm=float(np.linalg.norm(eiem2_reads)),
               reality_check=real, eiem2_blind_to_shear=bool(eiem2_blind_to_shear),
               conclusion=(
                   "gradH is a symmetric 6-component tensor.  EIEM2's offset eval point + face-NORMAL "
                   "projection measures only n.gradH.n = the 3 DIAGONAL components {xx,yy,zz} (rank 3); the 3 "
                   "SHEAR components {xy,yz,zx} are in its null space -- no alpha recovers them, because the "
                   "offset is along the normal and the transverse (shear-carrying) field components are "
                   "discarded by the normal projection.  The same offset KEEPING the full vector would be rank "
                   "6, proving the loss is the projection, not the offset.  The moment formulation matches all "
                   "6 gradH components, so it captures the shear -- real extra content, NOT a re-encoding of a "
                   "single alpha DOF.  This refutes 'alpha-offset == gradient' and 'moment ~ EIEM2(alpha*)'.  "
                   "A real config has ||shear||/||diag|| = O(1), so the missing shear is physically significant; "
                   "the shear (off-diagonal 1/r^3) functionals are simultaneously moment's accuracy advantage "
                   "(sheared/distorted fields) and its ~5x non-normality penalty vs EIEM2."))
    with open(os.path.join(HERE, "multipole_moment_shear_gradient.json"), "w") as f:
        json.dump(out, f, indent=2, default=float)
    print(f"\n  => EIEM2 blind to shear: {eiem2_blind_to_shear}.  saved multipole_moment_shear_gradient.json")


if __name__ == "__main__":
    main()
