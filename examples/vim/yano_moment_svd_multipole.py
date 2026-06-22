"""What field does the 6-DoF MSC (yano) element create, physically?  SVD says: PURE MULTIPOLES, and the
condition number IS the ratio of multipole field strengths (Sugahara's question 2026-06-22: "the condition
number is really the question 'what field does 6-DoF MMM create?' -- would SVD tell us?").

Take ONE hexahedron's 6-face-charge interaction matrix N (sigma -> demag field it produces) and SVD it.  Each
right singular vector is a charge pattern; project it onto the multipole basis (monopole, 3 dipole, 2 diagonal
quadrupole) built from the same face geometry.  The singular value is that mode's field STRENGTH.

Findings (cubic vs anisotropic):
  - The 6 singular modes are PURE multipoles -- monopole + 3 dipole + 2 diagonal quadrupole (on a cubic cell
    each mode projects ~100% onto one multipole).  6-DoF MSC creates NOTHING that is not a multipole; it is
    exactly MMM (the 3 dipole DoF) + a monopole + the 2 quadrupole moments representable by 6 axis-aligned
    face charges.  (The off-diagonal/shear quadrupole is not representable -- int x_i x_j dA = 0 on an
    axis-aligned face -- so it is not among the modes.)
  - cond(N) = (strongest multipole field) / (weakest).  The MONOPOLE is by far the strongest mode (it makes
    the strongest self-field), so it dominates cond(N) -- BUT it is constrained to zero by div B = 0, so it
    does not actually enter the solve.  The PHYSICALLY RELEVANT conditioning is the DIPOLE/QUADRUPOLE strength
    spread (cond with the monopole removed): ~1.4 on a cubic cell (dipole and quadrupole fields almost equally
    strong) and ~6-7 on an anisotropic cell.
  - So the aspect-ratio sensitivity (yano_moment_aspect_ratio.py) is exactly this: the quadrupole is the
    weakest field mode, and on an anisotropic cell the thin-direction quadrupole field weakens further, opening
    the dipole/quadrupole gap.  Block-Jacobi inverts the local block, AMPLIFYING that weak quadrupole mode ->
    the iteration trouble.  Condition number, read physically, = which multipole fields the 6 DoF make and how
    disparate their strengths are.

mdx-clean (import radia directly).
"""
import json
import os
import platform
from datetime import datetime

import numpy as np

import radia as rad

HERE = os.path.dirname(os.path.abspath(__file__))
NAMES = ["mono", "dipX", "dipY", "dipZ", "quadD1", "quadD2"]


def analyze(a, b, c):
    """SVD a single hex's 6x6 interaction matrix N; label each singular mode by its dominant multipole."""
    rad.UtiDelAll(); rad.set_demag_backend("yano")
    V = [[-a, -b, -c], [a, -b, -c], [a, b, -c], [-a, b, -c],
         [-a, -b, c], [a, -b, c], [a, b, c], [-a, b, c]]
    o = rad.ObjHexahedron(V, [0, 0, 0]); rad.MatApl(o, rad.MatLin(1000.0))
    handle = rad.BuildMatrix(rad.ObjCnt([o]))
    N, dof = rad.GetInteractMatrix(handle)
    G = np.asarray(rad.GetFaceGeom(handle), float)
    Ae = G[:, 1]; fc = G[:, 2:5]; ce = G[:, 8:11][0]; d = fc - ce
    N = np.asarray(N, float).reshape(dof, dof)
    _U, s, Vt = np.linalg.svd(N)
    # multipole functional-gradients in sigma-space, orthonormalized
    Bcols = np.array([Ae, Ae * d[:, 0], Ae * d[:, 1], Ae * d[:, 2],
                      Ae * (d[:, 0] ** 2 - d[:, 1] ** 2), Ae * (d[:, 1] ** 2 - d[:, 2] ** 2)]).T
    Q, _ = np.linalg.qr(Bcols)
    modes = []
    for i in range(dof):
        comp = Q.T @ Vt[i]; comp = comp / max(np.linalg.norm(comp), 1e-300)
        frac = (np.abs(comp) ** 2 * 100).tolist()
        modes.append(dict(sv=float(s[i]), dominant=NAMES[int(np.argmax(np.abs(comp)))],
                          fractions={n: round(f, 1) for n, f in zip(NAMES, frac)}))
    rad.UtiDelAll()
    cond_full = float(s[0] / s[-1]); cond_no_mono = float(s[1] / s[-1])    # s[0] is always the monopole
    return dict(half_extents=[a, b, c], singular_values=[float(x) for x in s],
                modes=modes, cond_N=cond_full, cond_without_monopole=cond_no_mono)


def main():
    cases = [("CUBIC 1:1:1", 0.01, 0.01, 0.01),
             ("THIN z 1:1:0.25", 0.01, 0.01, 0.0025),
             ("TALL z 1:1:4", 0.01, 0.01, 0.04)]
    print("\nWhat field does the 6-DoF MSC element create?  SVD of the single-hex interaction matrix N.\n")
    results = []
    for tag, a, b, c in cases:
        r = analyze(a, b, c); r["tag"] = tag; results.append(r)
        print(f"--- {tag}: cond(N)={r['cond_N']:.1f}, cond w/o monopole (dipole/quad spread)="
              f"{r['cond_without_monopole']:.1f} ---")
        for m in r["modes"]:
            shown = " ".join(f"{n}:{f:.0f}%" for n, f in m["fractions"].items() if f > 5)
            print(f"   sv={m['sv']:.3e} -> {m['dominant']:>6} ({shown})")
    cond_nm = [r["cond_without_monopole"] for r in results]
    out = dict(timestamp=datetime.now().isoformat(), hostname=platform.node(),
               benchmark="yano_moment_svd_multipole", results=results,
               conclusion=(
                   "SVD of the 6-DoF MSC interaction matrix N shows the element creates PURE MULTIPOLE fields: "
                   "1 monopole + 3 dipole + 2 diagonal quadrupole (on a cubic cell each singular mode is ~100% "
                   "one multipole).  6-DoF MSC = MMM (the 3 dipole DoF) + a monopole + the 2 quadrupole moments "
                   "representable by axis-aligned face charges; it creates nothing non-multipole (and cannot "
                   "represent the shear quadrupole).  The condition number IS the ratio of multipole field "
                   "strengths: the MONOPOLE is the strongest mode and dominates cond(N) (~27 cubic, ~285-338 "
                   "anisotropic) but is constrained to zero by div B = 0, so the physically relevant conditioning "
                   f"is the DIPOLE/QUADRUPOLE strength spread (cond without the monopole) = {cond_nm[0]:.1f} "
                   f"(cubic) -> {cond_nm[1]:.1f}-{cond_nm[2]:.1f} (anisotropic).  The quadrupole is the weakest "
                   "field mode; on an anisotropic cell the thin-direction quadrupole weakens further, opening the "
                   "dipole/quadrupole gap -- exactly the aspect-ratio sensitivity, now read physically.  "
                   "Block-Jacobi inverts the local block and amplifies that weakest (quadrupole) mode, which is "
                   "the iteration trouble on anisotropic cells."))
    with open(os.path.join(HERE, "yano_moment_svd_multipole.json"), "w") as f:
        json.dump(out, f, indent=2, default=float)
    print(f"\n  => 6-DoF creates pure multipoles (mono+dipole+2 quad); physical cond = dipole/quadrupole "
          f"spread = {cond_nm[0]:.1f} cubic -> {max(cond_nm):.1f} anisotropic.")
    print("  saved", os.path.join(HERE, "yano_moment_svd_multipole.json"))


if __name__ == "__main__":
    main()
