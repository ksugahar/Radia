"""Phase-1 verification (EIEM2 full-deletion track): the multipole-moment MMM path now supports IMA (image symmetry).
BuildCentroidFieldGrad adds the mirror-image contributions (field AND the rank-2 gradient, recomputed fresh
from the reflected geometry), so a HALF model solved with image= reproduces the EXPLICIT FULL model (the half
plus its hand-mirrored z<0 copy) solved WITHOUT image -- to machine precision, since IMA is just the
computational shortcut for the explicit mirror within the SAME moment formulation.

Covers the symmetric BC (applied field PARALLEL to the z=0 plane -> Hx -> image '+z') and the antisymmetric BC
(field PERPENDICULAR -> Hz -> image '-z').  This is the parity item that lets yano+IMA (e.g. the accelerator
MSC panel quarter/eighth models, calc_accel_msc.py --ima ...) run on moment instead of EIEM2.
"""
import numpy as np
import radia as rad

MU0 = 4e-7 * np.pi
MUR = 200.0


def _half_boxes():
    """2x2 layer of hexes entirely in z>0 (disjoint from its z<0 mirror -> no boundary elements on z=0)."""
    out = []
    for ix in range(2):
        for iy in range(2):
            x0, y0 = -0.02 + ix * 0.02, -0.02 + iy * 0.02
            out.append([[x0, y0, 0.006], [x0 + 0.02, y0, 0.006], [x0 + 0.02, y0 + 0.02, 0.006], [x0, y0 + 0.02, 0.006],
                        [x0, y0, 0.026], [x0 + 0.02, y0, 0.026], [x0 + 0.02, y0 + 0.02, 0.026], [x0, y0 + 0.02, 0.026]])
    return out


def _solve(boxes, Happ, image):
    rad.UtiDelAll(); rad.set_demag_backend("yano")
    objs = [rad.ObjHexahedron(b, [0, 0, 0]) for b in boxes]
    for h in objs:
        rad.MatApl(h, rad.MatLin(MUR))
    cont = rad.ObjCnt(objs + [rad.ObjBckg(lambda p: [MU0 * Happ[0], MU0 * Happ[1], MU0 * Happ[2]])])
    kw = {} if image is None else {"image": image}
    rad.Solve(cont, 1e-8, 500, 0, **kw)
    return np.asarray([m[1] for m in rad.ObjM(rad.ObjCnt(objs))], float)


def _case(name, Happ, image):
    half = _half_boxes()
    full = half + [[[p[0], p[1], -p[2]] for p in b] for b in half]
    M_ref = _solve(full, Happ, None)[:len(half)]
    M_ima = _solve(half, Happ, image)
    rel = np.linalg.norm(M_ima - M_ref) / max(np.linalg.norm(M_ref), 1e-30)
    print(f"{name}: Happ={Happ} image={image!r}  ||M||={np.linalg.norm(M_ima):.4e}  "
          f"rel_to_explicit_full={rel:.2e}  {'OK' if rel < 1e-6 else 'MISMATCH'}")
    return rel < 1e-6


def main():
    ok = True
    ok &= _case("symmetric     (Hx || z=0)", [1000.0, 0.0, 0.0], "+z")
    ok &= _case("antisymmetric (Hz _|_ z=0)", [0.0, 0.0, 1000.0], "-z")
    print("RESULT:", "PASS -- moment IMA reproduces the explicit full model" if ok else "FAIL")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
