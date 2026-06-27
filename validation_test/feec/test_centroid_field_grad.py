"""Golden test for the C++ surface-charge MSC centroid field+gradient accessor: rad.GetCentroidFieldGrad(handle).

This is the kernel of the parameter-free surface-charge MSC moment formulation (the fix for the eval-point alpha and
the finite-difference conditioning noise; see examples/vim/multipole_moment_analytic_selfterm.py).  For each hex
element it returns the demag field H and gradient gradH at the element CENTROID as linear functionals of
every source DOF charge:
  SELF face  -> bare charged-face field (interior centroid -> finite, no center charge);
  MUTUAL face -> collocation MMMM dipole layer = bare face - area*(point @ source center)  (finite -> singularity-free).
Field convention H = (1/4pi) int sigma (r-r')/|r-r'|^3 dA'.  Output is (nHex, 9, dof) with component order
(Hx,Hy,Hz, gxx,gyy,gzz,gxy,gxz,gyz); C[e,0:3,:] is the centroid field functional, C[e,3:9,:] the gradient.

This test locks:
  (1) shape (nHex, 9, dof) with dof == 6*nHex and rows aligned to GetFaceGeom/GetInteractMatrix;
  (2) the SELF term is patch-test exact -- a single cube with the uniform-M surface charge sigma_f = n_f.z
      reproduces the analytic demag tensor: H(center) = -(1/3) z_hat (N = 1/3) and gradH(center) = 0 (by the
      cube's symmetry).  This pins the 1/4pi factor, the sign, and that the self term needs no center charge;
  (3) the gradient is genuinely computed (non-trivially nonzero on a non-affine, non-centrally-symmetric hex).
"""
import sys

import numpy as np
import pytest

sys.path.insert(0, r"S:\Radia\01_GitHub\src\radia")
import radia as rad


def _cube(L):
    return [[-L, -L, -L], [L, -L, -L], [L, L, -L], [-L, L, -L],
            [-L, -L, L], [L, -L, L], [L, L, L], [-L, L, L]]


def _build_get(hexes, mu_r=1e3):
    rad.UtiDelAll(); rad.set_demag_backend("collocation_mmmm")
    objs = []
    for V in hexes:
        ob = rad.ObjHexahedron([list(v) for v in V], [0, 0, 0])
        rad.MatApl(ob, rad.MatLin(mu_r)); objs.append(ob)
    handle = rad.BuildMatrix(rad.ObjCnt(objs))
    G = np.asarray(rad.GetFaceGeom(handle), float)
    C = np.asarray(rad.GetCentroidFieldGrad(handle), float)
    rad.UtiDelAll()
    return G, C


def test_shape_and_alignment():
    # a 2x1x1 hex bar -> 2 hexes, 12 DOF; accessor is (2, 9, 12) and aligned to GetFaceGeom rows
    L = 0.02
    hexes = [[[x0, -L, -L], [x0 + 2 * L, -L, -L], [x0 + 2 * L, L, -L], [x0, L, -L],
              [x0, -L, L], [x0 + 2 * L, -L, L], [x0 + 2 * L, L, L], [x0, L, L]] for x0 in (-2 * L, 0.0)]
    G, C = _build_get(hexes)
    dof = G.shape[0]
    assert dof == 12
    assert C.shape == (2, 9, dof)
    # GetFaceGeom elem_local column identifies which hex each DOF belongs to (0 or 1), 6 faces each
    elem = G[:, 0].astype(int)
    assert sorted(np.bincount(elem).tolist()) == [6, 6]


def test_cube_demag_self_term():
    # single cube: the SELF centroid field of the uniform-M surface charge = the demag tensor N = 1/3.
    L = 0.02
    G, C = _build_get([_cube(L)])
    dof = G.shape[0]
    assert C.shape == (1, 9, dof)
    sigma = G[:, 5:8][:, 2]                          # M = z_hat -> sigma_f = n_f . z (top +1, bottom -1)
    H = C[0, 0:3, :] @ sigma
    gradH = C[0, 3:9, :] @ sigma
    assert abs(H[0]) < 1e-5 and abs(H[1]) < 1e-5, f"transverse field not ~0: {H}"
    assert abs(H[2] - (-1.0 / 3.0)) < 1e-4, f"H_z = {H[2]:.6f} != -1/3 (cube demag N=1/3)"
    assert np.abs(gradH).max() * L < 1e-6, f"gradH not ~0 at cube center (symmetry): {gradH}"


def test_gradient_is_computed_on_nonaffine_hex():
    # a frustum (top face shrunk) is NON-affine and NOT z-symmetric -> a real nonzero gradient at the
    # centroid (the cube test alone could pass with a hardcoded-zero gradient; this rules that out).
    L = 0.02
    fr = [[-L, -L, -L], [L, -L, -L], [L, L, -L], [-L, L, -L],
          [-L / 2, -L / 2, L], [L / 2, -L / 2, L], [L / 2, L / 2, L], [-L / 2, L / 2, L]]
    G, C = _build_get([fr])
    sigma = G[:, 5:8][:, 2]
    gradH = C[0, 3:9, :] @ sigma
    assert np.all(np.isfinite(gradH))
    # gzz (index 2) must be genuinely nonzero (top/bottom asymmetry), and everything bounded
    assert abs(gradH[2]) * L > 1e-3, f"expected a real gz gradient on the frustum, got {gradH}"
    assert np.abs(gradH).max() * L < 1e3


if __name__ == "__main__":
    test_shape_and_alignment(); print("OK shape/alignment")
    test_cube_demag_self_term(); print("OK cube demag self-term (N=1/3)")
    test_gradient_is_computed_on_nonaffine_hex(); print("OK gradient computed on non-affine hex")
