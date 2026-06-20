"""Golden test for the C++ yano-MSC loop (cell-graph cycle) basis: rad.GetLoopBasis(handle).

Stage 1 of wiring loop removal into the yano-MSC solve (see examples/vim/yano_loop_removal_scaling.py
for the exact-field + scaling evidence).  The loops are the field-null subspace of the surface-charge
collocation operator N (== the HDiv ker(B) cell/internal-face cycle space).  This test locks:

  (1) rad.GetLoopBasis returns L with shape (dof, nLoop);
  (2) nLoop == the topological cell-graph cycle count  n_internal_faces - (n_cells - 1)  (per component);
  (3) the columns of L span ker(N) on the EXACT C++ operator: ||N @ L|| / (||N|| ||L||) ~ 1e-12
      (this is the self-check the user asked for -- it confirms the geometry-only cycle basis, built in
      C++ from face-center matching, is aligned to the C++ DOF<->face order for ANY distortion);
  (4) rank(L) == nLoop (independent basis).

Self-contained: a structured + affine-sheared hex bar (affine keeps quad faces planar -> valid for
ObjHexahedron) is generated here; no external mesh.
"""
import sys
import numpy as np
import pytest

sys.path.insert(0, r"S:\Radia\01_GitHub\src\radia")
import radia as rad

# standard hex (VTK/Netgen MMB8T) corner order: bottom CCW 0-3, top CCW 4-7
_HEX_FACES = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]


def _build_bar(nx, ny, nz, h=1.0, distort=0.0):
    """structured hex bar; distort>0 -> graded spacing + affine shear (planar faces preserved)."""
    nxp, nyp, nzp = nx + 1, ny + 1, nz + 1

    def nid(i, j, k):
        return i + nxp * (j + nyp * k)

    def axis(n):
        w = [h * (1.0 + distort * 0.5 * ((((i * 7 + 3) % 5) / 4.0) - 0.5)) for i in range(n)]
        c = [0.0]
        for wi in w:
            c.append(c[-1] + wi)
        return c

    xs, ys, zs = axis(nx), axis(ny), axis(nz)
    S = np.eye(3)
    sh = distort * 0.4
    S[0, 1], S[0, 2], S[1, 2] = sh, sh * 0.7, sh * 0.5
    pos = {}
    for k in range(nzp):
        for j in range(nyp):
            for i in range(nxp):
                pos[nid(i, j, k)] = S @ np.array([xs[i], ys[j], zs[k]], float)
    hexes, elems = [], []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                c = [nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                     nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1)]
                elems.append(c)
                hexes.append(np.array([pos[n] for n in c]))
    return hexes, elems


def _cycle_count(elems):
    from collections import defaultdict
    occ = defaultdict(int)
    for c in elems:
        for f in _HEX_FACES:
            occ[tuple(sorted(c[i] for i in f))] += 1
    n_internal = sum(1 for v in occ.values() if v == 2)
    return n_internal - (len(elems) - 1)   # single connected component (a solid bar)


def _yano_matrix_and_loops(hexes, mu_r=1e3, scale=1e-2):
    rad.UtiDelAll()
    rad.set_demag_backend("yano")
    objs = []
    for V in hexes:
        ob = rad.ObjHexahedron([(V[i] * scale).tolist() for i in range(8)], [0, 0, 0])
        rad.MatApl(ob, rad.MatLin(mu_r))
        objs.append(ob)
    cont = rad.ObjCnt(objs)
    handle = rad.BuildMatrix(cont)
    N, dof = rad.GetInteractMatrix(handle)
    L, nLoop = rad.GetLoopBasis(handle)
    rad.UtiDelAll()
    return np.asarray(N, float), int(dof), np.asarray(L, float), int(nLoop)


@pytest.mark.parametrize("nx,ny,nz,distort", [(3, 3, 2, 0.0), (3, 3, 2, 0.6), (4, 3, 2, 0.4)])
def test_loop_basis_spans_kernel(nx, ny, nz, distort):
    hexes, elems = _build_bar(nx, ny, nz, distort=distort)
    cyc = _cycle_count(elems)
    N, dof, L, nLoop = _yano_matrix_and_loops(hexes)

    assert dof == 6 * len(elems)
    assert L.shape == (dof, nLoop), f"L shape {L.shape} != ({dof}, {nLoop})"
    # (2) topological cycle count
    assert nLoop == cyc, f"nLoop={nLoop} != cell-graph cycle count {cyc}"
    # (4) independent basis
    assert np.linalg.matrix_rank(L, tol=1e-9) == nLoop
    # (3) L spans ker(N) on the exact C++ operator -- the self-check
    rel = np.linalg.norm(N @ L) / (np.linalg.norm(N) * np.linalg.norm(L) + 1e-300)
    assert rel < 1e-9, f"||N L||/(||N|| ||L||) = {rel:.2e} (loop basis not in ker(N))"
    # cross-check: nLoop == numerical near-null dimension of N
    sv = np.linalg.svd(N, compute_uv=False)
    assert int(np.sum(sv < 1e-9 * sv[0])) == nLoop


def test_loop_basis_single_hex_has_no_loops():
    # one hex: no internal faces -> no loops
    hexes, _elems = _build_bar(1, 1, 1, distort=0.0)
    _N, dof, L, nLoop = _yano_matrix_and_loops(hexes)
    assert nLoop == 0
    assert L.shape == (dof, 0)


if __name__ == "__main__":
    for args in [(3, 3, 2, 0.0), (3, 3, 2, 0.6), (4, 3, 2, 0.4)]:
        test_loop_basis_spans_kernel(*args)
        print(f"OK {args}")
    test_loop_basis_single_hex_has_no_loops()
    print("OK single hex")
