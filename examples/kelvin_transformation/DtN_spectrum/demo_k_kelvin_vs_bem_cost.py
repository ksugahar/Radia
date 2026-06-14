# -*- coding: utf-8 -*-
# DEMO (k) (verified): Kelvin closure = SPARSE volume FE vs BEM exterior DtN = DENSE
# surface matrix. Both discretise the same exterior operator Lambda_ext on Gamma, but
# the matrix structure and assembly cost are opposite.
#
# Measured (this machine; BEM densifies V^{-1}(-1/2 M + K) column-by-column, O(ndof^2)):
#   Kelvin (sparse vol)  ndof  278 -> 1813 : nnz/row ~20 (fill->0), assemble 0.003 -> 0.020 s
#   BEM    (dense surf)  ndof  336 ->  456 : 100% dense (nnz=ndof^2), assemble 131 -> 191 s
# -> ~four orders of magnitude in assembly time at comparable ndof. Kelvin pays MORE DoF
# (the volume of the inverted ball) but each is cheap and sparse; BEM pays FEWER DoF
# (surface only) but a fully dense O(N^2) matrix that needs FMM/H-matrix to scale.
# Condensing the Kelvin ball to a DtN ON Gamma would re-create the dense N_Gamma^2 block,
# so the point is to KEEP the sparse volume rather than form the dense Schur complement.
#
# This script runs the Kelvin side fully and a SMALL BEM case (the dense assembly is slow
# on purpose -- that is the result); raise the BEM maxh down to reproduce the table above.
import os
import time
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import ngsolve as ng
from netgen.occ import Sphere, Pnt, OCCGeometry
from radia_mcp.radia_ngsolve.bem_integral import laplace_exterior_dtn

R = 1.0
x, y, z = ng.x, ng.y, ng.z
w = R * R / (x * x + y * y + z * z)


def kelvin_matrix(maxh, order=2):
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=maxh)).Curve(min(order + 1, 3))
    fes = ng.H1(mesh, order=order, dirichlet=".*")
    u, v = fes.TnT()
    t0 = time.time()
    a = ng.BilinearForm(w * ng.grad(u) * ng.grad(v) * ng.dx(bonus_intorder=6)); a.Assemble()
    t = time.time() - t0
    try:
        nze = int(a.mat.nze)
    except Exception:
        nze = len(a.mat.COO()[2])
    return fes.ndof, nze, t


def bem_matrix(maxh, order=1):
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=maxh)).Curve(4)
    fesS = ng.SurfaceL2(mesh, order=order, dual_mapping=True)
    t0 = time.time()
    Lam = laplace_exterior_dtn(fesS, intorder=10)
    t = time.time() - t0
    return fesS.ndof, int(np.count_nonzero(np.abs(Lam) > 1e-14)), t


print("Kelvin closure = SPARSE volume FE on the inverted ball:")
print("  {:>6} | {:>7} | {:>9} {:>9} {:>11}".format("maxh", "ndof", "nnz", "nnz/ndof", "assemble s"))
for mh in (0.5, 0.35, 0.25):
    nd, nnz, t = kelvin_matrix(mh)
    print("  {:>6.2f} | {:>7d} | {:>9d} {:>9.1f} {:>11.3f}".format(mh, nd, nnz, nnz / nd, t))

print("\nBEM exterior DtN = DENSE surface matrix (small case; finer = O(N^2) slower):")
print("  {:>6} | {:>7} | {:>9} {:>9} {:>11}".format("maxh", "ndof", "nnz", "fill %", "assemble s"))
for mh in (0.9,):
    nd, nnz, t = bem_matrix(mh)
    print("  {:>6.2f} | {:>7d} | {:>9d} {:>9.1f} {:>11.3f}".format(mh, nd, nnz, 100.0 * nnz / (nd * nd), t))
print("\n-> Kelvin: ~20 nnz/row, ms assembly. BEM: 100% dense, O(N^2) assembly (seconds->minutes).")
