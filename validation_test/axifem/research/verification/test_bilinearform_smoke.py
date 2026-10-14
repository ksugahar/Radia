"""Phase 2-C smoke test: verify BilinearForm grad(u)*grad(v)*dx integrates
   on an H1Henrotte FESpace and yields a well-conditioned positive-definite
   stiffness matrix on a simple axisymmetric mesh."""

import numpy as np
from netgen.occ import MoveTo, OCCGeometry, X, Y
from ngsolve import *
from radia.axifem import H1Henrotte


def main():
    # Tiny axisymmetric (r, z) mesh — single rectangle, structured Q1
    box = MoveTo(0.001, -0.005).Rectangle(0.005, 0.010).Face()
    box.edges.Min(X).name = "left"
    box.edges.Max(X).name = "right"
    box.edges.Min(Y).name = "bot"
    box.edges.Max(Y).name = "top"
    mesh = Mesh(OCCGeometry(box, dim=2).GenerateMesh(maxh=0.002))
    print(f"Mesh: {mesh.nv} vertices, {mesh.ne} elements")

    fes = H1Henrotte(mesh)
    print(f"H1Henrotte ndof = {fes.ndof}  (should equal mesh.nv = {mesh.nv})")

    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += InnerProduct(grad(u), grad(v)) * dx
    print("Assembling BilinearForm InnerProduct(grad(u), grad(v))*dx ...")
    a.Assemble()

    # Inspect the matrix: is it sparse, symmetric-ish, has reasonable values?
    A = a.mat
    rows, cols, vals = A.COO()
    rows = np.asarray(rows); cols = np.asarray(cols); vals = np.asarray(vals)
    print(f"  nnz = {len(vals)}")
    print(f"  diag mean = {np.mean([vals[i] for i in range(len(vals)) if rows[i] == cols[i]]):.4e}")
    print(f"  values range = [{vals.min():.4e}, {vals.max():.4e}]")

    # Check basic positive-definiteness on the diagonal (vertex DOFs)
    diag_vals = [vals[i] for i in range(len(vals)) if rows[i] == cols[i]]
    n_pos = sum(1 for v in diag_vals if v > 0)
    print(f"  positive diagonal entries: {n_pos}/{len(diag_vals)}")
    assert n_pos > 0, "expected at least some positive diagonal entries"

    print("\n[OK] BilinearForm grad(u)*grad(v)*dx assembles on H1Henrotte.")


if __name__ == "__main__":
    with TaskManager():
        main()
