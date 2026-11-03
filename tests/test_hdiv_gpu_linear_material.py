"""Device linear material solve against the CPU-assembled equation.

Covers conforming RT0 (order 0), BDM1 and BDM2, with ACA and the explicit
exact dense device Gram. Guards: ``_ResidentHDiv.solve_linear`` with a ``_ResidentMass`` solves
(B^T G B + inv_chi M) x = f with the device-built charge Gram, charge map and
reference mass, for every mass preconditioner: the solution matches a dense
solve of the same equation assembled from the CPU charge map, the CPU exact
dense normalized Gram and the NGSolve mass. Dropping or rescaling the
inv_chi M term, or a wrong device mass, moves the solution far beyond the
detector. Chebyshev settings without a Chebyshev mode are rejected.
Assumes: the CPU pieces are the reference on identical staged inputs; this is
parity of one discretization, not independent accuracy. LINEAR_REL is a
PROPOSED detector, not an adopted criterion. Optional GPU lane: skips visibly
unless RADIA_HDIV_NATIVE_TEST_DIR is set; a configured lane loads the binding
through the provenance-checked loader and fails if it cannot.
"""

import os

import numpy as np
import pytest

LINEAR_REL = 1e-8
INV_CHI = 0.1
PRECONDITIONERS = ("mass-diagonal", "mass-cholesky", "mass-polynomial", "mass-polynomial-graph",
                   "mass-chebyshev", "mass-chebyshev-graph")


@pytest.fixture(scope="module")
def cuda():
    if not os.environ.get("RADIA_HDIV_NATIVE_TEST_DIR"):
        pytest.skip("Optional GPU lane requires RADIA_HDIV_NATIVE_TEST_DIR")
    from radia.gpu import _native
    return _native("_hdiv_cuda")


def _case(tmp_path, order, cuda, compression="aca"):
    import ngsolve as ng
    from scipy.sparse import coo_matrix
    from netgen.csg import unit_cube
    import radia._radia_pybind as rp
    from radia import hdiv_gpu
    from radia.vim._vim import _prepare_gpu_charge_inputs
    path = tmp_path / f"cube-{order}.vol"
    unit_cube.GenerateMesh(maxh=.5).Save(str(path))
    with ng.TaskManager():
        mesh = ng.Mesh(str(path))
        fes = ng.HDiv(mesh, order=order)
        setup = hdiv_gpu._prepare_device_components(path, order, cuda, aca=None, compression=compression)
        B, inputs = _prepare_gpu_charge_inputs(fes, backend="gpu")
        gram = rp._ChargeGramHMatrix(**inputs, build=False)
        gram.build_exact_dense_normalized_gram(maximum_memory_mb=256)
        u, v = fes.TnT()
        form = ng.BilinearForm(fes, symmetric=True)
        form += ng.InnerProduct(u, v)*ng.dx
        form.Assemble()
        rows, columns, values = form.mat.COO()
        mass = coo_matrix((values, (rows, columns)), shape=(fes.ndof, fes.ndof)).tocsr()
    nd = fes.ndof
    columns = [B.T @ np.asarray(gram.matvec_sym(np.ascontiguousarray(B[:, k].toarray().ravel())))
               for k in range(nd)]
    system = np.column_stack(columns) + INV_CHI*mass.toarray()
    rhs = setup["rhs_operator"] @ np.array([0., 0., 1e5])
    return setup, np.linalg.solve(system, rhs), rhs


@pytest.mark.parametrize("order,compression", [(0, "aca"), (0, "dense"), (1, "aca"), (2, "aca"), (2, "dense")])
def test_device_linear_solve_matches_cpu_assembled_equation(tmp_path, request, cuda, order, compression):
    setup, expected, rhs = _case(tmp_path, order, cuda, compression)
    nd = setup["dofs"]
    mass = cuda._ResidentMass.from_flat_tables(**setup["evaluation_tables"], dofs=nd, backend="gpu")
    operator = cuda._ResidentHDiv.from_flat_charges(**setup["charge_gram"], dofs=nd, backend="gpu")
    try:
        worst = {}
        for preconditioner in PRECONDITIONERS:
            result = operator.solve_linear(mass, np.ascontiguousarray(rhs), INV_CHI, 1e-11, 2000, preconditioner)
            relative = np.linalg.norm(result["x"] - expected)/np.linalg.norm(expected)
            worst[preconditioner] = float(relative)
            assert result["relative_residual"] <= 1e-11 and result["preconditioner"] == preconditioner
            assert (result["preconditioner_graph_launches"] > 0) == preconditioner.endswith("-graph")
            assert relative <= LINEAR_REL, (
                f"order={order} {compression} preconditioner={preconditioner} relative={relative:.3e} limit={LINEAR_REL:.0e} "
                f"iterations={result['iterations']}")
        assert operator.assembly()["compression"] == compression
        request.node.user_properties += [("dofs", nd), ("max_relative", max(worst.values()))]
        with pytest.raises(ValueError, match="Chebyshev"):
            operator.solve_linear(mass, np.ascontiguousarray(rhs), INV_CHI, 1e-11, 2000, "mass-diagonal", degree=8)
        with pytest.raises(ValueError, match="preconditioner"):
            operator.solve_linear(mass, np.ascontiguousarray(rhs), INV_CHI, 1e-11, 2000, "reference-diagonal")
        mass.close()
        with pytest.raises((ValueError, RuntimeError), match="closed"):
            operator.solve_linear(mass, np.ascontiguousarray(rhs), INV_CHI, 1e-11, 2000, "mass-diagonal")
    finally:
        mass.close()
        operator.close()
