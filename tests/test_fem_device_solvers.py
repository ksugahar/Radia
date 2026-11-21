"""Device AMG, AMS and BDDC solves against sparse Cholesky on identical inputs.

Guards: H1 AMG with CG, real and complex curl-curl with AMS, the complex SIBC form with a boundary
term, and BDDC at p=2,3 for H1 and HCurl, with and without a boundary element group, each reach the
sparse Cholesky solution of the same form, with the Dirichlet lift kept and only permitted payload
families crossing. A wrong Galerkin product, gradient, interpolation, harmonic extension or Schur
complement breaks convergence or the agreement; a wrong true residual raises. The coarse kernel of a
gradient space or a pure-Neumann block is pinned on its stiffest row. The AMS gradient and interpolations,
and their transposes, equal the tables of the CPU route: a defect there barely moves the iterations.
Assumes: the CPU direct solve is the reference on identical NGSolve inputs; this is parity of one
discretization, not accuracy. SOLUTION_REL, ITERATION_LIMIT, H1_AMG_ITERATION_LIMIT and TRANSFER_REL are
PROPOSED detectors, not adopted criteria. Optional GPU lane: skips visibly unless RADIA_FEM_NATIVE_TEST_DIR is set.
"""
import os

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")

SOLUTION_REL = 1e-6
ITERATION_LIMIT = 120
H1_AMG_ITERATION_LIMIT = 25
TRANSFER_REL = 1e-13
TOLERANCE = 1e-10


@pytest.fixture(scope="module")
def route():
    if not os.environ.get("RADIA_FEM_NATIVE_TEST_DIR"):
        pytest.skip("Optional GPU lane requires RADIA_FEM_NATIVE_TEST_DIR")
    import radia.fem_gpu as module
    return module


def _cube(n, hexes=False):
    from ngsolve.meshes import MakeStructured3DMesh
    return MakeStructured3DMesh(hexes=hexes, nx=n, ny=n, nz=n,
                                mapping=lambda x, y, z: (x + 0.1 * y, y + 0.1 * z, z))


def _rhs(fes, kind):
    v = fes.TestFunction()
    form = ng.LinearForm(fes)
    if kind == "H1":
        form += (1 + ng.x * ng.y) * v * ng.dx
    else:
        form += ng.InnerProduct(ng.CF((1 + ng.y, ng.z, ng.x * ng.y)), v) * ng.dx
    with ng.TaskManager():
        form.Assemble()
    return np.array(form.vec.FV().NumPy(), copy=True)


def _curl_form_rhs(fes):
    # The production surface-impedance right-hand side is a curl form, free of gradients. With a gauge
    # mass far below the curl term, a gradient part would make the field huge and leave even sparse
    # Cholesky above TOLERANCE (true residual 2.9e-10 at mass 1e-3, 4.9e-8 at mass 1e-6).
    form = ng.LinearForm(fes)
    form += ng.InnerProduct(ng.CF((1 + ng.y, ng.z, ng.x * ng.y)), ng.curl(fes.TestFunction())) * ng.dx
    with ng.TaskManager():
        form.Assemble()
    return np.array(form.vec.FV().NumPy(), copy=True)


def _curl_difference(fes, actual, expected):
    # The gauge mass leaves the gradient part of the field weakly determined, so compare the curl.
    field, reference = ng.GridFunction(fes), ng.GridFunction(fes)
    field.vec.FV().NumPy()[:] = actual
    reference.vec.FV().NumPy()[:] = expected
    difference = ng.curl(field) - ng.curl(reference)
    denominator = ng.Integrate(ng.InnerProduct(ng.curl(reference), ng.Conj(ng.curl(reference))), fes.mesh).real
    return float(np.sqrt(ng.Integrate(ng.InnerProduct(difference, ng.Conj(difference)), fes.mesh).real / denominator))


def _check(route, solve, fes, coefficients, rhs, *, lift=None, bonus=(0, 0, 0), cycles=None, compare_curl=False, iteration_limit=ITERATION_LIMIT):
    with ng.TaskManager():
        expected = route.cpu_solve(fes, coefficients, rhs, lift=lift, preconditioner="direct", bonus=bonus)
    actual = solve(fes, coefficients, rhs, Backend="gpu", lift=lift, bonus=bonus, tolerance=TOLERANCE, cycles=cycles)
    assert actual["backend"] == "gpu"
    print("device iterations", actual["iterations"], "hierarchy", actual["hierarchy"])
    if compare_curl:
        difference = _curl_difference(fes, actual["solution"], expected["solution"])
    else:
        difference = float(np.linalg.norm(actual["solution"] - expected["solution"]) / np.linalg.norm(expected["solution"]))
    assert difference <= SOLUTION_REL, f"solution differs by {difference}, limit {SOLUTION_REL}; iterations {actual['iterations']}"
    assert 0 < actual["iterations"] <= iteration_limit, f"iterations {actual['iterations']}, limit {iteration_limit}"
    assert actual["true_residual"] <= 1e-6, f"true residual {actual['true_residual']}"
    names = {row["payload"] for row in actual["ledger"]}
    assert not any(name.startswith(("matrix", "operator")) for name in names), names
    assert not any(row["direction"] == "D2H" and row["payload"].startswith(("mesh.", "table.", "coef.", "rhs.", "lift."))
                   for row in actual["ledger"]), actual["ledger"]
    return actual


def test_h1_amg_with_lift(route):
    mesh = _cube(10)
    fes = ng.H1(mesh, order=1, dirichlet="left|right")
    coefficients = route.Coefficients(stiff=2.0, mass=0.5, stiff_aux=2.0, mass_aux=0.5)
    lift = ng.GridFunction(fes)
    lift.Set(1 + ng.x, definedon=mesh.Boundaries("left|right"))
    result = _check(route, route.solve_h1, fes, coefficients, _rhs(fes, "H1"), lift=np.array(lift.vec.FV().NumPy(), copy=True),
                    iteration_limit=H1_AMG_ITERATION_LIMIT)
    assert len(result["hierarchy"][0]["rows"]) >= 2, result["hierarchy"]


def _check_pure_neumann(route, fes, stiffness, min_coarse):
    # Without a Dirichlet boundary the constants are the kernel on every level, so sparse Cholesky has
    # no solution to offer: the reference is a dense solve with the stiffest row pinned, up to a constant.
    coefficients = route.Coefficients(stiff=stiffness, stiff_aux=stiffness)
    rhs = _rhs(fes, "H1")
    rhs -= rhs.mean()
    form = route.cpu_form(fes, coefficients)
    with ng.TaskManager():
        form.Assemble()
    rows, columns, values = form.mat.COO()
    dense = np.zeros((fes.ndof, fes.ndof))
    dense[np.array(rows), np.array(columns)] = np.array(values)
    kept = np.delete(np.arange(fes.ndof), np.argmax(np.diag(dense)))
    expected = np.zeros(fes.ndof)
    expected[kept] = np.linalg.solve(dense[np.ix_(kept, kept)], rhs[kept])
    system = route.DeviceSystem(fes, coefficients)
    hierarchy = system.setup("amg", min_coarse=min_coarse)
    assert hierarchy[0]["pinned"] == 1, hierarchy[0]
    assert hierarchy[0]["rows"][-1] > 1, hierarchy[0]
    actual = system.solve(rhs, tolerance=TOLERANCE, max_iterations=ITERATION_LIMIT)
    print("device iterations", actual["iterations"], "hierarchy", hierarchy[0])
    difference = actual["solution"] - actual["solution"].mean() - (expected - expected.mean())
    relative = float(np.linalg.norm(difference) / np.linalg.norm(expected - expected.mean()))
    assert relative <= SOLUTION_REL, f"solution differs by {relative} up to a constant, limit {SOLUTION_REL}; iterations {actual['iterations']}"


@pytest.mark.parametrize("min_coarse", [50, 10])
def test_h1_amg_pure_neumann_block_pins_the_kernel(route, min_coarse):
    # Coarsened to a single row, the coarse level would hold only a rounding value.
    _check_pure_neumann(route, ng.H1(_cube(4), order=1), 1.0, min_coarse)


@pytest.mark.parametrize("regions, maxh", [(("stiff", "soft"), 0.25), (("soft", "stiff"), 0.3)])
def test_h1_amg_kernel_seed_sits_on_the_stiffest_row(route, regions, maxh):
    # With the kernel left to the last pivot, the first numbering stalls CG; with the seed on the lowest
    # row, the second does. In each case the coarse inverse is the only part that differs.
    from netgen.occ import Box, Glue, OCCGeometry, Pnt

    boxes = []
    for position, name in enumerate(regions):
        box = Box(Pnt(position, 0, 0), Pnt(position + 1, 1, 1))
        box.mat(name)
        boxes.append(box)
    mesh = ng.Mesh(OCCGeometry(Glue(boxes)).GenerateMesh(maxh=maxh))
    _check_pure_neumann(route, ng.H1(mesh, order=1), {"stiff": 1e9, "soft": 1.0}, 50)


def test_real_curl_curl_ams(route):
    mesh = _cube(8)
    fes = ng.HCurl(mesh, order=1, nograds=True, dirichlet=".*")
    coefficients = route.Coefficients(stiff=1.0, mass=10.0, stiff_aux=1.0, mass_aux=10.0)
    result = _check(route, route.solve_hcurl, fes, coefficients, _rhs(fes, "HCurl"))
    assert len(result["hierarchy"]) == 4, result["hierarchy"]


def _csr(matrix):
    import scipy.sparse as sp

    rows, columns, values = matrix.COO()
    return sp.csr_matrix((np.array(values), (np.array(rows), np.array(columns))), shape=(matrix.height, matrix.width))


def test_ams_transfers_equal_the_cpu_tables_on_identical_inputs(route):
    # The oracles are the gradient and interpolation tables of the CPU route and the VectorH1 to HCurl table
    # of NGSolve. Only the probes cross. Spaces of different extents keep the three directions apart.
    import radia.sparsesolv_ngsolve as ssn
    import scipy.sparse as sp
    from netgen.occ import Box, OCCGeometry, Pnt

    mesh = ng.Mesh(OCCGeometry(Box(Pnt(0, 0, 0), Pnt(1, 0.7, 0.4))).GenerateMesh(maxh=0.2))
    fes = ng.HCurl(mesh, order=1, nograds=True, dirichlet=".*")
    coefficients = route.Coefficients(stiff=1.0, mass=10.0, stiff_aux=1.0, mass_aux=10.0)
    system = route.DeviceSystem(fes, coefficients)
    system.setup("ams")
    tables = route._cpu_ams(fes, coefficients, (0, 0, 0), ssn, fes.ndof)._accelerator_snapshot()["spaces"]
    assert len(tables) == 4, f"{len(tables)} CPU spaces, expected the gradient and three interpolations"
    gradient = _csr(ng.HCurl(mesh, order=1, nograds=True).CreateGradient()[0])
    vector_table = _csr(ng.ConvertOperator(ng.VectorH1(mesh, order=1), ng.HCurl(mesh, order=1, nograds=True)))
    columns = np.split(np.arange(vector_table.shape[1]), 3)
    independent = [gradient] + [vector_table[:, block] for block in columns]
    probes = np.random.default_rng(20261105)
    nodal, edge = probes.normal(size=mesh.nv), probes.normal(size=mesh.nedge)
    report = []
    for space, table in enumerate(tables):
        P, Pt = (sp.csr_matrix((t["values"], t["columns"], t["rowptr"]), shape=(t["rows"], t["cols"])) for t in (table["P"], table["Pt"]))
        for transpose, probe, cpu_table, other in ((False, nodal, P, independent[space]), (True, edge, Pt, independent[space].T)):
            actual = system.ams_transfer(space, probe, transpose=transpose)
            for name, expected in (("CPU route", cpu_table @ probe), ("NGSolve", other @ probe)):
                error = float(np.abs(actual - expected).max() / np.abs(expected).max())
                report.append(error)
                assert error <= TRANSFER_REL, (
                    f"space {space} transpose {transpose}: {error} against the {name} table, limit {TRANSFER_REL}; "
                    f"actual {actual[:4]}, expected {expected[:4]}")
    print("device transfer errors, largest", max(report))
    assert {row["payload"].split(".")[0] for row in system.ledger()} <= set(name.rstrip(".") for name in route.PERMITTED_PAYLOADS)


def test_complex_curl_curl_with_sibc_term(route):
    mesh = _cube(8)
    fes = ng.HCurl(mesh, order=1, nograds=True, complex=True, dirichlet="left|right")
    robin = 40.0
    coefficients = route.Coefficients(stiff=1.0, mass=1e-3, stiff_aux=1.0, mass_aux=1e-3, trace_im={"back": robin},
                                      trace_aux={"back": robin})
    _check(route, route.solve_hcurl, fes, coefficients, _curl_form_rhs(fes) + 0j, bonus=(4, 0, 0), compare_curl=True)


@pytest.mark.parametrize("order", [2, 3])
def test_h1_bddc(route, order):
    mesh = _cube(4)
    fes = ng.H1(mesh, order=order, dirichlet="left", wb_withedges=False)
    coefficients = route.Coefficients(stiff=3.0, stiff_aux=3.0)
    _check(route, route.solve_h1, fes, coefficients, _rhs(fes, "H1"))


@pytest.mark.parametrize("order", [2, 3])
def test_complex_hcurl_bddc(route, order):
    mesh = _cube(3)
    fes = ng.HCurl(mesh, order=order, complex=True, dirichlet=".*")
    coefficients = route.Coefficients(stiff=1.0, mass_im=30.0, stiff_aux=1.0, mass_aux=30.0)
    _check(route, route.solve_hcurl, fes, coefficients, _rhs(fes, "HCurl") + 0j, bonus=(2, 2, 2), cycles=3)


@pytest.mark.parametrize("order", [2, 3])
def test_complex_hcurl_bddc_with_sibc_term_and_gauge_mass(route, order):
    # The production surface-impedance form: curl-curl, a mass far below the curl term and a
    # complex boundary term, so the Schur complement carries a boundary element group.
    mesh = _cube(4)
    fes = ng.HCurl(mesh, order=order, nograds=True, complex=True, dirichlet="left|right")
    robin = 40.0 * (1 + 1j)
    coefficients = route.Coefficients(stiff=1.0, mass=1e-6, stiff_aux=1.0, mass_aux=1e-6, trace={"back": robin.real},
                                      trace_im={"back": robin.imag}, trace_aux={"back": abs(robin)})
    _check(route, route.solve_hcurl, fes, coefficients, _curl_form_rhs(fes) + 0j, bonus=(4, 0, 0), cycles=3, compare_curl=True)


def test_gradient_space_kernel_is_seeded_on_the_coarse_level(route):
    # The constants of the production surface-impedance form with a right-hand side free of gradients.
    # Left to the last pivot of the coarse elimination, the constant kernel of the gradient space
    # leaves a rounding value that stalls COCR.
    mesh = _cube(10)
    fes = ng.HCurl(mesh, order=2, nograds=True, complex=True, dirichlet="left|right")
    robin = 3.5355e8 * (1 + 1j)
    stiffness = 795774.7
    coefficients = route.Coefficients(stiff=stiffness, mass=0.796, stiff_aux=stiffness, mass_aux=0.796, trace={"back": robin.real},
                                      trace_im={"back": robin.imag}, trace_aux={"back": abs(robin)})
    rhs = _curl_form_rhs(fes)
    with ng.TaskManager():
        expected = route.cpu_solve(fes, coefficients, rhs, preconditioner="direct", bonus=(4, 0, 0))
    system = route.DeviceSystem(fes, coefficients, bonus=(4, 0, 0))
    hierarchy = system.setup("bddc", cycles=3, min_coarse=50)
    assert hierarchy[0]["pinned"] == 1, hierarchy[0]
    actual = system.solve(rhs, tolerance=TOLERANCE, max_iterations=ITERATION_LIMIT)
    print("device iterations", actual["iterations"], "hierarchy", hierarchy[0])
    relative = _curl_difference(fes, actual["solution"], expected["solution"])
    assert relative <= SOLUTION_REL, f"the curl differs by {relative}, limit {SOLUTION_REL}; iterations {actual['iterations']}"


def _other_mesh(name):
    if name == "hexes":
        return _cube(3, hexes=True)
    from validation_test.hmatrix_cuda.hdiv_mixed_meshes import mixed_grid, pyramid_cube
    return {"pyramids": lambda: pyramid_cube(0.4), "mixed": lambda: mixed_grid((2, 2, 2))}[name]()


@pytest.mark.parametrize("mesh_name", ["hexes", "pyramids", "mixed"])
@pytest.mark.parametrize("order", [1, 2])
def test_h1_on_other_elements(route, mesh_name, order):
    mesh = _other_mesh(mesh_name)
    fes = ng.H1(mesh, order=order, dirichlet=".*", wb_withedges=False)
    actual = _check(route, route.solve_h1, fes, route.Coefficients(stiff=3.0, stiff_aux=3.0), _rhs(fes, "H1"))
    assert actual["preconditioner"] == ("amg" if order == 1 else "bddc"), actual["preconditioner"]


@pytest.mark.parametrize("mesh_name", ["hexes", "pyramids", "mixed"])
@pytest.mark.parametrize("order", [1, 2])
def test_complex_hcurl_on_other_elements(route, mesh_name, order):
    mesh = _other_mesh(mesh_name)
    fes = ng.HCurl(mesh, order=order, nograds=True, complex=True, dirichlet=".*")
    coefficients = route.Coefficients(stiff=1.0, mass_im=30.0, stiff_aux=1.0, mass_aux=30.0)
    actual = _check(route, route.solve_hcurl, fes, coefficients, _rhs(fes, "HCurl") + 0j, bonus=(2, 2, 2), cycles=3)
    assert actual["preconditioner"] == "bddc", actual["preconditioner"]


def test_real_right_hand_side_is_complex_checked(route):
    mesh = _cube(2)
    fes = ng.H1(mesh, order=1, dirichlet="left")
    coefficients = route.Coefficients(stiff=1.0, stiff_aux=1.0)
    system = route.DeviceSystem(fes, coefficients)
    system.setup("amg")
    with pytest.raises(ValueError, match="complex but the system is real"):
        system.solve(np.ones(fes.ndof) * 1j)
    with pytest.raises(ValueError, match="must be a finite vector"):
        system.solve(np.ones(fes.ndof + 1))
    with pytest.raises(RuntimeError, match="right-hand side is zero"):
        system.solve(np.zeros(fes.ndof))


def test_true_residual_limit_raises(route, monkeypatch):
    mesh = _cube(4)
    fes = ng.H1(mesh, order=1, dirichlet="left")
    coefficients = route.Coefficients(stiff=1.0, stiff_aux=1.0)
    monkeypatch.setattr(route, "RELATIVE_LIMIT", 1e-30)
    with pytest.raises(RuntimeError, match="true relative residual"):
        route.solve_h1(fes, coefficients, _rhs(fes, "H1"), Backend="gpu", tolerance=1e-6)


def test_bddc_rejects_lowest_order(route):
    mesh = _cube(2)
    fes = ng.H1(mesh, order=1, dirichlet="left")
    system = route.DeviceSystem(fes, route.Coefficients(stiff=1.0, stiff_aux=1.0))
    with pytest.raises(ValueError, match="BDDC applies to p>1"):
        system.setup("bddc")
