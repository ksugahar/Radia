"""C++ PRIMA reduction of the fin volume-lattice PEEC port model."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

FIXTURES = Path(__file__).resolve().parent / "coil_from_cad" / "fixtures"
F_OP = 222e3


@pytest.fixture(scope="module")
def small_fin():
    """84-branch lattice; box_order 4 keeps setup short (the reduction, not
    the kernel, is under test)."""
    from radia.fin_volume import assemble_fin_volume

    _sweep, lat, _builder, topo, _energy = assemble_fin_volume(
        FIXTURES / "beak_fin_short.step", sigma=5.8e7, n_layers=1, n_lanes=8,
        n_stations=5, lane_grading="uniform", box_order=4)
    return lat, topo


class _DenseL:
    """apply_L shim so the Python PRIMA reference runs on the identical L."""

    def __init__(self, L):
        self.L = np.asarray(L, dtype=float)
        self.N = self.L.shape[0]

    def apply_L(self, x):
        return self.L @ x


def _python_prima(topo, q):
    import scipy.sparse as sp
    from radia.prima_hacapk import _descriptor_prima

    br = np.asarray(topo["segment_nodes"], dtype=np.int64)
    n, n_nodes = len(br), int(topo["n_nodes"])
    plus, minus = int(topo["ports"][0][0]), int(topo["ports"][0][1])
    A = sp.csr_matrix((np.r_[-np.ones(n), np.ones(n)],
                       (np.r_[br[:, 0], br[:, 1]], np.r_[np.arange(n), np.arange(n)])),
                      shape=(n_nodes, n))
    keep = np.array([k for k in range(n_nodes) if k != minus])
    e = np.zeros(n_nodes)
    e[plus], e[minus] = 1.0, -1.0
    R = np.asarray(topo["R"], dtype=float)
    return _descriptor_prima(_DenseL(topo["L"]), R, A[keep, :], e[keep], q, inv_R=1.0 / R)


def test_cpp_prima_matches_python_reference_on_identical_inputs(small_fin):
    """Guards: the C++ block Arnoldi spans the same Krylov space as the
    verified Python descriptor PRIMA and projects the same reduced
    E, G, B (compared through the subspace angle and the reduced port
    impedance over DC..1 MHz).
    Assumes: one port, expansion at s0 = 0 (the Python reference's point)."""
    from radia.fin_prima import reduce_topology

    _lat, topo = small_fin
    q = 8
    model, info = reduce_topology(topo, q)
    V_py, E_py, G_py, B_py, q_py = _python_prima(topo, q)
    assert info["q"] == q_py == q
    assert info["max_solve_residual"] <= 1e-10, info

    from radia.peec_matrices import prima_reduce
    raw = prima_reduce(np.asarray(topo["L"]), np.asarray(topo["R"]),
                       np.asarray(topo["segment_nodes"], dtype=np.int32), int(topo["n_nodes"]),
                       [topo["ports"][0][0]], [topo["ports"][0][1]], q)
    cosines = np.linalg.svd(raw["V"].T @ V_py, compute_uv=False)
    assert cosines.min() > 1 - 1e-8, cosines
    for f in (0.0, 1e3, F_OP, 1e6):
        s = 2j * math.pi * f
        z_py = complex(B_py @ np.linalg.solve(G_py + s * E_py, B_py.astype(complex)))
        z_cpp = complex(model.response(s)[0, 0])
        assert abs(z_cpp - z_py) <= 1e-8 * abs(z_py), (f, z_cpp, z_py)


def test_full_reference_matches_the_verified_circuit_solver(small_fin):
    """Guards: the numpy unreduced MNA used for order selection gives the
    same port impedance as the verified C++ circuit solver.
    Assumes: the circuit solver uses DC branch resistance (Zs=None)."""
    from radia.fin_prima import full_port_impedance
    from radia.peec_topology import PEECCircuitSolver

    _lat, topo = small_fin
    z = full_port_impedance(topo["L"], topo["R"], topo["segment_nodes"], topo["n_nodes"],
                            [topo["ports"][0][:2]], F_OP)[0, 0]
    z_ref = complex(PEECCircuitSolver(topo).compute_port_impedance(F_OP))
    assert abs(z - z_ref) <= 1e-9 * abs(z_ref), (z, z_ref)


def test_reduced_matches_unreduced_and_is_passive(small_fin):
    """Guards: the PRIMA port impedance converges to the unreduced one at
    222 kHz as the order grows (exact once the Krylov space is complete),
    and the congruence keeps E >= 0 and G + G^T >= 0 (passivity).
    Assumes: small lattice, dense full solve affordable."""
    from radia.fin_prima import order_error_table, reduce_topology, select_order

    _lat, topo = small_fin
    full = len(topo["R"]) + int(topo["n_nodes"]) - 1
    table = order_error_table(topo, [1, 2, 4, 8, 16, full], F_OP)
    err = [row["rel_error"] for row in table["rows"]]
    assert err[-1] <= 1e-9, table
    assert err[3] <= 5e-3 * err[0], err  # q=8; user-approved 2026-11-21 (was 1e-3)
    assert select_order(table, 1e-6)["q"] <= 16, err
    with pytest.raises(ValueError, match="no PRIMA order"):
        select_order({"frequency_Hz": F_OP, "rows": table["rows"][:1]}, 1e-12)
    with pytest.raises(ValueError, match="explicit tolerance"):
        select_order(table, None)

    model, _ = reduce_topology(topo, 8)
    G = -model.A
    for sym in (model.E, G + G.T):
        ev = np.linalg.eigvalsh(0.5 * (sym + sym.T))
        assert ev.min() >= -1e-12 * abs(ev).max(), ev


def _ladder(n_rungs=4):
    """Two rails with rungs: a small multi-port graph with exact R and a
    positive definite L (independent of the partial-element kernel)."""
    branches = []
    for k in range(n_rungs):
        branches += [(2 * k, 2 * k + 2), (2 * k + 1, 2 * k + 3), (2 * k + 2, 2 * k + 3)]
    branches = np.array(branches, dtype=np.int32)
    n = len(branches)
    rng = np.random.default_rng(7)
    X = rng.standard_normal((n, n))
    L = 1e-9 * (X @ X.T / n + np.eye(n))
    L = 0.5 * (L + L.T)
    R = 1e-3 * (1.0 + rng.random(n))
    return L, R, branches, 2 * n_rungs + 2


def test_block_arnoldi_two_ports_reaches_the_exact_impedance_matrix():
    """Guards: block (multi-port) start vectors, deflation and projection:
    a complete Krylov basis reproduces the full 2x2 impedance matrix.
    Assumes: ports (0, 1) and (2k+3, 1) share ground node 1."""
    from radia.fin_prima import full_port_impedance
    from radia.peec_matrices import prima_reduce

    L, R, br, n_nodes = _ladder()
    ports = [(0, 1), (n_nodes - 1, 1)]
    r = prima_reduce(L, R, br, n_nodes, [p[0] for p in ports], [p[1] for p in ports],
                     len(R) + n_nodes - 1)
    assert r["B"].shape[1] == 2 and r["q"] > 2, r["q"]
    s = 2j * math.pi * F_OP
    z_q = r["B"].T @ np.linalg.solve(r["G"] + s * r["E"], r["B"].astype(complex))
    z = full_port_impedance(L, R, br, n_nodes, ports, F_OP)
    assert np.linalg.norm(z_q - z) <= 1e-9 * np.linalg.norm(z), (z_q, z)
    r2 = prima_reduce(L, R, br, n_nodes, [p[0] for p in ports], [p[1] for p in ports], 2)
    assert r2["q"] == 2 and r2["B"].shape == (2, 2)
    with pytest.raises(ValueError, match="start-block rank"):
        prima_reduce(L, R, br, n_nodes, [p[0] for p in ports], [p[1] for p in ports], 1)


def test_prima_rejects_unsupported_input():
    """Guards: fail-loud input contract (no silent repair).
    Assumes: nothing beyond the C++ validation messages.  The parallel-edge
    case has a zero-incidence loop (1, -1) with negative energy that the
    reduced-E check cannot see; only the full-L admission rejects it."""
    from radia.peec_matrices import prima_reduce

    L, R, br, n_nodes = _ladder()
    ok = dict(L=L, R=R, branch_nodes=br, n_nodes=n_nodes, port_plus=[0], port_minus=[1], order=4)
    prima_reduce(**ok)
    bad_L = L.copy()
    bad_L[0, 1] += 1e-20
    cases = [
        (dict(L=bad_L), "symmetric"),
        (dict(R=np.where(np.arange(len(R)) == 2, 0.0, R)), "R\\[2\\]"),
        (dict(order=0), "order"),
        (dict(s0=-1.0), "s0"),
        (dict(port_plus=[1]), "port 0"),
        (dict(branch_nodes=np.vstack([br[:-1], [[0, 0]]]).astype(np.int32)), "self-loop"),
        (dict(n_nodes=n_nodes + 1), "Schur"),   # an isolated node: disconnected graph
    ]
    for change, message in cases:
        with pytest.raises((ValueError, RuntimeError), match=message):
            prima_reduce(**{**ok, **change})
    par = np.array([[1, 0], [1, 0]], dtype=np.int32)
    with pytest.raises(RuntimeError, match="L is not positive definite"):
        prima_reduce(np.array([[0.0, 1e-9], [1e-9, 0.0]]), np.full(2, 1e-3), par, 2,
                     [1], [0], 2)
    indefinite = L - 2.0 * np.linalg.eigvalsh(L).max() * np.eye(len(R))
    with pytest.raises(RuntimeError, match="positive definite|negative energy"):
        prima_reduce(**{**ok, "L": indefinite, "s0": 2 * math.pi * F_OP})


def test_solve_fin_volume_prima_matches_dense_and_gates_low_order():
    """Guards: the opt-in ``prima_order`` route of ``solve_fin_volume``
    recovers the dense branch currents from the PRIMA basis expanded at the
    operating frequency, records q and its residuals, and an order too low to
    pass the shared KCL/KVL true-residual gate raises instead of returning.
    Assumes: the 84-branch fixture lattice; at order 24 its Krylov space is
    complete enough to meet the 1e-6 residual limit."""
    from radia.fin_volume import solve_fin_volume

    kw = dict(frequency=F_OP, sigma=5.8e7, current_A=10.0, n_layers=1, n_lanes=8,
              n_stations=5, lane_grading="uniform", box_order=4)
    step = FIXTURES / "beak_fin_short.step"
    dense = solve_fin_volume(step, **kw)
    red = solve_fin_volume(step, prima_order=24, **kw)
    err = np.linalg.norm(red["current"] - dense["current"]) / np.linalg.norm(dense["current"])
    assert err < 1e-6, (err, red["prima"])
    assert red["prima"]["s0"] == pytest.approx(2 * math.pi * F_OP)
    assert red["coil_true_residual"]["kvl_relative"] <= 1e-6
    assert dense["prima"] is None
    with pytest.raises(RuntimeError, match="true relative residual"):
        solve_fin_volume(step, prima_order=1, **kw)
    with pytest.raises(ValueError, match="CPU-only"):
        solve_fin_volume(step, prima_order=4, backend="gpu", **kw)


def test_fin_prima_sweep_builds_once_and_gates_each_frequency(monkeypatch):
    """Guards: one assembly and basis, correct restored currents, and fail-loud gating."""
    from radia import fin_volume, peec_matrices

    step = FIXTURES / "beak_fin_short.step"
    points = (2 * F_OP, F_OP / 2, F_OP)
    kw = dict(sigma=5.8e7, current_A=10., n_layers=1, n_lanes=8,
              n_stations=5, lane_grading="uniform", box_order=4)
    dense = [fin_volume.solve_fin_volume(step, frequency=f, **kw) for f in points]
    assembled = fin_volume.assemble_fin_volume
    reduce = peec_matrices.prima_reduce
    calls = {"assembly": 0, "basis": 0}

    def count_assembly(*args, **kwargs):
        calls["assembly"] += 1
        return assembled(*args, **kwargs)

    def count_basis(*args, **kwargs):
        calls["basis"] += 1
        return reduce(*args, **kwargs)

    monkeypatch.setattr(fin_volume, "assemble_fin_volume", count_assembly)
    monkeypatch.setattr(peec_matrices, "prima_reduce", count_basis)
    sweep = fin_volume.solve_fin_volume_prima_sweep(
        step, frequencies=points, expansion_frequency=F_OP, prima_order=24, **kw)
    assert calls == {"assembly": 1, "basis": 1}
    assert len(sweep["points"]) == len(points)
    assert [p["frequency_hz"] for p in sweep["points"]] == list(points)
    for point, reference in zip(sweep["points"], dense):
        result = point["solution"]
        error = np.linalg.norm(result["I_fil"] - reference["I_fil"]) / np.linalg.norm(reference["I_fil"])
        assert error <= 1e-6, (point["frequency_hz"], error)
        assert max(result["coil_true_residual"][k] for k in ("kcl_relative", "kvl_relative")) <= 1e-6
        assert result["prima"]["s0"] == pytest.approx(2 * math.pi * F_OP)

    for invalid in (14.5, True, 0, -3, "14"):
        with pytest.raises(ValueError, match="positive integer"):
            fin_volume.solve_fin_volume_prima_sweep(
                step, frequencies=points, expansion_frequency=F_OP, prima_order=invalid, **kw)
    assert calls == {"assembly": 1, "basis": 1}
    with pytest.raises(RuntimeError, match="true relative residual"):
        fin_volume.solve_fin_volume_prima_sweep(
            step, frequencies=points, expansion_frequency=F_OP, prima_order=1, **kw)



def test_fin_prima_sweep_rejects_a_later_residual_failure(monkeypatch):
    """Guards: an invalid restored current at the second point aborts the sweep."""
    from radia import fin_volume

    points = (F_OP, 2 * F_OP, F_OP / 2)
    evaluate = fin_volume._evaluate_fin_prima
    evaluated = []

    def restore_with_second_point_failure(reduced, branch_count, frequency, current_A):
        current, info = evaluate(reduced, branch_count, frequency, current_A)
        evaluated.append(frequency)
        if frequency == points[1]:
            current = np.zeros_like(current)
        return current, info

    monkeypatch.setattr(fin_volume, "_evaluate_fin_prima", restore_with_second_point_failure)
    with pytest.raises(RuntimeError, match="true relative residual"):
        fin_volume.solve_fin_volume_prima_sweep(
            FIXTURES / "beak_fin_short.step", frequencies=points,
            expansion_frequency=F_OP, prima_order=24, sigma=5.8e7, current_A=10.0,
            n_layers=1, n_lanes=8, n_stations=5, lane_grading="uniform", box_order=4)
    assert evaluated == list(points[:2])



def test_uniform_field_emf_obeys_faraday_and_terminal_gauge():
    """Guards: directed flux, three field axes and collapsed-terminal closure."""
    from types import SimpleNamespace

    from radia.fin_prima import _uniform_field_emf_columns

    # Orthogonal loop edges of length 2 m and 3 m.
    corners = np.array([[0., 0., 0.], [1.2, 1.6, 0.], [1.2, 1.6, 3.], [0., 0., 3.]])
    branches = np.array([[0, 1], [1, 2], [2, 3], [3, 0]])
    bars = corners[branches]
    nodes = corners.copy()
    nodes[0] = [0.2, 0.3, 0.4]
    model = SimpleNamespace(nodes=nodes, branches=branches, branch_xyz=bars)
    emf = _uniform_field_emf_columns(model)
    expected_area = np.array([4.8, -3.6, 0.])
    np.testing.assert_allclose(emf.sum(axis=0), -4e-7 * np.pi * expected_area,
                               rtol=1e-13, atol=1e-19)
    shift = np.array([2., -4., 8.])
    shifted = SimpleNamespace(nodes=nodes + shift, branches=branches,
                              branch_xyz=bars + shift)
    shifted_emf = _uniform_field_emf_columns(shifted)
    node_step = nodes[branches[:, 1]] - nodes[branches[:, 0]]
    np.testing.assert_allclose(shifted_emf - emf,
                               -2e-7 * np.pi * np.cross(shift, node_step),
                               rtol=1e-13, atol=1e-19)
    omega = 2 * np.pi * F_OP
    impedance = 2.0 + 1j * omega * 8e-6
    drives = np.column_stack((np.eye(3), [0.4, -0.9, 0.6]))
    expected = -1j * omega * 4e-7 * np.pi * (expected_area @ drives) / impedance
    for columns in (emf, shifted_emf):
        current = 1j * omega * (columns.sum(axis=0) @ drives) / impedance
        np.testing.assert_allclose(current, expected, rtol=1e-13, atol=1e-14)



def test_uniform_field_sweep_matches_dense_mna(monkeypatch):
    """Guards: a single basis handles the four inputs and gates lifted currents."""
    from types import SimpleNamespace

    from radia import fin_volume

    nodes = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    branches = np.array([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]])
    lat = SimpleNamespace(nodes=nodes, branches=branches, branch_xyz=nodes[branches])
    rng = np.random.default_rng(3)
    factor = rng.normal(size=(6, 6))
    L = 1e-6 * (factor.T @ factor + np.eye(6))
    R = np.linspace(.1, .2, 6)
    topology = {"L": L, "R": R, "ports": [(1, 0)]}
    calls = []

    def assemble(*args, **kwargs):
        calls.append(kwargs["max_branches"])
        return None, lat, None, topology, {}

    monkeypatch.setattr(fin_volume, "assemble_fin_volume", assemble)
    field = np.array([.2 + .4j, -.9, .6])
    drive = 1.5
    frequencies = (F_OP / 2, F_OP, 2 * F_OP)
    trial = fin_volume.solve_fin_volume_prima_sweep(
        None, frequencies=frequencies, expansion_frequency=F_OP, prima_order=8,
        sigma=5.8e7, current_A=drive, n_layers=1, n_lanes=8, n_stations=5,
        max_branches=13000, uniform_field_A_per_m=field)
    assert calls == [12000]
    assert tuple(p["frequency_hz"] for p in trial["points"]) == frequencies

    # Independently assemble grounded MNA and rectangular-loop Faraday sources.
    incidence = (np.eye(4)[branches[:, 1]] - np.eye(4)[branches[:, 0]]).T[1:]
    flux_rhs = np.r_[np.zeros(3), -2e-7 * np.pi * field[[2, 1, 0]] * [1, -1, 1]]
    for point in trial["points"]:
        s = 2j * np.pi * point["frequency_hz"]
        matrix = np.block([[np.diag(R) + s * L, -incidence.T],
                           [incidence, np.zeros((3, 3))]])
        reference = np.linalg.solve(matrix, np.r_[s * flux_rhs, drive, 0., 0.])[:6]
        result = point["solution"]
        np.testing.assert_allclose(result["I_fil"], reference, rtol=1e-6, atol=1e-10)
        assert max(result["prima"]["kcl_relative"]) <= 1e-6
        assert max(result["prima"]["kvl_relative"]) <= 1e-6
        assert result["coil_power_W"] == pytest.approx(.5 * np.vdot(reference, R * reference).real)
