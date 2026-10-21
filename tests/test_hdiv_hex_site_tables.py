"""Rule-key isolation and live-instance lifetime for shared reference tables."""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
from ngsolve.meshes import MakeStructured3DMesh
import radia.vim._vim as vm


@pytest.fixture
def site_table_inputs(monkeypatch):
    mesh = MakeStructured3DMesh(hexes=True, nx=1, ny=1, nz=1)
    original = vm._rp._ChargeGramHMatrix
    captured = {}
    def capture(**kwargs):
        captured.update(kwargs)
        return original(**kwargs)
    with ng.TaskManager():
        with monkeypatch.context() as patch:
            patch.setattr(vm._rp, "_ChargeGramHMatrix", capture)
            fes = ng.HDiv(mesh, order=1)
            _, first, _ = vm.build_charge_gram(fes, _materialize_mass=False, _build_hmatrix=False)
        yield original, captured, first


def test_site_tables_key_weights_and_preserve_evicted_live_instances(site_table_inputs):
    original, captured, first = site_table_inputs
    initial = first.hex_state_check()
    assert initial["ctor"] == initial["now"]
    first_tables = first.hex_state_breakdown()
    changed = dict(captured, gw_in=np.asarray(captured["gw_in"]) * 1.001)
    # Same length/nodes, different weights must not alias the reference data.
    second = original(**changed)
    assert second.hex_state_breakdown()["cellSiteRad"] != first_tables["cellSiteRad"]
    third = original(**dict(captured, gl_in=np.asarray(captured["gl_in"])[:-1],
                           gw_in=np.asarray(captured["gw_in"])[:-1]))
    for owner in (first, second, third):
        check = owner.hex_state_check()
        assert check["ctor"] == check["now"]
    assert first.hex_state_check() == initial
    assert first.hex_state_breakdown() == first_tables


def test_concurrent_cold_rule_constructors_preserve_each_instance(site_table_inputs):
    from concurrent.futures import ThreadPoolExecutor
    import threading
    original, captured, first = site_table_inputs
    first_state = first.hex_state_check()
    # The binding releases the GIL during native construction. Two synchronized
    # workers exercise cold different-rule publication under the shared mutex.
    # This is a correctness canary, with no scheduling or throughput assertion.
    with ThreadPoolExecutor(max_workers=2) as pool:
        for phase in range(4):
            configs = [dict(captured, gw_in=np.asarray(captured["gw_in"]) *
                            (1. + .001 * (2 * phase + worker + 1))) for worker in range(2)]
            expected = []
            for config in configs:
                serial = original(**config)
                expected.append(serial.hex_state_breakdown())
                del serial
            sentinel = original(**captured)  # Evict both candidate rule keys.
            barrier = threading.Barrier(2, timeout=15)
            def construct(config):
                worker = threading.get_ident()
                barrier.wait()
                return worker, original(**config)
            pending = [pool.submit(construct, config) for config in configs]
            results = [future.result(timeout=30) for future in pending]
            assert len({worker for worker, _ in results}) == 2
            for (_, owner), reference in zip(results, expected):
                state = owner.hex_state_check()
                assert state["ctor"] == state["now"]
                assert owner.hex_state_breakdown() == reference
            assert first.hex_state_check() == first_state
            assert sentinel.hex_state_check()["ctor"] == sentinel.hex_state_check()["now"]
