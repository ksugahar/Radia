"""Fail-closed contract of the fin P2 adoption entry.

Guards: validation_records/induction_heating/20261116_fin_peec_axisym_p2/entry.py writes a passing adoption
result only from a physically passing, clean, fully certified verdict, and takes the reference from the terminal
checkpoint of the converged order 2, exterior 3 sequence.
Assumes: the verdict layout of validation_test/induction_heating/axisym_volume_peec_fin.py; the numbers below are
synthetic and carry no physical meaning.
"""
import copy
import importlib.util
import json
from pathlib import Path

import pytest

ENTRY = (Path(__file__).resolve().parents[1] / "validation_records" / "induction_heating"
         / "20261116_fin_peec_axisym_p2" / "entry.py")
QUANTITIES = dict(R_ohm=8.0e-3, L_H=1.06e-6, root_heat_W_per_A2=2.07e-4, tip_heat_W_per_A2=8.1e-4)


def load_entry():
    spec = importlib.util.spec_from_file_location("fin_p2_entry", ENTRY)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def reference_rows():
    rows = []
    for order in (1, 2):
        for exterior in (1.5, 3.0):
            for level in range(3):
                rows.append(dict(order=order, exterior_factor=exterior, level=level, ndof=1000*(level + 1)*order,
                                 checkpoint=level in (0, 2), within_budget=True,
                                 **{k: v*(1.0 + 1e-3*(2 - level)) for k, v in QUANTITIES.items()}))
    return rows


def verdict():
    peec = [dict(boundary_h_mm=h, cells=c, true_relative_residual=1e-15, assembly_s=1.0, solve_s=1.0,
                 **{k: v*(1.0 + f) for k, v in QUANTITIES.items()})
            for h, c, f in ((0.06, 100, 4e-3), (0.04, 200, 2e-3), (0.028, 400, 1e-3))]
    sequences = [dict(order=o, exterior_factor=e, converged=True, stop="converged") for o in (1, 2) for e in (1.5, 3.0)]
    return dict(schema="synthetic", git_head="0"*40, git_dirty=False, script_sha256="0"*64,
                peec_runtime=dict(peec_matrices_sha256="1"*64), reference_runtime=dict(axifem_sha256="2"*64),
                radius_mm=200.0, frequency_hz=1.5e5, sigma_S_per_m=5.8e7, skin_depth_mm=0.17, cuts_mm=[1.6, 3.5],
                normalization="1 A peak", reference_levels=reference_rows(), reference_sequences=sequences,
                reference_spread_max={}, reference_resolved=True, reference_converged=True, reference_swept=True,
                reference_loss_consistent=True, reference_sequences_without_checkpoint=[],
                reference_certification={"o1_f1.5": [], "o1_f3.0": [], "o2_f1.5": [], "o2_f3.0": []},
                peec_levels=peec, verdict=dict(status="PASS"), not_claimed=[])


def run(tmp_path, data):
    vp, hp, op = tmp_path / "verdict.json", tmp_path / "hashes.json", tmp_path / "out.json"
    vp.write_text(json.dumps(data))
    hp.write_text(json.dumps({str(r["boundary_h_mm"]): dict(cells=r["cells"], sha256="3"*64)
                              for r in data["peec_levels"]}))
    load_entry().main(vp, hp, op)
    return json.loads(op.read_text())


def test_passing_verdict_takes_the_final_checkpoint(tmp_path):
    out = run(tmp_path, verdict())
    assert out["adoption_result"]["verdict"] == "PASS", out["adoption_result"]["verdict"]
    assert out["reference"]["level"] == 2 and out["reference"]["ndof"] == 6000, out["reference"]
    shuffled = verdict()
    shuffled["reference_levels"] = shuffled["reference_levels"][::-1]
    assert run(tmp_path, shuffled)["reference"]["level"] == 2


def _mutate(path, value):
    def apply(data):
        target = data
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value
    return apply


def _o2f3(data):
    return [r for r in data["reference_levels"] if r["order"] == 2 and r["exterior_factor"] == 3.0]


NEGATIVE = {
    "physical FAIL": _mutate(("verdict", "status"), "FAIL"),
    "dirty source": _mutate(("git_dirty",), True),
    "reference unresolved": _mutate(("reference_resolved",), False),
    "reference not converged": _mutate(("reference_converged",), False),
    "open certification item": _mutate(("reference_certification", "o2_f3.0"), ["tip_heat"]),
    "sequence without checkpoint": _mutate(("reference_sequences_without_checkpoint",), ["o1_f1.5"]),
    "final level not a checkpoint": lambda d: _o2f3(d)[-1].update(checkpoint=False),
    "final level over budget": lambda d: _o2f3(d)[-1].update(within_budget=False),
    "duplicate level": lambda d: _o2f3(d)[1].update(level=2),
    "sequence stopped by budget": lambda d: d["reference_sequences"][3].update(converged=False, stop="budget"),
}


@pytest.mark.parametrize("name", sorted(NEGATIVE))
def test_entry_fails_closed(tmp_path, name):
    data = copy.deepcopy(verdict())
    NEGATIVE[name](data)
    with pytest.raises(SystemExit):
        run(tmp_path, data)
    assert not (tmp_path / "out.json").exists(), f"{name}: an adoption result was written"
