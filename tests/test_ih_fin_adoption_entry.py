"""Fail-closed contract of the fin P2 adoption entry.

Guards: validation_records/induction_heating/20261116_fin_peec_axisym_p2/entry.py writes a passing adoption
result only when every reference sequence of the verdict re-certifies under the approved settings, the stored
checks, certification and physical verdict equal their re-derivation by the runner's own functions, and the
source is clean; the reference is the final checkpoint of order 2 at exterior 3.
Assumes: the verdict layout and the resolve_reference/certify_sequence functions of
validation_test/induction_heating/axisym_volume_peec_fin.py; the numbers below are synthetic.
"""
import copy
import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
ENTRY = ROOT / "validation_records" / "induction_heating" / "20261116_fin_peec_axisym_p2" / "entry.py"
BASE = dict(R_ohm=8.0e-3, L_H=1.06e-6, root_heat_W_per_A2=2.07e-4, tip_heat_W_per_A2=8.1e-4)


def load_entry():
    spec = importlib.util.spec_from_file_location("fin_p2_entry", ENTRY)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sequence(order, exterior):
    """Three levels, checkpoints at 0 and 2 (four times apart in DOFs), converged within a tenth of the limits."""
    levels = []
    for level in range(3):
        scale = 1.0 + 2e-4*(2 - level)
        levels.append(dict(level=level, exterior_factor=exterior, ndof=1000*4**level*order, elapsed_s=10.0*(level + 1),
                           within_budget=True, checkpoint=level != 1, terminal_vs_energy_R=1e-9,
                           **{k: v*scale for k, v in BASE.items()}))
    return dict(order=order, exterior_factor=exterior, start_maxh_mm=0.05, air_maxh_mm=20.0, air_divisor=20.0,
                theta=0.5, checkpoint_dof_growth=2.0, max_ndof=6e6, max_levels=16, budget_s=7200.0,
                elapsed_s=40.0, converged=True, stop="converged", levels=levels)


def verdict():
    entry = load_entry()
    runner, targets, physical_verdict = entry._runner()
    tol = {k: 0.1*targets[k] for k in ("R", "partial_L", "root_heat", "tip_heat")}
    sequences = [sequence(o, f) for o in entry.ORDERS for f in entry.EXTERIORS]
    for s in sequences:
        s["tolerance"] = tol
    res = runner.resolve_reference(sequences, entry.ORDERS, entry.EXTERIORS,
                                   lambda o: {k: v for k, v in sequences[0].items()
                                              if k not in ("order", "exterior_factor", "elapsed_s", "converged",
                                                           "stop", "levels")}, targets)
    assert res["resolved"], res
    peec = [dict(boundary_h_mm=h, cells=c, true_relative_residual=1e-15, assembly_s=1.0, solve_s=1.0,
                 **{k: v*(1.0 + f) for k, v in BASE.items()})
            for h, c, f in ((0.06, 100, 4e-3), (0.04, 200, 2e-3), (0.028, 400, 1e-3))]
    return dict(schema="synthetic", git_head="0"*40, git_dirty=False, script_sha256="0"*64,
                peec_runtime=dict(peec_matrices_sha256="1"*64), reference_runtime=dict(axifem_sha256="2"*64),
                radius_mm=200.0, frequency_hz=1.5e5, sigma_S_per_m=5.8e7, skin_depth_mm=0.17, cuts_mm=[1.6, 3.5],
                normalization="1 A peak", reference_settings=dict(entry.SETTINGS, tolerance=tol),
                reference_levels=[dict(level, order=s["order"]) for s in sequences for level in s["levels"]],
                reference_sequences=[{k: v for k, v in s.items() if k != "levels"} for s in sequences],
                reference_checks=res["checks"], reference_spread_max=res["spread"],
                reference_resolved=res["resolved"], reference_converged=res["converged"],
                reference_swept=res["swept"], reference_loss_consistent=res["consistent"],
                reference_sequences_without_checkpoint=res["missing"], reference_certification=res["certification"],
                peec_levels=peec, verdict=physical_verdict(runner.quantities(peec[-1]), res["reference"], True),
                not_claimed=[])


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
    assert (out["reference"]["level"], out["reference"]["ndof"]) == (2, 32000), out["reference"]


def _o2f3(data):
    return [r for r in data["reference_levels"] if r["order"] == 2 and r["exterior_factor"] == 3.0]


NEGATIVE = {
    "physical FAIL": lambda d: d["peec_levels"][-1].update(R_ohm=d["peec_levels"][-1]["R_ohm"]*1.2),
    "stored verdict altered": lambda d: d["verdict"].update(status="FAIL"),
    "dirty source": lambda d: d.update(git_dirty=True),
    "certification keys missing": lambda d: d.update(reference_certification={"o1_f1.5": []}),
    "stored checks altered": lambda d: d["reference_checks"]["order_o1"].update(R=1e-3),
    "final level over budget": lambda d: _o2f3(d)[-1].update(elapsed_s=7201.0),
    "final level not a checkpoint": lambda d: _o2f3(d)[-1].update(checkpoint=False),
    "duplicate level": lambda d: _o2f3(d)[1].update(level=2),
    "reference value altered": lambda d: _o2f3(d)[-1].update(tip_heat_W_per_A2=_o2f3(d)[-1]["tip_heat_W_per_A2"]*1.01),
    "sequence stopped by budget": lambda d: d["reference_sequences"][3].update(converged=False, stop="budget"),
    "settings changed": lambda d: d["reference_settings"].update(budget_s=1e9),
    "sequence missing": lambda d: d["reference_sequences"].pop(0),
}


@pytest.mark.parametrize("name", sorted(NEGATIVE))
def test_entry_fails_closed(tmp_path, name):
    data = copy.deepcopy(verdict())
    NEGATIVE[name](data)
    with pytest.raises(SystemExit):
        run(tmp_path, data)
    assert not (tmp_path / "out.json").exists(), f"{name}: an adoption result was written"
