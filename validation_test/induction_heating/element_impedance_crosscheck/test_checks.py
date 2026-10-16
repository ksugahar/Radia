"""Lightweight protections against invalid spatial/phase comparisons."""
import json
import hashlib
import zipfile
from pathlib import Path

import numpy as np
import pytest

from boundary_fem import digest, read_case, read_samples
from compare import spatial_error
from reaction_gate import evaluate
import reaction_gate

CASE = Path(__file__).with_name("case.json")


def reaction_record():
    return dict(phasor_convention="peak amplitude, exp(+i omega t)",
                case_kind="nonmagnetic-reactive", supplied_coil_impedance_ohm=[0, 0],
                material_contract="constant-mu0-positive-resistivity-no-remanence-no-bh",
                accepted_updates=[dict(update_index=1, voltage_update_index=1,
                                       heat_update_index=1, current_update_index=1,
                                       state="accepted", current_A=[0, 2],
                                       port_voltage_V=[4, 2], face_heat_W=[.5, 1.5],
                                       linear_true_relative_residual=1e-10)])


def test_port_work_uses_conjugate_current_and_does_not_claim_qualification():
    plan = json.loads(CASE.with_name("reaction_gate_plan.json").read_text())
    result = evaluate(reaction_record(), plan)
    assert result["recorded_updates_pass"]
    assert result["updates"][0]["reactive_power_var"] == -4
    assert not result["b1_closed"] and not result["release_accepted"]


def test_reaction_gate_rejects_wrong_sign_and_noncoherent_updates():
    plan = json.loads(CASE.with_name("reaction_gate_plan.json").read_text())
    data = reaction_record()
    data["accepted_updates"][0]["port_voltage_V"] = [-4, 2]
    assert not evaluate(data, plan)["recorded_updates_pass"]
    data["accepted_updates"][0]["heat_update_index"] = 0
    with pytest.raises(ValueError, match="same accepted update"):
        evaluate(data, plan)


def test_each_update_must_balance_and_large_coil_subtraction_is_rejected():
    plan = json.loads(CASE.with_name("reaction_gate_plan.json").read_text())
    data = reaction_record()
    data["accepted_updates"][0]["face_heat_W"] = [.1, .1]
    assert not evaluate(data, plan)["recorded_updates_pass"]
    data["supplied_coil_impedance_ohm"] = [.1, 0]
    with pytest.raises(ValueError, match="cancellation"):
        evaluate(data, plan)


def test_cli_binds_acquisition_to_wheel_bytes_and_rejects_module_drift(tmp_path, monkeypatch):
    mex = tmp_path/"radia_cuda_ih.mexw64"
    mex.write_bytes(b"synthetic module bytes; never executed")
    wheel = tmp_path/"synthetic.whl"
    with zipfile.ZipFile(wheel, "w") as archive:
        archive.write(mex, "radia_gpu/matlab/radia_cuda_ih.mexw64")
    sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    plan = CASE.with_name("reaction_gate_plan.json")
    record = reaction_record()
    record["artifact_identity"] = dict(source_commit="a"*40, source_dirty=False,
        wheel_sha256=sha(wheel), loaded_mex_sha256=sha(mex),
        loaded_mex_path=str(mex.resolve()), plan_sha256=sha(plan))
    source = tmp_path/"records.json"
    source.write_text(json.dumps(record))
    output = tmp_path/"result.json"
    monkeypatch.setattr("sys.argv", ["reaction_gate.py", str(source), "--output", str(output),
        "--wheel", str(wheel), "--loaded-mex", str(mex), "--source-commit", "a"*40])
    reaction_gate.main()
    result = json.loads(output.read_text())
    assert result["artifact_identity"]["wheel_sha256"] == sha(wheel)
    assert "loaded_mex_path" not in result["artifact_identity"]
    assert not result["b1_closed"]
    mex.write_bytes(b"different module")
    with pytest.raises(ValueError, match="match the unique installed-wheel"):
        reaction_gate.main()


def test_complex_vector_phase_is_not_removed():
    reference = np.array([[1., 0., 0.], [0., 2., 0.]], dtype=complex)
    assert spatial_error(1j*reference, reference, [1., 3.], 1e-12) == pytest.approx(np.sqrt(2))


def test_spatial_norm_uses_physical_area_and_rejects_negative_weight():
    reference = np.array([1., 1.])
    assert spatial_error([2., 1.], reference, [1., 3.], 1e-18) == pytest.approx(.5)
    with pytest.raises(ValueError, match="positive area"):
        spatial_error([2., 1.], reference, [-1., 3.], 1e-18)


def test_fractional_geometry_indices_are_not_truncated(tmp_path):
    data = json.loads(CASE.read_text())
    data["triangles_one_based"][0][0] = 1.5
    path = tmp_path/"case.json"
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="Integer"):
        read_case(path)


def test_material_patch_permutation_cannot_pass_sample_identity(tmp_path):
    _, vertices, triangles, _ = read_case(CASE)
    points = [vertices[t].mean(axis=0).tolist() for t in triangles]
    areas = [np.linalg.norm(np.cross(vertices[t[1]]-vertices[t[0]], vertices[t[2]]-vertices[t[0]]))/2 for t in triangles]
    data = dict(case_file_sha256=digest(CASE), points_m=points,
                patch_one_based=[1, 2, 3, 4], area_weights_m2=areas)
    path = tmp_path/"samples.json"
    path.write_text(json.dumps(data))
    assert read_samples(path, CASE, vertices, triangles)["patch_one_based"] == [1, 2, 3, 4]
    data["patch_one_based"] = [2, 1, 3, 4]
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="patch"):
        read_samples(path, CASE, vertices, triangles)
