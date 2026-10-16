"""Fail-closed aggregation of the predeclared same-boundary-law experiment."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def read(path):
    def reject(value):
        raise ValueError("Nonfinite JSON constant: " + value)
    record = json.loads(Path(path).read_text(encoding="utf-8-sig"), parse_constant=reject)
    if record.get("status") != "computed-not-yet-crosschecked":
        raise ValueError("Only completed independent computation records are comparable")
    return record


def scalar_error(value, reference, floor):
    if not np.isfinite([value, reference]).all():
        raise ValueError("Nonfinite comparison scalar")
    return float(abs(value-reference)/max(abs(reference), floor))


def spatial_error(value, reference, weights, floor):
    value, reference, weights = np.asarray(value), np.asarray(reference), np.asarray(weights)
    if value.shape != reference.shape or value.shape[0] != len(weights):
        raise ValueError("Spatial comparisons require identical physical sample shapes")
    if not np.isfinite(value).all() or not np.isfinite(reference).all() or not np.isfinite(weights).all() or np.any(weights <= 0):
        raise ValueError("Finite spatial fields and positive area weights required")
    def norm(x):
        squared = np.abs(x)**2
        if squared.ndim == 2:
            squared = squared.sum(axis=1)
        return float(np.sqrt(np.dot(weights, squared)))
    return norm(value-reference)/max(norm(reference), floor)


def reaction(record):
    return complex(record["reaction_impedance_ohm_real"], record["reaction_impedance_ohm_imag"])


def fields(record):
    return np.asarray(record["H_t_real"]) + 1j*np.asarray(record["H_t_imag"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--artifact-receipt", type=Path, required=True)
    for role in ("bem-coarse", "bem-fine", "fem-coarse", "fem-fine", "fem-exterior", "fem-quadrature", "fem-offset"):
        parser.add_argument("--"+role, type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError("Refusing to overwrite previous comparison evidence")
    case = json.loads(args.case.read_text())
    case_hash = hashlib.sha256(args.case.read_bytes()).hexdigest()
    roles = ("bem_coarse", "bem_fine", "fem_coarse", "fem_fine", "fem_exterior", "fem_quadrature", "fem_offset")
    records = {role: read(getattr(args, role)) for role in roles}
    receipt = json.loads(args.artifact_receipt.read_text())
    if receipt.get("source_dirty") is not False or not receipt.get("binary_sha256"):
        raise ValueError("Clean exact-source artifact receipt required")
    if receipt.get("accepted") is not True and receipt.get("schema") != "radia.gpu.ih-research-payload.v1":
        raise ValueError("Unaccepted artifact build receipt")
    for role in ("bem_coarse", "bem_fine"):
        record = records[role]
        if record["case_file_sha256"] != case_hash or record["mex_sha256"] != receipt["binary_sha256"]:
            raise ValueError("BEM case or actual MEX artifact identity mismatch")
        if record["frequency_Hz"] != case["frequency_Hz"] or record["current_peak_A"] != case["current_peak_A"]:
            raise ValueError("BEM excitation differs from the declared FEM case")
    if records["bem_coarse"]["refinement_level"] != 1 or records["bem_fine"]["refinement_level"] != 2:
        raise ValueError("Both predeclared BEM refinements are required")
    fem = [records[x] for x in roles if x.startswith("fem_")]
    if any(r["plan"]["case_sha256"] != case_hash for r in fem):
        raise ValueError("Independent FEM case mismatch")
    if any(r["plan"].get("source_dirty") is not False or r.get("label_check_passed") is not True for r in fem):
        raise ValueError("Clean committed FEM sources and checked mesh labels are required")
    if len({r["plan"]["source_sha256"] for r in fem}) != 1 or len({r["ngsolve"] for r in fem}) != 1:
        raise ValueError("FEM controls must use one implementation and runtime")
    for role, index in (("fem_coarse", 0), ("fem_fine", 1), ("fem_exterior", 2), ("fem_quadrature", 2), ("fem_offset", 2)):
        if records[role]["run"] != case["fem_runs"][index] or records[role]["plan"]["run_index"] != index:
            raise ValueError("Predeclared FEM run identity mismatch")
    baseline = records["fem_exterior"]
    if any(records[role]["mesh_sha256"] != baseline["mesh_sha256"] for role in ("fem_quadrature", "fem_offset")):
        raise ValueError("Quadrature and trace controls must reuse the exact exterior reference mesh")
    if baseline["plan"]["filament_gauss_order"] != 16 or records["fem_quadrature"]["plan"]["filament_gauss_order"] != 24:
        raise ValueError("Both predeclared filament quadrature orders are required")
    if baseline["plan"]["surface_probe_offset_m"] != 1e-9 or records["fem_offset"]["plan"]["surface_probe_offset_m"] != 5e-10:
        raise ValueError("Both predeclared air-side offsets are required")
    if records["fem_fine"]["ndof"] <= records["fem_coarse"]["ndof"]:
        raise ValueError("FEM refinement did not increase resolution")
    b = records["bem_fine"]
    weights = np.asarray(b["area_weights_m2"]).reshape(-1)
    for r in fem:
        for key in ("points_m", "patch_one_based", "area_weights_m2"):
            if not np.array_equal(np.asarray(r[key]), np.asarray(b[key])):
                raise ValueError("Field comparison lacks identical physical quadrature")
    patches = np.asarray(b["patch_one_based"], dtype=int).reshape(-1)-1
    zs_real = np.asarray(case["patch_zs_ohm_real"])[patches]
    q = lambda r: .5*zs_real*np.sum(np.abs(fields(r))**2, axis=1)
    gates = {}
    def gate(name, error, limit):
        if not np.isfinite(error):
            raise ValueError("Nonfinite gate value: " + name)
        gates[name] = dict(error=float(error), tolerance=limit, passed=bool(error <= limit))
    gate("body_power", scalar_error(b["total_power_W"], baseline["total_power_W"], 1e-18), .02)
    gate("body_reaction", scalar_error(reaction(b), reaction(baseline), 1e-12), .03)
    gate("complex_tangential_H", spatial_error(fields(b), fields(baseline), weights, 1e-12), .05)
    gate("q_surf", spatial_error(np.asarray(b["q_surf_W_per_m2"]).reshape(-1), q(baseline), weights, 1e-18), .10)
    for label, first, second, limit in (
        ("fem_mesh", "fem_coarse", "fem_fine", .005),
        ("fem_exterior", "fem_fine", "fem_exterior", .005),
        ("bem_refinement", "bem_coarse", "bem_fine", .01),
        ("filament_quadrature", "fem_exterior", "fem_quadrature", .001)):
        a, r = records[first], records[second]
        gate(label+"_power", scalar_error(a["total_power_W"], r["total_power_W"], 1e-18), limit)
        gate(label+"_reaction", scalar_error(reaction(a), reaction(r), 1e-12), limit)
    offset = records["fem_offset"]
    gate("probe_offset_H", spatial_error(fields(baseline), fields(offset), weights, 1e-12), .001)
    gate("probe_offset_q", spatial_error(q(baseline), q(offset), weights, 1e-18), .001)
    for role, r in records.items():
        if r["total_power_W"] <= 0 or r["linear_true_relative_residual"] < 0 or r["power_identity_relative"] < 0:
            raise ValueError("Positive physical power required: " + role)
        # The identity's declared denominator is dissipated power, not source work.
        recomputed = abs(r["total_power_W"]-.5*reaction(r).real*case["current_peak_A"]**2)/max(r["total_power_W"], 1e-18)
        if abs(recomputed-r["power_identity_relative"]) > 1e-12*max(1, recomputed):
            raise ValueError("Stored power identity disagrees with recorded work: " + role)
        gate(role+"_true_residual", r["linear_true_relative_residual"], 1e-6)
        gate(role+"_power_identity", r["power_identity_relative"], .01)
    out = dict(schema="radia.ih-element-impedance-comparison.v1",
               same_boundary_law_pass=all(g["passed"] for g in gates.values()),
               nonlinear_volumetric_validation="not-performed", experimental_validation="not-performed",
               release_accepted=False, source_commit=receipt["source_commit"], native_sha256=receipt["binary_sha256"],
               case_file_sha256=case_hash, gates=gates,
               input_sha256={role: hashlib.sha256(getattr(args, role).read_bytes()).hexdigest() for role in roles})
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"same_boundary_law_pass": out["same_boundary_law_pass"], "failed_gates": [k for k,g in gates.items() if not g["passed"]]}))
    raise SystemExit(0 if out["same_boundary_law_pass"] else 1)


if __name__ == "__main__":
    main()
