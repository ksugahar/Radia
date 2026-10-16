"""Evaluate recorded port work without reusing a BEM assembly matrix.

This is a diagnostic evaluator, not an artifact or physical qualification gate.
The preregistered plan is still awaiting independent formula review.
"""
import argparse
import hashlib
import json
import math
import re
import zipfile
from pathlib import Path


def phasor(value):
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError("Phasors require [real, imaginary] SI components")
    if any(isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x)
           for x in value):
        raise ValueError("Finite numeric phasor components required")
    return complex(*value)


def evaluate(data, plan):
    if data.get("phasor_convention") != plan["phasor_convention"]:
        raise ValueError("Peak phasor/sign convention mismatch")
    mode = data.get("case_kind")
    if mode not in ("temperature-bh-zero-resistance", "nonmagnetic-reactive"):
        raise ValueError("Explicit dedicated case kind required")
    coil = phasor(data["supplied_coil_impedance_ohm"])
    if coil.real != 0:
        raise ValueError("Zero coil resistance required to avoid cancellation")
    if mode == "nonmagnetic-reactive":
        if coil != 0 or data.get("material_contract") != "constant-mu0-positive-resistivity-no-remanence-no-bh":
            raise ValueError("Reactive sign requires the separate nonmagnetic case")
    updates = data.get("accepted_updates")
    if not isinstance(updates, list) or not updates:
        raise ValueError("Nonempty accepted updates required")
    balance = plan["power_balance"]
    rows = []
    previous = 0
    for row in updates:
        step = row.get("update_index")
        if type(step) is not int or step <= previous:
            raise ValueError("Strictly increasing positive accepted update indices required; omit placeholders")
        previous = step
        if row.get("state") != "accepted" or row.get("voltage_update_index") != step or row.get("heat_update_index") != step or row.get("current_update_index") != step:
            raise ValueError("Voltage, current and heat must identify the same accepted update")
        current = phasor(row["current_A"])
        voltage = phasor(row["port_voltage_V"])
        if abs(current) == 0:
            raise ValueError("Nonzero current required")
        heat = row["face_heat_W"]
        if not isinstance(heat, list) or not heat or any(
            isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) or x < 0 for x in heat
        ):
            raise ValueError("Finite nonnegative face-integrated heat required")
        total = math.fsum(heat)
        work = .5 * voltage * current.conjugate()
        reaction = voltage/current - coil
        if not all(math.isfinite(x) for x in (total, work.real, work.imag, reaction.real, reaction.imag)):
            raise ValueError("Nonfinite derived work or reaction")
        error = abs(work.real - total)/max(abs(total), balance["normalization_floor_W"])
        residual = row["linear_true_relative_residual"]
        if isinstance(residual, bool) or not isinstance(residual, (int, float)) or not math.isfinite(residual) or residual < 0:
            raise ValueError("Finite nonnegative true relative residual required")
        power_pass = total > balance["minimum_body_heat_W"] and error <= balance["relative_limit"]
        sign_pass = (reaction.imag < plan["reactive_sign"]["maximum_imaginary_reaction_impedance_ohm"]
                     if mode == "nonmagnetic-reactive" else None)
        rows.append(dict(update_index=step, port_power_W=work.real, body_heat_W=total,
                         coil_loss_W=0., power_balance_relative_error=error,
                         power_balance_pass=power_pass, reactive_power_var=work.imag,
                         reaction_impedance_ohm=[reaction.real, reaction.imag],
                         reactive_sign_pass=sign_pass,
                         true_residual_pass=residual <= plan["linear_true_relative_residual_limit"]))
    return dict(schema_version=1, status="diagnostic-evaluation", case_kind=mode,
                recorded_updates_pass=all(r["power_balance_pass"] and r["true_residual_pass"]
                                         and r["reactive_sign_pass"] is not False for r in rows),
                updates=rows, b1_closed=False, release_accepted=False,
                limitations=["Input provenance and material declarations are not independently certified by this evaluator",
                             "Independent FEM and two-refinement convergence must be reviewed separately",
                             "CPU/GPU parity and installed-wheel/Toolkit acceptance are separate gates"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("records", type=Path)
    parser.add_argument("--plan", type=Path, default=Path(__file__).with_name("reaction_gate_plan.json"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--loaded-mex", type=Path, required=True,
                        help="Actual loaded module path recorded by the MATLAB acquisition driver")
    parser.add_argument("--source-commit", required=True, help="Expected clean candidate full SHA")
    args = parser.parse_args()
    raw, plan_raw = args.records.read_bytes(), args.plan.read_bytes()
    data = json.loads(raw)
    if not re.fullmatch(r"[0-9a-f]{40}", args.source_commit):
        raise ValueError("Full candidate source SHA required")
    wheel_hash = hashlib.sha256(args.wheel.read_bytes()).hexdigest()
    mex_hash = hashlib.sha256(args.loaded_mex.read_bytes()).hexdigest()
    with zipfile.ZipFile(args.wheel) as archive:
        modules = [name for name in archive.namelist() if name.endswith("/radia_cuda_ih.mexw64")]
        if len(modules) != 1 or hashlib.sha256(archive.read(modules[0])).hexdigest() != mex_hash:
            raise ValueError("Loaded IH MEX must match the unique installed-wheel module")
    expected = dict(source_commit=args.source_commit, wheel_sha256=wheel_hash,
                    loaded_mex_sha256=mex_hash,
                    loaded_mex_path=str(args.loaded_mex.resolve()),
                    plan_sha256=hashlib.sha256(plan_raw).hexdigest(), source_dirty=False)
    identity = data.get("artifact_identity", {})
    for key, value in expected.items():
        if key not in identity or type(identity[key]) is not type(value) or identity[key] != value:
            raise ValueError(f"Acquisition artifact identity mismatch: {key}")
    result = evaluate(data, json.loads(plan_raw))
    result.update(input_sha256=hashlib.sha256(raw).hexdigest(),
                  plan_sha256=hashlib.sha256(plan_raw).hexdigest(),
                  evaluator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    result["artifact_identity"] = {key: value for key, value in expected.items() if key != "loaded_mex_path"}
    result["limitations"].append("Loaded-path declaration and source-to-wheel build receipt require external review; matching bytes alone is insufficient")
    # Exclusive creation protects previous evidence even if the requested path
    # accidentally equals the input or plan path.
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    if not result["recorded_updates_pass"]:
        raise SystemExit("Recorded reaction diagnostic gates failed; see output JSON")


if __name__ == "__main__":
    main()
