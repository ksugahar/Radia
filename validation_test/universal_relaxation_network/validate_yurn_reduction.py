"""Private SA/RM checkpoint replay and uncertainty sensitivity, without raw exports.

Run only with caller-trusted research checkpoints. Outputs contain aggregate
metrics and hashes, never measured curves or learned parameter arrays.
"""

import argparse
import hashlib
import importlib
import json
import platform
import sys
import time
import types
from pathlib import Path

import numpy as np
import torch


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--epochs", type=int, default=1000)
    parser.add_argument("--restarts", type=int, default=1)
    parser.add_argument("--trusted-checkpoints", action="store_true")
    args = parser.parse_args()
    if not args.trusted_checkpoints:
        parser.error(
            "explicit --trusted-checkpoints is required for local research pickle files"
        )
    torch.set_num_threads(1)
    # URN is pure Python. Load its actual package-relative modules without the
    # unrelated Radia native root; no numerical implementation is mocked.
    package = types.ModuleType("_urn_validation")
    package.__path__ = [str(args.source)]
    sys.modules[package.__name__] = package
    core = importlib.import_module("_urn_validation.y_admittance_urn")
    reduction = importlib.import_module("_urn_validation.reduction")
    report = {
        "schema": "radia.yurn-private-reduction-validation.v1",
        "host": platform.node(),
        "python": sys.version,
        "torch": torch.__version__,
        "numpy": np.__version__,
        "threads": 1,
        "epochs": args.epochs,
        "restarts": args.restarts,
        "source_sha256": {
            p.name: digest(p)
            for p in [args.source / "y_admittance_urn.py", args.source / "reduction.py"]
        },
        "measurement_uncertainty_established": False,
        "scope": "scenario sensitivity on all training samples; no held-out accuracy claim",
        "checkpoints": [],
        "reductions": [],
    }
    for name in ("pcb", "nl87"):
        for stem in ("yurn34", "yurn_compact"):
            prefix = f"{stem}_{name}_log_components"
            pt, curve = args.inputs / f"{prefix}.pt", args.inputs / f"{prefix}.npz"
            payload = torch.load(pt, map_location="cpu", weights_only=False)
            values = dict(payload["config"])
            source_device = values.pop("device", "cpu")
            config = core.YAdmittanceURNConfig(**values)
            freq = np.asarray(payload["frequency_hz"], dtype=float)
            with np.load(curve, allow_pickle=False) as saved:
                target = saved["measured_impedance"]
                expected = saved["fitted_impedance"]
                np.testing.assert_array_equal(freq, saved["frequency_hz"])
            model = core.YAdmittanceURN(freq, config, z_data=target)
            model.load_state_dict(payload["state_dict"], strict=True)
            prediction = model.predict(freq)
            replay_error = float(
                np.max(
                    np.abs(prediction - expected) / np.maximum(np.abs(expected), 1e-30)
                )
            )
            if replay_error > 1e-7:
                raise RuntimeError(
                    f"{prefix}: checkpoint replay mismatch {replay_error}"
                )
            record = {
                "name": prefix,
                "bases": config.total_basis_functions,
                "samples": len(freq),
                "hashes": {pt.name: digest(pt), curve.name: digest(curve)},
                "source_device": source_device,
                "replay_device": "cpu",
                "replay_max_relative_error": replay_error,
                "s_domain_rmse": core.s_domain_rmse(prediction, target),
                "log_component_rmse": core.log_component_rmse(prediction, target),
                "redundant_pairs": reduction.redundant_y_bases(model, freq),
            }
            report["checkpoints"].append(record)
            print(json.dumps(record), flush=True)
            if stem != "yurn_compact":
                continue
            started = time.perf_counter()
            _, trace = reduction.reduce_y_admittance_urn(
                model,
                freq,
                target,
                uncertainty_ohm=0.01 * np.abs(target),
                refit_epochs=args.epochs,
                refit_restarts=args.restarts,
            )
            scenarios = []
            for percent in (0.5, 1, 2, 5, 10, 20):
                eligible = [
                    t for t in trace["trials"] if t["max_uncertainty_ratio"] <= percent
                ]
                scenarios.append(
                    {
                        "assumed_relative_budget_percent": percent,
                        "selected_basis_count": min(
                            (t["basis_count"] for t in eligible), default=None
                        ),
                    }
                )
            report["reductions"].append(
                {
                    "dataset": name,
                    "elapsed_seconds": time.perf_counter() - started,
                    "scenarios": scenarios,
                    "trace": trace,
                }
            )
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(
                json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
            )
            print(json.dumps(report["reductions"][-1]), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
