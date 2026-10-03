import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
RESULT = HERE / "team28_arnoldi_height_sweep.json"


def test_arnoldi_height_sweep_is_fresh_and_meets_the_reduced_gates():
    payload = json.loads(RESULT.read_text(encoding="utf-8"))

    assert payload["schema"] == "radia.team28-arnoldi-height-sweep.v1"
    assert payload["pass"] is True
    assert all(payload["checks"].values())
    for relative, digest in payload["source_sha256"].items():
        data = (REPO_ROOT / relative).read_bytes().replace(b"\r\n", b"\n")
        assert hashlib.sha256(data).hexdigest() == digest, relative
    assert payload["dZ_mm"] == list(range(-7, 18))
    # Gates of the retired reduced sweep, on the supported reduction.
    assert payload["max_abs_reduced_minus_lab_N"] < 1.0e-3
    assert payload["max_abs_reduced_minus_full_fem_sweep_N"] < 1.0e-3
    assert abs(payload["equilibrium_abs_height_mm"]
               - payload["published_stationary_height_mm"]) < 0.6
    # The order is the smallest one meeting the measured selection error.
    errors = payload["max_abs_reduced_minus_direct_N_by_order"]
    selected = payload["selected_order"]
    assert errors[selected - 1] < payload["order_selection_N"]
    assert all(error >= payload["order_selection_N"] for error in errors[:selected - 1])
