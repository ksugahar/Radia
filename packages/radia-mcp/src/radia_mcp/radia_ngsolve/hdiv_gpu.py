"""Thin artifact-producing HDiv GPU workflow; no domain math in MCP."""
from __future__ import annotations

import json
import hashlib
import os
import tempfile
from pathlib import Path


def hdiv_nonlinear_gpu_run(vol_path: str, label_contract: str, law_path: str,
                           applied_fields: list[list[float]], result_path: str,
                           backend: str, order: int = 1,
                           preconditioner: str = "reference-diagonal",
                           tolerance: float = 5e-14, linear_tolerance: float = 5e-14) -> dict:
    """Reuse one prepared shared GPU owner for up to eight uniform H cases.

    Input is an inverse-BH polynomial payload, or a radia.bh-table.v1 object
    with ``table`` rows [H in A/m, B in tesla], using Radia's existing law fit.
    Output is private execution evidence, not independent CPU certification.
    Optional numerical dependencies load only when this workflow is called.
    """
    if backend != "gpu":
        raise ValueError("Select backend='gpu' explicitly; CPU uses the radia.vim workflow")
    if not isinstance(applied_fields, list) or not 1 <= len(applied_fields) <= 8:
        raise ValueError("Provide one to eight applied fields in A/m")
    if preconditioner not in ("reference-diagonal", "reference-chebyshev", "reference-chebyshev-graph"):
        raise ValueError("Unsupported explicit shared Newton preconditioner")
    if isinstance(order, bool) or order not in (1, 2):
        raise ValueError("HDiv order must be 1 or 2")
    import math
    for value in (tolerance, linear_tolerance):
        if isinstance(value, bool) or not math.isfinite(value) or not 0 < value <= 1e-6:
            raise ValueError("Explicit solve tolerances must be finite in (0,1e-6]")
    law_file = Path(law_path).resolve(strict=True)
    if law_file.stat().st_size > 2**20:
        raise ValueError("Inverse BH law exceeds input budget")
    law_bytes = law_file.read_bytes()
    law = json.loads(law_bytes.decode("utf-8"))
    if not isinstance(law, dict):
        raise ValueError("Expected an explicit BH material JSON object")
    output = Path(result_path).resolve()
    # Reserve the result exclusively; never overwrite another run's evidence.
    with output.open("x", encoding="utf-8"):
        pass
    record = dict(schema="radia.hdiv-gpu-workflow.v1", accepted=False,
                  backend="gpu", independent_cpu_certified=False, cases=[],
                  order=order, preconditioner=preconditioner,
                  tolerance=tolerance, linear_tolerance=linear_tolerance,
                  material_sha256=hashlib.sha256(law_bytes).hexdigest())

    def save():
        serialized = json.dumps(record, indent=2, allow_nan=False)
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8",
                    dir=output.parent, prefix=output.name + ".", suffix=".tmp",
                    delete=False) as stream:
                temporary = Path(stream.name)
                stream.write(serialized)
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(temporary, output)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)

    save()
    try:
        import ngsolve as ng
        from radia.hdiv_gpu import PreparedHDivNonlinear
        with ng.TaskManager():
            factory = PreparedHDivNonlinear
            material = law
            if law.get("schema") == "radia.bh-table.v1":
                factory = PreparedHDivNonlinear.from_bh_table
                material = law["table"]
            elif "schema" in law:
                raise ValueError("Unsupported explicit BH material schema")
            with factory(vol_path, material, backend=backend,
                    label_contract=label_contract, order=order) as prepared:
                for field in applied_fields:
                    result = prepared.solve(field, preconditioner=preconditioner,
                                            tolerance=tolerance, linear_tolerance=linear_tolerance)
                    # Normalize NumPy scalars/arrays without importing numerical
                    # packages on the cold MCP discovery path. Serialize before
                    # mutating the durable partial record.
                    case = json.loads(json.dumps(
                        dict(applied_field_A_per_m=field, **result),
                        default=lambda value: value.tolist(), allow_nan=False))
                    record["cases"].append(case)
                    save()
            record["closed"] = prepared.stats()["closed"]
        record["accepted"] = True
        save()
    except BaseException as exc:
        record["error"] = str(exc)
        save()
        raise
    return dict(accepted=True, backend="gpu", case_count=len(record["cases"]),
                result_path=str(output), closed=record["closed"], independent_cpu_certified=False)
