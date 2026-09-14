"""Reidentify the committed analytical Potter-Schmulian fixture, without replacing it."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import platform
from pathlib import Path

import numpy as np


def verify(output):
    here = Path(__file__).resolve().parent
    root = here.parents[1]
    implementation = root / "src/radia/hysteresis_io.py"
    spec = importlib.util.spec_from_file_location("lineage_hysteresis_io", implementation)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    source = here / "fixtures/potter_schmulian/B_input.mat"
    generator = source.with_name("CASE_02.m")
    fixture = here / "binput_play_fixture.npz"
    loops, db, bmax = module.load_mat(source)
    eta, tables, bplay = module.build_shape_functions(loops, db)
    reconstructed = {"K": len(eta), "eta": eta, "chi": eta,
                     "dB": db, "Bplay": bplay,
                     "K_sub": len(eta[::2]), "eta_sub": eta[::2], "chi_sub": eta[::2]}
    for k, (r, f) in enumerate(tables):
        reconstructed.update({f"r_{k}": r, f"f_{k}": f})
    for k, (r, f) in enumerate(tables[::2]):
        reconstructed.update({f"r_sub_{k}": r, f"f_sub_{k}": f})
    with np.load(fixture, allow_pickle=False) as saved:
        assert set(saved.files) == set(reconstructed), "fixture field set changed"
        errors = {}
        for key, actual in reconstructed.items():
            actual = np.asarray(actual)
            assert np.isfinite(actual).all(), key
            np.testing.assert_allclose(actual, saved[key], rtol=1e-12, atol=1e-12,
                                       err_msg=key)
            errors[key] = float(np.max(np.abs(actual - saved[key])))
    report = {
        "schema": "radia.hysteresis-fixture-lineage.v1",
        "host": platform.node(), "python": platform.python_version(), "numpy": np.__version__,
        "source_kind": "analytical Potter-Schmulian branch data; not a measured specimen",
        "source_generator": "CASE_02.m uses S=0.9, Ms=1.91, Hc=500; fminsearch and makima interpolation",
        "verification_scope": "Reidentification from preserved MAT branches, not a rerun of MATLAB generation or the coupled field solve",
        "loops": len(loops), "dB_T": db, "Bmax_T": float(np.max(bmax)),
        "K": len(eta), "verified_fields": len(errors), "max_absolute_difference": max(errors.values()),
        "sha256": {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in [source, generator, fixture, implementation]},
        "passed": True,
    }
    Path(output).write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    print(json.dumps(verify(parser.parse_args().output), indent=2))
