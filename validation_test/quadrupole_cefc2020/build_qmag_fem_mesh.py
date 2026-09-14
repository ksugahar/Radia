"""Full-domain FEM mesh (iron + physical air + Kelvin exterior) of the CEFC 2020 quadrupole.

Usage:
    python build_qmag_fem_mesh.py --output-dir <dir> [--cubit <coreform_cubit.com>]

Reuses the ESRF coil-yoke Kelvin builder (validation_test/esrf_three_engine/
build_esrf_coil_yoke_kelvin_mesh.py): the iron STEP is webcut at the z = 0 mesh seam, the positive
half is tet-meshed inside a physical air sphere with a bore refinement, reflected, and wrapped in the
periodic Kelvin exterior; ``validate_fem_mesh`` then checks the Kelvin identification and the label
contract.  The mesh is the reduced-A / mixed total-reduced Omega side of the Q-mag three-engine
comparison (``run_qmag_three_engine.py``).  Parameters follow ESRF example 6 (the same yoke size class).
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "src" / "radia").exists())
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(HERE))

import qmag_case as Q  # noqa: E402

ESRF_BUILDER = REPO / "validation_test" / "esrf_three_engine" / "build_esrf_coil_yoke_kelvin_mesh.py"
FEM_MESH_PARAMETERS = {
    "kelvin_radius_m": 0.16,          # yoke corner at 0.128 m
    "iron_size_m": 0.006,
    "air_size_m": 0.016,
    "kelvin_size_m": 0.035,
    "gap_size_m": 0.003,
    "gap_refinement_radius_m": 0.026,     # pole faces at 0.020 m
    "gap_refinement_half_length_m": 0.032,   # iron length +-0.030 m
    "beam_axis": 2,
    "curve_order": 2,
}


def _load_builder():
    spec = importlib.util.spec_from_file_location("_esrf_coil_yoke_builder", ESRF_BUILDER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {ESRF_BUILDER}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cubit", type=Path, default=None)
    options = parser.parse_args(argv)
    builder = _load_builder()
    cubit = (options.cubit or builder.DEFAULT_CUBIT).resolve()
    if not cubit.is_file():
        raise FileNotFoundError(cubit)
    out = options.output_dir.resolve()
    assets = out / "assets"
    assets.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(Q.STEP_PATH, assets / "iron.step")     # the builder expects <assets>/iron.step
    journal, vol = builder.write_journal(
        out, assets_dir=assets, kelvin_helper=builder.DEFAULT_KELVIN_HELPER.resolve(), **FEM_MESH_PARAMETERS)
    completed = subprocess.run([str(cubit), "-nographics", "-batch", "-nojournal", "-input", str(journal)],
                               cwd=str(out), check=False, capture_output=True, text=True)
    (out / "cubit.log").write_text(completed.stdout + completed.stderr, encoding="utf-8", errors="replace")
    if not vol.is_file():
        raise RuntimeError(f"Cubit did not create {vol}; returncode={completed.returncode}; see {out / 'cubit.log'}")
    report = builder.validate_fem_mesh(vol)
    report.update({
        "schema": "radia.qmag-cefc2020-fem-mesh.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(), "host": platform.node(),
        "iron_step_sha256": _sha256(Q.STEP_PATH), "vol": str(vol), "vol_sha256": _sha256(vol),
        "parameters": FEM_MESH_PARAMETERS, "cubit_returncode": int(completed.returncode),
    })
    report_path = out / "qmag_fem_kelvin.mesh.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({k: report[k] for k in report if k not in ("parameters",)}, indent=1, sort_keys=True))
    print("wrote", report_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
