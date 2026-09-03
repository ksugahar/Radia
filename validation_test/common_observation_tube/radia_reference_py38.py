"""Legacy Radia BEM oracle for the common finite-yoke observation tube.

This script is intentionally Python-3.8-only.  It loads the unmodified 2023
ESRF Radia extension instead of the repository's mesh-backed ``rad.Solve``
dispatcher, so the reference remains independent of HDiv-MMM.
"""

import argparse
import hashlib
import importlib.util
import json
import os

DEFAULT_LEGACY_RADIA = (
    r"S:\Radia\00_installer\2023_09_27_Radia-master\env\radia_python"
)
LEGACY_EXTENSION = os.environ.get(
    "RADIA_LEGACY_EXTENSION",
    os.path.join(DEFAULT_LEGACY_RADIA, "radia_py3_8_x64.pyd"),
)
if not os.path.isfile(LEGACY_EXTENSION):
    raise ImportError(
        "legacy Radia extension not found; set RADIA_LEGACY_EXTENSION to "
        "radia_py3_8_x64.pyd"
    )
_spec = importlib.util.spec_from_file_location("radia", LEGACY_EXTENSION)
if _spec is None or _spec.loader is None:
    raise ImportError("cannot load the legacy non-MPI Radia extension")
rad = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rad)


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _bar(center_mm, size_mm, divisions, material):
    bar = rad.ObjRecMag(center_mm, size_mm)
    rad.ObjDivMag(bar, divisions)
    rad.MatApl(bar, material)
    return bar


def solve_reference(points_m, mu_r, applied_b_t, segmentation):
    """Solve the four-bar finite rectangular yoke with classic Radia."""
    rad.UtiDelAll()
    chi = float(mu_r) - 1.0
    # The legacy API exposes the isotropic law through its anisotropic-shaped
    # signature.  An effectively zero remanence vector supplies the otherwise
    # mandatory easy axis without affecting this 0.1 T benchmark.
    material = rad.MatLin([chi, chi], [1.0e-30, 0.0, 0.0])
    transverse = int(segmentation)
    longitudinal = 2 * transverse
    bars = [
        _bar(
            [0, 15, 0],
            [40, 10, 60],
            [2 * transverse, transverse, longitudinal],
            material,
        ),
        _bar(
            [0, -15, 0],
            [40, 10, 60],
            [2 * transverse, transverse, longitudinal],
            material,
        ),
        _bar(
            [15, 0, 0], [10, 20, 60], [transverse, transverse, longitudinal], material
        ),
        _bar(
            [-15, 0, 0], [10, 20, 60], [transverse, transverse, longitudinal], material
        ),
    ]
    iron = rad.ObjCnt(bars)
    source = rad.ObjBckg([0.0, float(applied_b_t), 0.0])
    model = rad.ObjCnt([iron, source])
    degrees_of_freedom = int(rad.ObjDegFre(model))
    solve_result = list(rad.Solve(model, 1.0e-8, 5000))
    points_mm = [[1000.0 * value for value in point] for point in points_m]
    field = rad.Fld(model, "b", points_mm)
    if points_mm and not isinstance(field[0], (list, tuple)):
        field = [field]
    result = {
        "engine": "radia_reference_2023",
        "radia_extension_file": os.path.basename(LEGACY_EXTENSION),
        "radia_extension_sha256": _sha256(LEGACY_EXTENSION),
        "mu_r": float(mu_r),
        "applied_b_t": float(applied_b_t),
        "segmentation": transverse,
        "degrees_of_freedom": degrees_of_freedom,
        "solve_result": solve_result,
        "b_t": [[float(value) for value in row] for row in field],
    }
    rad.UtiDelAll()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-json", required=True)
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--segmentation", type=int, default=3)
    options = parser.parse_args()
    with open(options.input_json, "r") as stream:
        request = json.load(stream)
    result = solve_reference(
        request["points_m"],
        request["mu_r"],
        request["applied_b_t"],
        options.segmentation,
    )
    with open(options.output_json, "w") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(
        json.dumps(
            {
                "output_json": options.output_json,
                "degrees_of_freedom": result["degrees_of_freedom"],
                "solve_result": result["solve_result"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
