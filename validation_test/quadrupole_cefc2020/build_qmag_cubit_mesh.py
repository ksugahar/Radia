"""Build conforming z-swept HEX meshes of the CEFC 2020 quadrupole iron with headless Cubit.

Usage:
    python build_qmag_cubit_mesh.py --output-dir <dir> [--sizes 0.010 0.006 0.004] [--order 2]
                                    [--cubit "C:/Program Files/Coreform Cubit 2025.12/bin/coreform_cubit.com"]

Every size writes ``qmag_h<mm>.vol`` plus its journal and Cubit log, runs the conformity survey
(no hanging facets, no duplicated boundary faces) and ``check-vol``, and records everything in
``mesh_manifest.json``.  A non-conforming or incomplete mesh is an error, never a warning.
"""
from __future__ import annotations

import argparse
import hashlib
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

import ngsolve as ng  # noqa: E402

import qmag_case as Q  # noqa: E402
from radia.vim import mesh_conformity_report  # noqa: E402

DEFAULT_CUBIT = Path(r"C:\Program Files\Coreform Cubit 2025.12\bin\coreform_cubit.com")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--sizes", type=float, nargs="+", default=(0.010, 0.006, 0.004), help="mesh sizes in metres")
    parser.add_argument("--order", type=int, default=2, help="curved-element order of the .vol export")
    parser.add_argument("--cubit", type=Path, default=DEFAULT_CUBIT)
    options = parser.parse_args(argv)
    cubit = options.cubit.resolve()
    if not cubit.is_file():
        raise FileNotFoundError(cubit)
    out = options.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    meshes = []
    for size in options.sizes:
        tag = f"qmag_h{size * 1000.0:g}"
        vol = out / f"{tag}.vol"
        journal = out / f"{tag}.jou"
        log = out / f"{tag}.cubit.log"
        Q.write_mesh_journal(journal, vol, size_m=size, order=options.order)
        Q.run_cubit(cubit, journal, log)
        mesh = ng.Mesh(str(vol))
        conformity = mesh_conformity_report(mesh)
        if not conformity["conforming"]:
            raise RuntimeError(f"{vol}: non-conforming mesh {conformity}")
        check_vol = shutil.which("check-vol")
        if check_vol is None:
            raise RuntimeError("check-vol (cubit-mesh-export) is not on PATH")
        check = subprocess.run([check_vol, str(vol), "--report-json", str(out / f"{tag}.vol-check.json")],
                               capture_output=True, text=True, check=False)
        if check.returncode != 0:
            raise RuntimeError(f"check-vol failed for {vol}:\n{check.stdout}\n{check.stderr}")
        families = {}
        for el in mesh.Elements(ng.VOL):
            families[str(el.type).rsplit(".", 1)[-1]] = families.get(str(el.type).rsplit(".", 1)[-1], 0) + 1
        entry = {"size_m": float(size), "vol": str(vol), "sha256": _sha256(vol), "journal": str(journal),
                 "ne": int(mesh.ne), "nbnd": int(mesh.GetNE(ng.BND)), "nv": int(mesh.nv), "families": families,
                 "materials": list(mesh.GetMaterials()), "boundaries": sorted(set(mesh.GetBoundaries())),
                 "conformity": conformity, "order": int(options.order)}
        print(f"{tag}: {entry['ne']} elements ({families}), {entry['nbnd']} boundary faces, conforming", flush=True)
        meshes.append(entry)
    manifest = {
        "schema": "radia.qmag-cefc2020-mesh-manifest.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(), "host": platform.node(),
        "cubit": str(cubit), "step": str(Q.STEP_PATH), "step_sha256": _sha256(Q.STEP_PATH),
        "iron_volume_m3_cad": Q.IRON_VOLUME_M3, "meshes": meshes,
    }
    (out / "mesh_manifest.json").write_text(json.dumps(manifest, indent=1), encoding="utf-8")
    print("wrote", out / "mesh_manifest.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
