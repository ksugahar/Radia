"""calc_equivalence_source -- Stage 2 CLI for the future NFS panel.

Layer 4 (headless computation) per the Cubit panel 4-layer architecture
(CLAUDE.md).  Extracts an equivalence-theorem near-field source from an
NGSolve FEM solution on a named BND surface, saves the NFS artifact,
and optionally evaluates the reconstructed external field at a list
of observation points.

This script does NOT import cubit (Layer 4 rule).  It accepts a .vol
mesh + a .sol GridFunction file (NGSolve native), plus a surface BND
label that must enclose all primary sources.

Usage:
    python calc_equivalence_source.py \\
        --vol model.vol \\
        --sol fem_H.sol \\
        --surface nfs_surface \\
        --omega 0 \\
        --output model.nfs.json \\
        [--obs-csv obs_points.csv --probe-out probed_field.csv]

Output (stdout, JSON):
    {
      "status": "OK" | "error",
      "vol": "model.vol",
      "surface_label": "nfs_surface",
      "n_faces": 12345,
      "artifact": "model.nfs.json",
      "artifact_bytes": 1234567,
      "omega": 0,
      "probe_points": 100,            # if --obs-csv provided
      "probe_out": "probed_field.csv" # if --obs-csv provided
    }

Stage promotion:
    Stage 1 (variables enumerated):
        --vol (.vol mesh, must have a CLOSED BND surface enclosing
               all sources, named via --surface)
        --sol (.sol GridFunction file, NGSolve native format)
        --surface (BND label, must enclose all sources)
        --omega (rad/s; 0 = static magnetostatic)
        --field-component (h | b | a -- which field the .sol stores)
        --mu-r (relative permeability OUTSIDE the source region; default 1)
        --output (NFS .json artifact path)
        --obs-csv (optional: CSV of observation points "x,y,z")
        --probe-out (optional: CSV of probed field at obs)

    Stage 2 (this CLI):
        argparse + JSON-on-stdout.  Validated by phase1_static_coil.py
        + phase2_wpt_harmonic.py golden tests in
        validation_test/equivalence_source/.

    Stage 3 (application integration):
        Wrap as an Electromagnet Simulink-block mode or a future NFS block.
        The mask selects the .vol and sideset, configures omega and the
        observation grid, and delegates execution to this CLI.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(
        description="Extract / evaluate an equivalence-theorem near-field source"
    )
    parser.add_argument("--vol", required=True, type=Path,
                         help="NGSolve .vol mesh file")
    parser.add_argument("--sol", required=True, type=Path,
                         help="NGSolve .sol GridFunction file")
    parser.add_argument("--surface", required=True,
                         help="BND label of the closed extraction surface "
                              "(must enclose all primary sources)")
    parser.add_argument("--omega", type=float, default=0.0,
                         help="Angular frequency (rad/s); 0 = static")
    parser.add_argument("--field-component", choices=["h", "b", "a"], default="b",
                         help="Which field the .sol stores: 'h' = H "
                              "(A/m), 'b' = B (T, will be divided by "
                              "mu_0 mu_r), 'a' = A (T*m, will be "
                              "curled to get B then H)")
    parser.add_argument("--mu-r", type=float, default=1.0,
                         help="Relative permeability outside source "
                              "(used when --field-component=b or a)")
    parser.add_argument("--output", type=Path, required=True,
                         help="Output NFS artifact (.nfs.json)")
    parser.add_argument("--obs-csv", type=Path, default=None,
                         help="CSV of observation points (x,y,z per row)")
    parser.add_argument("--probe-out", type=Path, default=None,
                         help="CSV output: x,y,z,Hx_re,Hx_im,Hy_re,...")
    args = parser.parse_args()

    # Late imports so --help works without NGSolve
    try:
        from ngsolve import Mesh, GridFunction, HCurl, H1, VectorH1, curl
        from ngsolve import TaskManager
    except ImportError as e:
        print(json.dumps({"status": "error",
                          "error": f"NGSolve not importable: {e}"}))
        return 2

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from equivalence_source import NearFieldSource, MU_0

    if not args.vol.exists():
        print(json.dumps({"status": "error",
                          "error": f"vol file not found: {args.vol}"}))
        return 1
    if not args.sol.exists():
        print(json.dumps({"status": "error",
                          "error": f"sol file not found: {args.sol}"}))
        return 1

    with TaskManager():
        mesh = Mesh(str(args.vol))
        # Validate surface label exists
        bnds = list(mesh.GetBoundaries())
        if args.surface not in bnds:
            print(json.dumps({
                "status": "error",
                "error": f"surface label '{args.surface}' not found",
                "available_boundaries": bnds,
            }))
            return 1

        # Load the .sol.  The .sol's FE space must MATCH what the user
        # solved with: HCurl for A, VectorH1 for B or H (3-component
        # vector), H1 for scalar (e.g. A_z in 2D).  We pick the FE space
        # by file size / field-component:
        #   --field-component=a : HCurl order >= 1
        #   --field-component=b | h : VectorH1 (3 components)
        #     (legacy scalar H1 path is reserved for explicit 2D scalar
        #     A_z .sol files and is not used in this v1 release)
        if args.field_component == "a":
            fes = HCurl(mesh, order=1)
            gf_A = GridFunction(fes)
            gf_A.Load(str(args.sol))
            # H = (1/(mu_0 mu_r)) curl(A)
            H_cf = (1.0 / (MU_0 * args.mu_r)) * curl(gf_A)
            E_cf = None  # static omega = 0 assumed when field is A only
        else:
            fes = VectorH1(mesh, order=1)
            gf = GridFunction(fes)
            gf.Load(str(args.sol))
            if args.field_component == "h":
                H_cf = gf
            else:
                # B / (mu_0 mu_r) = H
                H_cf = gf * (1.0 / (MU_0 * args.mu_r))
            E_cf = None

        # Extract NFS
        nfs = NearFieldSource.extract_ngsolve(
            mesh, gf_E=E_cf, gf_H=H_cf,
            surface_label=args.surface,
            omega=args.omega,
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        nfs.save(args.output)

        out = {
            "status": "OK",
            "vol": str(args.vol),
            "sol": str(args.sol),
            "surface_label": args.surface,
            "n_faces": int(nfs.n_faces),
            "omega": float(args.omega),
            "artifact": str(args.output),
            "artifact_bytes": int(args.output.stat().st_size),
        }

        # Optional: probe at observation points
        if args.obs_csv is not None:
            if not args.obs_csv.exists():
                print(json.dumps({"status": "error",
                                  "error": f"obs-csv not found: {args.obs_csv}"}))
                return 1
            obs = []
            with args.obs_csv.open("r") as fh:
                reader = csv.reader(fh)
                for row in reader:
                    if not row or row[0].startswith("#"):
                        continue
                    try:
                        obs.append([float(row[0]), float(row[1]), float(row[2])])
                    except (ValueError, IndexError):
                        continue
            obs = np.asarray(obs, dtype=np.float64)
            if obs.size == 0:
                print(json.dumps({"status": "error",
                                  "error": "obs-csv has no numeric rows"}))
                return 1

            E_rec, H_rec = nfs.evaluate(obs)

            out["probe_points"] = int(obs.shape[0])
            if args.probe_out is not None:
                args.probe_out.parent.mkdir(parents=True, exist_ok=True)
                with args.probe_out.open("w", newline="") as fh:
                    w = csv.writer(fh)
                    w.writerow(["x", "y", "z",
                                  "Hx_re", "Hx_im", "Hy_re", "Hy_im", "Hz_re", "Hz_im",
                                  "Ex_re", "Ex_im", "Ey_re", "Ey_im", "Ez_re", "Ez_im"])
                    for p, e_row, h_row in zip(obs, E_rec, H_rec):
                        w.writerow([p[0], p[1], p[2],
                                      h_row[0].real, h_row[0].imag,
                                      h_row[1].real, h_row[1].imag,
                                      h_row[2].real, h_row[2].imag,
                                      e_row[0].real, e_row[0].imag,
                                      e_row[1].real, e_row[1].imag,
                                      e_row[2].real, e_row[2].imag])
                out["probe_out"] = str(args.probe_out)

        print(json.dumps(out))
        return 0


if __name__ == "__main__":
    sys.exit(main())
