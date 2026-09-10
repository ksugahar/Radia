"""Validate the high-order TET charge Gram with mirror images on flat and curved meshes.

Matrix: geometry {flat_box, curved_sphere} x order {1 (BDM1), 2 (BDM2)}.  Each reduced
quarter model (x > 0, y > 0) carries ``image="-x-y"`` -- the quadrupole parity of the ESRF
quadrupole cases -- and is compared with its full model, which carries no image.

Checks per case
  symmetry : sampled |G_ab - G_ba| / max |G_aa| from the folded entry oracle (symmetric by
             construction, so this must be roundoff)
  on_plane : folded self entries of the cut-face charges lying on the mirror planes, which the
             antisymmetric image annihilates, relative to the median self entry -- roundoff
  psd      : the eigenvalues of the sigma-normalized dense Gram assembled from matvec_sym; the raw
             O(n^2) quadratic form on the minimizing vector confirms the sign without the H-matrix
  physics  : quarter+image demag field == full-model demag field at eight points spread over all
             four quadrants.  flat_box: the full mesh is the exact mirrored union of the quarter mesh,
             so the agreement is roundoff + ACA level; curved_sphere: independent meshes, so the
             agreement is discretization level
  legacy   : the same Gram built by the legacy per-entry scalar fold (subprocess with
             RADIA_HDIV_DISABLE_HO_IMAGE_BLOCK=1 RADIA_HDIV_DISABLE_HO_IMAGE_FAR=1): entry agreement on
             the sampled pairs, field agreement, and the build-time ratio

Build times are a relative, same-host smoke.  Decision-grade timing belongs on an idle mdx / hibino
(Benchmark Policy); the JSON records the host so nobody mistakes a LAB number for one.

Usage
  python validate_hdiv_vim_tet_image_dispatch.py            # full matrix, writes the JSON next to this file
  python validate_hdiv_vim_tet_image_dispatch.py --quick    # coarser meshes, no speed case
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "validate_hdiv_vim_tet_image_dispatch_results.json"

MU_R = 100.0
GRADIENT = 1.0e4          # A/m per m: H_ext = GRADIENT * (y, x, 0) = -grad(GRADIENT x y), odd in x and y
IMAGE = "-x-y"
GRAM_EPS = 1.0e-12
CG_TOL = 1.0e-12
PAIR_SAMPLES = 4000
SEED = 20260905
OBSERVATION_POINTS = np.array([
    [1.5, 0.5, 0.2], [-1.5, 0.5, 0.2],
    [0.5, -1.5, 0.2], [-0.5, -1.5, -0.2],
    [1.2, 1.2, 0.5], [-1.2, -1.2, -0.5],
    [0.3, 0.3, 1.6], [-0.3, 0.3, -1.6],
])
LEGACY_ENV = {
    "RADIA_HDIV_DISABLE_HO_IMAGE_BLOCK": "1",
    "RADIA_HDIV_DISABLE_HO_IMAGE_FAR": "1",
}
GATES = {
    "symmetry": 1.0e-12,
    "on_plane_residue": 1.0e-10,
    "normalized_min_eigenvalue": -1.0e-9,
    "physics_flat_box": 1.0e-8,
    "physics_curved_sphere": 3.0e-2,
    "legacy_field_agreement": 2.0e-3,
}


# --------------------------------------------------------------------------- geometry
def quarter_box_mesh(maxh):
    import ngsolve as ng
    from netgen.occ import Box, OCCGeometry, Pnt

    with ng.TaskManager():
        return ng.Mesh(OCCGeometry(Box(Pnt(0, 0, -1), Pnt(1, 1, 1))).GenerateMesh(maxh=maxh))


def mirrored_union(quarter, axes=(0, 1), plane_tol=1.0e-12):
    """Reflect a reduced netgen mesh into its full model, sharing the plane vertices.

    Every non-empty subset of ``axes`` reflects the quarter once; an odd number of reflections
    flips the orientation, which a vertex swap restores.  Boundary faces lying on a mirror plane
    become internal and are dropped.
    """
    import ngsolve as ng
    from netgen.meshing import Element2D, Element3D, FaceDescriptor, Mesh as NetgenMesh, MeshPoint, Pnt

    source = quarter.ngmesh
    coordinates = [tuple(float(v) for v in point.p) for point in source.Points()]
    full = NetgenMesh(dim=3)
    full.Add(FaceDescriptor(surfnr=1, domin=1, domout=0, bc=1))
    index = {}
    subsets = [()]
    for mask in range(1, 1 << len(axes)):
        subsets.append(tuple(axes[k] for k in range(len(axes)) if mask & (1 << k)))
    new_id = {}
    for subset in subsets:
        for old, xyz in enumerate(coordinates, start=1):
            reflected = list(xyz)
            for axis in subset:
                reflected[axis] = -reflected[axis]
            key = tuple(round(v, 12) + 0.0 for v in reflected)
            if key not in index:
                index[key] = full.Add(MeshPoint(Pnt(*reflected)))
            new_id[(subset, old)] = index[key]
    for subset in subsets:
        flip = len(subset) % 2 == 1
        for element in source.Elements3D():
            vertices = [new_id[(subset, v.nr)] for v in element.vertices]
            if flip:
                vertices[0], vertices[1] = vertices[1], vertices[0]
            full.Add(Element3D(1, vertices))
        for element in source.Elements2D():
            on_plane = any(
                all(abs(coordinates[v.nr - 1][axis]) <= plane_tol for v in element.vertices)
                for axis in axes)
            if on_plane:
                continue
            vertices = [new_id[(subset, v.nr)] for v in element.vertices]
            if flip:
                vertices[0], vertices[1] = vertices[1], vertices[0]
            full.Add(Element2D(1, vertices))
    full.SetMaterial(1, "iron")
    return ng.Mesh(full)


def quarter_sphere_mesh(maxh):
    import ngsolve as ng
    from netgen.occ import Box, OCCGeometry, Pnt, Sphere

    shape = Sphere(Pnt(0, 0, 0), 1.0) * Box(Pnt(0, 0, -2), Pnt(2, 2, 2))
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(shape).GenerateMesh(maxh=maxh))
        mesh.Curve(2)
    return mesh


def full_sphere_mesh(maxh):
    import ngsolve as ng
    from netgen.occ import OCCGeometry, Pnt, Sphere

    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0)).GenerateMesh(maxh=maxh))
        mesh.Curve(2)
    return mesh


def quadrupole_field():
    import ngsolve as ng

    return ng.CoefficientFunction((GRADIENT * ng.y, GRADIENT * ng.x, 0.0))


# --------------------------------------------------------------------------- worker pieces
def solve_case(mesh, order, curved, image):
    import ngsolve as ng
    from radia import vim

    started = time.perf_counter()
    with ng.TaskManager():
        result = vim.Solve(
            mesh, mu_r=MU_R, H_ext=quadrupole_field(), order=order,
            curve_order=2 if curved else None, gram_eps=GRAM_EPS, tol=CG_TOL, image=image)
        field = vim.FieldFromSolution(result, OBSERVATION_POINTS, algorithm="direct")
    return {
        "field_H_Am": np.asarray(field, float).tolist(),
        "demag": float(result["demag"]),
        "iters": int(result["iters"]),
        "ndof": int(result["ndof"]),
        "n_el": int(result["n_el"]),
        "charge_gram_wall_s": float(result.get("charge_gram_wall_s", float("nan"))),
        "hmat_stats": {k: (float(v) if isinstance(v, (int, float, np.floating)) else v)
                       for k, v in dict(result.get("hmat_stats") or {}).items()},
        "wall_s": time.perf_counter() - started,
    }


def plane_face_dofs(mesh, fes, axes=(0, 1), plane_tol=1.0e-12):
    import ngsolve as ng

    dofs = set()
    for element in mesh.Elements(ng.BND):
        coordinates = np.asarray([mesh.vertices[v.nr].point for v in element.vertices], float)
        if any(np.max(np.abs(coordinates[:, axis])) <= plane_tol for axis in axes):
            dofs.update(int(d) for d in fes.GetDofNrs(element) if int(d) >= 0)
    return dofs


def gram_case(mesh, order, curved, image):
    """Gram-level metrics on the reduced model: symmetry, on-plane residue, PSD, sampled entries."""
    import ngsolve as ng
    from radia import vim
    from radia.vim._image import image_group, parse_image_string

    planes = parse_image_string(image)
    masks, signs = [], []
    for axes, sign in image_group(planes):
        masks.append(int(sum(1 << axis for axis in axes)))
        signs.append(float(sign))
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=order)
        started = time.perf_counter()
        B, G, _ = vim.ChargeGram(
            fes, eps=GRAM_EPS, curve_order=2 if curved else None,
            image_masks=masks, image_signs=signs)
        build_wall_s = time.perf_counter() - started
        n = int(G.ndof())
        stats = {k: (float(v) if isinstance(v, (int, float, np.floating)) else v)
                 for k, v in dict(G.stats()).items()}

        # --- symmetry of the folded entry oracle on sampled pairs
        rng = np.random.default_rng(SEED)
        pairs = rng.integers(0, n, size=(PAIR_SAMPLES, 2))
        diagonal = np.array([G.entry(int(p), int(p)) for p in range(n)])
        scale = float(np.max(np.abs(diagonal)))
        entries_ab = np.array([G.entry(int(a), int(b)) for a, b in pairs])
        entries_ba = np.array([G.entry(int(b), int(a)) for a, b in pairs])
        symmetry = float(np.max(np.abs(entries_ab - entries_ba)) / scale)

        # --- on-plane cut-face charges: rows of B supported only on plane-face DOFs
        Bc = B.tocsr()
        plane_dofs = plane_face_dofs(mesh, fes)
        on_plane = []
        for row in range(n):
            columns = Bc.indices[Bc.indptr[row]:Bc.indptr[row + 1]]
            if len(columns) and all(int(c) in plane_dofs for c in columns):
                on_plane.append(row)
        off_plane = np.setdiff1d(np.arange(n), np.array(on_plane, int))
        median_diag = float(np.median(diagonal[off_plane]))
        on_plane_residue = (float(np.max(np.abs(diagonal[on_plane])) / median_diag)
                            if on_plane else 0.0)

        # --- dense sigma-normalized Gram from the symmetric H-matrix apply
        sigma = np.asarray(G.charge_sigma(), float)
        if sigma.size != n:
            raise RuntimeError("charge_sigma is not populated after the H-matrix build")
        dense = np.empty((n, n))
        unit = np.zeros(n)
        for j in range(n):
            unit[:] = 0.0
            unit[j] = 1.0
            dense[:, j] = np.asarray(G.matvec_sym(unit), float)
        normalized = dense / sigma[:, None] / sigma[None, :]
        asymmetry = float(np.max(np.abs(normalized - normalized.T)))
        eigenvalues, vectors = np.linalg.eigh(0.5 * (normalized + normalized.T))
        v_min = vectors[:, 0]
        raw_form = float(G.raw_symmetric_quadratic_form(v_min / sigma))
        hmat_form = float(v_min @ (normalized @ v_min))
    return {
        "n_charge": n,
        "ndof": int(fes.ndof),
        "n_el": int(mesh.ne),
        "build_wall_s": build_wall_s,
        "hmat_stats": stats,
        "symmetry": symmetry,
        "on_plane_charge_count": len(on_plane),
        "on_plane_residue": on_plane_residue,
        "normalized_min_eigenvalue": float(eigenvalues[0]),
        "normalized_max_eigenvalue": float(eigenvalues[-1]),
        "normalized_dense_asymmetry": asymmetry,
        "raw_quadratic_form_at_min_vector": raw_form,
        "hmat_quadratic_form_at_min_vector": hmat_form,
        "sampled_pairs": pairs.tolist(),
        "sampled_entries": entries_ab.tolist(),
        "entry_scale": scale,
    }


def build_only(mesh, order, curved, image):
    """Speed case: the reduced Gram build alone."""
    import ngsolve as ng
    from radia import vim
    from radia.vim._image import image_group, parse_image_string

    masks, signs = [], []
    for axes, sign in image_group(parse_image_string(image)):
        masks.append(int(sum(1 << axis for axis in axes)))
        signs.append(float(sign))
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=order)
        started = time.perf_counter()
        _, G, _ = vim.ChargeGram(
            fes, eps=GRAM_EPS, curve_order=2 if curved else None,
            image_masks=masks, image_signs=signs)
        build_wall_s = time.perf_counter() - started
        stats = {k: (float(v) if isinstance(v, (int, float, np.floating)) else v)
                 for k, v in dict(G.stats()).items()}
    return {"n_charge": int(G.ndof()), "n_el": int(mesh.ne), "build_wall_s": build_wall_s,
            "hmat_stats": stats}


def worker(args):
    os.environ.setdefault("RADIA_HDIV_BLOCK_CACHE_STATS", "1")
    import ngsolve as ng
    import radia

    if args.threads > 0:
        ng.SetNumThreads(args.threads)
    curved = args.geometry in ("curved_sphere", "speed_sphere")
    payload = {"geometry": args.geometry, "order": args.order, "maxh": args.maxh,
               "legacy": bool(args.legacy), "radia_file": radia.__file__}
    if args.geometry == "speed_sphere":
        payload["reduced_build"] = build_only(quarter_sphere_mesh(args.maxh), args.order, True, IMAGE)
    else:
        if args.geometry == "flat_box":
            quarter = quarter_box_mesh(args.maxh)
            full = None if args.skip_full else mirrored_union(quarter)
        else:
            quarter = quarter_sphere_mesh(args.maxh)
            full = None if args.skip_full else full_sphere_mesh(args.maxh)
        payload["reduced"] = solve_case(quarter, args.order, curved, IMAGE)
        payload["gram"] = gram_case(quarter, args.order, curved, IMAGE)
        if full is not None:
            payload["full"] = solve_case(full, args.order, curved, None)
    Path(args.out).write_text(json.dumps(payload), encoding="utf-8")
    return 0


# --------------------------------------------------------------------------- driver
def run_worker(geometry, order, maxh, out, *, legacy, skip_full, threads):
    command = [sys.executable, str(Path(__file__).resolve()), "--worker", "--geometry", geometry,
               "--order", str(order), "--maxh", str(maxh), "--out", str(out),
               "--threads", str(threads)]
    if legacy:
        command.append("--legacy")
    if skip_full:
        command.append("--skip-full")
    env = dict(os.environ)
    env["RADIA_HDIV_BLOCK_CACHE_STATS"] = "1"
    if legacy:
        env.update(LEGACY_ENV)
    completed = subprocess.run(command, env=env, capture_output=True, text=True)
    if completed.returncode != 0:
        raise RuntimeError(
            f"worker failed ({geometry}, order {order}, legacy={legacy}):\n"
            f"{completed.stdout[-4000:]}\n{completed.stderr[-4000:]}")
    return json.loads(Path(out).read_text(encoding="utf-8"))


def relative_field_error(reference, candidate):
    reference = np.asarray(reference, float)
    candidate = np.asarray(candidate, float)
    scale = np.maximum(np.linalg.norm(reference, axis=1),
                       GRADIENT * np.linalg.norm(OBSERVATION_POINTS[:, :2], axis=1))
    return float(np.max(np.linalg.norm(candidate - reference, axis=1) / scale))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--geometry", choices=("flat_box", "curved_sphere", "speed_sphere"))
    parser.add_argument("--order", type=int, default=1)
    parser.add_argument("--maxh", type=float, default=0.5,
                        help="matrix mesh size; the legacy curved BDM2 fold bounds it from below")
    parser.add_argument("--out", type=Path)
    parser.add_argument("--legacy", action="store_true")
    parser.add_argument("--skip-full", action="store_true")
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--speed-maxh", type=float, default=0.3,
                        help="curved quarter-sphere size for the BDM1 build-time A/B")
    parser.add_argument("--results", type=Path, default=RESULTS)
    args = parser.parse_args(argv)
    if args.worker:
        return worker(args)

    maxh = 0.6 if args.quick else args.maxh
    scratch = Path(os.environ.get("TEMP", ".")) / "radia_tet_image_dispatch"
    scratch.mkdir(parents=True, exist_ok=True)
    cases = []
    failures = []
    for geometry in ("flat_box", "curved_sphere"):
        for order in (1, 2):
            tag = f"{geometry}_bdm{order}"
            new = run_worker(geometry, order, maxh, scratch / f"{tag}_new.json",
                             legacy=False, skip_full=False, threads=args.threads)
            legacy = run_worker(geometry, order, maxh, scratch / f"{tag}_legacy.json",
                                legacy=True, skip_full=True, threads=args.threads)
            gram_new, gram_legacy = new["gram"], legacy["gram"]
            if gram_new["sampled_pairs"] != gram_legacy["sampled_pairs"]:
                raise RuntimeError("sampled pairs differ between the new and legacy workers")
            entries_new = np.asarray(gram_new["sampled_entries"])
            entries_legacy = np.asarray(gram_legacy["sampled_entries"])
            entry_delta = np.abs(entries_new - entries_legacy) / gram_new["entry_scale"]
            physics = relative_field_error(new["full"]["field_H_Am"], new["reduced"]["field_H_Am"])
            legacy_field = relative_field_error(new["reduced"]["field_H_Am"],
                                                legacy["reduced"]["field_H_Am"])
            counters = {k: gram_new["hmat_stats"].get(k, 0.0)
                        for k in ("ho_image_far_entries", "ho_image_block_entries",
                                  "ho_image_scalar_entries")}
            case = {
                "geometry": geometry, "order": order, "maxh": maxh,
                "n_el_reduced": gram_new["n_el"], "n_el_full": new["full"]["n_el"],
                "n_charge": gram_new["n_charge"], "ndof": gram_new["ndof"],
                "symmetry": gram_new["symmetry"],
                "on_plane_charge_count": gram_new["on_plane_charge_count"],
                "on_plane_residue": gram_new["on_plane_residue"],
                "normalized_min_eigenvalue": gram_new["normalized_min_eigenvalue"],
                "normalized_max_eigenvalue": gram_new["normalized_max_eigenvalue"],
                "raw_quadratic_form_at_min_vector": gram_new["raw_quadratic_form_at_min_vector"],
                "hmat_quadratic_form_at_min_vector": gram_new["hmat_quadratic_form_at_min_vector"],
                "legacy_normalized_min_eigenvalue": gram_legacy["normalized_min_eigenvalue"],
                "physics_reduced_vs_full": physics,
                "demag_reduced": new["reduced"]["demag"], "demag_full": new["full"]["demag"],
                "legacy_entry_delta_max": float(np.max(entry_delta)),
                "legacy_entry_delta_p99": float(np.percentile(entry_delta, 99.0)),
                "legacy_field_agreement": legacy_field,
                "build_wall_s_new": gram_new["build_wall_s"],
                "build_wall_s_legacy": gram_legacy["build_wall_s"],
                "build_speedup": gram_legacy["build_wall_s"] / max(gram_new["build_wall_s"], 1e-9),
                "image_dispatch_counters": counters,
                "solve_iters_reduced": new["reduced"]["iters"],
                "solve_iters_full": new["full"]["iters"],
            }
            gate_key = "physics_flat_box" if geometry == "flat_box" else "physics_curved_sphere"
            checks = {
                "symmetry": case["symmetry"] <= GATES["symmetry"],
                "on_plane_residue": case["on_plane_residue"] <= GATES["on_plane_residue"],
                "normalized_min_eigenvalue":
                    case["normalized_min_eigenvalue"] >= GATES["normalized_min_eigenvalue"],
                "raw_form_sign": case["raw_quadratic_form_at_min_vector"] >= GATES["normalized_min_eigenvalue"],
                "physics": physics <= GATES[gate_key],
                "legacy_field_agreement": legacy_field <= GATES["legacy_field_agreement"],
                "block_or_far_used": (counters["ho_image_far_entries"]
                                      + counters["ho_image_block_entries"]) > 0.0,
            }
            case["checks"] = checks
            case["passed"] = all(checks.values())
            if not case["passed"]:
                failures.append(tag)
            cases.append(case)
            print(f"{tag:22s} n_charge={case['n_charge']:6d} sym={case['symmetry']:.1e} "
                  f"on_plane={case['on_plane_residue']:.1e} lam_min={case['normalized_min_eigenvalue']:+.2e} "
                  f"physics={physics:.2e} legacy_dE={case['legacy_entry_delta_max']:.1e} "
                  f"legacy_dF={legacy_field:.1e} speedup={case['build_speedup']:.1f}x "
                  f"{'PASS' if case['passed'] else 'FAIL'}", flush=True)

    # The BDM2 build-time A/B lives in the matrix above: on a curved mesh the legacy
    # per-entry fold of six quadratic face modes is so slow that even the 0.5 quarter
    # sphere takes it a quarter of an hour, so a larger BDM2 speed case is pointless.
    speed = []
    speed_cases = [] if args.quick else [(1, args.speed_maxh)]
    for order, speed_maxh in speed_cases:
        new = run_worker("speed_sphere", order, speed_maxh, scratch / f"speed_bdm{order}_new.json",
                         legacy=False, skip_full=True, threads=args.threads)
        legacy = run_worker("speed_sphere", order, speed_maxh,
                            scratch / f"speed_bdm{order}_legacy.json",
                            legacy=True, skip_full=True, threads=args.threads)
        entry = {
            "order": order, "maxh": speed_maxh,
            "n_el": new["reduced_build"]["n_el"], "n_charge": new["reduced_build"]["n_charge"],
            "build_wall_s_new": new["reduced_build"]["build_wall_s"],
            "build_wall_s_legacy": legacy["reduced_build"]["build_wall_s"],
            "hmat_build_time_new": new["reduced_build"]["hmat_stats"].get("build_time"),
            "hmat_build_time_legacy": legacy["reduced_build"]["hmat_stats"].get("build_time"),
            "compression": new["reduced_build"]["hmat_stats"].get("compression"),
            "image_dispatch_counters": {
                k: new["reduced_build"]["hmat_stats"].get(k, 0.0)
                for k in ("ho_image_far_entries", "ho_image_block_entries", "ho_image_scalar_entries")},
        }
        entry["build_speedup"] = entry["build_wall_s_legacy"] / max(entry["build_wall_s_new"], 1e-9)
        speed.append(entry)
        print(f"speed_sphere_bdm{order}     n_el={entry['n_el']:5d} n_charge={entry['n_charge']:6d} "
              f"new={entry['build_wall_s_new']:.1f}s legacy={entry['build_wall_s_legacy']:.1f}s "
              f"speedup={entry['build_speedup']:.1f}x", flush=True)

    import ngsolve
    import radia

    payload = {
        "schema": "radia.validation.hdiv-vim-tet-image-dispatch.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "hostname": platform.node(),
        "timing_role": "relative same-host smoke; decision-grade timing belongs on idle mdx/hibino",
        "radia_version": getattr(radia, "__version__", None),
        "ngsolve_version": getattr(ngsolve, "__version__", None),
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
        "image": IMAGE,
        "mu_r": MU_R,
        "quadrupole_gradient_A_per_m2": GRADIENT,
        "gram_eps": GRAM_EPS,
        "gates": GATES,
        "legacy_environment": LEGACY_ENV,
        "observation_points_m": OBSERVATION_POINTS.tolist(),
        "cases": cases,
        "speed": speed,
        "passed": not failures,
    }
    args.results.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"results: {args.results}")
    if failures:
        raise SystemExit(f"FAILED cases: {', '.join(failures)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
