"""ESRF #6 HEX BDM1 charge-Gram definiteness gate (the regression that reproduced the iteration-86 CG breakdown).

Background (2026-09-05): ``vim.Solve`` on the ESRF example-6 quadrupole iron mesh (1648 Cubit HEX cells,
592 trilinear-distorted, BDM1) failed in the chi0 warmstart of the energy-Newton path -- the production
Jacobi CG hit ``p^T A p < 0`` at iteration 86 on hibino and LAB alike -- because the charge Gram
``N = B^T G B`` itself was indefinite: ``lambda_min(M^-1 N)`` between -5e-3 and -2e-3 against the physical
band [0, 1], while the table's initial permeability (mu_r 2001) gives a mass floor of only 5e-4.  The
defect is the HEX near family on distorted cells (touching pairs classified far, static-site radial inner,
graded outer rule shared with the far tensor product); the raw O(n^2) Gram is negative along the same
direction, so ACA is not the cause.

This lane rebuilds the production Gram on the real asset and checks, in order:

1. the exact 13-host element clusters of ``inv_chi0*M + N`` are positive semidefinite (M-metric);
2. the PRODUCTION Jacobi CG (``solve_configured_linear_material_auto_prec``) converges at the case-6
   chi0 floor and at 2x, 4x and 8x larger initial permeability (mu_r ~ 2001, 4001, 8001, 16001);
3. LOBPCG on the generalized problem ``N v = lambda M v`` reports ``lambda_min`` above ``--lambda-floor``
   and ``lambda_max`` below 1 + 1e-3;
4. along the LOBPCG minimizing direction the raw O(n^2) quadratic form and the H-matrix quadratic form
   agree to ``--compression-tolerance`` of the M-scale (the separate compression-error evaluation).

Any failed check exits non-zero.  Timings are recorded with the host name and are relative only unless
the host is an idle mdx/hibino.  Run::

    python validate_hex_gram_definiteness.py --assets-dir <dir with model.vol> --output results/hex_gram_definiteness_<host>.json
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "src" / "radia").exists())
if str(REPO / "src") not in sys.path:
    sys.path.insert(0, str(REPO / "src"))
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import ngsolve as ng  # noqa: E402
import scipy.linalg as sla  # noqa: E402
import scipy.sparse as sp  # noqa: E402
import scipy.sparse.linalg as spla  # noqa: E402

import radia  # noqa: E402
from radia import vim  # noqa: E402
from radia.esrf_examples import get_esrf_bh_table  # noqa: E402
from radia.vim import _vim as V  # noqa: E402
from radia.vim._nonlinear import _bh_inverse_funcs  # noqa: E402

SCHEMA = "radia.validation.esrf-hex-gram-definiteness.v1"
CASE = 6


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _element_dof_sets(mesh, fes):
    sets = []
    for el in mesh.Elements(ng.VOL):
        sets.append(sorted(set(int(d) for d in fes.GetDofNrs(el) if int(d) >= 0)))
    return sets


def _colour_classes(dof_sets):
    """Greedy colouring so that every class has disjoint DOF sets (the block extractor needs unique DOFs)."""
    owner: dict[int, set[int]] = {}
    classes: list[list[int]] = []
    for k, ds in enumerate(dof_sets):
        used = set().union(*[owner.get(d, set()) for d in ds]) if ds else set()
        c = 0
        while c in used:
            c += 1
        for d in ds:
            owner.setdefault(d, set()).add(c)
        while len(classes) <= c:
            classes.append([])
        classes[c].append(k)
    return classes


def cluster_check(G, mesh, fes, inv_chi0):
    dof_sets = _element_dof_sets(mesh, fes)
    gen_min = np.zeros(mesh.ne)
    gen_max = np.zeros(mesh.ne)
    for cls in _colour_classes(dof_sets):
        cand, offs = [], [0]
        for k in cls:
            cand.extend(dof_sets[k])
            offs.append(len(cand))
        cand = np.asarray(cand, np.int32)
        offs = np.asarray(offs, np.int32)
        blk_a = np.asarray(G.configured_linear_material_element_blocks(1.0, cand, offs), float)
        blk_n = np.asarray(G.configured_linear_material_element_blocks(0.0, cand, offs), float)
        pos = 0
        for j, k in enumerate(cls):
            w = offs[j + 1] - offs[j]
            a_e = blk_a[pos:pos + w * w].reshape(w, w)
            n_e = blk_n[pos:pos + w * w].reshape(w, w)
            pos += w * w
            n_e = 0.5 * (n_e + n_e.T)
            m_e = 0.5 * ((a_e - n_e) + (a_e - n_e).T) / inv_chi0
            vals = sla.eigh(n_e, m_e, eigvals_only=True)
            gen_min[k] = vals[0]
            gen_max[k] = vals[-1]
    return {"lambda_min": float(gen_min.min()), "lambda_max": float(gen_max.max()),
            "negative_clusters": int(np.sum(gen_min < -1.0e-10))}


def floor_scan(G, rhs, chi0, factors, tol, maxit):
    import re
    scan = []
    for s in factors:
        started = time.perf_counter()
        entry = {"floor_factor": float(s), "mu_r_equivalent": float(1.0 + chi0 / s)}
        try:
            res = G.solve_configured_linear_material_auto_prec(float(s), rhs, tol, maxit)
            entry.update(converged=True, iterations=int(res["iters"]), prec_min=float(res["prec_min"]),
                         prec_max=float(res["prec_max"]))
        except RuntimeError as exc:
            message = str(exc)
            match = re.search(r"iteration (\d+) -- p\^T A p = ([-+0-9.eE]+)", message)
            entry.update(converged=False, error=message[:200],
                         breakdown_iteration=int(match.group(1)) if match else None,
                         pAp=float(match.group(2)) if match else None)
        entry["wall_s"] = time.perf_counter() - started
        scan.append(entry)
        print(f"  floor x{s:g} (mu_r ~ {entry['mu_r_equivalent']:.0f}): "
              + (f"converged in {entry['iterations']} iterations" if entry["converged"]
                 else f"BREAKDOWN {entry.get('error', '')[:110]}"), flush=True)
        if not entry["converged"]:
            break
    return scan


def lobpcg_generalized(G, M_geom, n_face, k, maxiter, seed=6):
    """Edges of the generalized spectrum N v = lambda M v, preconditioned by the geometry mass Riesz map
    (the persistent PARDISO factor of M that the production CG also uses), which makes LOBPCG converge on
    the clustered low end where the unpreconditioned iteration stalls."""
    def apply_n(X):
        X = np.asarray(X, float)
        if X.ndim == 1:
            return np.asarray(G.apply_configured_demag(np.ascontiguousarray(X), True), float)
        return np.column_stack([np.asarray(G.apply_configured_demag(np.ascontiguousarray(X[:, j]), True), float)
                                for j in range(X.shape[1])])

    def apply_minv(X):
        X = np.asarray(X, float)
        if X.ndim == 1:
            return np.asarray(G.apply_configured_mass_riesz(np.ascontiguousarray(X)), float)
        return np.column_stack([np.asarray(G.apply_configured_mass_riesz(np.ascontiguousarray(X[:, j])), float)
                                for j in range(X.shape[1])])

    op_n = spla.LinearOperator((n_face, n_face), matvec=apply_n, matmat=apply_n, dtype=float)
    op_minv = spla.LinearOperator((n_face, n_face), matvec=apply_minv, matmat=apply_minv, dtype=float)
    rng = np.random.default_rng(seed)
    X0 = rng.standard_normal((n_face, k))
    started = time.perf_counter()
    vals, vecs = spla.lobpcg(op_n, X0, B=M_geom, M=op_minv, largest=False, maxiter=maxiter, tol=1.0e-8,
                             verbosityLevel=0)
    order = np.argsort(vals)
    vals, vecs = vals[order], vecs[:, order]
    vals_hi, _ = spla.lobpcg(op_n, rng.standard_normal((n_face, 2)), B=M_geom, M=op_minv, largest=True,
                             maxiter=maxiter, tol=1.0e-8, verbosityLevel=0)
    return vals, vecs, float(np.max(vals_hi)), time.perf_counter() - started


def mesh_conformity(mesh) -> dict:
    """Hanging facets / duplicated boundary faces (see ``radia.vim.mesh_conformity_report``)."""
    return vim.mesh_conformity_report(mesh)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--assets-dir", type=Path, required=True, help="directory holding the case-6 model.vol")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--gram-eps", type=float, default=1.0e-10)
    parser.add_argument("--cg-tol", type=float, default=1.0e-6, help="production warmstart tolerance")
    parser.add_argument("--cg-maxit", type=int, default=12000)
    parser.add_argument("--floor-factors", type=float, nargs="+", default=(1.0, 0.5, 0.25, 0.125))
    parser.add_argument("--lambda-floor", type=float, default=-1.0e-8)
    parser.add_argument("--lobpcg-vectors", type=int, default=4)
    parser.add_argument("--lobpcg-maxiter", type=int, default=300)
    parser.add_argument("--compression-tolerance", type=float, default=1.0e-6,
                        help="allowed |raw - hmatrix| quadratic-form gap along the minimizing direction, M-scale")
    parser.add_argument("--skip-raw", action="store_true", help="skip the O(n^2) raw quadratic form")
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument("--expect-iron-sha256", default=None)
    parser.add_argument("--allow-nonconforming", action="store_true",
                        help="diagnose a mesh with hanging nodes / duplicated faces instead of refusing it")
    options = parser.parse_args(argv)
    if options.threads > 0:
        ng.SetNumThreads(options.threads)

    mesh_path = (options.assets_dir / "model.vol").resolve()
    if not mesh_path.is_file():
        raise FileNotFoundError(mesh_path)
    sha = _sha256(mesh_path)
    if options.expect_iron_sha256 and sha != options.expect_iron_sha256:
        raise RuntimeError(f"iron mesh sha256 {sha} does not match the expected {options.expect_iron_sha256}")
    bh = np.asarray(get_esrf_bh_table(CASE), float)
    fields, _, _ = _bh_inverse_funcs(bh[:, 0], bh[:, 1])
    _, nd0 = fields(np.array([1.0e-12]))
    chi0 = 1.0 / max(float(nd0[0]), 1.0e-30)
    inv_chi0 = 1.0 / max(chi0, 1.0)
    report = {
        "schema": SCHEMA, "generated_at_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "host": platform.node(), "radia_version": getattr(radia, "__version__", None), "radia_file": radia.__file__,
        "case": CASE, "iron_mesh": str(mesh_path), "iron_mesh_sha256": sha,
        "chi0": chi0, "inv_chi0": inv_chi0, "gram_eps": options.gram_eps, "checks": {},
    }
    failures = []
    mesh = ng.Mesh(str(mesh_path))
    conformity = mesh_conformity(mesh)
    print(f"mesh: ne {mesh.ne} bnd {mesh.GetNE(ng.BND)} hanging facets {conformity['hanging_facets']} "
          f"duplicated boundary face pairs {conformity['duplicated_boundary_face_pairs']}", flush=True)
    if not conformity["conforming"] and not options.allow_nonconforming:
        raise RuntimeError("the iron mesh is non-conforming (%d hanging facets, %d duplicated boundary face "
                           "pairs): regenerate it with imprint/merge (export_esrf_cubit_assets) or pass "
                           "--allow-nonconforming to diagnose it" % (conformity["hanging_facets"],
                                                                     conformity["duplicated_boundary_face_pairs"]))
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=1)
        n_face = int(fes.ndof)
        affinity = V._hex_mapping_affinity_report(mesh)
        report["mesh"] = {"ne": int(mesh.ne), "ndof": n_face, "nonaffine_cells": int(affinity["nonaffine_cell_count"]),
                          "conformity": conformity}
        started = time.perf_counter()
        B, G, M_mass = vim.ChargeGram(fes, eps=options.gram_eps)
        report["gram_build_s"] = time.perf_counter() - started
        stats = dict(G.stats())
        report["gram_stats"] = {k: (float(v) if isinstance(v, (int, float, np.floating)) else str(v))
                                for k, v in stats.items() if str(k).startswith("hex_") or k in ("compression", "max_rank")}
        print(f"gram built in {report['gram_build_s']:.0f} s on {platform.node()} (n_face {n_face}, "
              f"near_inner_exact {stats.get('hex_near_inner_exact')}, glnear {stats.get('hex_glnear_n')})", flush=True)
        M_geom = sp.csr_matrix(M_mass)
        l2 = ng.L2(mesh, order=0)
        invchi = ng.GridFunction(l2)
        invchi.vec.FV().NumPy()[:] = inv_chi0
        u, v = fes.TnT()
        W = ng.BilinearForm(fes)
        W += invchi * u * v * ng.dx
        W.Assemble()
        G.configure_mass_matrix_ngsolve(W.mat)

        # 1. exact element clusters
        started = time.perf_counter()
        clusters = cluster_check(G, mesh, fes, inv_chi0)
        clusters["wall_s"] = time.perf_counter() - started
        report["checks"]["clusters"] = clusters
        print(f"clusters: lambda_min {clusters['lambda_min']:+.3e} lambda_max {clusters['lambda_max']:.4f} "
              f"negative {clusters['negative_clusters']}", flush=True)
        if clusters["negative_clusters"]:
            failures.append("element clusters lost positive semidefiniteness")

        # 2. production CG at the chi0 floor and beyond
        rng = np.random.default_rng(66)
        rhs = np.ascontiguousarray(M_geom @ rng.standard_normal(n_face))
        scan = floor_scan(G, rhs, chi0, options.floor_factors, options.cg_tol, options.cg_maxit)
        report["checks"]["floor_scan"] = scan
        if not all(entry["converged"] for entry in scan) or len(scan) != len(options.floor_factors):
            failures.append("production CG breakdown inside the floor scan")

        # 3. generalized spectrum edges
        vals, vecs, lam_max, lob_s = lobpcg_generalized(G, M_geom, n_face, options.lobpcg_vectors, options.lobpcg_maxiter)
        spectrum = {"lambda_min_ritz": [float(x) for x in vals], "lambda_max_ritz": lam_max, "wall_s": lob_s}
        report["checks"]["spectrum"] = spectrum
        print(f"lobpcg: lambda_min Ritz {vals[0]:+.3e} lambda_max Ritz {lam_max:.5f} ({lob_s:.0f} s)", flush=True)
        if vals[0] < options.lambda_floor:
            failures.append(f"lambda_min {vals[0]:.3e} below the floor {options.lambda_floor:.1e}")
        if lam_max > 1.0 + 1.0e-3:
            failures.append(f"lambda_max {lam_max:.5f} above the physical band")

        # 4. compression: raw O(n^2) versus H-matrix along the minimizing direction
        if not options.skip_raw:
            p = np.ascontiguousarray(vecs[:, 0])
            q = np.ascontiguousarray(np.asarray(sp.csr_matrix(B) @ p, float))
            m_scale = float(p @ (M_geom @ p))
            started = time.perf_counter()
            raw = float(G.raw_symmetric_quadratic_form(q))
            hm = float(q @ np.asarray(G.matvec_sym(q), float))
            gap = abs(raw - hm) / m_scale
            compression = {"raw_form_M_metric": raw / m_scale, "hmatrix_form_M_metric": hm / m_scale,
                           "gap_M_metric": gap, "wall_s": time.perf_counter() - started}
            report["checks"]["compression"] = compression
            print(f"compression along the minimizing direction: raw {raw / m_scale:+.3e} hmatrix {hm / m_scale:+.3e} "
                  f"gap {gap:.2e} (M-metric)", flush=True)
            if gap > options.compression_tolerance:
                failures.append(f"raw/H-matrix quadratic-form gap {gap:.2e} exceeds {options.compression_tolerance:.1e}")
    report["failures"] = failures
    report["passed"] = not failures
    if options.output is not None:
        options.output.parent.mkdir(parents=True, exist_ok=True)
        options.output.write_text(json.dumps(report, indent=1), encoding="utf-8")
        print("wrote", options.output, flush=True)
    print("PASSED" if not failures else "FAILED: " + "; ".join(failures), flush=True)
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
