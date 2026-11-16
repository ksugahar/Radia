"""Axisymmetric volume PEEC of the revolved beak fin against the ring FEM.

The 5% fin adjudication (user decision 2026-11-01, limits in
``fin_return_principle.PHYSICAL_TARGETS``) asks whether a volume PEEC, which
resolves the conductor interior instead of imposing a surface impedance,
reproduces an independent interior-resolved reference.  The reference is the
A_phi ring FEM with an exact Dirichlet-to-Neumann exterior
(``radia.eddy_axisym_ring``), whose drive is the loop voltage.

Under that axisymmetric drive the three-dimensional volume PEEC of the
revolved solid has an azimuthally uniform current, so each meridian cell is
one ring loop and the circuit is exact with ring partial elements
(native kernel ``radia.peec_axisym``, C++ ``rad_peec_axisym``):

* P0 azimuthal current density on each meridian triangle ``k``;
* ``R_k = 2 pi r_c / (sigma A_k)`` (r is linear on the triangle);
* ``M_kl`` = the coaxial-filament mutual inductance averaged over both
  triangles; self and near pairs subtract ``mu0 r (ln(8 r) - 2 - ln|x-y|)``
  and integrate its logarithm in closed form;
* every loop sees the same loop voltage, ``(R + j omega M) I = V``.

Shared with the reference: the section, sigma, frequency, the loop-voltage
port and the region cuts.  Independent: meshes, current space (P0 integral
equation vs H1 A_phi FEM), assembly, exterior (free-space kernel vs DtN) and
solver.  This settles the volume-PEEC physics at 150 kHz on this section; it
does not exercise the three-dimensional Cartesian partial elements of
``radia.volume_peec``.

Quantities per 1 A peak: R, loop L, root heat (1.6 <= x < 3.5 mm) and tip heat
(x >= 3.5 mm), x measured from the ring radius.  Mapping the user's "root" and
"tip" onto these cuts is a proposal, not a criterion change.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))

MU0 = 4e-7*math.pi
ROOT_X_MM, TIP_X_MM = 1.6, 3.5

def section_mesh(radius_mm, boundary_h_mm, interior_h_mm, grading):
    """Meridian triangles of the tracked beak section, split at the region cuts."""
    from build123d import Plane, export_step
    from make_beak_fin_step import _beak_face
    from netgen.occ import Box, Glue, OCCGeometry, Pnt
    import ngsolve as ng
    import tempfile
    face = _beak_face(Plane.XY, radial_offset=radius_mm)
    with tempfile.TemporaryDirectory(prefix="peec_section_") as tmp:   # keep the checkout clean
        step = Path(tmp) / "section.step"
        export_step(face, str(step))
        section = OCCGeometry(str(step), dim=2).shape.Scale(Pnt(0, 0, 0), 1e-3)
    cuts = [(radius_mm + ROOT_X_MM)*1e-3, (radius_mm + TIP_X_MM)*1e-3]
    big = 10.0

    def slab(x0, x1):
        return Box(Pnt(x0, -big, -big), Pnt(x1, big, big))
    parts = [section*slab(-big, cuts[0]), section*slab(cuts[0], cuts[1]), section*slab(cuts[1], big)]
    for part, name in zip(parts, ("body", "root", "tip")):
        part.faces.name = name
        part.faces.maxh = interior_h_mm*1e-3
    shape = Glue(parts)
    shape.edges.maxh = boundary_h_mm*1e-3
    mesh = ng.Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=interior_h_mm*1e-3, grading=grading))
    points = np.array([mesh[v].point for v in mesh.vertices])
    tris = np.array([[v.nr for v in el.vertices] for el in mesh.Elements(ng.VOL)])
    region = np.array([el.mat for el in mesh.Elements(ng.VOL)])
    if set(region) != {"body", "root", "tip"}:
        raise RuntimeError(f"unexpected section regions {sorted(set(region))}")
    return points, tris, region


def _peec_runtime():
    """Identity of the native PEEC kernel that produced the candidate."""
    from radia import peec_matrices
    path = Path(peec_matrices.__file__)
    return dict(peec_matrices=str(path), peec_matrices_sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def solve_peec(points, tris, region, frequency, sigma, threads):
    """The native axisymmetric volume PEEC (radia.peec_axisym) on one meridian mesh."""
    import ngsolve as ng
    from radia.peec_axisym import solve_axisym_ring_peec
    ng.SetNumThreads(threads)
    with ng.TaskManager():
        out = solve_axisym_ring_peec(points, tris, frequency, sigma, region=region)
    regional = out["region_loss_W_per_A2"]
    total = float(out["cell_loss_W_per_A2"].sum())
    return dict(cells=out["cells"], R_ohm=out["R_ohm"], L_H=out["L_H"], power_balance=out["power_balance"],
                root_heat_W_per_A2=regional["root"], tip_heat_W_per_A2=regional["tip"],
                beak_loss_fraction=(regional["root"] + regional["tip"])/total,
                tip_loss_fraction=regional["tip"]/total,
                true_relative_residual=out["true_relative_residual"], assembly_s=out["assembly_s"],
                solve_s=out["solve_s"], min_cell_mm=float(np.sqrt(2*out["area_m2"].min())*1e3), threads=threads)


def _branch_ring_solver():
    """The corrected ring solver of this checkout on the installed radia natives.

    The installed package predates the A_phi drive fix, and the checkout has no
    built natives, so only ``eddy_axisym_ring`` is taken from the checkout.
    """
    import importlib.util
    import radia
    name = "radia.eddy_axisym_ring"
    path = ROOT/"src"/"radia"/"eddy_axisym_ring.py"
    if name not in sys.modules or Path(sys.modules[name].__file__).resolve() != path.resolve():
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    from radia import axifem
    axifem_path = Path(axifem.__file__)
    return sys.modules[name].solve_axisym_ring, dict(
        natives=str(Path(radia.__file__).parent), version=radia.__version__,
        axifem=str(axifem_path), axifem_sha256=hashlib.sha256(axifem_path.read_bytes()).hexdigest(),
        ring_solver=str(path.relative_to(ROOT)), ring_solver_sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def solve_level(mesh, ext_radius, ext_z, radius_mm, frequency, sigma, order, modes=30):
    """One reference solve on ``mesh``; returns the report row and the ring result."""
    solve_axisym_ring, _ = _branch_ring_solver()
    t0 = time.perf_counter()
    res = solve_axisym_ring(mesh, conductor="conductor", frequency_hz=frequency, sigma=sigma, order=order,
                            exterior_radius=ext_radius, exterior_center_z=ext_z, exterior_modes=modes)
    seconds = time.perf_counter() - t0
    total = res.total_loss()
    beak = res.loss_beyond((radius_mm + ROOT_X_MM)*1e-3)
    tip = res.loss_beyond((radius_mm + TIP_X_MM)*1e-3)
    scale = 0.5*res.resistance_ohm/total
    row = dict(order=order, ndof=int(res.ndof), elements=int(mesh.ne), seconds=seconds,
               R_ohm=res.resistance_ohm, L_H=res.inductance_H,
               root_heat_W_per_A2=float((beak - tip)*scale), tip_heat_W_per_A2=float(tip*scale),
               beak_loss_fraction=float(beak/total), tip_loss_fraction=float(tip/total),
               terminal_vs_energy_R=float(abs(res.terminal_resistance_ohm - res.resistance_ohm)
                                          /res.resistance_ohm))
    return row, res


def zz_indicator(res, order):
    """Per-element ZZ recovery indicator of grad A_phi in the revolved volume measure.

    The conductor loss is |J|^2/(2 sigma) with J = sigma (1/(2 pi r) - j omega A_phi), so its
    discretization error follows the error in A_phi, which is driven by the unresolved gradient.
    The gradient is recovered by averaging into a continuous vector H1 space of the same order.
    """
    import ngsolve as ng
    mesh = res.mesh
    d = ng.grad(res.gf_a)
    recovered = ng.GridFunction(ng.VectorH1(mesh, order=order, complex=True))
    recovered.Set(d)
    e = d - recovered
    density = 2*math.pi*ng.x/MU0*(e[0]*ng.Conj(e[0]) + e[1]*ng.Conj(e[1])).real
    eta2 = ng.Integrate(density, mesh, element_wise=True, order=2*order + 2)
    return np.maximum(np.asarray(eta2, dtype=float), 0.0)


def mark_whole_mesh(mesh, eta2, theta):
    """Doerfler marking over every element, conductor and air alike.

    Conductor-only marking stalled order 1 at a biased value because the error left in the
    unrefined air was never marked (user decision 2026-11-16: mark the whole mesh).
    Returns the flags and the marked conductor and air counts.
    """
    import ngsolve as ng
    cond = np.array([el.mat == "conductor" for el in mesh.Elements(ng.VOL)])
    if eta2.sum() <= 0:
        raise RuntimeError("no error to mark")
    ranked = np.argsort(eta2)[::-1]
    keep = np.searchsorted(np.cumsum(eta2[ranked]), theta*eta2.sum()) + 1
    flags = np.zeros(len(eta2), dtype=bool)
    flags[ranked[:keep]] = True
    return flags, int((flags & cond).sum()), int((flags & ~cond).sum())


def level_ndof(mesh, order):
    """DOFs of the reference space on ``mesh``, known before the solve."""
    from radia.axifem import H1Henrotte
    return int(H1Henrotte(mesh, order=order, complex=True, dirichlet="axis").ndof)


def adaptive_sequence(radius_mm, frequency, sigma, order, factor, start_maxh, air_mm, air_divisor,
                      theta, tol, max_ndof, max_levels, budget_s, clock=time.monotonic):
    """Refine one (order, exterior) reference adaptively until two successive levels agree within tol.

    The budget is checked on actual quantities: the deadline before every solve, the DOF count of
    the next level before it is solved, and the deadline again before a checkpoint may be accepted.
    A level that finishes past the deadline is kept in the record but cannot converge the sequence.
    """
    import ngsolve as ng
    from aphi_beak_section_sweep import build_section_ring
    if not (isinstance(max_levels, int) and not isinstance(max_levels, bool) and max_levels >= 1):
        raise ValueError(f"max_levels must be a positive integer, got {max_levels!r}")
    if not (math.isfinite(max_ndof) and max_ndof > 0 and math.isfinite(budget_s) and budget_s > 0):
        raise ValueError(f"budget must be positive and finite, got max_ndof={max_ndof!r} budget_s={budget_s!r}")
    started = clock()
    deadline = started + budget_s
    mesh, ext_radius, ext_z = build_section_ring(radius_mm, maxh_conductor_mm=start_maxh, maxh_air_mm=air_mm,
                                                exterior_factor=factor, air_divisor=air_divisor)
    # A Doerfler step can touch a few percent of the elements, and two such neighbours agree trivially.
    # Successive levels are therefore compared as checkpoints at least `growth` times apart in DOFs.
    growth_min = 2.0
    rows, checkpoints, stop = [], [], None
    while stop is None:
        if clock() >= deadline:
            stop = "budget: deadline"
            break
        ndof = level_ndof(mesh, order)
        if ndof > max_ndof:
            stop = f"budget: next level {ndof} DOFs above {max_ndof:.6g}"
            break
        if clock() >= deadline:                     # counting the DOFs may itself take the budget
            stop = "budget: deadline"
            break
        row, res = solve_level(mesh, ext_radius, ext_z, radius_mm, frequency, sigma, order)
        if row["ndof"] != ndof:
            raise RuntimeError(f"solved {row['ndof']} DOFs, admitted {ndof}")
        eta2 = zz_indicator(res, order)
        del res
        now = clock()                               # one sample for the record and the decision
        in_budget = now <= deadline
        row.update(level=len(rows), exterior_factor=factor, indicator=float(math.sqrt(eta2.sum())),
                   elapsed_s=now - started, within_budget=in_budget)
        is_check = in_budget and (not checkpoints or row["ndof"] >= growth_min*rows[checkpoints[-1]]["ndof"])
        row["checkpoint"] = is_check
        if is_check and checkpoints:
            row["change"] = relative_change(row, rows[checkpoints[-1]], tol)
        rows.append(row)
        print(json.dumps({"phase": "reference", **row}), flush=True)
        if not in_budget:
            stop = "budget: deadline"
            break
        if is_check:
            checkpoints.append(len(rows) - 1)
            if len(checkpoints) >= 2 and all(row["change"][k] <= tol[k] for k in tol):
                stop = "converged"
                break
        if len(rows) >= max_levels:
            stop = "budget: level limit"
            break
        flags, n_cond, n_air = mark_whole_mesh(mesh, eta2, theta)
        row.update(marked_conductor=n_cond, marked_air=n_air)
        for el in mesh.Elements(ng.VOL):
            mesh.SetRefinementFlag(el, bool(flags[el.nr]))
        mesh.Refine()
    end = clock()                                   # one final sample for the record and the decision
    if stop == "converged" and end > deadline:
        stop = "budget: deadline after the last checkpoint"
    return dict(order=order, exterior_factor=factor, start_maxh_mm=start_maxh, air_maxh_mm=air_mm,
                air_divisor=air_divisor, theta=theta, checkpoint_dof_growth=growth_min,
                tolerance=tol, max_ndof=max_ndof, max_levels=max_levels, budget_s=budget_s,
                elapsed_s=end - started, converged=stop == "converged", stop=stop, levels=rows)


def relative_change(new, old, keys):
    """Change of each quantity between two levels, relative to the newer (finer) level."""
    return {k: abs(quantities(new)[k] - quantities(old)[k])/abs(quantities(new)[k]) for k in keys}


def certify_sequence(seq, expected):
    """Re-derive convergence of a recorded sequence from its levels and the expected settings.

    Returns the reasons the sequence does not qualify (empty when it does).
    """
    bad = []
    for key, want in expected.items():
        if seq.get(key) != want:
            bad.append(f"{key} {seq.get(key)!r} != {want!r}")
    levels = seq.get("levels")
    if not isinstance(levels, list) or not levels:
        return bad + ["no recorded levels"]
    # Structure first: an unordered NaN or a non-integer would slip through every comparison below.
    is_count = lambda v: isinstance(v, int) and not isinstance(v, bool)
    is_time = lambda v: isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) and v >= 0
    structural = []
    for i, r in enumerate(levels):
        if not (is_count(r.get("ndof")) and r["ndof"] > 0):
            structural.append(f"level {i} ndof {r.get('ndof')!r} is not a positive integer")
        if not is_count(r.get("level")):
            structural.append(f"level {i} number {r.get('level')!r} is not an integer")
        if not is_time(r.get("elapsed_s")):
            structural.append(f"level {i} elapsed_s {r.get('elapsed_s')!r} is not finite and nonnegative")
        q = {k: r.get(v) for k, v in (("R", "R_ohm"), ("partial_L", "L_H"), ("root_heat", "root_heat_W_per_A2"),
                                       ("tip_heat", "tip_heat_W_per_A2"))}
        if not all(isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) and v != 0
                   for v in q.values()):
            structural.append(f"level {i} quantities {q} are not finite and nonzero")
    if not is_time(seq.get("elapsed_s")):
        structural.append(f"sequence elapsed_s {seq.get('elapsed_s')!r} is not finite and nonnegative")
    if structural:
        return bad + structural
    times = [r["elapsed_s"] for r in levels]
    if any(b < a for a, b in zip(times, times[1:])):
        bad.append(f"level elapsed times {times} decrease")
    if seq["elapsed_s"] < times[-1]:
        bad.append(f"sequence elapsed_s {seq['elapsed_s']} before its last level {times[-1]}")
    if seq["elapsed_s"] > expected["budget_s"]:
        bad.append(f"sequence elapsed_s {seq['elapsed_s']} above the budget {expected['budget_s']}")
    tol = expected["tolerance"]
    if seq.get("stop") != "converged" or seq.get("converged") is not True:
        bad.append(f"stop {seq.get('stop')!r}")
    if [r.get("level") for r in levels] != list(range(len(levels))):
        bad.append("levels are not numbered 0..n-1")
    if len(levels) > expected["max_levels"]:
        bad.append(f"{len(levels)} levels above {expected['max_levels']}")
    for r in levels:
        if r["ndof"] > expected["max_ndof"]:
            bad.append(f"level {r.get('level')} has {r['ndof']} DOFs above {expected['max_ndof']:.6g}")
        if r["elapsed_s"] > expected["budget_s"] or r.get("within_budget") is not True:
            bad.append(f"level {r.get('level')} past the wall-clock budget")
    cps = [r for r in levels if r.get("checkpoint")]
    for a, b in zip(cps, cps[1:]):
        if b["ndof"] < expected["checkpoint_dof_growth"]*a["ndof"]:
            bad.append(f"checkpoints {a.get('level')}->{b.get('level')} less than "
                       f"{expected['checkpoint_dof_growth']}x apart in DOFs")
    if len(cps) < 2:
        bad.append("fewer than two checkpoints")
    else:
        change = relative_change(cps[-1], cps[-2], tol)
        if not all(change[k] <= tol[k] for k in tol):
            bad.append(f"last checkpoints differ by {change}")
        if levels[-1] is not cps[-1]:
            bad.append("the sequence does not end on its converged checkpoint")
    return bad


def quantities(row):
    return dict(R=row["R_ohm"], partial_L=row["L_H"], root_heat=row["root_heat_W_per_A2"],
                tip_heat=row["tip_heat_W_per_A2"])


SEQUENCE_SCHEMA = "radia.ih.axisym-ring-adaptive-sequence.v1"


def script_sha256():
    """Hash of this script with normalized line endings, so CRLF and LF checkouts agree."""
    return hashlib.sha256(Path(__file__).read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def git_state():
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain"], cwd=ROOT, capture_output=True,
                                text=True).stdout.strip())
    return head, dirty


def expected_settings(args, settings, order, tol):
    """The settings a sequence of this order must carry, for saved and live sequences alike."""
    s, a, d = settings[order]
    return dict(start_maxh_mm=s, air_maxh_mm=a, air_divisor=d, theta=args.ref_theta,
                checkpoint_dof_growth=2.0, tolerance=tol, max_ndof=args.ref_max_ndof,
                max_levels=args.ref_max_levels, budget_s=3600*args.ref_hours)


def load_sequences(paths, args, provenance, settings, tol):
    """Admit saved sequences only from this source, script and runtime, with the configured settings.

    Every configured (order, exterior) must appear exactly once, and each sequence's convergence is
    re-derived from its levels (certify_sequence).  Any failure raises, naming every reason.
    """
    head, _ = git_state()
    script = script_sha256()
    want = {(o, f) for o in args.ref_order for f in args.ref_exterior_factor}
    seen, sequences, problems = set(), [], []
    for p in paths:
        doc = json.loads(Path(p).read_text(encoding="utf-8"))
        bad = []
        for key, value in (("schema", SEQUENCE_SCHEMA), ("git_head", head), ("git_dirty", False),
                           ("script_sha256", script), ("reference_runtime", provenance),
                           ("radius_mm", args.radius_mm), ("frequency_hz", args.frequency_hz),
                           ("sigma_S_per_m", args.sigma)):
            if doc.get(key) != value:
                bad.append(f"{key} {doc.get(key)!r} != {value!r}")
        seq = doc.get("sequence") or {}
        key = (seq.get("order"), seq.get("exterior_factor"))
        if key not in want:
            bad.append(f"(order, exterior) {key} not configured")
        elif key in seen:
            bad.append(f"(order, exterior) {key} given twice")
        else:
            seen.add(key)
            bad += certify_sequence(seq, expected_settings(args, settings, key[0], tol))
        if bad:
            problems.append(f"{p}: " + "; ".join(bad))
        sequences.append(seq)
    if seen != want and not problems:
        problems.append(f"sequences {sorted(seen)} do not cover orders x exteriors {sorted(want)}")
    if problems:
        raise ValueError("saved sequences refused:\n" + "\n".join(problems))
    return sequences


def resolve_reference(sequences, orders, exteriors, expected_of, targets):
    """Certify every sequence and resolve the reference from the final checkpoints.

    The reference is the final checkpoint of the highest order at the widest exterior. Its checks compare the last two
    checkpoints of every sequence (mesh), every exterior at the highest order and every order at the widest exterior,
    each relative to the finer or reference level; the reference is resolved when every sequence is certified, both
    orders and exteriors are swept, the terminal Re(V/I) agrees with the dissipation within 1e-3 and every check is
    within one tenth of its target. Saved verdicts are re-derived by the same function.
    """
    by = {(s["order"], s["exterior_factor"]): s for s in sequences}
    want = {(o, f) for o in orders for f in exteriors}
    if set(by) != want or len(by) != len(sequences):
        raise ValueError(f"sequences {sorted(by)} do not cover orders x exteriors {sorted(want)} once each")
    top, wide = max(orders), max(exteriors)
    certification = {f"o{s['order']}_f{s['exterior_factor']}": certify_sequence(s, expected_of(s["order"]))
                     for s in sequences}
    converged = all(s["converged"] for s in sequences) and not any(certification.values())
    swept = len(set(orders)) >= 2 and len(set(exteriors)) >= 2
    # A sequence stopped by its budget before any checkpoint leaves no reference to compare with.
    missing = sorted(f"o{o}_f{f}" for (o, f), s in by.items() if not any(r["checkpoint"] for r in s["levels"]))
    out = dict(certification=certification, converged=converged, swept=swept, missing=missing, reference=None,
               reference_level=None, checks=None, spread=None, consistent=False, resolved=False)
    if missing:
        return out
    final = {key: [r for r in s["levels"] if r["checkpoint"]][-1] for key, s in by.items()}
    best = final[(top, wide)]
    ref_q = quantities(best)
    rel = lambda a, b: {k: abs(quantities(a)[k] - quantities(b)[k])/abs(quantities(b)[k]) for k in ref_q}
    checks = {}
    for (o, f), s in sorted(by.items()):            # mesh: the last two adaptive checkpoints
        lv = [r for r in s["levels"] if r["checkpoint"]]
        # relative to the finer level, as in the stopping test
        checks[f"mesh_o{o}_f{f}"] = rel(lv[-2], lv[-1]) if len(lv) >= 2 else {k: math.inf for k in ref_q}
    for f in exteriors:                              # exterior, at the highest order
        checks[f"exterior_f{f}"] = rel(final[(top, f)], best)
    for o in orders:                                 # order, at the widest exterior
        checks[f"order_o{o}"] = rel(final[(o, wide)], best)
    spread = {k: max(c[k] for c in checks.values()) for k in ref_q}
    consistent = best["terminal_vs_energy_R"] <= 1e-3
    resolved = swept and consistent and converged and all(v <= 0.1*targets[k] for k, v in spread.items())
    out.update(reference=ref_q, reference_level=best, checks=checks, spread=spread, consistent=consistent,
               resolved=resolved)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radius-mm", type=float, default=200.0)
    parser.add_argument("--frequency-hz", type=float, default=150e3)
    parser.add_argument("--sigma", type=float, default=5.8e7)
    parser.add_argument("--peec-boundary-h-mm", type=float, nargs="+", default=[0.06, 0.04, 0.028])
    parser.add_argument("--peec-interior-h-mm", type=float, default=0.5)
    parser.add_argument("--grading", type=float, default=0.25)
    # Adaptive reference (user decision 2026-11-16): each (order, exterior) sequence starts from a coarse
    # mesh and refines the elements (conductor and air) carrying half of the ZZ indicator,
    # until two successive checkpoints agree within one tenth of every limit, under a DOF,
    # level and wall-clock budget.  A sequence that runs out of budget leaves the reference unresolved.
    parser.add_argument("--ref-order", type=int, nargs="+", default=[1, 2])
    parser.add_argument("--ref-exterior-factor", type=float, nargs="+", default=[1.5, 3.0])
    parser.add_argument("--ref-start-maxh-mm", type=float, nargs="+", default=None,
                        help="starting conductor mesh size per --ref-order entry")
    parser.add_argument("--ref-air-mm", type=float, nargs="+", default=None,
                        help="far-air mesh size per --ref-order entry")
    parser.add_argument("--ref-air-divisor", type=float, nargs="+", default=None,
                        help="far air no finer than exterior radius / divisor, per --ref-order entry")
    parser.add_argument("--ref-theta", type=float, default=0.5, help="Doerfler fraction of the indicator")
    parser.add_argument("--ref-max-ndof", type=float, default=6e6)
    parser.add_argument("--ref-max-levels", type=int, default=16)
    parser.add_argument("--ref-hours", type=float, default=2.0, help="wall-clock budget per sequence")
    parser.add_argument("--ref-sequence", type=float, nargs=2, metavar=("ORDER", "FACTOR"),
                        help="run only this adaptive sequence and write it to --output")
    parser.add_argument("--ref-sequence-json", type=Path, nargs="+",
                        help="use these finished sequences instead of running them")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    _, provenance = _branch_ring_solver()      # before any module puts the checkout's src first
    from fin_return_principle import PHYSICAL_TARGETS, physical_verdict

    n = len(args.ref_order)
    starts = args.ref_start_maxh_mm or [0.05]*n
    airs = args.ref_air_mm or [20.0]*n
    divisors = args.ref_air_divisor or [20.0]*n
    if not (len(starts) == len(airs) == len(divisors) == n):
        raise ValueError("one starting mesh size, air size and divisor per reference order")
    settings = {o: (s, a, d) for o, s, a, d in zip(args.ref_order, starts, airs, divisors)}
    tol = {k: 0.1*PHYSICAL_TARGETS[k] for k in ("R", "partial_L", "root_heat", "tip_heat")}

    def run(order, factor):
        s, a, d = settings[order]
        return adaptive_sequence(args.radius_mm, args.frequency_hz, args.sigma, order, factor, s, a, d,
                                 args.ref_theta, tol, args.ref_max_ndof, args.ref_max_levels,
                                 3600*args.ref_hours)

    if args.ref_sequence:
        order, factor = int(args.ref_sequence[0]), float(args.ref_sequence[1])
        seq = run(order, factor)
        head, dirty = git_state()
        args.output.write_text(json.dumps(dict(
            schema=SEQUENCE_SCHEMA, git_head=head, git_dirty=dirty,
            script_sha256=script_sha256(),
            reference_runtime=provenance, radius_mm=args.radius_mm, frequency_hz=args.frequency_hz,
            sigma_S_per_m=args.sigma, sequence=seq), indent=1) + "\n", encoding="utf-8")
        print(json.dumps({"phase": "sequence", "order": order, "exterior_factor": factor,
                          "converged": seq["converged"], "stop": seq["stop"],
                          "certification": certify_sequence(seq, expected_settings(args, settings, order, tol))}),
              flush=True)
        return
    if args.ref_sequence_json:
        sequences = load_sequences(args.ref_sequence_json, args, provenance, settings, tol)
    else:
        sequences = [run(o, f) for o in args.ref_order for f in args.ref_exterior_factor]
    res = resolve_reference(sequences, args.ref_order, args.ref_exterior_factor,
                            lambda order: expected_settings(args, settings, order, tol), PHYSICAL_TARGETS)
    certification, converged, swept, missing = res["certification"], res["converged"], res["swept"], res["missing"]
    ref_q, checks, spread = res["reference"], res["checks"], res["spread"]
    consistent, resolved = res["consistent"], res["resolved"]
    reference = [dict(level, order=s["order"]) for s in sequences for level in s["levels"]]

    peec = []
    for h in args.peec_boundary_h_mm:
        points, tris, region = section_mesh(args.radius_mm, h, args.peec_interior_h_mm, args.grading)
        row = dict(boundary_h_mm=h, interior_h_mm=args.peec_interior_h_mm, grading=args.grading,
                   **solve_peec(points, tris, region, args.frequency_hz, args.sigma, args.threads))
        peec.append(row)
        print(json.dumps({"phase": "peec", **row}), flush=True)
    candidate = quantities(peec[-1])
    if missing:
        verdict = dict(status="NOT_ADOPTED: reference sequences without a checkpoint", within_raw_limits=None,
                       differences=None, sequences_without_checkpoint=missing)
    else:
        verdict = physical_verdict(candidate, ref_q, reference_resolved=resolved)
    head, dirty = git_state()
    report = dict(
        schema="radia.ih.axisym-volume-peec-fin.v1", generated_at_utc=datetime.now(timezone.utc).isoformat(),
        git_head=head, git_dirty=dirty, script_sha256=script_sha256(),
        reference_runtime=provenance, peec_runtime=_peec_runtime(), platform_class=platform.system(),
        radius_mm=args.radius_mm, frequency_hz=args.frequency_hz, sigma_S_per_m=args.sigma,
        skin_depth_mm=1e3/math.sqrt(math.pi*args.frequency_hz*MU0*args.sigma),
        cuts_mm=dict(root=[ROOT_X_MM, TIP_X_MM], tip=[TIP_X_MM, None]),
        normalization="per 1 A peak terminal current; loop-voltage port; exp(+j omega t)",
        reference_settings=dict(theta=args.ref_theta, tolerance=tol, max_ndof=args.ref_max_ndof,
                                max_levels=args.ref_max_levels, budget_s=3600*args.ref_hours,
                                per_order={str(o): dict(start_maxh_mm=v[0], air_maxh_mm=v[1], air_divisor=v[2])
                                           for o, v in settings.items()}),
        reference_levels=reference, reference_sequences=[{k: v for k, v in s.items() if k != "levels"}
                                                         for s in sequences],
        reference_checks=checks, reference_spread_max=spread,
        reference_resolved_rule=("adaptive (user decision 2026-11-16): two orders and two exteriors, each an "
                                 "adaptive sequence that refines the elements of the whole mesh (conductor and "
                                 "air) holding the Doerfler fraction of the ZZ grad A_phi indicator, "
                                 "until two successive checkpoints (levels at least twice apart in DOFs) agree "
                                 "within one tenth of every limit, inside a DOF, "
                                 "level and wall-clock budget; then every quantity within one tenth of its limit for "
                                 "(b) each exterior against the widest at the highest order, (c) each order's final "
                                 "level against the highest order's at the widest exterior; the highest order's "
                                 "terminal Re(V/I) within 1e-3 of the dissipation resistance; a sequence out of "
                                 "budget leaves the reference unresolved"),
        reference_certification=certification, reference_sequences_without_checkpoint=missing,
        reference_swept=swept, reference_loss_consistent=consistent, reference_converged=converged,
        reference_resolved=resolved, peec_levels=peec,
        peec_change_last_two={k: abs(quantities(peec[-1])[k] - quantities(peec[-2])[k])/abs(quantities(peec[-1])[k])
                              for k in candidate} if len(peec) > 1 else None,
        verdict=verdict,
        not_claimed=("an axisymmetric drive only; the three-dimensional Cartesian partial elements of "
                     "radia.volume_peec and any non-axisymmetric port are not exercised; root/tip mapping "
                     "onto the 1.6/3.5 mm cuts is a proposal"))
    args.output.write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({"phase": "verdict", **verdict}), flush=True)


if __name__ == "__main__":
    import ngsolve as ng

    with ng.TaskManager():
        main()
