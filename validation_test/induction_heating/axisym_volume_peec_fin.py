"""Axisymmetric volume PEEC of the revolved beak fin against the ring FEM.

The 5% fin adjudication (user decision 2026-11-01, limits in
``fin_return_principle.PHYSICAL_TARGETS``) asks whether a volume PEEC, which
resolves the conductor interior instead of imposing a surface impedance,
reproduces an independent interior-resolved reference.  The reference is the
A_phi ring FEM with an exact Dirichlet-to-Neumann exterior
(``radia.eddy_axisym_ring``), whose drive is the loop voltage.

Under that axisymmetric drive the three-dimensional volume PEEC of the
revolved solid has an azimuthally uniform current, so each meridian cell is
one ring loop and the circuit is exact with ring partial elements:

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
import os
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

# Dunavant rules (barycentric, weights sum to 1).
_R1 = (np.array([[1/3, 1/3, 1/3]]), np.array([1.0]))
_R3 = (np.array([[2/3, 1/6, 1/6], [1/6, 2/3, 1/6], [1/6, 1/6, 2/3]]), np.full(3, 1/3))
_a, _b = 0.445948490915965, 0.091576213509771
_R6 = (np.array([[1-2*_a, _a, _a], [_a, 1-2*_a, _a], [_a, _a, 1-2*_a],
                 [1-2*_b, _b, _b], [_b, 1-2*_b, _b], [_b, _b, 1-2*_b]]),
       np.array([0.223381589678011]*3 + [0.109951743655322]*3))
_c, _d = 0.470142064105115, 0.101286507323456
_R7 = (np.array([[1/3, 1/3, 1/3], [1-2*_c, _c, _c], [_c, 1-2*_c, _c], [_c, _c, 1-2*_c],
                 [1-2*_d, _d, _d], [_d, 1-2*_d, _d], [_d, _d, 1-2*_d]]),
       np.array([0.225] + [0.132394152788506]*3 + [0.125939180544827]*3))


def filament_mutual(r1, z1, r2, z2):
    """Mutual inductance of two coaxial circular filaments (Maxwell)."""
    from scipy.special import ellipe, ellipkm1
    s = (r1+r2)**2 + (z1-z2)**2
    m1 = ((r1-r2)**2 + (z1-z2)**2)/s
    m = 1.0 - m1
    k = np.sqrt(m)
    return MU0*np.sqrt(r1*r2)*((2/k - k)*ellipkm1(m1) - (2/k)*ellipe(m))


def log_potential(x, tri):
    """int_T ln|x-y| dA_y in closed form; x (..., q, 2) points, tri (..., 3, 2) triangles."""
    x = np.asarray(x, dtype=float)
    tri = np.asarray(tri, dtype=float)
    out = np.zeros(x.shape[:-1])
    for i in range(3):
        p, q = tri[..., i, None, :], tri[..., (i+1) % 3, None, :]
        e = q - p
        length = np.linalg.norm(e, axis=-1)
        u = e/length[..., None]
        cross = (p[..., 0]-x[..., 0])*(q[..., 1]-x[..., 1]) - (p[..., 1]-x[..., 1])*(q[..., 0]-x[..., 0])
        h = np.abs(cross)/length
        t1 = np.sum((p - x)*u, axis=-1)
        t2 = t1 + length
        safe = np.where(h > 0, h, 1.0)

        def F(t):
            return (h/4)*(t*np.log(h*h + t*t) - 3*t + 2*h*np.arctan(t/safe))
        out += np.where(h > 0, np.sign(cross)*(F(t2) - F(t1)), 0.0)
    return out


def section_mesh(radius_mm, boundary_h_mm, interior_h_mm, grading):
    """Meridian triangles of the tracked beak section, split at the region cuts."""
    from build123d import Plane, export_step
    from make_beak_fin_step import _beak_face
    from netgen.occ import Box, Glue, OCCGeometry, Pnt
    import ngsolve as ng
    face = _beak_face(Plane.XY, radial_offset=radius_mm)
    step = HERE / f"_peec_section_{os.getpid()}.step"
    export_step(face, str(step))
    try:
        section = OCCGeometry(str(step), dim=2).shape.Scale(Pnt(0, 0, 0), 1e-3)
    finally:
        step.unlink(missing_ok=True)
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


def assemble(points, tris, near_factor=3.0, moderate_factor=10.0, chunk=2048):
    """Galerkin P0 ring partial inductance matrix (symmetric)."""
    from scipy.spatial import cKDTree
    V = points[tris]                                  # (n,3,2)
    n = len(tris)
    area = 0.5*np.abs((V[:, 1, 0]-V[:, 0, 0])*(V[:, 2, 1]-V[:, 0, 1])
                      - (V[:, 2, 0]-V[:, 0, 0])*(V[:, 1, 1]-V[:, 0, 1]))
    if np.any(area <= 0):
        raise RuntimeError("degenerate meridian triangle")
    c = V.mean(axis=1)
    size = np.sqrt(2*area)
    M = np.empty((n, n))
    for s in range(0, n, chunk):          # centroid rule everywhere first
        e = min(n, s+chunk)
        M[s:e] = filament_mutual(c[s:e, None, 0], c[s:e, None, 1], c[None, :, 0], c[None, :, 1])
    np.fill_diagonal(M, 0.0)
    tree = cKDTree(c)

    def pts(rule, idx):
        bary, w = rule
        return np.einsum("qj,njd->nqd", bary, V[idx]), w

    def pair_list(factor):
        """Pairs a < b with centroid distance below factor * (size_a + size_b)."""
        found = tree.query_ball_point(c, r=factor*4*size)     # graded mesh: neighbour size <= 3x
        a = np.repeat(np.arange(n), [len(f) for f in found])
        b = np.concatenate([np.asarray(f, dtype=np.int64) for f in found])
        keep = (a < b) & (np.linalg.norm(c[a] - c[b], axis=1) < factor*(size[a] + size[b]))
        return np.stack([a[keep], b[keep]], axis=1)
    moderate = pair_list(moderate_factor)
    if len(moderate):
        a, b = moderate[:, 0], moderate[:, 1]
        xa, wa = pts(_R3, a)
        xb, wb = pts(_R3, b)
        val = filament_mutual(xa[:, :, None, 0], xa[:, :, None, 1], xb[:, None, :, 0], xb[:, None, :, 1])
        val = np.einsum("pij,i,j->p", val, wa, wb)
        M[a, b] = val
        M[b, a] = val
    near = pair_list(near_factor)
    pairs = np.concatenate([np.stack([np.arange(n), np.arange(n)], axis=1), near])
    xo, wo = _R7
    xi, wi = _R6

    def directed(t, s_, chunk=20000):
        """Galerkin entry with target triangles t and source triangles s_, singular part in closed form."""
        out = np.empty(len(t))
        for k in range(0, len(t), chunk):
            tt, ss = t[k:k+chunk], s_[k:k+chunk]
            X = np.einsum("qj,pjd->pqd", xo, V[tt])
            Y = np.einsum("qj,pjd->pqd", xi, V[ss])
            r = X[..., 0]
            logpart = MU0*r*((np.log(8*r) - 2.0)*area[ss, None] - log_potential(X, V[ss]))
            dist = np.linalg.norm(X[:, :, None, :] - Y[:, None, :, :], axis=-1)
            subtr = MU0*r[..., None]*(np.log(8*r[..., None]) - 2.0 - np.log(dist))
            rem = filament_mutual(X[:, :, None, 0], X[:, :, None, 1], Y[:, None, :, 0], Y[:, None, :, 1]) - subtr
            total = area[tt]*(logpart @ wo) + area[tt]*area[ss]*np.einsum("pij,i,j->p", rem, wo, wi)
            out[k:k+chunk] = total/(area[tt]*area[ss])
        return out
    a, b = pairs[:, 0], pairs[:, 1]
    val = 0.5*(directed(a, b) + directed(b, a))
    M[a, b] = val
    M[b, a] = val
    return M, area, c


def solve_peec(points, tris, region, frequency, sigma, threads):
    t0 = time.perf_counter()
    M, area, c = assemble(points, tris)
    t1 = time.perf_counter()
    omega = 2*math.pi*frequency
    R = 2*math.pi*c[:, 0]/(sigma*area)
    Z = 1j*omega*M
    Z[np.diag_indices_from(Z)] += R
    rhs = np.ones(len(R), dtype=complex)
    I = np.linalg.solve(Z, rhs)
    t2 = time.perf_counter()
    residual = float(np.linalg.norm(Z @ I - rhs)/np.linalg.norm(rhs))
    if not residual <= 1e-6:
        raise RuntimeError(f"true relative residual {residual:.3e} above 1e-6")
    current = complex(I.sum())
    impedance = 1.0/current
    loss = 0.5*R*np.abs(I)**2/abs(current)**2         # per 1 A peak
    regional = {name: float(loss[region == name].sum()) for name in ("body", "root", "tip")}
    total = float(loss.sum())
    return dict(cells=int(len(R)), R_ohm=float(impedance.real), L_H=float(impedance.imag/omega),
                power_balance=abs(2*total - impedance.real)/impedance.real,
                root_heat_W_per_A2=regional["root"], tip_heat_W_per_A2=regional["tip"],
                beak_loss_fraction=(regional["root"]+regional["tip"])/total,
                tip_loss_fraction=regional["tip"]/total,
                true_relative_residual=residual, assembly_s=t1-t0, solve_s=t2-t1,
                min_cell_mm=float(np.sqrt(2*area.min())*1e3), threads=threads)


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
    return sys.modules[name].solve_axisym_ring, dict(
        natives=str(Path(radia.__file__).parent), version=radia.__version__,
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


def adaptive_sequence(radius_mm, frequency, sigma, order, factor, start_maxh, air_mm, air_divisor,
                      theta, tol, max_ndof, max_levels, deadline):
    """Refine one (order, exterior) reference adaptively until two successive levels agree within tol."""
    import ngsolve as ng
    from aphi_beak_section_sweep import build_section_ring
    mesh, ext_radius, ext_z = build_section_ring(radius_mm, maxh_conductor_mm=start_maxh, maxh_air_mm=air_mm,
                                                exterior_factor=factor, air_divisor=air_divisor)
    # A Doerfler step can touch a few percent of the elements, and two such neighbours agree trivially.
    # Successive levels are therefore compared as checkpoints at least `growth` times apart in DOFs.
    growth_min = 2.0
    rows, checkpoints, stop = [], [], None
    while stop is None:
        row, res = solve_level(mesh, ext_radius, ext_z, radius_mm, frequency, sigma, order)
        eta2 = zz_indicator(res, order)
        del res
        row.update(level=len(rows), exterior_factor=factor, indicator=float(math.sqrt(eta2.sum())))
        is_check = not checkpoints or row["ndof"] >= growth_min*rows[checkpoints[-1]]["ndof"]
        row["checkpoint"] = is_check
        if is_check and checkpoints:
            prev = rows[checkpoints[-1]]
            row["change"] = {k: abs(quantities(row)[k] - quantities(prev)[k])/abs(quantities(row)[k])
                             for k in tol}
        rows.append(row)
        print(json.dumps({"phase": "reference", **row}), flush=True)
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
        ne = mesh.ne
        for el in mesh.Elements(ng.VOL):
            mesh.SetRefinementFlag(el, bool(flags[el.nr]))
        mesh.Refine()
        growth = mesh.ne/ne
        if row["ndof"]*growth > max_ndof:
            stop = f"budget: next level ~{row['ndof']*growth:.3g} DOFs above {max_ndof:.3g}"
        elif time.time() + 1.5*growth*row["seconds"] > deadline:
            stop = "budget: deadline"
    return dict(order=order, exterior_factor=factor, start_maxh_mm=start_maxh, air_maxh_mm=air_mm,
                air_divisor=air_divisor, theta=theta, checkpoint_dof_growth=growth_min,
                tolerance=tol, max_ndof=max_ndof,
                max_levels=max_levels, converged=stop == "converged", stop=stop, levels=rows)


def quantities(row):
    return dict(R=row["R_ohm"], partial_L=row["L_H"], root_heat=row["root_heat_W_per_A2"],
                tip_heat=row["tip_heat_W_per_A2"])


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
                                 time.time() + 3600*args.ref_hours)

    if args.ref_sequence:
        order, factor = int(args.ref_sequence[0]), float(args.ref_sequence[1])
        seq = run(order, factor)
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
        args.output.write_text(json.dumps(dict(
            schema="radia.ih.axisym-ring-adaptive-sequence.v1", git_head=head,
            script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            reference_runtime=provenance, radius_mm=args.radius_mm, frequency_hz=args.frequency_hz,
            sigma_S_per_m=args.sigma, sequence=seq), indent=1) + "\n", encoding="utf-8")
        print(json.dumps({"phase": "sequence", "order": order, "exterior_factor": factor,
                          "converged": seq["converged"], "stop": seq["stop"]}), flush=True)
        return
    if args.ref_sequence_json:
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
        sequences = []
        for p in args.ref_sequence_json:
            doc = json.loads(p.read_text(encoding="utf-8"))
            if doc["git_head"] != head or doc["radius_mm"] != args.radius_mm \
                    or doc["frequency_hz"] != args.frequency_hz or doc["sigma_S_per_m"] != args.sigma:
                raise ValueError(f"{p} was run at {doc['git_head']} or on another operating point")
            sequences.append(doc["sequence"])
    else:
        sequences = [run(o, f) for o in args.ref_order for f in args.ref_exterior_factor]
    by = {(s["order"], s["exterior_factor"]): s for s in sequences}
    want = {(o, f) for o in args.ref_order for f in args.ref_exterior_factor}
    if set(by) != want:
        raise ValueError(f"sequences {sorted(by)} do not cover orders x exteriors {sorted(want)}")
    top = max(args.ref_order)
    wide = max(args.ref_exterior_factor)
    final = {key: [r for r in s["levels"] if r["checkpoint"]][-1] for key, s in by.items()}
    best = final[(top, wide)]
    ref_q = quantities(best)
    rel = lambda a, b: {k: abs(quantities(a)[k] - quantities(b)[k])/abs(quantities(b)[k]) for k in ref_q}
    checks = {}
    for (o, f), s in sorted(by.items()):                # mesh: the last two adaptive checkpoints
        lv = [r for r in s["levels"] if r["checkpoint"]]
        checks[f"mesh_o{o}_f{f}"] = rel(lv[-1], lv[-2]) if len(lv) >= 2 else {k: math.inf for k in ref_q}
    for f in args.ref_exterior_factor:                   # exterior, at the highest order
        checks[f"exterior_f{f}"] = rel(final[(top, f)], best)
    for o in args.ref_order:                             # order, at the widest exterior
        checks[f"order_o{o}"] = rel(final[(o, wide)], best)
    spread = {k: max(c[k] for c in checks.values()) for k in ref_q}
    swept = len(set(args.ref_order)) >= 2 and len(set(args.ref_exterior_factor)) >= 2
    consistent = best["terminal_vs_energy_R"] <= 1e-3
    converged = all(s["converged"] for s in sequences)
    resolved = (swept and consistent and converged
                and all(v <= 0.1*PHYSICAL_TARGETS[k] for k, v in spread.items()))
    reference = [dict(level, order=s["order"]) for s in sequences for level in s["levels"]]

    peec = []
    for h in args.peec_boundary_h_mm:
        points, tris, region = section_mesh(args.radius_mm, h, args.peec_interior_h_mm, args.grading)
        row = dict(boundary_h_mm=h, interior_h_mm=args.peec_interior_h_mm, grading=args.grading,
                   **solve_peec(points, tris, region, args.frequency_hz, args.sigma, args.threads))
        peec.append(row)
        print(json.dumps({"phase": "peec", **row}), flush=True)
    candidate = quantities(peec[-1])
    verdict = physical_verdict(candidate, ref_q, reference_resolved=resolved)
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain"], cwd=ROOT, capture_output=True, text=True).stdout.strip())
    report = dict(
        schema="radia.ih.axisym-volume-peec-fin.v1", generated_at_utc=datetime.now(timezone.utc).isoformat(),
        git_head=head, git_dirty=dirty, script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        reference_runtime=provenance, platform_class=platform.system(),
        radius_mm=args.radius_mm, frequency_hz=args.frequency_hz, sigma_S_per_m=args.sigma,
        skin_depth_mm=1e3/math.sqrt(math.pi*args.frequency_hz*MU0*args.sigma),
        cuts_mm=dict(root=[ROOT_X_MM, TIP_X_MM], tip=[TIP_X_MM, None]),
        normalization="per 1 A peak terminal current; loop-voltage port; exp(+j omega t)",
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
