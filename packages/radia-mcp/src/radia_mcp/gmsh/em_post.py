"""Electromagnetics-native post verbs: flux, circulation, force, harmonics.

The generic lane (``post_process``) exposes gmsh's own plugins; this
module adds the four reductions an EM engineer actually asks a field
file for, all built on the same primitive -- sample a GEOMETRIC LOCUS
with ``gmsh.view.probe``, then reduce:

* ``flux_integral``   -- int B.n dA over a rectangle or a disc (Wb)
* ``line_integral``   -- int H.dl around a circle or a polyline (A)
* ``maxwell_force``   -- surface integral of the Maxwell stress over an
  axis-aligned box (N), with the optional torque about a point (N.m)
* ``gap_harmonics``   -- space-harmonic spectrum on an air-gap circle

Measured semantics that shape the implementation (gmsh 4.15.2):

- ``gmsh.view.probe`` returns ``(values, distance)`` and a LIST-BASED
  VECTOR view answers with the NEAREST value at ANY distance, so a
  point outside the data is detected by ``distance > 0``, never by the
  values being empty.  On this build ``values`` is a numpy array --
  gate on ``len(values)``, never on truthiness.
- A locus that pokes outside the mesh would silently integrate over a
  PARTIAL surface, so every verb counts ``n_outside`` and refuses to
  return a number when it is nonzero.

Quadrature choices (all measured, see the test module):

- rectangle: midpoint tensor grid; disc: POLAR midpoint grid with
  ``dA = r dr dtheta`` so the rim is not over-weighted (the naive
  cartesian mask over-counts the boundary ring).
- circulation: the integral runs in PARAMETER space against the
  ANALYTIC tangent.  Integrating in chord space instead costs 1.0e-4
  relative on a 256-chord circle -- polygon-vs-arc geometry, not
  quadrature error.
- harmonics: an ENDPOINT-EXCLUDED circle (theta = 2 pi k / N,
  k = 0..N-1).  Including both endpoints double-counts theta = 0 and
  smears every bin.
"""

from __future__ import annotations

import json
import math
import tempfile
from pathlib import Path
from typing import Any

from ._gmsh_subprocess import run_gmsh_json_subprocess
from .post_process import _check_input, _finite_float, _finite_vector

# Vacuum permeability (CODATA 2018), the default for the Maxwell stress.
MU0 = 1.25663706212e-6

_SURFACE_KINDS = ("circle", "rect")
_PATH_KINDS = ("circle", "polyline")
_COMPONENTS = ("auto", "radial", "tangential", "axial", "magnitude",
               "scalar")

_EM_SCRIPT = r"""
import json
import math
import sys

cfg_path, out_path = sys.argv[1], sys.argv[2]
with open(cfg_path, encoding="utf-8") as f:
    cfg = json.load(f)
result = {"ok": False, "ran": False}


def _view_tag_by_selector(tags, selector):
    # selector: None -> index 0; int -> view index; str -> view name.
    import gmsh
    if selector is None:
        return tags[0]
    if isinstance(selector, int):
        if not 0 <= selector < len(tags):
            raise RuntimeError(
                f"view index {selector} out of range (0..{len(tags)-1})")
        return tags[selector]
    names = {}
    for tag in tags:
        idx = gmsh.view.getIndex(tag)
        names[gmsh.option.getString(f"View[{idx}].Name")] = tag
    if selector not in names:
        raise RuntimeError(
            f"view {selector!r} not found; available: {sorted(names)}")
    return names[selector]


def _view_name(tag):
    import gmsh
    idx = gmsh.view.getIndex(tag)
    return gmsh.option.getString(f"View[{idx}].Name")


def _view_ncomp(tag):
    # Storage-kind dispatch: list-based views encode the component count
    # in the data-type letter (S/V/T); model-based views report it via
    # getHomogeneousModelData.
    import gmsh
    try:
        dtypes, _ne, _data = gmsh.view.getListData(tag)
    except Exception:
        dtypes = ()
    if dtypes:
        return {"S": 1, "V": 3, "T": 9}[str(dtypes[0])[0]]
    _dt, _tags, _data, _time, ncomp = gmsh.view.getHomogeneousModelData(
        tag, 0)
    return int(ncomp)


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a, b):
    return [a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]]


def _unit(v):
    n = math.sqrt(_dot(v, v))
    return [c / n for c in v]


def _frame(n_hat):
    # Deterministic right-handed (e1, e2, n_hat) frame: seed with the
    # cartesian axis n_hat leans on least, then Gram-Schmidt.  Fixed
    # tie-breaking keeps the harmonic phases reproducible run to run.
    k = min(range(3), key=lambda i: abs(n_hat[i]))
    seed = [1.0 if i == k else 0.0 for i in range(3)]
    proj = _dot(seed, n_hat)
    e1 = _unit([seed[i] - proj * n_hat[i] for i in range(3)])
    e2 = _cross(n_hat, e1)
    return e1, e2


try:
    import gmsh
    gmsh.initialize(["-noconfig"])
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(cfg["path"])
        tags = list(gmsh.view.getTags())
        if not tags:
            raise RuntimeError("no post-processing views in the file")
        tag = _view_tag_by_selector(tags, cfg.get("view"))
        ncomp = _view_ncomp(tag)
        step = int(cfg.get("time_step", 0))
        op = cfg["op"]

        def _probe(p):
            # None = the point is outside the data.  A list-based
            # VECTOR view hands back the nearest value at ANY distance
            # (measured gmsh 4.15.2), so the distance IS the test;
            # values is a numpy array, so gate on len().
            values, distance = gmsh.view.probe(
                tag, p[0], p[1], p[2], step=step)
            if len(values) < ncomp or float(distance) > 0.0:
                return None
            return [float(v) for v in values[:ncomp]]

        def _need_vector(verb):
            if ncomp != 3:
                raise RuntimeError(
                    f"{verb} needs a VECTOR view; view "
                    f"{_view_name(tag)!r} has {ncomp} component(s)")

        if op == "flux":
            _need_vector("flux_integral")
            n_grid = int(cfg["n_grid"])
            samples = []          # (point, weight_dA)
            if cfg["kind"] == "rect":
                c = cfg["center"]
                u = cfg["u_vec"]
                v = cfg["v_vec"]
                nrm = _cross(u, v)
                area_total = math.sqrt(_dot(nrm, nrm))
                n_hat = [c_ / area_total for c_ in nrm]
                d_a = area_total / (n_grid * n_grid)
                for i in range(n_grid):
                    s = (i + 0.5) / n_grid - 0.5
                    for j in range(n_grid):
                        t = (j + 0.5) / n_grid - 0.5
                        samples.append(
                            ([c[k] + s * u[k] + t * v[k] for k in range(3)],
                             d_a))
            else:
                c = cfg["center"]
                n_hat = _unit(cfg["normal"])
                e1, e2 = _frame(n_hat)
                radius = float(cfg["radius"])
                # Polar midpoint: dA = r dr dtheta reproduces pi r^2
                # exactly, unlike an unweighted (r, theta) grid.
                d_r = radius / n_grid
                d_th = 2.0 * math.pi / n_grid
                for i in range(n_grid):
                    r_i = (i + 0.5) * d_r
                    for j in range(n_grid):
                        th = (j + 0.5) * d_th
                        cs, sn = math.cos(th), math.sin(th)
                        samples.append(
                            ([c[k] + r_i * (cs * e1[k] + sn * e2[k])
                              for k in range(3)],
                             r_i * d_r * d_th))
            flux = 0.0
            area = 0.0
            n_outside = 0
            for point, d_a in samples:
                field = _probe(point)
                if field is None:
                    n_outside += 1
                    continue
                flux += _dot(field, n_hat) * d_a
                area += d_a
            result.update({
                "ok": n_outside == 0, "ran": True,
                "view": _view_name(tag),
                "flux": flux, "area": area,
                "n_points": len(samples), "n_outside": n_outside,
                "normal": n_hat, "kind": cfg["kind"],
            })
            if n_outside:
                result["error"] = (
                    f"{n_outside} of {len(samples)} surface samples fall "
                    "outside the field data: the patch pokes out of the "
                    "mesh and the flux would cover only part of it")

        elif op == "circulation":
            _need_vector("line_integral")
            samples = []          # (point, tangent * weight)
            if cfg["kind"] == "circle":
                c = cfg["center"]
                n_hat = _unit(cfg["normal"])
                e1, e2 = _frame(n_hat)
                radius = float(cfg["radius"])
                n_int = int(cfg["n"])
                # Trapezoid in the PARAMETER u against the analytic
                # tangent dp/du = r(-e1 sin u + e2 cos u).  Chord-space
                # integration biases the answer by the polygon deficit.
                d_u = 2.0 * math.pi / n_int
                for k in range(n_int + 1):
                    u = k * d_u
                    w = d_u * (0.5 if k in (0, n_int) else 1.0)
                    cs, sn = math.cos(u), math.sin(u)
                    point = [c[i] + radius * (cs * e1[i] + sn * e2[i])
                             for i in range(3)]
                    tang = [radius * (-sn * e1[i] + cs * e2[i])
                            for i in range(3)]
                    samples.append((point, [t * w for t in tang]))
            else:
                pts = [list(p) for p in cfg["points"]]
                if cfg["closed"]:
                    pts = pts + [list(pts[0])]
                n_seg = len(pts) - 1
                per_seg = max(2, int(round(int(cfg["n"]) / n_seg)))
                for s in range(n_seg):
                    a, b = pts[s], pts[s + 1]
                    d_vec = [b[i] - a[i] for i in range(3)]
                    # Segment parameter t in [0, 1]: dl = d_vec dt, so
                    # the integrand is F.d_vec -- LINEAR in t for a
                    # linear field, hence trapezoid-exact.
                    d_t = 1.0 / (per_seg - 1)
                    for k in range(per_seg):
                        t = k * d_t
                        w = d_t * (0.5 if k in (0, per_seg - 1) else 1.0)
                        point = [a[i] + t * d_vec[i] for i in range(3)]
                        samples.append((point, [c_ * w for c_ in d_vec]))
            integral = 0.0
            n_outside = 0
            for point, weighted_tangent in samples:
                field = _probe(point)
                if field is None:
                    n_outside += 1
                    continue
                integral += _dot(field, weighted_tangent)
            result.update({
                "ok": n_outside == 0, "ran": True,
                "view": _view_name(tag),
                "integral": integral,
                "n_samples": len(samples), "n_outside": n_outside,
                "kind": cfg["kind"],
            })
            if n_outside:
                result["error"] = (
                    f"{n_outside} of {len(samples)} path samples fall "
                    "outside the field data: the loop leaves the mesh "
                    "and the circulation would cover only part of it")

        elif op == "maxwell_force":
            _need_vector("maxwell_force")
            c = cfg["center"]
            half = cfg["half"]
            n_grid = int(cfg["n_grid"])
            mu0 = float(cfg["mu0"])
            about = cfg.get("torque_about")
            force = [0.0, 0.0, 0.0]
            torque = [0.0, 0.0, 0.0]
            per_face = []
            n_points = 0
            n_outside = 0
            for axis in range(3):
                b_ax, c_ax = (axis + 1) % 3, (axis + 2) % 3
                for sign in (1.0, -1.0):
                    n_hat = [0.0, 0.0, 0.0]
                    n_hat[axis] = sign
                    d_a = (4.0 * half[b_ax] * half[c_ax]
                           / (n_grid * n_grid))
                    face_force = [0.0, 0.0, 0.0]
                    face_area = 0.0
                    for i in range(n_grid):
                        s = ((i + 0.5) / n_grid - 0.5) * 2.0 * half[b_ax]
                        for j in range(n_grid):
                            t = ((j + 0.5) / n_grid - 0.5) * 2.0 * half[c_ax]
                            point = list(c)
                            point[axis] += sign * half[axis]
                            point[b_ax] += s
                            point[c_ax] += t
                            n_points += 1
                            field = _probe(point)
                            if field is None:
                                n_outside += 1
                                continue
                            # T.n = (B (B.n) - |B|^2 n / 2) / mu0
                            bn = _dot(field, n_hat)
                            b2 = _dot(field, field)
                            tn = [(field[k] * bn - 0.5 * b2 * n_hat[k]) / mu0
                                  for k in range(3)]
                            face_area += d_a
                            for k in range(3):
                                face_force[k] += tn[k] * d_a
                            if about is not None:
                                arm = [point[k] - about[k] for k in range(3)]
                                mom = _cross(arm, tn)
                                for k in range(3):
                                    torque[k] += mom[k] * d_a
                    for k in range(3):
                        force[k] += face_force[k]
                    per_face.append({
                        "face": ("+" if sign > 0 else "-") + "xyz"[axis],
                        "force": face_force, "area": face_area})
            result.update({
                "ok": n_outside == 0, "ran": True,
                "view": _view_name(tag),
                "force_n": force,
                "torque_nm": None if about is None else torque,
                "per_face": per_face,
                "n_points": n_points, "n_outside": n_outside,
                "mu0": mu0,
            })
            if n_outside:
                result["error"] = (
                    f"{n_outside} of {n_points} box-face samples fall "
                    "outside the field data: the integration box pokes "
                    "out of the mesh and the stress integral would cover "
                    "only part of it")

        elif op == "gap_harmonics":
            component = cfg["component"]
            if component == "auto":
                component = "radial" if ncomp == 3 else "scalar"
            if component == "scalar" and ncomp != 1:
                raise RuntimeError(
                    f"component 'scalar' needs a 1-component view; view "
                    f"{_view_name(tag)!r} has {ncomp} component(s)")
            if component != "scalar" and ncomp != 3:
                raise RuntimeError(
                    f"component {component!r} needs a VECTOR view; view "
                    f"{_view_name(tag)!r} has {ncomp} component(s)")
            c = cfg["center"]
            n_hat = _unit(cfg["axis"])
            e1, e2 = _frame(n_hat)
            radius = float(cfg["radius"])
            n_samples = int(cfg["n_samples"])
            # Endpoint-EXCLUDED sampling: theta = 2 pi k / N covers the
            # circle exactly once.  Including theta = 2 pi repeats the
            # first sample and smears every bin.
            values = []
            thetas = []
            n_outside = 0
            for k in range(n_samples):
                th = 2.0 * math.pi * k / n_samples
                cs, sn = math.cos(th), math.sin(th)
                e_r = [cs * e1[i] + sn * e2[i] for i in range(3)]
                point = [c[i] + radius * e_r[i] for i in range(3)]
                field = _probe(point)
                if field is None:
                    n_outside += 1
                    values.append(0.0)
                    thetas.append(th)
                    continue
                if component == "scalar":
                    val = field[0]
                elif component == "radial":
                    val = _dot(field, e_r)
                elif component == "tangential":
                    e_t = [-sn * e1[i] + cs * e2[i] for i in range(3)]
                    val = _dot(field, e_t)
                elif component == "axial":
                    val = _dot(field, n_hat)
                else:
                    val = math.sqrt(_dot(field, field))
                values.append(val)
                thetas.append(th)
            result.update({
                "ok": n_outside == 0, "ran": True,
                "view": _view_name(tag),
                "component": component,
                "samples": values, "thetas": thetas,
                "n_samples": n_samples, "n_outside": n_outside,
            })
            if n_outside:
                result["error"] = (
                    f"{n_outside} of {n_samples} gap samples fall outside "
                    "the field data: the circle leaves the mesh and the "
                    "spectrum would be built from a partial revolution")

        else:
            raise RuntimeError(f"unknown op {op!r}")
    finally:
        gmsh.finalize()
except Exception as exc:
    result["error"] = f"{type(exc).__name__}: {exc}"
with open(out_path, "w", encoding="utf-8") as f:
    json.dump(result, f)
"""


def _run_em(cfg: dict[str, Any], timeout_s: float) -> dict[str, Any]:
    with tempfile.TemporaryDirectory(prefix="radia_mcp_gmsh_em_") as work:
        cfg_path = Path(work) / "em.json"
        cfg_path.write_text(json.dumps(cfg), encoding="utf-8")
        result = run_gmsh_json_subprocess(
            _EM_SCRIPT, [str(cfg_path)],
            timeout_s=timeout_s, prefix="radia_mcp_gmsh_em_")
    result.setdefault("input", cfg.get("path"))
    return result


def _one_kind(spec: Any, kinds: tuple[str, ...], label: str) -> str:
    """Return the single selected variant key, or raise naming the set."""
    if not isinstance(spec, dict):
        raise ValueError(
            f"{label} must be a dict naming exactly one of: "
            f"{', '.join(kinds)}")
    unknown = sorted(set(spec) - set(kinds))
    if unknown:
        raise ValueError(
            f"unknown {label} kind(s) {', '.join(unknown)}; valid: "
            f"{', '.join(kinds)}")
    present = [k for k in kinds if k in spec]
    if len(present) != 1:
        got = ", ".join(present) if present else "nothing"
        raise ValueError(
            f"{label} must name exactly ONE of: {', '.join(kinds)} "
            f"(got {got})")
    if not isinstance(spec[present[0]], dict):
        raise ValueError(f"{label}.{present[0]} must be a dict")
    return present[0]


def _need_keys(block: dict[str, Any], keys: tuple[str, ...],
               label: str) -> None:
    missing = [k for k in keys if k not in block]
    if missing:
        raise ValueError(
            f"{label} needs keys: {', '.join(keys)} "
            f"(missing: {', '.join(missing)})")


def _positive_int(value: Any, name: str, minimum: int = 1) -> int:
    try:
        clean = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be an integer >= {minimum}") from exc
    if clean < minimum:
        raise ValueError(f"{name} must be >= {minimum}")
    return clean


def _nonzero_vector(value: Any, name: str) -> list[float]:
    vec = _finite_vector(value, name, 3)
    if math.sqrt(sum(c * c for c in vec)) == 0.0:
        raise ValueError(f"{name} must be nonzero")
    return vec


def harmonic_series(samples: list[float], *,
                    max_harmonic: int | None = None) -> dict[str, Any]:
    """Reduce ONE revolution of uniform samples to space harmonics.

    ``samples`` are N values on theta = 2 pi k / N (endpoint EXCLUDED).
    With ``s(theta) = a_0 + sum_n [a_n cos(n theta) + b_n sin(n theta)]``
    the rfft convention is ``a_n = 2 Re(F_n) / N``, ``b_n = -2 Im(F_n)
    / N`` and ``a_0 = Re(F_0) / N``.

    Bins run to ``n <= (N - 1) // 2``, i.e. STRICTLY below Nyquist: for
    even N the Nyquist bin carries a single cosine whose coefficient is
    ``Re(F_N/2) / N``, so applying the factor-2 rule there would report
    double the truth.  ``max_harmonic`` only trims the returned list;
    the THD always uses every resolved bin so it cannot be flattered by
    a narrow report window.

    ``thd`` is referenced to n = 1 and is ``None`` when the signal has
    no n = 1 content (MEASURED: a uniform axial field gives
    ``fundamental`` 1.4e-16 against a DC of 2.0, and dividing by it
    reports a THD of 2.86 -- pure roundoff dressed as a result).  A
    machine whose fundamental sits at the pole-pair order should be
    sampled over ONE pole pair, or read from ``harmonics`` directly.
    """
    import numpy as np

    n = len(samples)
    if n < 3:
        raise ValueError("harmonic_series needs at least 3 samples")
    spectrum = np.fft.rfft(np.asarray(samples, dtype=float))
    n_max = (n - 1) // 2
    a_0 = float(spectrum[0].real) / n
    rows: list[dict[str, Any]] = [
        {"n": 0, "a_n": a_0, "b_n": 0.0, "amplitude": abs(a_0),
         "phase_deg": 0.0}]
    for order in range(1, n_max + 1):
        a_n = 2.0 * float(spectrum[order].real) / n
        b_n = -2.0 * float(spectrum[order].imag) / n
        rows.append({
            "n": order, "a_n": a_n, "b_n": b_n,
            "amplitude": math.hypot(a_n, b_n),
            "phase_deg": math.degrees(math.atan2(b_n, a_n))})
    fundamental = rows[1]["amplitude"] if n_max >= 1 else 0.0
    rest = math.sqrt(sum(r["amplitude"] ** 2 for r in rows[2:]))
    # A fundamental at the roundoff floor of the spectrum is NOT a
    # denominator: 1e-12 x the largest coefficient present separates
    # "no n = 1 content" from a small but real fundamental.
    scale = max([abs(a_0)] + [r["amplitude"] for r in rows[1:]])
    thd = (rest / fundamental
           if fundamental > 1e-12 * scale and scale > 0.0 else None)
    if max_harmonic is not None:
        rows = [r for r in rows if r["n"] <= max_harmonic]
    return {"harmonics": rows, "fundamental": fundamental, "thd": thd,
            "n_bins": n_max}


def flux_integral(path: str | Path, surface: dict[str, Any], *,
                  view: str | int | None = None,
                  n_grid: int = 32, time_step: int = 0,
                  timeout_s: float = 300.0) -> dict[str, Any]:
    """Integrate ``B.n`` over a rectangle or a disc -- the flux (Wb).

    ``surface`` names exactly ONE of::

        {"rect":   {"center": [3], "u_vec": [3], "v_vec": [3]}}
        {"circle": {"center": [3], "normal": [3], "radius": r}}

    The rectangle spans ``center +- u_vec/2 +- v_vec/2`` and its normal
    is ``normalize(u x v)``; the disc's is ``normalize(normal)``.  Both
    are sampled with a MIDPOINT rule (``n_grid`` per direction, polar
    with ``dA = r dr dtheta`` for the disc), which is exact for a
    uniform field and second-order otherwise.

    A patch that pokes out of the mesh returns ``ok: False`` with
    ``n_outside`` -- a partial surface silently reporting a "flux" is
    exactly the number nobody can audit later.
    """
    err = _check_input(path)
    if err:
        return err
    try:
        kind = _one_kind(surface, _SURFACE_KINDS, "surface")
        block = surface[kind]
        cfg: dict[str, Any] = {"op": "flux", "path": str(path),
                               "view": view, "kind": kind,
                               "n_grid": _positive_int(n_grid, "n_grid"),
                               "time_step": _positive_int(time_step,
                                                          "time_step", 0)}
        if kind == "rect":
            _need_keys(block, ("center", "u_vec", "v_vec"), "surface.rect")
            u_vec = _nonzero_vector(block["u_vec"], "u_vec")
            v_vec = _nonzero_vector(block["v_vec"], "v_vec")
            cross = [u_vec[1] * v_vec[2] - u_vec[2] * v_vec[1],
                     u_vec[2] * v_vec[0] - u_vec[0] * v_vec[2],
                     u_vec[0] * v_vec[1] - u_vec[1] * v_vec[0]]
            if math.sqrt(sum(c * c for c in cross)) == 0.0:
                raise ValueError(
                    "u_vec and v_vec must not be parallel (u x v = 0 "
                    "leaves the patch normal undefined)")
            cfg.update({"center": _finite_vector(block["center"],
                                                 "center", 3),
                        "u_vec": u_vec, "v_vec": v_vec})
        else:
            _need_keys(block, ("center", "normal", "radius"),
                       "surface.circle")
            radius = _finite_float(block["radius"], "radius")
            if radius <= 0.0:
                raise ValueError("radius must be positive")
            cfg.update({"center": _finite_vector(block["center"],
                                                 "center", 3),
                        "normal": _nonzero_vector(block["normal"],
                                                  "normal"),
                        "radius": radius})
    except (TypeError, ValueError) as exc:
        return {"ok": False, "error": str(exc)}
    return _run_em(cfg, timeout_s)


def line_integral(path: str | Path, path_spec: dict[str, Any], *,
                  view: str | int | None = None, n: int = 512,
                  expected_ni: float | None = None,
                  time_step: int = 0,
                  timeout_s: float = 300.0) -> dict[str, Any]:
    """Circulation ``int H.dl`` around a loop -- Ampere's law (A).

    ``path_spec`` names exactly ONE of::

        {"circle":   {"center": [3], "normal": [3], "radius": r}}
        {"polyline": {"points": [[x, y, z], ...], "closed": bool}}

    The circle is integrated in PARAMETER space against the analytic
    tangent ``dp/du``; the polyline against each segment's direction.
    Both are trapezoid rules, so a field that is linear along the locus
    is integrated exactly.  A circle is positively oriented by the
    right-hand rule about ``normal``.

    ``expected_ni`` (the enclosed ampere-turns) is compared against the
    result and reported as ``rel_err_vs_expected`` -- the usual way to
    check that a cut loop really links the conductor it should.
    """
    err = _check_input(path)
    if err:
        return err
    try:
        kind = _one_kind(path_spec, _PATH_KINDS, "path_spec")
        block = path_spec[kind]
        cfg: dict[str, Any] = {"op": "circulation", "path": str(path),
                               "view": view, "kind": kind,
                               "n": _positive_int(n, "n", 2),
                               "time_step": _positive_int(time_step,
                                                          "time_step", 0)}
        if kind == "circle":
            _need_keys(block, ("center", "normal", "radius"),
                       "path_spec.circle")
            radius = _finite_float(block["radius"], "radius")
            if radius <= 0.0:
                raise ValueError("radius must be positive")
            cfg.update({"center": _finite_vector(block["center"],
                                                 "center", 3),
                        "normal": _nonzero_vector(block["normal"],
                                                  "normal"),
                        "radius": radius})
        else:
            _need_keys(block, ("points",), "path_spec.polyline")
            points = [_finite_vector(p, "points", 3)
                      for p in block["points"]]
            closed = bool(block.get("closed", False))
            if len(points) < 2:
                raise ValueError("path_spec.polyline needs >= 2 points")
            if closed and len(points) < 3:
                raise ValueError(
                    "a closed path_spec.polyline needs >= 3 points")
            cfg.update({"points": points, "closed": closed})
        clean_ni = (None if expected_ni is None
                    else _finite_float(expected_ni, "expected_ni"))
        if clean_ni == 0.0:
            # A relative error against zero has no meaning; an enclosed
            # current of zero is checked against the raw integral.
            raise ValueError(
                "expected_ni must be nonzero (compare a zero-linkage "
                "loop against the raw integral instead)")
    except (TypeError, ValueError) as exc:
        return {"ok": False, "error": str(exc)}
    result = _run_em(cfg, timeout_s)
    result["expected_ni"] = clean_ni
    result["rel_err_vs_expected"] = (
        (result["integral"] - clean_ni) / clean_ni
        if result.get("ok") and clean_ni is not None else None)
    return result


def maxwell_force(path: str | Path, *, box: dict[str, Any] | None = None,
                  view: str | int | None = None, n_grid: int = 24,
                  torque_about: list[float] | None = None,
                  mu0: float = MU0, time_step: int = 0,
                  timeout_s: float = 300.0) -> dict[str, Any]:
    """Maxwell stress force (N) on an axis-aligned box, and its torque.

    ``box = {"center": [3], "half": [3]}``.  Each of the 6 faces is
    sampled with a midpoint grid and the stress

        T.n = (B (B.n) - |B|^2 n / 2) / mu0

    is integrated with the ANALYTIC outward normal; ``torque_about``
    adds ``int (r - r0) x (T.n) dA``.

    ``per_face`` is reported alongside the total precisely because the
    total is often a small difference of large opposing face terms --
    seeing 1e6 N per face cancel to 0 is what distinguishes a real
    cancellation from a zero integrand.
    """
    err = _check_input(path)
    if err:
        return err
    try:
        if box is None:
            raise ValueError(
                "box is required: {'center': [3], 'half': [3]}")
        if not isinstance(box, dict):
            raise ValueError("box must be a dict with keys: center, half")
        _need_keys(box, ("center", "half"), "box")
        half = _finite_vector(box["half"], "half", 3)
        if any(h <= 0.0 for h in half):
            raise ValueError("every box half-extent must be positive")
        clean_mu0 = _finite_float(mu0, "mu0")
        if clean_mu0 <= 0.0:
            raise ValueError("mu0 must be positive")
        cfg: dict[str, Any] = {
            "op": "maxwell_force", "path": str(path), "view": view,
            "center": _finite_vector(box["center"], "center", 3),
            "half": half, "mu0": clean_mu0,
            "n_grid": _positive_int(n_grid, "n_grid"),
            "time_step": _positive_int(time_step, "time_step", 0),
            "torque_about": (None if torque_about is None
                             else _finite_vector(torque_about,
                                                 "torque_about", 3))}
    except (TypeError, ValueError) as exc:
        return {"ok": False, "error": str(exc)}
    return _run_em(cfg, timeout_s)


def gap_harmonics(path: str | Path, center: list[float],
                  axis: list[float], radius: float, *,
                  view: str | int | None = None, n_samples: int = 360,
                  component: str = "auto",
                  max_harmonic: int | None = None,
                  time_step: int = 0,
                  timeout_s: float = 300.0) -> dict[str, Any]:
    """Space-harmonic spectrum on an air-gap circle.

    Samples ``n_samples`` points at theta = 2 pi k / N (endpoint
    EXCLUDED) on the circle of ``radius`` about ``center`` in the plane
    normal to ``axis``, projects the field onto ``component``
    (``radial`` | ``tangential`` | ``axial`` | ``magnitude`` for vector
    views, ``scalar`` for 1-component views, ``auto`` = radial /
    scalar), and reduces to ``a_n``/``b_n``/amplitude/phase plus the
    THD -- the motor air-gap reading of a field file.

    ``n_samples`` bounds the resolvable order at ``(N - 1) // 2``; pick
    at least 4 samples per pole pair of the highest harmonic of
    interest.
    """
    err = _check_input(path)
    if err:
        return err
    try:
        if component not in _COMPONENTS:
            raise ValueError(
                f"component must be one of: {', '.join(_COMPONENTS)} "
                f"(got {component!r})")
        clean_radius = _finite_float(radius, "radius")
        if clean_radius <= 0.0:
            raise ValueError("radius must be positive")
        clean_max = (None if max_harmonic is None
                     else _positive_int(max_harmonic, "max_harmonic"))
        cfg = {"op": "gap_harmonics", "path": str(path), "view": view,
               "center": _finite_vector(center, "center", 3),
               "axis": _nonzero_vector(axis, "axis"),
               "radius": clean_radius,
               "component": component,
               "n_samples": _positive_int(n_samples, "n_samples", 3),
               "time_step": _positive_int(time_step, "time_step", 0)}
    except (TypeError, ValueError) as exc:
        return {"ok": False, "error": str(exc)}
    result = _run_em(cfg, timeout_s)
    if not result.get("ok"):
        return result
    result.update(harmonic_series(result.pop("samples"),
                                  max_harmonic=clean_max))
    result.pop("thetas", None)
    return result
