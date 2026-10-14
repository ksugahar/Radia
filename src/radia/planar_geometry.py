"""Shared 2D planar mesh helpers.

This module contains only method-neutral NGSolve mesh plumbing used by the
planar HDiv-VIM, anisotropic, hysteresis, and eddy-coupling layers.
"""
from __future__ import annotations

import numpy as np
import ngsolve as ng

from radia.planar_materials import (
    per_region_chi as _region_chi_for,
    per_region_law as _per_region_law_mats,
    region_ids as _region_ids,
)


def extract_geometry(mesh):
    """NGSolve 2D mesh -> ``(verts, offsets, centroids, areas)``."""
    if mesh.dim != 2:
        raise ValueError("planar_geometry: mesh.dim must be 2 (got %d)" % mesh.dim)
    pts = np.array([list(mesh[v].point)[:2] for v in mesh.vertices])
    verts, offsets, centroids, areas = [], [0], [], []
    for el in mesh.Elements(ng.VOL):
        V = pts[[v.nr for v in el.vertices]]
        x, y = V[:, 0], V[:, 1]
        a2 = np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)
        cx = np.sum((x + np.roll(x, -1)) * (x * np.roll(y, -1) - np.roll(x, -1) * y)) / (3 * a2)
        cy = np.sum((y + np.roll(y, -1)) * (x * np.roll(y, -1) - np.roll(x, -1) * y)) / (3 * a2)
        verts.append(V)
        offsets.append(offsets[-1] + len(V))
        centroids.append((cx, cy))
        areas.append(0.5 * abs(a2))
    return (np.vstack(verts).astype(np.float64), np.asarray(offsets, dtype=np.int32),
            np.asarray(centroids, float), np.asarray(areas, float))


def eval_Hext(H_ext, centroids, mesh):
    """Applied field at element centroids -> ``(nElem, 2)`` array."""
    arr = np.asarray(H_ext, dtype=object) if not isinstance(H_ext, ng.CoefficientFunction) else None
    if arr is not None:
        a = np.asarray(H_ext, float)
        if a.shape == (2,):
            return np.tile(a, (len(centroids), 1))
        if a.shape == (len(centroids), 2):
            return a
        raise ValueError("planar_geometry: H_ext array must be (2,) uniform or (nElem, 2)")
    out = np.zeros((len(centroids), 2))
    for i, (px, py) in enumerate(centroids):
        val = H_ext(mesh(px, py))
        out[i] = (val[0], val[1])
    return out


def M_avg(M, areas):
    """Area-weighted average of per-element planar magnetization."""
    w = areas / areas.sum()
    return np.array([float(w @ M[:, 0]), float(w @ M[:, 1])])


def element_materials(mesh):
    """Per-element material name in the same order as ``extract_geometry``."""
    return [el.mat for el in mesh.Elements(ng.VOL)]


def per_region_chi(mesh, mu_r_dict):
    """Per-element chi from a ``{region_name: mu_r}`` dict."""
    return _region_chi_for(element_materials(mesh), mu_r_dict)


def sub_geometry(verts, offsets, ids):
    """Extract the geometry arrays for an element subset."""
    blocks = [verts[offsets[i]:offsets[i + 1]] for i in ids]
    vs = np.vstack(blocks) if blocks else np.zeros((0, 2))
    offs = np.zeros(len(ids) + 1, np.int32)
    for j, i in enumerate(ids):
        offs[j + 1] = offs[j] + (offsets[i + 1] - offsets[i])
    return np.ascontiguousarray(vs, np.float64), offs


def pm_hard_M(mesh, pm, mats, nElem):
    """Full-length fixed magnetization array for permanent-magnet regions."""
    M = np.zeros((nElem, 2))
    rid = _region_ids(mats)
    for name, mv in pm.items():
        if name not in rid:
            raise ValueError("planar_geometry: pm region %r not in mesh; regions: %s"
                             % (name, sorted(rid)))
        mv = np.asarray(mv, float)
        if mv.shape == (2,) or mv.shape == (len(rid[name]), 2):
            M[rid[name]] = mv
        else:
            raise ValueError("planar_geometry: pm[%r] must be [Mx,My] or (nRegionElem,2)" % name)
    return M


def per_region_law(mesh, bh_dict):
    """Per-element ``(M_of_h, chi_sec, chi0)`` from a ``{region: B-H table}`` dict."""
    return _per_region_law_mats(element_materials(mesh), bh_dict)


# Private compatibility names for modules that share these helpers internally.
_extract_geometry = extract_geometry
_eval_Hext = eval_Hext
_M_avg = M_avg
_element_materials = element_materials
_per_region_chi = per_region_chi
_sub_geometry = sub_geometry
_pm_hard_M = pm_hard_M
_per_region_law = per_region_law
