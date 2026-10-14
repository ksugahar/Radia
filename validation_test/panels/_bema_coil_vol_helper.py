"""Test helper: pre-mesh a coil STEP into a Netgen .vol with the
'source' / 'sink' boundary labels BEM-A expects on its --coil-vol input.

After the BEM-A redesign (radia 4.33.0+) the production calc_inductance.py
no longer meshes coil STEP files on the fly; it expects a pre-meshed
surface .vol from Cubit / Netgen.  Existing test fixtures only ship the
STEP, so this helper does the one-shot OCC -> Netgen mesh conversion
inside the test runner (cached per pytest session).

The labelling logic mirrors the retired _build_bema_coil_mesh() in
calc_inductance.py: face.name == 'source'/'sink' if present, else the
two smallest PLANE faces are auto-detected by |y-centroid|.
"""
from __future__ import annotations

import os

from ngsolve import TaskManager


def step_to_coil_vol(step_path, vol_path, maxh=0.012,
                       source_name="source", sink_name="sink",
                       volume_mesh=False):
    """Mesh ``step_path`` into ``vol_path`` with source/sink boundary labels.

    Parameters
    ----------
    step_path : str
        Input CAD STEP file.  Must exist.
    vol_path : str
        Output .vol path.  Returned for convenience; written next to the
        STEP if relative.  Will be overwritten if it already exists.
    maxh : float
        Surface mesh maxh in meters (default 0.012 m -- the converged
        value used in the retired BEM-A panel default).
    source_name, sink_name : str
        Boundary labels written to the .vol.
    volume_mesh : bool
        False (default): surface-only mesh (Glue of the labelled faces),
        the classic helper output.  True: mesh the SOLID with volume
        tetrahedra -- the common Cubit-export shape.  Used by the
        regression test for the volume-coil-vol failure class (a raw
        ``Mesh(coil.vol)`` handed to HDivSurface picks up extra null
        modes and the EFIE saddle LU goes singular).

    Returns
    -------
    str
        Absolute path to the written .vol.
    """
    from netgen.occ import OCCGeometry, Glue
    from netgen.meshing import MeshingParameters

    if not os.path.isfile(step_path):
        raise FileNotFoundError(f"STEP not found: {step_path}")

    geo = OCCGeometry(step_path)
    shape = geo.shape
    faces = list(shape.faces)
    src_idx = snk_idx = None
    for i, f in enumerate(faces):
        nm = getattr(f, "name", None)
        if nm == source_name and src_idx is None:
            src_idx = i
        elif nm == sink_name and snk_idx is None:
            snk_idx = i

    src_ids, snk_ids = (None if src_idx is None else [src_idx],
                        None if snk_idx is None else [snk_idx])
    if src_ids is None or snk_ids is None:
        # Auto-detect the two current-injection caps.
        #
        # A cap may be SPLIT into several sub-faces: the ih_fem_kelvin demo
        # coil's caps are each split into a z<0 / z>0 pair.  The old
        # "two smallest faces by mass" rule then labelled the two HALVES OF
        # THE SAME cap as source/sink, so the EFIE drove current through the
        # 2.6 mm thickness instead of around the 66 mm ring and returned
        # L = 0.28 nH instead of ~90 nH (PEEC on the same coil: 104.9 nH) --
        # a 373x error that passed silently (2026-07-15).
        #
        # Sub-faces of ONE cap share their VERTEX centroid exactly, while the
        # two caps have distinct ones -- and unlike face centres / pairwise
        # distances (all ~2.6 mm on that coil, hence useless) this separates
        # them cleanly.  So: group every face by vertex centroid, then take
        # the two groups of least total area as the caps and label ALL
        # sub-faces of each.
        import numpy as np

        def _vertex_centroid(f):
            P = np.array([[v.p[0], v.p[1], v.p[2]] for v in f.vertices],
                         dtype=float)
            return P.mean(axis=0)

        try:
            groups = {}
            for i, f in enumerate(faces):
                key = tuple(np.round(_vertex_centroid(f), 9))
                groups.setdefault(key, []).append(i)
            ranked = sorted(
                groups.items(),
                key=lambda kv: sum(faces[i].mass for i in kv[1]))
        except Exception as exc:
            raise ValueError(
                f"could not auto-detect source/sink faces in {step_path}: "
                f"{exc}.  Add face.name = {source_name!r}/{sink_name!r} "
                f"in build123d before export.") from exc
        if len(ranked) < 2:
            raise ValueError(
                f"could not auto-detect source/sink faces in {step_path}: "
                f"found {len(ranked)} distinct face group(s), need 2 caps.  "
                f"Add face.name = {source_name!r}/{sink_name!r} in build123d "
                f"before export.")
        (ca, ga), (cb, gb) = ranked[0], ranked[1]
        # Keep the historical tie-break: the cap nearer the y=0 plane is
        # the source.
        if abs(ca[1]) <= abs(cb[1]):
            src_ids, snk_ids = ga, gb
        else:
            src_ids, snk_ids = gb, ga

    src_set, snk_set = set(src_ids), set(snk_ids)
    for i, f in enumerate(faces):
        if i in src_set:
            f.name = source_name
        elif i in snk_set:
            f.name = sink_name
        else:
            f.name = "body"

    with TaskManager():
        # volume_mesh: mesh the original solid (face names set above are
        # carried on the TopoDS faces, so boundary labels survive) -> the
        # .vol contains volume tets like a Cubit export.  Default: Glue the
        # faces for the classic surface-only mesh.
        mesh_shape = shape if volume_mesh else Glue(faces)
        ngmesh = OCCGeometry(mesh_shape).GenerateMesh(
            mp=MeshingParameters(maxh=maxh, curvaturesafety=1.0,
                                  segmentsperedge=1))
        vol_path = os.path.abspath(vol_path)
        ngmesh.Save(vol_path)

    # Self-validate: the cap grouping above is a heuristic, so verify the
    # labels actually give a current path AROUND the coil rather than a
    # shortcut across it.  Fail loud here instead of writing a .vol that
    # silently poisons every BEM-A solve consuming it (the 0.28 nH vs
    # 90 nH incident, 2026-07-15).
    from ngsolve import Mesh
    from calc_inductance import _arrays_from_bema_coil_mesh
    from radia.bem.coil_inductance_ngsolve import (
        check_source_sink_current_path)

    _m = Mesh(vol_path)
    _v, _t, _sm, _km = _arrays_from_bema_coil_mesh(
        _m, source_name=source_name, sink_name=sink_name)
    check_source_sink_current_path(
        _v, _t, _sm, _km, source_label=source_name, sink_label=sink_name)
    return vol_path


def coil_vol_for(step_path, *, maxh=0.012, cache_dir=None):
    """Cached wrapper: mesh ``step_path`` -> ``<step_path>.bema.vol``.

    Re-uses an existing .vol if it is newer than the STEP (avoids re-meshing
    on every test run).  If ``cache_dir`` is given, the .vol lives there
    rather than alongside the STEP (useful when the STEP fixture lives in
    a read-only / shared location).
    """
    if cache_dir:
        os.makedirs(cache_dir, exist_ok=True)
        base = os.path.basename(step_path).rsplit(".", 1)[0]
        vol_path = os.path.join(cache_dir, f"{base}_bema.vol")
    else:
        vol_path = step_path + ".bema.vol"
    if (os.path.isfile(vol_path)
            and os.path.getmtime(vol_path) >= os.path.getmtime(step_path)):
        return os.path.abspath(vol_path)
    return step_to_coil_vol(step_path, vol_path, maxh=maxh)
