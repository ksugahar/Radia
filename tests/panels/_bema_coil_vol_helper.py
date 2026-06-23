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
                       source_name="source", sink_name="sink"):
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

    if src_idx is None or snk_idx is None:
        try:
            fs_sorted = sorted(enumerate(faces), key=lambda x: x[1].mass)
        except Exception as exc:
            raise ValueError(
                f"could not auto-detect source/sink faces in {step_path}: "
                f"{exc}.  Add face.name = {source_name!r}/{sink_name!r} "
                f"in build123d before export.") from exc
        ca, fa = fs_sorted[0]
        cb, fb = fs_sorted[1]
        if abs(fa.center[1]) < abs(fb.center[1]):
            src_idx, snk_idx = ca, cb
        else:
            src_idx, snk_idx = cb, ca

    for i, f in enumerate(faces):
        if i == src_idx:
            f.name = source_name
        elif i == snk_idx:
            f.name = sink_name
        else:
            f.name = "body"

    with TaskManager():
        ngmesh = OCCGeometry(Glue(faces)).GenerateMesh(
            mp=MeshingParameters(maxh=maxh, curvaturesafety=1.0,
                                  segmentsperedge=1))
        vol_path = os.path.abspath(vol_path)
        ngmesh.Save(vol_path)
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
