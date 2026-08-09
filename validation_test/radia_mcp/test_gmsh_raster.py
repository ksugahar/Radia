"""Ray-cast volume rendering and LIC: the two intrinsic-limit closers.

Both used to be "honest gaps": the slice-stack is per-slice (not
per-ray) and the dense-streamline texture leaves pixels empty.  These
are the real algorithms, run in numpy over gmsh-probed fields, and the
tests hold them to properties only the real algorithm has:

* ray casting -- for a uniform medium the minimum transmittance equals
  ``(1 - alpha)^n`` with n the deepest ray's inside-sample count,
  EXACTLY (both numbers come out of the same compositing loop, so this
  pins the front-to-back math, not an approximation of it); and looking
  at a gradient from opposite sides gives opposite dominant colours,
  which only per-ray occlusion produces;
* LIC -- noise convolved along a uniform field must be smeared ALONG
  the field: the finite-difference energy along the flow drops far
  below the transverse energy, and the anisotropy rotates with the
  field.
"""

from __future__ import annotations

import numpy as np
import pytest
from radia_mcp.gmsh.raster import lic, volume_raycast

pytest.importorskip("gmsh", reason="gmsh not installed")
Image = pytest.importorskip("PIL.Image", reason="Pillow not installed")


def _box_msh(path, value_at, *, ncomp=1, n=6, lo=-1.0, hi=1.0):
    """Filled unit box (5-tet cells) with one NodeData view ``f``."""
    nodes, idx, tag = {}, {}, 1
    for i in range(n):
        for j in range(n):
            for k in range(n):
                nodes[tag] = [lo + (hi - lo) * i / (n - 1),
                              lo + (hi - lo) * j / (n - 1),
                              lo + (hi - lo) * k / (n - 1)]
                idx[(i, j, k)] = tag
                tag += 1
    elements, et = {}, 1
    for i in range(n - 1):
        for j in range(n - 1):
            for k in range(n - 1):
                c = [idx[(i, j, k)], idx[(i + 1, j, k)],
                     idx[(i + 1, j + 1, k)], idx[(i, j + 1, k)],
                     idx[(i, j, k + 1)], idx[(i + 1, j, k + 1)],
                     idx[(i + 1, j + 1, k + 1)], idx[(i, j + 1, k + 1)]]
                for a, b, cc, d in ((0, 1, 3, 4), (1, 2, 3, 6), (1, 3, 4, 6),
                                    (1, 4, 5, 6), (3, 4, 6, 7)):
                    elements[et] = [c[a], c[b], c[cc], c[d]]
                    et += 1
    rows = {t: value_at(p) for t, p in nodes.items()}
    out = ["$MeshFormat", "4.1 0 8", "$EndMeshFormat",
           "$Entities", "0 0 0 1",
           f"1 {lo:g} {lo:g} {lo:g} {hi:g} {hi:g} {hi:g} 0 0",
           "$EndEntities",
           "$Nodes", f"1 {len(nodes)} 1 {len(nodes)}",
           f"3 1 0 {len(nodes)}"]
    out += [str(t) for t in nodes]
    out += [" ".join(f"{value:.15e}" for value in p)
            for p in nodes.values()]
    out += ["$EndNodes", "$Elements",
            f"1 {len(elements)} 1 {len(elements)}",
            f"3 1 4 {len(elements)}"]
    out += [f"{t} " + " ".join(str(x) for x in ns)
            for t, ns in elements.items()]
    out += ["$EndElements", "$NodeData", "1", '"f"', "1", "0",
            "3", "0", str(ncomp), str(len(rows))]
    out += [f"{t} " + " ".join(f"{v:.9e}" for v in vals)
            for t, vals in rows.items()]
    out += ["$EndNodeData"]
    path.write_text("\n".join(out) + "\n", encoding="utf-8")
    return path


def _gray(png):
    return np.asarray(Image.open(png).convert("L")).astype(float) / 255.0


def _rgb(png):
    return np.asarray(Image.open(png).convert("RGB")).astype(float) / 255.0


# --------------------------------------------------------------------
# volume_raycast
# --------------------------------------------------------------------

def test_raycast_transmittance_matches_the_closed_form(tmp_path):
    """Uniform medium: T_min == (1 - alpha)^(deepest inside count),
    exactly -- the front-to-back product, not an approximation."""
    f = _box_msh(tmp_path / "u.msh", lambda p: [5.0])
    a = 0.11
    res = volume_raycast(f, tmp_path / "u.png", view="f", grid=12,
                         view_dir="+z", image_size=96, n_steps=24,
                         value_range=[0.0, 5.0], alpha=a, alpha_power=0.0)
    assert res["ok"], res
    n = res["max_inside_samples"]
    assert n > 0
    assert res["transmittance_min"] == pytest.approx((1.0 - a) ** n,
                                                     rel=1e-9)


def test_raycast_occlusion_depends_on_the_viewing_side(tmp_path):
    """x-gradient medium at high opacity: the first thing each ray hits
    dominates, so +x and -x views have opposite dominant colours (jet:
    high = red, low = blue).  Only per-ray occlusion does this; a
    per-slice composite shows the same mix from both sides."""
    f = _box_msh(tmp_path / "g.msh",
                 lambda p: [max(0.0, min(1.0, 0.5 * (p[0] + 1.0)))])
    kw = {"view": "f", "grid": 14, "image_size": 96, "n_steps": 28,
          "value_range": [0.0, 1.0], "alpha": 0.9, "alpha_power": 1.0,
          "colorbar": False}
    hi = volume_raycast(f, tmp_path / "px.png", view_dir="+x", **kw)
    lo = volume_raycast(f, tmp_path / "mx.png", view_dir="-x", **kw)
    assert hi["ok"] and lo["ok"], (hi, lo)

    def centre_mean(png):
        img = _rgb(png)
        h, w, _ = img.shape
        c = img[h // 3:2 * h // 3, w // 3:2 * w // 3]
        return c[..., 0].mean(), c[..., 2].mean()      # red, blue

    r_hi, b_hi = centre_mean(tmp_path / "px.png")
    r_lo, b_lo = centre_mean(tmp_path / "mx.png")
    assert r_hi > b_hi, "+x view must be red-dominant (sees t=1 first)"
    assert b_lo > r_lo, "-x view must be blue-dominant (sees t=0 first)"


def test_raycast_outside_is_transparent(tmp_path):
    f = _box_msh(tmp_path / "u.msh", lambda p: [1.0])
    res = volume_raycast(f, tmp_path / "iso.png", view="f", grid=12,
                         view_dir="iso", image_size=96, n_steps=24,
                         alpha=0.4, colorbar=False)
    assert res["ok"], res
    assert res["constant_auto_range"] is True
    assert res["transmittance_min"] < 1.0
    # an oblique view leaves image corners outside the cube: white
    img = _rgb(tmp_path / "iso.png")
    assert res["found_fraction"] == pytest.approx(1.0)
    assert float(img[3, 3].min()) > 0.95


def test_raycast_rejects_bad_arguments(tmp_path):
    f = _box_msh(tmp_path / "u.msh", lambda p: [1.0])
    with pytest.raises(ValueError, match="grid"):
        volume_raycast(f, tmp_path / "x.png", grid=4)
    with pytest.raises(ValueError, match="alpha"):
        volume_raycast(f, tmp_path / "x.png", alpha=0.0)
    with pytest.raises(ValueError, match="view_dir"):
        volume_raycast(f, tmp_path / "x.png", view_dir="+w")
    bad = volume_raycast(f, tmp_path / "x.png", view="nope", grid=8,
                         image_size=64)
    assert not bad["ok"] and "no view named" in bad["error"]


# --------------------------------------------------------------------
# lic
# --------------------------------------------------------------------

def _aniso(gray, mask=None):
    """(along-x, along-y) mean squared finite differences."""
    g = gray if mask is None else np.where(mask, gray, np.nan)
    dx = np.nanmean(np.diff(g, axis=1) ** 2)
    dy = np.nanmean(np.diff(g, axis=0) ** 2)
    return dx, dy


def test_lic_smears_along_a_uniform_field(tmp_path):
    f = _box_msh(tmp_path / "vx.msh", lambda p: [1.0, 0.0, 0.0], ncomp=3)
    res = lic(f, tmp_path / "vx.png", view="f", plane="xy",
              resolution=96, kernel=12, color_by_magnitude=False)
    assert res["ok"], res
    assert res["found_fraction"] == pytest.approx(1.0)
    g = _gray(tmp_path / "vx.png")
    h, w = g.shape
    core = g[h // 4:3 * h // 4, w // 4:3 * w // 4]     # inside the axes box
    dx, dy = _aniso(core)
    assert dx < 0.35 * dy, (dx, dy)


def test_lic_anisotropy_rotates_with_the_field(tmp_path):
    f = _box_msh(tmp_path / "vy.msh", lambda p: [0.0, 1.0, 0.0], ncomp=3)
    res = lic(f, tmp_path / "vy.png", view="f", plane="xy",
              resolution=96, kernel=12, color_by_magnitude=False)
    assert res["ok"], res
    g = _gray(tmp_path / "vy.png")
    h, w = g.shape
    core = g[h // 4:3 * h // 4, w // 4:3 * w // 4]
    dx, dy = _aniso(core)
    assert dy < 0.35 * dx, (dx, dy)


def test_lic_is_deterministic_for_a_seed(tmp_path):
    f = _box_msh(tmp_path / "v.msh", lambda p: [1.0, 0.5, 0.0], ncomp=3)
    a = lic(f, tmp_path / "a.png", view="f", plane="xy", resolution=80,
            kernel=8, seed=7)
    b = lic(f, tmp_path / "b.png", view="f", plane="xy", resolution=80,
            kernel=8, seed=7)
    assert a["ok"] and b["ok"]
    assert np.array_equal(np.asarray(Image.open(tmp_path / "a.png")),
                          np.asarray(Image.open(tmp_path / "b.png")))


def test_lic_guards(tmp_path):
    vec = _box_msh(tmp_path / "v.msh", lambda p: [1.0, 0.0, 0.0], ncomp=3)
    with pytest.raises(ValueError, match="unknown plane"):
        lic(vec, plane="ab")
    with pytest.raises(ValueError, match="resolution"):
        lic(vec, resolution=32)
    with pytest.raises(ValueError, match="kernel"):
        lic(vec, kernel=1)
    off = lic(vec, tmp_path / "x.png", view="f", plane="xy", offset=5.0,
              resolution=64)
    assert not off["ok"] and "outside the mesh" in off["error"]
    scal = _box_msh(tmp_path / "s.msh", lambda p: [1.0])
    res = lic(scal, tmp_path / "y.png", view="f", plane="xy",
              resolution=64)
    assert not res["ok"] and "vector view" in res["error"]


# --------------------------------------------------------------------
# STEP overlay: depth-composited CAD (raycast) + section outline (lic)
# --------------------------------------------------------------------

def test_raycast_step_depth_ordering(tmp_path):
    """A CAD plate on the NEAR side hides the volume (achromatic gray
    centre); the same plate on the FAR side is hidden BY a nearly
    opaque volume (red-dominant centre).  Only genuine per-ray depth
    compositing distinguishes the two -- pasting the CAD on top would
    make both gray."""
    occ = pytest.importorskip("netgen.occ", reason="netgen.occ not installed")

    f = _box_msh(tmp_path / "u.msh", lambda p: [1.0])
    near = tmp_path / "near.step"
    far = tmp_path / "far.step"
    occ.Box(occ.Pnt(-2, -2, 1.2), occ.Pnt(2, 2, 1.4)).WriteStep(str(near))
    occ.Box(occ.Pnt(-2, -2, -1.4), occ.Pnt(2, 2, -1.2)).WriteStep(str(far))
    kw = {"view": "f", "grid": 12, "view_dir": "+z", "image_size": 96,
          "n_steps": 24, "value_range": [0.0, 1.0], "alpha": 0.9,
          "alpha_power": 0.0, "colorbar": False,
          "step_color": (0.6, 0.6, 0.6)}

    a = volume_raycast(f, tmp_path / "near.png", step_files=[near], **kw)
    b = volume_raycast(f, tmp_path / "far.png", step_files=[far], **kw)
    assert a["ok"] and b["ok"], (a, b)
    assert a["cad_triangles"] > 0
    assert a["cad_covered_fraction"] > 0.5

    def centre(png):
        img = _rgb(png)
        h, w, _ = img.shape
        return img[h // 3:2 * h // 3, w // 3:2 * w // 3]

    c_near = centre(tmp_path / "near.png")
    c_far = centre(tmp_path / "far.png")
    # near plate: gray CAD (channels equal); far plate: red volume wins
    assert abs(float(c_near[..., 0].mean() - c_near[..., 2].mean())) < 0.06
    assert float(c_far[..., 0].mean() - c_far[..., 2].mean()) > 0.3


def test_lic_step_outline_draws_on_the_texture(tmp_path):
    occ = pytest.importorskip("netgen.occ", reason="netgen.occ not installed")

    f = _box_msh(tmp_path / "vx.msh", lambda p: [1.0, 0.0, 0.0], ncomp=3)
    step = tmp_path / "square.step"
    occ.Box(occ.Pnt(-0.5, -0.5, -0.5), occ.Pnt(0.5, 0.5, 0.5)
            ).WriteStep(str(step))
    plain = lic(f, tmp_path / "plain.png", view="f", plane="xy",
                resolution=96, kernel=8, seed=3)
    lined = lic(f, tmp_path / "lined.png", view="f", plane="xy",
                resolution=96, kernel=8, seed=3, step_files=[step])
    assert plain["ok"] and lined["ok"], (plain, lined)
    assert lined["step_outline_segments"] > 0

    a = np.asarray(Image.open(tmp_path / "plain.png").convert("L"),
                   dtype=float)
    b = np.asarray(Image.open(tmp_path / "lined.png").convert("L"),
                   dtype=float)
    changed = np.abs(a - b) > 8
    # the outline changes SOME pixels, darkens them, and stays thin
    assert changed.any()
    assert changed.mean() < 0.15
    assert float(b[changed].mean()) < float(a[changed].mean())


def test_raster_step_missing_file_fails_loudly(tmp_path):
    f = _box_msh(tmp_path / "u.msh", lambda p: [1.0])
    r = volume_raycast(f, tmp_path / "x.png", view="f", grid=8,
                       image_size=64, step_files=[tmp_path / "no.step"])
    assert not r["ok"] and "not found" in r["error"]
    v = _box_msh(tmp_path / "v.msh", lambda p: [1.0, 0.0, 0.0], ncomp=3)
    l2 = lic(v, tmp_path / "y.png", view="f", plane="xy", resolution=64,
             step_files=[tmp_path / "no.step"])
    assert not l2["ok"] and "not found" in l2["error"]
