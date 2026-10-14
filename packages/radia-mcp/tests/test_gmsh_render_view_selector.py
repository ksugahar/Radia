"""View selector + fail-loud result contract for the gmsh render lane.

Fixture: the house Kuhn 6-tet cube [-0.5, 0.5]^3 (see
test_gmsh_particle_trace.py) carrying TWO scalar NodeData views, "A" =
the node x coordinate and "B" = the node z coordinate.  P1 interpolation
is exact for both, and the two fields have orthogonal gradients, so
isolating one or the other is visible in the rendered pixels.

The view-selector goldens are EQUIVALENCE goldens: ``view="A"`` must
produce the same picture, pixel for pixel, as the hand-written
``options={"View[0].Visible": 1, "View[1].Visible": 0}`` route it
replaces.  Renders are compared on DECODED PIXELS, never on file bytes:
MEASURED on gmsh 4.15.2, two renders of identical settings give
identical pixels but different PNG bytes.
"""

import hashlib
import importlib.util

import pytest
from radia_mcp.gmsh.render import export_animation, render_png

_GMSH_AVAILABLE = importlib.util.find_spec("gmsh") is not None
_PIL_AVAILABLE = importlib.util.find_spec("PIL") is not None

# Kuhn decomposition of the cube [-0.5, 0.5]^3: 6 positively oriented
# tets sharing the diagonal node1->node8.
_CUBE = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 8 1 8
3 1 0 8
1
2
3
4
5
6
7
8
-0.5 -0.5 -0.5
0.5 -0.5 -0.5
-0.5 0.5 -0.5
0.5 0.5 -0.5
-0.5 -0.5 0.5
0.5 -0.5 0.5
-0.5 0.5 0.5
0.5 0.5 0.5
$EndNodes
$Elements
1 6 1 6
3 1 4 6
1 1 2 4 8
2 1 2 8 6
3 1 3 8 4
4 1 3 7 8
5 1 5 6 8
6 1 5 8 7
$EndElements
"""

_P = {1: (-0.5, -0.5, -0.5), 2: (0.5, -0.5, -0.5), 3: (-0.5, 0.5, -0.5),
      4: (0.5, 0.5, -0.5), 5: (-0.5, -0.5, 0.5), 6: (0.5, -0.5, 0.5),
      7: (-0.5, 0.5, 0.5), 8: (0.5, 0.5, 0.5)}

# Wide enough that the FLTK sidebar still leaves a usable canvas: the
# exported width is NOT the requested width on Windows builds
# (MEASURED gmsh 4.15.2: width=300 exported 100 px).
_W, _H = 700, 400


def _nodedata(name, ncomp, rows, time=0.0, step=0):
    lines = ["$NodeData", "1", f'"{name}"', "1", str(time), "3",
             str(step), str(ncomp), str(len(rows))]
    for tag, vals in rows:
        lines.append(str(tag) + " " + " ".join(f"{v:.16g}" for v in vals))
    lines.append("$EndNodeData")
    return "\n".join(lines) + "\n"


@pytest.fixture
def two_view_msh(tmp_path):
    """Cube with view "A" = x coordinate and view "B" = z coordinate."""
    path = tmp_path / "two_views.msh"
    path.write_text(
        _CUBE
        + _nodedata("A", 1, [(t, [p[0]]) for t, p in _P.items()])
        + _nodedata("B", 1, [(t, [p[2]]) for t, p in _P.items()]),
        encoding="utf-8")
    return path


def _skip_if_no_graphics(result):
    if not result.get("ran"):
        pytest.skip(f"no gmsh graphics context: {result.get('error')}")


def _pixels(path):
    """Decoded-pixel digest.  PNG bytes are not reproducible across two
    identical gmsh renders (MEASURED); the pixels are."""
    from PIL import Image

    with Image.open(path) as im:
        return hashlib.sha256(im.convert("RGB").tobytes()).hexdigest()


# ======================================================================
# num_steps guard + "ok" can never survive a later failure
# ======================================================================

@pytest.mark.parametrize("bad", [0, -2])
def test_export_animation_rejects_nonpositive_num_steps(tmp_path,
                                                        two_view_msh, bad):
    """MEASURED before the guard: num_steps=0 (and -2) returned
    {"ok": True, "error": "IndexError: list index out of range",
    "num_steps": 0} and wrote no GIF -- a contradictory result the caller
    could not act on.  Reject the value where it is still in hand."""
    gif = tmp_path / "zero.gif"
    result = export_animation(two_view_msh, gif, num_steps=bad)
    assert result["ok"] is False, result
    assert "num_steps" in result["error"]
    assert str(bad) in result["error"]
    assert ">= 1" in result["error"]        # names the valid alternative
    assert "None" in result["error"]
    assert not gif.exists()


@pytest.mark.parametrize("bad", [1.5, True, "3"])
def test_export_animation_rejects_non_integer_num_steps(tmp_path,
                                                        two_view_msh, bad):
    result = export_animation(two_view_msh, tmp_path / "x.gif",
                              num_steps=bad)
    assert result["ok"] is False, result
    assert "integer" in result["error"]
    assert repr(bad) in result["error"]


@pytest.mark.skipif(not _GMSH_AVAILABLE, reason="gmsh package not installed")
def test_worker_never_reports_ok_after_a_late_failure(tmp_path,
                                                      two_view_msh):
    """The worker sets ok=True BEFORE its last work item (GIF assembly),
    so a failure there used to be reported as success.  Drive the worker
    directly with the rejected config to prove the second half of the
    fix: ok must be False even though the animation stage already
    succeeded."""
    from radia_mcp.gmsh.render import _run_render

    gif = tmp_path / "raw.gif"
    cfg = {
        "mode": "animation", "path": str(two_view_msh), "merge_files": [],
        "geometry_display": False, "gif_out": str(gif),
        "frames_dir": str(tmp_path).replace("\\", "/"),
        "width": 200, "height": 160, "numsubedges": 4, "rotation": None,
        "view_select": None, "view_indices": None,
        "num_steps": 0,                     # -> zero frames
        "delay_ms": 40, "color": None, "glyphs": None, "clip_planes": [],
        "axes": None, "annotations": [], "options": {},
        "string_options": {}, "auto_mesh_display": False,
        "adapt_views": True, "smooth_normals": True, "link_views": True,
        "clip_views": False, "orbit": None, "time_step": None,
    }
    result = _run_render(cfg, 300.0)
    _skip_if_no_graphics(result)
    assert result["ok"] is False, result
    assert result["error"]
    assert "GIF" in result["error"]         # names what actually failed
    assert not gif.exists()


# ======================================================================
# view selector: shape guards (no gmsh subprocess involved)
# ======================================================================

@pytest.mark.parametrize("bad, exc", [(True, TypeError), (-1, ValueError),
                                      ([], ValueError), (1.5, TypeError),
                                      ("", ValueError), ([None], TypeError)])
def test_view_selector_rejects_bad_shapes(two_view_msh, tmp_path, bad, exc):
    with pytest.raises(exc):
        render_png(two_view_msh, tmp_path / "x.png", view=bad)


# ======================================================================
# view selector: pixel goldens
# ======================================================================

@pytest.mark.skipif(not _GMSH_AVAILABLE or not _PIL_AVAILABLE,
                    reason="gmsh or Pillow not installed")
def test_view_selector_matches_the_raw_option_route(tmp_path, two_view_msh):
    """EQUIVALENCE golden: view="A" is exactly the 2-option hand-written
    isolation it replaces, and isolating A, isolating B, and showing both
    are three different pictures."""
    both = tmp_path / "both.png"
    only_a = tmp_path / "a.png"
    only_b = tmp_path / "b.png"
    raw_a = tmp_path / "raw_a.png"

    r_both = render_png(two_view_msh, both, width=_W, height=_H)
    _skip_if_no_graphics(r_both)
    assert r_both["ok"] is True, r_both
    assert r_both["n_views"] == 2

    r_a = render_png(two_view_msh, only_a, width=_W, height=_H, view="A")
    assert r_a["ok"] is True, r_a
    assert r_a["view_names"] == ["A", "B"]
    assert r_a["view_selected"] == [0]
    assert r_a["blank_check"]["looks_blank"] is False, r_a

    r_b = render_png(two_view_msh, only_b, width=_W, height=_H, view="B")
    assert r_b["ok"] is True, r_b
    assert r_b["view_selected"] == [1]

    r_raw = render_png(two_view_msh, raw_a, width=_W, height=_H,
                       options={"View[0].Visible": 1, "View[1].Visible": 0})
    assert r_raw["ok"] is True, r_raw

    assert _pixels(only_a) == _pixels(raw_a)
    assert _pixels(only_a) != _pixels(only_b)
    assert _pixels(only_a) != _pixels(both)
    assert _pixels(only_b) != _pixels(both)


@pytest.mark.skipif(not _GMSH_AVAILABLE or not _PIL_AVAILABLE,
                    reason="gmsh or Pillow not installed")
def test_view_selector_list_shows_every_named_view(tmp_path, two_view_msh):
    """view=[0, 1] hides nothing, so it must reproduce the no-selector
    render pixel for pixel."""
    both = tmp_path / "both.png"
    listed = tmp_path / "listed.png"
    by_name = tmp_path / "by_name.png"

    r_both = render_png(two_view_msh, both, width=_W, height=_H)
    _skip_if_no_graphics(r_both)
    r_listed = render_png(two_view_msh, listed, width=_W, height=_H,
                          view=[0, 1])
    assert r_listed["ok"] is True, r_listed
    assert r_listed["view_selected"] == [0, 1]
    r_named = render_png(two_view_msh, by_name, width=_W, height=_H,
                         view=["B", "A"])
    assert r_named["view_selected"] == [0, 1]     # normalized, deduped

    assert _pixels(listed) == _pixels(both)
    assert _pixels(by_name) == _pixels(both)


@pytest.mark.skipif(not _GMSH_AVAILABLE or not _PIL_AVAILABLE,
                    reason="gmsh or Pillow not installed")
def test_explicit_visible_option_overrides_the_selector(tmp_path,
                                                        two_view_msh):
    """The selector is applied with the structured controls, BEFORE the
    raw option passthrough, so a hand-set View[i].Visible still wins."""
    both = tmp_path / "both.png"
    override = tmp_path / "override.png"

    r_both = render_png(two_view_msh, both, width=_W, height=_H)
    _skip_if_no_graphics(r_both)
    r_over = render_png(two_view_msh, override, width=_W, height=_H,
                        view="A", options={"View[1].Visible": 1})
    assert r_over["ok"] is True, r_over
    assert r_over["view_selected"] == [0]      # the selector still ran
    assert _pixels(override) == _pixels(both)  # the raw option won


@pytest.mark.skipif(not _GMSH_AVAILABLE, reason="gmsh package not installed")
def test_unknown_view_name_fails_loud_with_the_available_names(
        tmp_path, two_view_msh):
    # No graphics guard: name resolution runs before fltk.initialize, so
    # this path is reachable even on a machine with no OpenGL context.
    result = render_png(two_view_msh, tmp_path / "x.png", width=300,
                        height=240, view="B_magnitude")
    assert result["ok"] is False, result
    assert "B_magnitude" in result["error"]
    assert "'A'" in result["error"] and "'B'" in result["error"]


@pytest.mark.skipif(not _GMSH_AVAILABLE, reason="gmsh package not installed")
def test_out_of_range_view_index_fails_loud(tmp_path, two_view_msh):
    result = render_png(two_view_msh, tmp_path / "x.png", width=300,
                        height=240, view=7)
    assert result["ok"] is False, result
    assert "7" in result["error"]
    assert "[0, 2)" in result["error"]
    assert "'A'" in result["error"] and "'B'" in result["error"]


# ======================================================================
# animation: the selector doubles as the animation target
# ======================================================================

@pytest.mark.skipif(not _GMSH_AVAILABLE or not _PIL_AVAILABLE,
                    reason="gmsh or Pillow not installed")
def test_animation_view_defaults_view_indices(tmp_path, two_view_msh):
    """One argument instead of two: view="B" both hides view A and makes
    B the stepped target (the notebook passed view_indices=[4] PLUS six
    View[i].Visible entries for exactly this)."""
    gif = tmp_path / "anim.gif"
    result = export_animation(two_view_msh, gif, view="B",
                              width=300, height=240)
    _skip_if_no_graphics(result)
    assert result["ok"] is True, result
    assert result["view_selected"] == [1]
    assert result["view_indices"] == [1]
    assert gif.is_file()


@pytest.mark.skipif(not _GMSH_AVAILABLE or not _PIL_AVAILABLE,
                    reason="gmsh or Pillow not installed")
def test_explicit_view_indices_beats_the_selector(tmp_path, two_view_msh):
    result = export_animation(two_view_msh, tmp_path / "anim.gif",
                              view="B", view_indices=[0],
                              width=300, height=240)
    _skip_if_no_graphics(result)
    assert result["ok"] is True, result
    assert result["view_selected"] == [1]      # visibility still from view
    assert result["view_indices"] == [0]       # stepping still explicit
