"""Fail-fast guards for the gmsh post-processing verbs.

Every test here pins a defect that used to be SILENT or FATAL:

* ``gmsh.view.probe(step=N)`` with N >= NbTimeStep is a native
  out-of-bounds read that KILLS the child process ("access violation
  reading 0xFFFFFFFFFFFFFFFF", MEASURED gmsh 4.15.2) -- every probing
  verb must validate the step first.
* ``Plugin(ModulusPhase)`` / ``Plugin(HarmonicToTime)`` on a
  single-step view fail with gmsh's useless "Unknown plugin or plugin
  action", naming neither the view nor the step count.
* ``field_histogram`` used to rewrite a reversed user ``value_range``
  into ``[lo, lo + 1]`` and return all-zero counts.
* List-data views exported into a ``.msh`` target come back as legacy
  MSH v2.2, which this lane's own v4.1-only readers reject.

Fixtures follow the house pattern: a hand-written Kuhn 6-tet cube
[-0.5, 0.5]^3 with P1-exact NodeData, so every expected number is a
closed form.
"""

import importlib.util

import pytest
from radia_mcp.gmsh import post_process
from radia_mcp.gmsh.post_process import (
    field_histogram,
    harmonic_to_time,
    modulus_phase,
    probe_field,
    smooth_to_nodes,
    streamlines,
    streamlines_2d,
    warp_view,
)

_GMSH_AVAILABLE = importlib.util.find_spec("gmsh") is not None
pytestmark_gmsh = pytest.mark.skipif(not _GMSH_AVAILABLE,
                                     reason="gmsh package not installed")

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


def _nodedata(name, ncomp, rows, time=0.0, step=0):
    lines = ["$NodeData", "1", f'"{name}"', "1", str(time), "3",
             str(step), str(ncomp), str(len(rows))]
    for tag, vals in rows:
        lines.append(str(tag) + " " + " ".join(f"{v:.16g}" for v in vals))
    lines.append("$EndNodeData")
    return "\n".join(lines) + "\n"


def _write(tmp_path, text, name):
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return p


@pytest.fixture
def one_step_vector(tmp_path):
    """Uniform B = (0, 0, 1) T -- ONE time step."""
    return _write(tmp_path,
                  _CUBE + _nodedata("B", 3,
                                    [(i, [0.0, 0.0, 1.0]) for i in _P]),
                  "one_step.msh")


@pytest.fixture
def two_step_scalar(tmp_path):
    """A re/im pair: step 0 = 3, step 1 = 4 -> |Z| = 5, phase = 0.9273."""
    text = (_CUBE
            + _nodedata("Z", 1, [(i, [3.0]) for i in _P], step=0)
            + _nodedata("Z", 1, [(i, [4.0]) for i in _P], step=1,
                        time=1.0))
    return _write(tmp_path, text, "two_step.msh")


@pytest.fixture
def three_step_scalar(tmp_path):
    text = _CUBE + "".join(
        _nodedata("Z", 1, [(i, [float(s + 3)]) for i in _P], step=s,
                  time=float(s))
        for s in range(3))
    return _write(tmp_path, text, "three_step.msh")


# ----------------------------------------------------------------------
# E1: an out-of-range step must not fault the child process
# ----------------------------------------------------------------------

@pytestmark_gmsh
def test_probe_field_rejects_a_step_past_the_last_one(one_step_vector):
    result = probe_field(one_step_vector, [[0.0, 0.0, 0.0]], step=7)
    assert result["ok"] is False
    assert "access violation" not in result["error"]
    assert "step 7 out of range" in result["error"]
    assert "view 'B' has 1 step(s) (0..0)" in result["error"]


@pytestmark_gmsh
def test_probe_field_still_accepts_the_all_steps_sentinel(one_step_vector):
    """step=-1 means "every step" and must survive the new gate."""
    result = probe_field(one_step_vector, [[0.0, 0.0, 0.0]], step=-1)
    assert result["ok"] is True, result.get("error")
    entry = result["views"][0]["points"][0]
    assert entry["found"] is True
    assert entry["values"] == pytest.approx([0.0, 0.0, 1.0], abs=1e-12)


@pytestmark_gmsh
def test_probe_field_rejects_a_step_below_minus_one(one_step_vector):
    result = probe_field(one_step_vector, [[0.0, 0.0, 0.0]], step=-5)
    assert result["ok"] is False
    assert "must be -1 (all steps) or non-negative" in result["error"]


@pytestmark_gmsh
def test_streamlines_rejects_a_step_past_the_last_one(one_step_vector,
                                                     tmp_path):
    result = streamlines(one_step_vector, [0.0, 0.0, 0.0],
                         [0.1, 0.0, 0.0], n_seeds=1, time_step=5,
                         max_steps=4, out_file=tmp_path / "s.pos")
    assert result["ok"] is False
    assert "access violation" not in result["error"]
    assert "time_step 5 out of range" in result["error"]


@pytestmark_gmsh
def test_streamlines_2d_rejects_a_step_past_the_last_one(one_step_vector,
                                                        tmp_path):
    result = streamlines_2d(one_step_vector, [-0.4, -0.4, 0.0],
                            [0.4, -0.4, 0.0], [-0.4, 0.4, 0.0],
                            time_step=4, max_lines=2,
                            out_file=tmp_path / "s2.pos")
    assert result["ok"] is False
    assert "access violation" not in result["error"]
    assert "time_step 4 out of range" in result["error"]


# ----------------------------------------------------------------------
# LENSA-04/05: the complex-pair verbs need a real two-step view
# ----------------------------------------------------------------------

@pytestmark_gmsh
def test_modulus_phase_rejects_a_single_step_view(one_step_vector,
                                                  tmp_path):
    result = modulus_phase(one_step_vector, out_file=tmp_path / "mp.pos")
    assert result["ok"] is False
    assert "Unknown plugin" not in result["error"]
    assert "view 'B' has 1 time step" in result["error"]
    assert "modulus/phase needs a two-step re/im view" in result["error"]


@pytestmark_gmsh
def test_harmonic_to_time_rejects_a_single_step_view(one_step_vector,
                                                     tmp_path):
    result = harmonic_to_time(one_step_vector, out_file=tmp_path / "h.pos")
    assert result["ok"] is False
    assert "Unknown plugin" not in result["error"]
    assert "view 'B' has 1 time step" in result["error"]
    assert "harmonic_to_time needs a two-step re/im view" in result["error"]


@pytestmark_gmsh
def test_modulus_phase_on_a_true_pair_is_unchanged_and_unwarned(
        two_step_scalar, tmp_path):
    """3 + 4i -> |Z| = 5, phase = atan2(4, 3) = 0.9273 (the documented
    measurement); a genuine 2-step pair must not gain a warning."""
    out = tmp_path / "mp.pos"
    result = modulus_phase(two_step_scalar, out_file=out)
    assert result["ok"] is True, result.get("error")
    assert "warning" not in result
    assert result["n_steps"] == 2
    got = probe_field(out, [[0.0, 0.0, 0.0]], step=-1)
    assert got["ok"] is True, got.get("error")
    steps = got["views"][0]["points"][0]["steps"]
    assert steps[0][0] == pytest.approx(5.0, abs=1e-9)
    assert steps[1][0] == pytest.approx(0.9272952180016122, abs=1e-9)


@pytestmark_gmsh
def test_multi_step_view_with_default_pair_steps_warns(three_step_scalar,
                                                       tmp_path):
    """LENSA-04: >2 steps + untouched real_step/imag_step is a guess.

    It must not FAIL (a legitimate multi-step file may still carry re/im
    at 0/1) but the assumption has to be stated -- the beam-animation
    views this lane mints are many-step and are NOT a re/im pair.
    """
    result = modulus_phase(three_step_scalar, out_file=tmp_path / "m.pos")
    assert result["ok"] is True, result.get("error")
    assert "3 time steps" in result["warning"]
    assert "0 / 1 defaults" in result["warning"]

    hto = harmonic_to_time(three_step_scalar, n_steps=4,
                           out_file=tmp_path / "h.pos")
    assert hto["ok"] is True, hto.get("error")
    assert "3 time steps" in hto["warning"]


@pytestmark_gmsh
def test_explicit_pair_steps_silence_the_warning(three_step_scalar,
                                                 tmp_path):
    result = modulus_phase(three_step_scalar, real_step=0, imag_step=2,
                           out_file=tmp_path / "m2.pos")
    assert result["ok"] is True, result.get("error")
    assert "warning" not in result


@pytestmark_gmsh
def test_pair_step_index_out_of_range_is_named(three_step_scalar,
                                               tmp_path):
    result = modulus_phase(three_step_scalar, real_step=5,
                           out_file=tmp_path / "m3.pos")
    assert result["ok"] is False
    assert "real_step 5 out of range" in result["error"]
    assert "3 step(s) (0..2)" in result["error"]


# ----------------------------------------------------------------------
# LENSA-06: field_histogram must not invent a value_range
# ----------------------------------------------------------------------

@pytest.mark.parametrize(
    ("value_range", "message"),
    [
        ([5.0, 1.0], "value_range must be increasing, got [5.0, 1.0]"),
        ([1.0, 1.0], "value_range must be increasing, got [1.0, 1.0]"),
        ([0.0, float("inf")], "value_range must be finite"),
        ([1.0], "value_range must be [lo, hi]"),
    ],
)
def test_field_histogram_rejects_a_bad_value_range(tmp_path, value_range,
                                                   message):
    msh = _write(tmp_path,
                 _CUBE + _nodedata("f", 1, [(i, [float(i)]) for i in _P]),
                 "f.msh")
    result = field_histogram(msh, value_range=value_range)
    assert result["ok"] is False
    assert message in result["error"]
    # the old code answered ok=true with an invented range + zero counts
    assert "counts" not in result


def test_field_histogram_keeps_a_valid_range(tmp_path):
    """Values 1..8 on the 8 nodes: [0.5, 8.5] in 8 unit bins centres
    each value in its own bin, so every count is exactly 1."""
    msh = _write(tmp_path,
                 _CUBE + _nodedata("f", 1, [(i, [float(i)]) for i in _P]),
                 "f.msh")
    result = field_histogram(msh, bins=8, value_range=[0.5, 8.5])
    assert result["ok"] is True, result.get("error")
    assert result["counts"] == [1] * 8
    assert result["n_samples"] == 8
    assert result["bin_edges"][0] == pytest.approx(0.5)
    assert result["bin_edges"][-1] == pytest.approx(8.5)


def test_field_histogram_reports_a_constant_field_widening(tmp_path):
    """A DERIVED degenerate range is data, not user input: widen but say so."""
    msh = _write(tmp_path,
                 _CUBE + _nodedata("f", 1, [(i, [2.5]) for i in _P]),
                 "const.msh")
    result = field_histogram(msh, bins=4)
    assert result["ok"] is True, result.get("error")
    assert result["bin_edges"][0] == pytest.approx(2.5)
    assert result["bin_edges"][-1] == pytest.approx(3.5)
    assert sum(result["counts"]) == 8
    assert "every sample equals" in result["note"]


# ----------------------------------------------------------------------
# LENSA-03: list-data views export as .pos, never as .msh
# ----------------------------------------------------------------------

@pytest.fixture
def forbid_gmsh_launch(monkeypatch):
    def unexpected(*_args, **_kwargs):
        raise AssertionError("invalid input reached the Gmsh subprocess")

    monkeypatch.setattr(post_process, "run_gmsh_json_subprocess",
                        unexpected)


@pytest.mark.parametrize("suffix", [".msh", ".vtk", ""])
def test_streamlines_rejects_a_non_pos_out_file(tmp_path, suffix,
                                                forbid_gmsh_launch):
    msh = _write(tmp_path,
                 _CUBE + _nodedata("B", 3,
                                   [(i, [0.0, 0.0, 1.0]) for i in _P]),
                 "b.msh")
    result = streamlines(msh, [0.0, 0.0, 0.0], [0.1, 0.0, 0.0],
                         n_seeds=1, max_steps=4,
                         out_file=tmp_path / f"lines{suffix}")
    assert result["ok"] is False
    assert "must end in '.pos'" in result["error"]
    assert "v2.2" in result["error"]


@pytest.mark.parametrize(
    ("verb", "call"),
    [
        ("isosurface", lambda p, o: post_process.isosurface(p, 4.5,
                                                            out_file=o)),
        ("cut_plane", lambda p, o: post_process.cut_plane_extract(
            p, [0.0, 0.0, 1.0], 0.0, out_file=o)),
        ("math_eval", lambda p, o: post_process.math_eval(p, ["v0"],
                                                          out_file=o)),
        ("threshold", lambda p, o: post_process.threshold(p, 0.0, 9.0,
                                                          out_file=o)),
        ("flux_lines", lambda p, o: post_process.flux_lines(p, n_levels=3,
                                                            out_file=o)),
        ("extract_skin", lambda p, o: post_process.extract_skin(
            p, out_file=o)),
    ],
)
def test_every_list_data_verb_rejects_a_msh_out_file(tmp_path, verb, call,
                                                     forbid_gmsh_launch):
    """MEASURED: each of these wrote "$MeshFormat / 2.2 0 8" into .msh."""
    msh = _write(tmp_path,
                 _CUBE + _nodedata("f", 1, [(i, [float(i)]) for i in _P]),
                 "f.msh")
    result = call(msh, tmp_path / "out.msh")
    assert result["ok"] is False, verb
    assert "must end in '.pos'" in result["error"], verb


@pytestmark_gmsh
def test_in_place_verbs_still_accept_a_msh_out_file(one_step_vector,
                                                    two_step_scalar,
                                                    tmp_path):
    """warp / smooth / modulus_phase rewrite the MODEL-based view.

    MEASURED: those three produced "$MeshFormat / 4.1 0 8" in a .msh
    target while every list-data verb produced 2.2 -- so the guard must
    NOT fire for them.
    """
    cases = [
        ("warp", tmp_path / "w.msh",
         lambda p: warp_view(one_step_vector, 1.0, out_file=p)),
        ("smooth", tmp_path / "s.msh",
         lambda p: smooth_to_nodes(one_step_vector, out_file=p)),
        ("modulus_phase", tmp_path / "m.msh",
         lambda p: modulus_phase(two_step_scalar, out_file=p)),
    ]
    for name, target, run in cases:
        result = run(target)
        assert result["ok"] is True, (name, result.get("error"))
        lines = target.read_text(encoding="utf-8").splitlines()
        assert lines[0] == "$MeshFormat", name
        assert lines[1].startswith("4.1"), (name, lines[1])
