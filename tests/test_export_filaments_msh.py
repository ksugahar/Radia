"""Tests for export_filaments_msh (PEEC filament -> GMSH view export).

Locks the base |I| ElementData contract and the opt-in extensions:
"I direction" (unit tangent x |I| arrows) and "I_complex" (re/im two
time steps, the harmonic-to-time contract for AC current animation).
"""

import importlib.util
import math
import re
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from radia.gmsh_post_export import export_filaments_msh  # noqa: E402

_GMSH_AVAILABLE = importlib.util.find_spec("gmsh") is not None

# two straight filaments along +x at y = 0 / 0.01, two segments each
_PATHS = [
    [((0.0, 0.0, 0.0), (0.1, 0.0, 0.0)),
     ((0.1, 0.0, 0.0), (0.2, 0.0, 0.0))],
    [((0.0, 0.01, 0.0), (0.1, 0.01, 0.0)),
     ((0.1, 0.01, 0.0), (0.2, 0.01, 0.0))],
]
_CURRENTS = np.array([3.0 + 4.0j, 1.0 - 1.0j])


def _sections(path, kind="ElementData"):
    text = Path(path).read_text(encoding="utf-8")
    out = []
    for m in re.finditer(rf"\${kind}\n(.*?)\$End{kind}", text, re.S):
        body = m.group(1).strip().splitlines()
        name = body[1].strip('"')
        n_str = int(body[0])
        idx = 1 + n_str
        n_real = int(body[idx])
        idx += 1 + n_real
        n_int = int(body[idx])
        ints = [int(v) for v in body[idx + 1:idx + 1 + n_int]]
        rows = [ln.split() for ln in body[idx + 1 + n_int:]]
        out.append({"name": name, "step": ints[0], "ncomp": ints[1],
                    "rows": rows})
    return out


@pytest.mark.parametrize(
    ("paths", "currents", "kwargs", "message"),
    [
        ([], None, {}, "at least one filament"),
        ([[]], None, {}, "has no segments"),
        (_PATHS, [1.0], {}, "shape"),
        (_PATHS, [1.0, np.nan], {}, "non-finite"),
        (_PATHS, None, {"direction_view": True}, "currents is required"),
        (_PATHS, _CURRENTS, {"viz_subdivide_n": 0}, "positive integer"),
    ],
)
def test_invalid_input_fails_before_writing(
        tmp_path, paths, currents, kwargs, message):
    out = tmp_path / "invalid.msh"
    with pytest.raises(ValueError, match=message):
        export_filaments_msh(paths, str(out), currents=currents, **kwargs)
    assert not out.exists()


def test_base_current_magnitude_elementdata(tmp_path):
    out = tmp_path / "fil.msh"
    export_filaments_msh(_PATHS, str(out), currents=_CURRENTS)
    secs = _sections(out)
    assert [s["name"] for s in secs] == ["|I| [A]"]
    sec = secs[0]
    assert sec["ncomp"] == 1
    assert len(sec["rows"]) == 4  # 2 filaments x 2 segments
    values = sorted({float(r[1]) for r in sec["rows"]})
    assert values == [pytest.approx(math.sqrt(2.0)), pytest.approx(5.0)]


def test_direction_view_unit_tangent_times_current(tmp_path):
    out = tmp_path / "fil_dir.msh"
    export_filaments_msh(_PATHS, str(out), currents=_CURRENTS,
                         direction_view=True)
    secs = _sections(out)
    assert [s["name"] for s in secs] == ["|I| [A]", "I direction [A]"]
    dir_sec = secs[1]
    assert dir_sec["ncomp"] == 3
    for row in dir_sec["rows"]:
        vx, vy, vz = (float(row[1]), float(row[2]), float(row[3]))
        assert vy == pytest.approx(0.0, abs=1e-12)
        assert vz == pytest.approx(0.0, abs=1e-12)
        # +x tangent scaled by |I| of the owning filament
        assert vx in (pytest.approx(5.0), pytest.approx(math.sqrt(2.0)))


def test_complex_steps_write_re_im_pair(tmp_path):
    out = tmp_path / "fil_cplx.msh"
    export_filaments_msh(_PATHS, str(out), currents=_CURRENTS,
                         complex_steps=True)
    secs = _sections(out)
    names = [s["name"] for s in secs]
    assert names == ["|I| [A]", "I_complex [A]", "I_complex [A]"]
    re_sec = next(s for s in secs if s["name"] == "I_complex [A]"
                  and s["step"] == 0)
    im_sec = next(s for s in secs if s["name"] == "I_complex [A]"
                  and s["step"] == 1)
    re_vals = sorted({float(r[1]) for r in re_sec["rows"]})
    im_vals = sorted({float(r[1]) for r in im_sec["rows"]})
    assert re_vals == [pytest.approx(1.0), pytest.approx(3.0)]
    assert im_vals == [pytest.approx(-1.0), pytest.approx(4.0)]


def test_complex_steps_accept_real_currents_with_zero_imaginary_step(tmp_path):
    out = tmp_path / "fil_real.msh"
    export_filaments_msh(_PATHS, str(out), currents=[2.0, 3.0],
                         complex_steps=True)
    secs = [s for s in _sections(out) if s["name"] == "I_complex [A]"]
    assert [s["step"] for s in secs] == [0, 1]
    assert {float(row[1]) for row in secs[0]["rows"]} == {2.0, 3.0}
    assert {float(row[1]) for row in secs[1]["rows"]} == {0.0}


@pytest.mark.skipif(not _GMSH_AVAILABLE
                    or importlib.util.find_spec("netgen") is None,
                    reason="gmsh or netgen not installed")
def test_step_merge_with_shaded_faces_does_not_crash(tmp_path):
    # Regression: discrete-curve entity tags 1..N collide with merged
    # OCC curve tags; gmsh binds the filament mesh to the CAD curve
    # and hard-crashes (0xC000041D) once shaded geometry faces are on.
    # The entity-tag offset in export_filaments_msh prevents this.
    from netgen.occ import Box, Pnt

    step = tmp_path / "box.step"
    Box(Pnt(0, 0, -0.05), Pnt(0.2, 0.02, 0.05)).WriteStep(str(step))
    fil = tmp_path / "fil.msh"
    export_filaments_msh(_PATHS, str(fil), currents=_CURRENTS)

    import subprocess
    code = (
        "import gmsh\n"
        "gmsh.initialize(['-noconfig'])\n"
        "gmsh.option.setNumber('General.Terminal', 0)\n"
        f"gmsh.open({str(fil)!r})\n"
        f"gmsh.merge({str(step)!r})\n"
        "gmsh.option.setNumber('Geometry.Surfaces', 1)\n"
        "gmsh.option.setNumber('Geometry.SurfaceType', 2)\n"
        "gmsh.option.setNumber('General.GraphicsWidth', 320)\n"
        "gmsh.option.setNumber('General.GraphicsHeight', 280)\n"
        "gmsh.option.setNumber('General.GraphicsPositionX', 100)\n"
        "gmsh.option.setNumber('General.GraphicsPositionY', 100)\n"
        "gmsh.fltk.initialize()\n"
        "gmsh.fltk.update()\n"
        f"gmsh.write({str(tmp_path / 'overlay.png')!r})\n"
        "gmsh.fltk.finalize()\n"
        "gmsh.finalize()\n"
        "print('OK')\n"
    )
    proc = subprocess.run([sys.executable, "-u", "-c", code], text=True,
                          capture_output=True, timeout=180,
                          stdin=subprocess.DEVNULL)
    if proc.returncode != 0 and "OK" not in proc.stdout:
        if "fltk" in (proc.stdout + proc.stderr).lower():
            pytest.skip("no gmsh graphics context")
    assert proc.returncode == 0, proc.stdout + proc.stderr


@pytest.mark.skipif(not _GMSH_AVAILABLE, reason="gmsh not installed")
def test_gmsh_roundtrip_views(tmp_path):
    out = tmp_path / "fil_all.msh"
    export_filaments_msh(_PATHS, str(out), currents=_CURRENTS,
                         direction_view=True, complex_steps=True)
    import subprocess
    code = (
        "import json, sys, gmsh\n"
        "gmsh.initialize(['-noconfig'])\n"
        "gmsh.option.setNumber('General.Terminal', 0)\n"
        f"gmsh.open({str(out)!r})\n"
        "tags = list(gmsh.view.getTags())\n"
        "names = [gmsh.option.getString(f'View[{gmsh.view.getIndex(t)}]"
        ".Name') for t in tags]\n"
        "steps = [int(gmsh.option.getNumber(f'View[{gmsh.view.getIndex(t)}]"
        ".NbTimeStep')) for t in tags]\n"
        "gmsh.finalize()\n"
        "print(json.dumps({'names': names, 'steps': steps}))\n"
    )
    proc = subprocess.run([sys.executable, "-c", code], text=True,
                          capture_output=True, timeout=120,
                          stdin=subprocess.DEVNULL)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    import json
    info = json.loads(proc.stdout.strip().splitlines()[-1])
    assert info["names"] == ["|I| [A]", "I direction [A]", "I_complex [A]"]
    assert info["steps"] == [1, 1, 2]  # complex view carries re/im steps
