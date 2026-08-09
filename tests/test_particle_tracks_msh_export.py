"""export_particle_tracks_msh: Lorentz tracks -> GMSH line elements.

The fixture track is the CLOSED-FORM uniform-B orbit from
radia.particle_tracking.uniform_magnetic_trajectory (numpy only, no
scipy), wrapped in the radia-time-domain-lorentz-track/v1 schema that
track_lorentz_ivp emits -- so the exporter is tested against exact
circle geometry, and the schema contract is locked without the 'beam'
extra installed.
"""

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from radia.gmsh_post_export import export_particle_tracks_msh  # noqa: E402
from radia.particle_tracking import (  # noqa: E402
    ParticleSpecies,
    uniform_magnetic_trajectory,
)

QE = 1.602176634e-19
ME = 9.1093837015e-31
ELECTRON = ParticleSpecies(charge_c=-QE, mass_kg=ME)


def _circle_track(n=65, voltage=100.0, bz=1.0e-3):
    """One full nonrelativistic gyration as a v1 track dict."""
    omega = abs(ELECTRON.charge_c) * bz / ELECTRON.mass_kg
    period = 2.0 * math.pi / omega
    times = np.linspace(0.0, period, n)
    closed = uniform_magnetic_trajectory(
        ELECTRON, voltage, times, magnetic_flux_density_t=bz)
    speed2 = np.sum(closed["velocity_m_s"] ** 2, axis=1)
    return {
        "schema": "radia-time-domain-lorentz-track/v1",
        "time_s": times.tolist(),
        "position_m": closed["position_m"].tolist(),
        "velocity_m_s": closed["velocity_m_s"].tolist(),
        "kinetic_energy_j": (0.5 * ELECTRON.mass_kg * speed2).tolist(),
    }


def test_exports_circle_track_with_all_views(tmp_path):
    track = _circle_track()
    out = tmp_path / "orbit.msh"
    summary = export_particle_tracks_msh(track, out)
    assert out.is_file()
    assert summary["n_tracks"] == 1
    assert summary["per_track_segments"] == [64]
    assert summary["n_elements"] == 64
    assert summary["n_degenerate_skipped"] == 0
    assert summary["views"] == ["time [s]", "speed [m/s]",
                                "kinetic energy [eV]", "v [m/s]"]

    text = out.read_text()
    assert text.count("$ElementData") == 4
    assert '"track_0"' in text
    # kinetic energy view is in eV: 100 V acceleration -> 100 eV
    assert "$MeshFormat\n4.1 0 8" in text


def test_two_tracks_become_two_groups_and_stride_decimates(tmp_path):
    t0 = _circle_track(n=101)
    t1 = _circle_track(n=11, bz=2.0e-3)
    out = tmp_path / "pair.msh"
    summary = export_particle_tracks_msh([t0, t1], out, stride=10)
    # 101 samples / stride 10 -> keep 0,10,...,100 = 11 points = 10 segs
    assert summary["per_track_segments"] == [10, 1]
    assert summary["n_tracks"] == 2
    text = out.read_text()
    assert '"track_0"' in text and '"track_1"' in text


@pytest.mark.skipif(importlib.util.find_spec("gmsh") is None,
                    reason="gmsh package not installed")
def test_gmsh_round_trip_geometry_and_views(tmp_path):
    import gmsh

    voltage, bz = 100.0, 1.0e-3
    track = _circle_track(voltage=voltage, bz=bz)
    out = tmp_path / "orbit.msh"
    export_particle_tracks_msh(track, out)

    gmsh.initialize(["-noconfig"])
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(str(out))
        etypes, etags, _ = gmsh.model.mesh.getElements(dim=1)
        assert [int(t) for t in etypes] == [1]
        assert len(etags[0]) == 64
        names = []
        for tag in gmsh.view.getTags():
            idx = gmsh.view.getIndex(tag)
            names.append(gmsh.option.getString(f"View[{idx}].Name"))
        assert names == ["time [s]", "speed [m/s]",
                         "kinetic energy [eV]", "v [m/s]"]
        # the node cloud spans the analytic gyro diameter in y
        _tags, coords, _ = gmsh.model.mesh.getNodes()
        pts = np.asarray(coords, dtype=float).reshape(-1, 3)
        speed = math.sqrt(2.0 * QE * voltage / ME)
        diameter = 2.0 * ME * speed / (QE * bz)
        assert pts[:, 1].max() - pts[:, 1].min() == pytest.approx(
            diameter, rel=2e-3)
    finally:
        gmsh.finalize()


def test_rejects_non_track_inputs(tmp_path):
    out = tmp_path / "bad.msh"
    with pytest.raises(ValueError, match="track dict"):
        export_particle_tracks_msh(42, out)
    with pytest.raises(ValueError, match="at least one track"):
        export_particle_tracks_msh([], out)
    with pytest.raises(ValueError, match="schema"):
        export_particle_tracks_msh({"schema": "something-else/v1"}, out)
    track = _circle_track()
    del track["velocity_m_s"]
    with pytest.raises(ValueError, match="missing keys"):
        export_particle_tracks_msh(track, out)


def test_rejects_inconsistent_and_degenerate_tracks(tmp_path):
    out = tmp_path / "bad.msh"
    track = _circle_track()
    track["kinetic_energy_j"] = track["kinetic_energy_j"][:-1]
    with pytest.raises(ValueError, match="inconsistent sample counts"):
        export_particle_tracks_msh(track, out)

    frozen = _circle_track(n=5)
    frozen["position_m"] = [[0.0, 0.0, 0.0]] * 5
    with pytest.raises(ValueError, match="never moved"):
        export_particle_tracks_msh(frozen, out)

    short = _circle_track(n=3)
    with pytest.raises(ValueError, match="stride"):
        export_particle_tracks_msh(short, out, stride=0)
