import asyncio
import json
import os
from pathlib import Path
import subprocess
import sys

from radia_mcp.radia_acoustic import acoustic_capabilities, acoustic_usage, cq_grid_gate, fsi_preflight_gate

def test_capabilities_are_ngsolve_first_and_not_matlab():
    c=acoustic_capabilities()
    assert c["owner"]=="radia.acoustics" and c["matlab_runtime"] is False
    assert "ngsolve.bem" in c["numerical_backends"] and "fsi" in c["apis"] and "cq" in c["apis"]

def test_usage_covers_bem_fsi_and_cq():
    assert "HelmholtzSL" in acoustic_usage("ngsolve_bem")
    assert "VectorH1" in acoustic_usage("fsi")
    assert "s=delta(zeta)/dt" in acoustic_usage("cq")

def test_fsi_preflight_rejects_bad_geometry_and_material():
    assert fsi_preflight_gate(wavenumber=2.0)["ok"]
    bad=fsi_preflight_gate(wavenumber=2.0,c_longitudinal=1.0,c_transverse=1.0,radius_deviation=.1)
    assert not bad["ok"] and not bad["checks"]["positive_lame_lambda"] and not bad["checks"]["spherical_dtn_geometry"]

def test_cq_grid_pins_convention_and_symmetry():
    for method, num_time in (("BDF1", 15), ("BDF2", 16)):
        g=cq_grid_gate(num_time=num_time,time_step=.1,method=method)
        assert g["ok"] and g["convention"]=="s=delta(zeta)/dt; kappa=i*s/c"
    assert not cq_grid_gate(num_time=2,time_step=.1)["ok"]

def test_cq_grid_import_is_minimal_dependency_safe(tmp_path):
    (tmp_path / "numpy.py").write_text(
        "raise RuntimeError('NumPy must not be imported by the acoustic gate')\n",
        encoding="utf-8",
    )
    (tmp_path / "radia.py").write_text(
        "raise RuntimeError('Radia must not be imported by the acoustic gate')\n",
        encoding="utf-8",
    )
    package_source = Path(__file__).resolve().parents[1] / "src"
    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join(
        [str(tmp_path), str(package_source), environment.get("PYTHONPATH", "")]
    )
    program = (
        "from radia_mcp.radia_acoustic import cq_grid_gate; "
        "assert cq_grid_gate(num_time=15,time_step=.1,method='BDF1')['ok']; "
        "assert cq_grid_gate(num_time=16,time_step=.1,method='BDF2')['ok']"
    )
    completed = subprocess.run(
        [sys.executable, "-c", program],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr

def test_server_registers_four_domain_tools():
    from radia_mcp.radia_acoustic.server import mcp, radia_acoustic_capabilities
    names={t.name for t in asyncio.run(mcp.list_tools())}
    assert {"radia_acoustic_usage","radia_acoustic_capabilities","radia_acoustic_fsi_preflight","radia_acoustic_cq_grid"} <= names
    assert json.loads(radia_acoustic_capabilities())["owner"]=="radia.acoustics"
