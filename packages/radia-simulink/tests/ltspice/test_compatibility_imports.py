"""The independent circuit implementation must not import the CPU solver."""
import os
from pathlib import Path
import subprocess
import sys


def test_ltspice_import_is_independent_of_cpu_solver():
    source = Path(__file__).resolve().parents[2] / "src"
    code = "import sys; from radia_simulink.ltspice import topology_signature; assert topology_signature('R1 a 0 1k\\n.end\\n'); assert 'radia' not in sys.modules"
    env = dict(os.environ, PYTHONPATH=str(source))
    subprocess.run([sys.executable, "-c", code], env=env, check=True, capture_output=True, text=True)
