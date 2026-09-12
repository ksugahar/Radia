"""ESRF #6 factory. Do not execute the large case before candidate approval."""
import importlib.util
from pathlib import Path
import sys

import ngsolve as ng
import radia as rad
from radia.kelvin_identify_ngsolve import detect_kelvin_offset


def create_case(mesh_path):
    source = Path(__file__).resolve().parents[1] / 'esrf_three_engine/esrf_coil_yoke.py'
    spec = importlib.util.spec_from_file_location('_esrf_omega_source', source)
    model = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = model
    spec.loader.exec_module(model)
    mesh = ng.Mesh(str(mesh_path))
    rad.UtiDelAll()
    coil, _ = model.build_radia_coil_source(6)
    center = tuple(float(v) for v in detect_kelvin_offset(mesh))
    centres, samples = model.observation_volume_quadrature(6)
    import hashlib
    return {'mesh': mesh, 'H_s': rad.RadiaField(coil, 'h'),
            'observation_centres': centres, 'observation_samples': samples,
            'H_ext': rad.KelvinRadiaFieldStrength(coil, center, 0.16, (0., 0., 0.)),
            'radius': 0.16, 'center': center, 'mu_r': 1000.0,
            'controls': {'case': 'esrf6', 'mu_r': 1000.0, 'radius': 0.16,
                         'kelvin_center': center, 'physical_center': [0., 0., 0.],
                         'coil_source_file': str(source),
                         'coil_source_sha256': hashlib.sha256(source.read_bytes()).hexdigest()}}
