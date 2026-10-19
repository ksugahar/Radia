"""Independent NumPy fixture for the real thermal shadow-vector regression."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np


def fixture():
    a = np.array([[4+1j, .2-.1j, 0], [.1+.2j, 5-.5j, .3], [.2, 0, 6+.3j]])
    prescribed = np.array([1+.2j, -.3+.4j, .7-.1j])
    rhs = a @ prescribed
    electromagnetic = np.linalg.solve(a, rhs)
    p = np.array([[1, .2, 0], [0, 1, .1], [.3, 0, 1],
                  [1, 0, .4], [.2, 1, 0], [0, .3, 1]])
    k = np.array([[2., -1, 0], [-1, 2, -1], [0, -1, 2]])
    c = .1*np.eye(3)
    projection = np.array([[1., 0], [.5, .5], [0, 1]])
    system = np.eye(3)+.01*k+.002*np.eye(3)
    heat = (.5*np.abs(p @ electromagnetic)**2).reshape(2, 3).sum(axis=1)
    temperature = [np.array([300., 301., 302.])]
    for _ in range(5):
        old = temperature[-1]
        load = .01*(projection @ heat-k @ old+2*c @ (300*np.ones(3)-old))
        temperature.append(old+np.linalg.solve(system, load))
    return dict(schema='radia.ih-real-thermal-fixture.v1', numpy=np.__version__,
                matrix_real=a.real.tolist(), matrix_imag=a.imag.tolist(),
                rhs_real=rhs.real.tolist(), rhs_imag=rhs.imag.tolist(),
                expected_electromagnetic_real=electromagnetic.real.tolist(),
                expected_electromagnetic_imag=electromagnetic.imag.tolist(),
                thermal_system=system.tolist(), stiffness=k.tolist(),
                convection=c.tolist(), heat_projection=projection.tolist(),
                temperature=np.array(temperature).tolist(), heat=heat.tolist(),
                temperature_relative_limit=1e-8, heat_relative_limit=1e-9,
                true_residual_limit=1e-10)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    result = fixture()
    result['generator_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
