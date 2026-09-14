"""Axisymmetric source attribution; diagnostic changes, not proposed physics fixes."""
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import ngsolve as ng
from netgen.geom2d import SplineGeometry
import argparse
import csv
import platform
import sys
import radia
import radia.panels.calc_heat_axisym as calc_heat_axisym

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, required=True, help="Received TKE08 verification directory")
parser.add_argument("--out", type=Path, required=True)
args_cli = parser.parse_args()
data = args_cli.input / "data"

def contour_points():
    with (data / "meridian.csv").open(encoding="utf-8-sig") as stream:
        raw = [(float(row["r_mm"])*1e-3, float(row["z_mm"])*1e-3)
               for row in csv.DictReader(stream)]
    if len(raw) != 300 or raw[12] != (.0085, .065) or raw[-1] != raw[12]:
        raise ValueError("Unexpected received contour layout")
    return raw[12:-1]

base = SimpleNamespace(
    contour_points=contour_points,
    VOL3=data / "work_occ_h0.0018.vol",
    QFILES={"A": data / "qsurf_caseA_9048A.sol", "B": data / "qsurf_caseB_8293A.sol"},
    calc_heat_axisym=calc_heat_axisym,
)

OUT = args_cli.out


def mesh_split():
    pts = base.contour_points()
    geo = SplineGeometry()
    ids = [geo.AddPoint(*p) for p in pts]
    for i, p in enumerate(pts):
        end = (i + 1) % len(pts)
        inner = abs(p[0] - .0075) < 1e-9 and abs(pts[end][0] - .0075) < 1e-9
        geo.Append(['line', ids[i], ids[end]], leftdomain=1, rightdomain=0,
                   bc='inner' if inner else 'outer')
    geo.SetMaterial(1, 'default')
    return ng.Mesh(geo.GenerateMesh(maxh=.0018))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    payload = {'method': 'axisymmetric standard H1 p2, 2*pi*r, dt=0.1s, 5s',
               'runtime': {'host': platform.node(), 'python': sys.version,
                           'radia_version': radia.__version__, 'ngsolve_version': ng.__version__,
                           'solver_file': str(Path(calc_heat_axisym.__file__).resolve()),
                           'solver_sha256': hashlib.sha256(Path(calc_heat_axisym.__file__).read_bytes()).hexdigest()},
               'scope': 'source attribution only; deleting inner loss is not a validated correction',
               'input_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in [base.VOL3, *base.QFILES.values()]}, 'cases': {}}
    ng.SetNumThreads(4)
    with ng.TaskManager():
        for case in ('A', 'B'):
            mesh = mesh_split()
            fes = ng.H1(mesh, order=2)
            u, v = fes.TnT()
            w = 2*np.pi*ng.x
            volume = ng.Integrate(w, mesh)
            args = SimpleNamespace(q_uniform=None, qsurf_sol=str(base.QFILES[case]),
                                   em_vol=str(base.VOL3), qsurf_order=1, n_phi_samples=128)
            qgf, q, transfer = base.calc_heat_axisym._build_axisym_qsurf_gf(
                mesh, {'inner', 'outer'}, args)
            qgf.Save(str(OUT / f'{case}_q.sol'))
            power = {name: float(ng.Integrate(q*w, mesh, ng.BND,
                         definedon=mesh.Boundaries(name), order=8))
                     for name in ('inner', 'outer')}
            area = float(ng.Integrate(w, mesh, ng.BND, definedon=mesh.Boundaries('inner')))
            a = ng.BilinearForm(fes, symmetric=True)
            a += 46.6*ng.grad(u)*ng.grad(v)*w*ng.dx + 10*u*v*w*ng.ds
            a.Assemble()
            m = ng.BilinearForm(fes, symmetric=True)
            m += 7800*467*u*v*w*ng.dx
            m.Assemble()
            star = m.mat.CreateMatrix()
            star.AsVector().data = m.mat.AsVector() + .1*a.mat.AsVector()
            inv = star.Inverse(fes.FreeDofs(), inverse='sparsecholesky')
            fields, rows = {}, {}
            for role, selector, scale in [('all', 'inner|outer', 1),
                                         ('outer', 'outer', 1), ('inner', 'inner', 1),
                                         ('outer_equal_power', 'outer', sum(power.values())/power['outer'])]:
                f = ng.LinearForm(fes)
                f += scale*q*v*w*ng.ds(selector) + 10*25*v*w*ng.ds
                f.Assemble()
                g = ng.GridFunction(fes)
                g.Set(ng.CF(25))
                residual = g.vec.CreateVector()
                conv_energy = 0.
                for _ in range(50):
                    residual.data = f.vec - a.mat*g.vec
                    g.vec.data += .1*(inv*residual)
                    conv_energy += .1*float(ng.Integrate(10*(g-25)*w, mesh, ng.BND))
                input_power = scale*sum(power[n] for n in selector.split('|'))
                stored = float(ng.Integrate(7800*467*(g-25)*w, mesh))
                closure = (stored + conv_energy - 5*input_power)/(5*input_power)
                assert abs(closure) < 1e-9, closure
                rows[role] = {'power_W': input_power,
                    'Tmean_C': float(ng.Integrate(g*w, mesh))/volume,
                    'Tinner_mean_C': float(ng.Integrate(g*w, mesh, ng.BND,
                                        definedon=mesh.Boundaries('inner')))/area,
                    'energy_relative_residual': closure}
                fields[role] = g
                g.Save(str(OUT / f'{case}_{role}_T.sol'))
                print(case, role, rows[role], flush=True)
            rel = float(ng.sqrt(ng.Integrate((fields['all']-fields['outer']-fields['inner']+25)**2*w, mesh)
                       / ng.Integrate((fields['all']-25)**2*w, mesh)))
            assert rel < 1e-10, rel
            payload['cases'][case] = {'power_W': power,
                'inner_power_fraction': power['inner']/sum(power.values()),
                'transfer': transfer, 'runs': rows, 'superposition_relative_L2': rel}
    payload['status'] = 'pass'
    (OUT/'result.json').write_text(json.dumps(payload, indent=2), encoding='utf-8')
    print('PASS', flush=True)


if __name__ == '__main__':
    main()
