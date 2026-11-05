"""Uniform applied field / permeable sphere with a star-shaped Kelvin exterior, on the mixed Omega route.

The analytic solution of a permeable sphere is independent of the truncation shape, so the same exact field judges a
spherical and a rounded-disk Kelvin boundary. The unbounded applied field stays analytic and only the decaying reaction
enters the Kelvin ball. Reports the probes inside the physical region and the physical reaction recovered from the
exterior through the pullback.
"""
import argparse
import json
from pathlib import Path
import sys
import time
import types


def exact_h(point, mu_r, radius, applied=1.0):
    import numpy as np
    r = np.asarray(point, dtype=float)
    distance = np.linalg.norm(r)
    if distance < radius:
        return np.array([0.0, 0.0, 3.0/(mu_r + 2.0)])*applied
    return exact_reaction(point, mu_r, radius, applied) + np.array([0.0, 0.0, applied])


def exact_reaction(point, mu_r, radius, applied=1.0):
    import numpy as np
    r = np.asarray(point, dtype=float)
    distance = np.linalg.norm(r)
    alpha = (mu_r - 1.0)/(mu_r + 2.0)
    return applied*alpha*radius**3*(3.0*r*r[2]/distance**5 - np.array([0.0, 0.0, 1.0])/distance**3)


def body(boundary):
    """The solid inside the Kelvin boundary, about the origin."""
    from netgen.occ import Axis, ArcOfCircle, Face, Pnt, Segment, Sphere, Wire, Z
    from radia.kelvin_star import RoundedDiskBoundary
    if not isinstance(boundary, RoundedDiskBoundary):
        return Sphere(Pnt(0, 0, 0), boundary.radius)
    flat, edge = boundary.flat_radius, boundary.edge_radius
    bottom, rim_low, rim, rim_high, top = Pnt(0, 0, -edge), Pnt(flat, 0, -edge), Pnt(flat + edge, 0, 0), Pnt(flat, 0, edge), Pnt(0, 0, edge)
    section = Face(Wire([Segment(bottom, rim_low), ArcOfCircle(rim_low, rim, rim_high), Segment(rim_high, top), Segment(top, bottom)]))
    return section.Revolve(Axis(Pnt(0, 0, 0), Z), 360)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--maxh', type=float, default=.06)
    parser.add_argument('--iron-maxh', type=float, default=.02)
    parser.add_argument('--center-maxh', type=float, default=.01)
    parser.add_argument('--order', type=int, choices=(1, 2, 3), default=2)
    parser.add_argument('--geometry-order', type=int, choices=(1, 2, 3), default=2)
    parser.add_argument('--shape', choices=('sphere', 'disk', 'thin-disk'), default='disk')
    args = parser.parse_args()
    import numpy as np
    import ngsolve as ng
    from netgen.meshing import IdentificationType, MeshingParameters
    from netgen.occ import Glue, OCCGeometry, Pnt, Sphere, Vec, gp_Trsf
    root = Path(__file__).resolve().parents[2]/'src'/'radia'
    namespace = types.ModuleType('radia')
    namespace.__path__ = [str(root)]
    sys.modules['radia'] = namespace
    from radia.kelvin_material import MU_0
    from radia.kelvin_solver import solve_magnetostatic_mixed_total_reduced_omega_kelvin as solve
    from radia.kelvin_star import RoundedDiskBoundary, SphereBoundary, star_kelvin_map, star_kelvin_mu_cf, star_kelvin_physical_h
    boundary = {'sphere': SphereBoundary(.3), 'disk': RoundedDiskBoundary(.2, .1), 'thin-disk': RoundedDiskBoundary(.25, .05)}[args.shape]
    mu_r, iron_radius, center = 100.0, 0.07, (0.9, 0.0, 0.0)
    start = time.perf_counter()
    record = dict(schema='mixed-omega.star-kelvin-sphere.v1', passed=False, shape=args.shape, mu_r=mu_r, iron_radius=iron_radius, maxh=args.maxh, center_maxh=args.center_maxh,
                  order=args.order, geometry_order=args.geometry_order, ngsolve=ng.__version__)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    try:
        iron = Sphere(Pnt(0, 0, 0), iron_radius)
        iron.mat('iron')
        iron.faces.name = 'interface'
        iron.maxh = args.iron_maxh
        inner = body(boundary)
        inner.faces.name = 'kelvin_int'
        air = inner - iron
        air.mat('air')
        exterior = body(boundary).Move(center)
        exterior.mat('kelvin')
        exterior.faces.name = 'kelvin_ext'
        for vertex in exterior.vertices:
            vertex.name = 'GND'
            break
        geo = Glue([iron, air, exterior])
        inner_faces = [face for solid in geo.solids for face in solid.faces if face.name == 'kelvin_int']
        outer_faces = [face for solid in geo.solids for face in solid.faces if face.name == 'kelvin_ext']
        record['identified_faces'] = len(inner_faces)
        for face in outer_faces:
            c = face.center
            partner = min(inner_faces, key=lambda other: (other.center - Pnt(c.x - center[0], c.y, c.z)).Norm())
            partner.Identify(face, 'kelvin', IdentificationType.PERIODIC, gp_Trsf.Translation(Vec(*center)))
        with ng.TaskManager():
            # A flat boundary maps distant points close to the centre: |T x| = rho^2/|x|.
            parameters = MeshingParameters(maxh=args.maxh, grading=.3)
            parameters.RestrictH(center[0], center[1], center[2], args.center_maxh)
            mesh = ng.Mesh(OCCGeometry(geo).GenerateMesh(mp=parameters))
            mesh.Curve(args.geometry_order)
            space = ng.H1(mesh, order=args.order)
            source = ng.GridFunction(space)
            source.Set(-ng.z)
            zero = ng.GridFunction(space)
            mu = star_kelvin_mu_cf(mesh, boundary, center, mu_r_by_material={'iron': mu_r})
            result = solve(mesh, ng.CF((0, 0, 1)), source, 1.0, center, mu_cf=mu, reduced_materials=('air',), total_materials=('iron', 'kelvin'),
                           interface_boundary='interface', kelvin_interface_boundary='kelvin_int', reduced_normal_flux=MU_0*ng.specialcf.normal(3)[2],
                           reduced_flux_boundary='kelvin_int', kelvin_source_potential=zero, order=args.order)
        shell = np.array([(0.25, 0.0, 0.02), (0.15, 0.15, 0.03), (0.0, 0.1, 0.08), (0.05, 0.2, 0.0)])
        shell_actual = np.array([np.array(result['H_cf'](mesh(*point))) for point in shell]) - np.array([0.0, 0.0, 1.0])
        shell_exact = np.array([exact_reaction(point, mu_r, iron_radius) for point in shell])
        points = [(0.0, 0.0, 0.0), (0.02, 0.0, 0.01), (0.03, 0.03, 0.03), (0.05, 0.05, 0.05), (0.12, 0.03, 0.02), (0.02, 0.15, 0.02)]
        actual = np.array([np.array(result['H_cf'](mesh(*point))) for point in points])
        exact = np.array([exact_h(point, mu_r, iron_radius) for point in points])
        far = np.array([(0.6, 0.0, 0.1), (0.2, 0.5, 0.2), (0.0, 0.1, 0.45), (1.2, 0.5, 0.8)])
        computational = star_kelvin_map(far, boundary, (0.0, 0.0, 0.0)) + np.asarray(center)
        comp_h = np.array([np.array(result['H_cf'](mesh(*point))) for point in computational])
        recovered = star_kelvin_physical_h(comp_h, computational, boundary, center)
        reaction = np.array([exact_reaction(point, mu_r, iron_radius) for point in far])
        record.update(elements=mesh.ne, dofs=int(result.get('ndof', 0)) if isinstance(result.get('ndof', 0), int) else None, points=points, actual_H=actual.tolist(), exact_H=exact.tolist(),
                      relative_rms=float(np.linalg.norm(actual - exact)/np.linalg.norm(exact)),
                      reaction_points=far.tolist(), recovered_reaction=recovered.tolist(), exact_reaction=reaction.tolist(),
                      reaction_relative_rms=float(np.linalg.norm(recovered - reaction)/np.linalg.norm(reaction)),
                      shell_points=shell.tolist(), shell_reaction=shell_actual.tolist(), shell_exact=shell_exact.tolist(),
                      shell_relative_rms=float(np.linalg.norm(shell_actual - shell_exact)/np.linalg.norm(shell_exact)),
                      residual=result.get('linear_residual'), elapsed_seconds=time.perf_counter() - start)
        record['passed'] = True
    except Exception as error:
        record['failure'] = repr(error)
        raise
    finally:
        args.out.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps({key: record.get(key) for key in ('shape', 'order', 'maxh', 'elements', 'relative_rms', 'shell_relative_rms', 'reaction_relative_rms', 'elapsed_seconds')}))


if __name__ == '__main__':
    main()
