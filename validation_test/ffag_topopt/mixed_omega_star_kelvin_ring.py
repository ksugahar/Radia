"""A polygonal ring yoke in a uniform axial field on the mixed Omega route with a star-shaped Kelvin exterior.

A second analysis route for the HDiv-MMM ring: the same polygon, material and field, reporting the mean axial
magnetization. ``--reduction`` selects the full model, the ``z >= 0`` half (the potential vanishes on the midplane because
it is odd in z), one rotational sector, or the half sector, so the reductions can be compared on one route before the two
routes are compared. No limit is set here: acceptance belongs to the user.
"""
import argparse
import json
import math
from pathlib import Path
import sys
import time
import types


def polygon_prism(radius, stations, low, high):
    from netgen.occ import Face, Pnt, Segment, Wire, Vec
    corners = [Pnt(radius*math.cos(2.0*math.pi*k/stations), radius*math.sin(2.0*math.pi*k/stations), low) for k in range(stations)]
    wire = Wire([Segment(corners[k], corners[(k + 1) % stations]) for k in range(stations)])
    return Face(wire).Extrude(Vec(0, 0, high - low))


def disk_body(flat, edge):
    from netgen.occ import ArcOfCircle, Axis, Face, Pnt, Segment, Wire, Z
    bottom, rim_low, rim, rim_high, top = Pnt(0, 0, -edge), Pnt(flat, 0, -edge), Pnt(flat + edge, 0, 0), Pnt(flat, 0, edge), Pnt(0, 0, edge)
    section = Face(Wire([Segment(bottom, rim_low), ArcOfCircle(rim_low, rim, rim_high), Segment(rim_high, top), Segment(top, bottom)]))
    return section.Revolve(Axis(Pnt(0, 0, 0), Z), 360)


def wedge(fold, reach, half_height, half_only):
    """The prism over the sector 0 <= angle <= 2 pi/fold, and z >= 0 when `half_only`."""
    from netgen.occ import Face, Pnt, Segment, Wire, Vec
    angle = 2.0*math.pi/fold
    low = 0.0 if half_only else -half_height
    a, b, c = Pnt(0, 0, low), Pnt(reach, 0, low), Pnt(reach*math.cos(angle), reach*math.sin(angle), low)
    return Face(Wire([Segment(a, b), Segment(b, c), Segment(c, a)])).Extrude(Vec(0, 0, half_height - low))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--reduction', choices=('full', 'half', 'sector', 'half-sector'), default='full')
    parser.add_argument('--fold', type=int, default=8)
    parser.add_argument('--stations', type=int, default=32)
    parser.add_argument('--maxh', type=float, default=.03)
    parser.add_argument('--iron-maxh', type=float, default=.006)
    parser.add_argument('--center-maxh', type=float, default=.005)
    parser.add_argument('--order', type=int, choices=(1, 2, 3), default=2)
    parser.add_argument('--flat', type=float, default=.06)
    parser.add_argument('--edge', type=float, default=.03)
    args = parser.parse_args()
    import numpy as np
    import ngsolve as ng
    from netgen.meshing import IdentificationType, MeshingParameters
    from netgen.occ import Axis, Glue, OCCGeometry, Pnt, Vec, Z, gp_Trsf
    root = Path(__file__).resolve().parents[2]/'src'/'radia'
    namespace = types.ModuleType('radia')
    namespace.__path__ = [str(root)]
    sys.modules['radia'] = namespace
    from radia.kelvin_material import MU_0
    from radia.kelvin_solver import solve_magnetostatic_mixed_total_reduced_omega_kelvin as solve
    from radia.kelvin_star import RoundedDiskBoundary, star_kelvin_mu_cf
    mu_r, applied, inner_radius, outer_radius, half_height = 10.0, 1e5, 0.03, 0.05, 0.004
    half_only = args.reduction in ('half', 'half-sector')
    sector = args.reduction in ('sector', 'half-sector')
    boundary = RoundedDiskBoundary(args.flat, args.edge)
    shift = 3.0*(args.flat + args.edge)
    start = time.perf_counter()
    record = dict(schema='mixed-omega.star-kelvin-ring.v1', passed=False, reduction=args.reduction, fold=args.fold, stations=args.stations, maxh=args.maxh,
                  iron_maxh=args.iron_maxh, center_maxh=args.center_maxh, order=args.order, flat=args.flat, edge=args.edge, mu_r=mu_r, applied=applied,
                  inner_radius=inner_radius, outer_radius=outer_radius, half_height=half_height, ngsolve=ng.__version__)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    try:
        iron = polygon_prism(outer_radius, args.stations, -half_height, half_height) - polygon_prism(inner_radius, args.stations, -2*half_height, 2*half_height)
        iron.faces.name = 'interface'
        inner = disk_body(args.flat, args.edge)
        inner.faces.name = 'kelvin_int'
        exterior = disk_body(args.flat, args.edge)
        exterior.faces.name = 'kelvin_ext'
        reach = 4.0*(args.flat + args.edge)
        if sector or half_only:
            cut = wedge(args.fold, reach, 2.0*args.edge, half_only) if sector else None
            if cut is None:
                from netgen.occ import Box
                cut = Box(Pnt(-reach, -reach, 0.0), Pnt(reach, reach, 2.0*args.edge))
            iron, inner, exterior = iron*cut, inner*cut, exterior*cut
        air = inner - iron
        iron.mat('iron')
        iron.maxh = args.iron_maxh
        air.mat('air')
        exterior = exterior.Move((shift, 0.0, 0.0))
        exterior.mat('kelvin')
        for solid, kind in ((air, 'air'), (iron, 'iron'), (exterior, 'ext')):
            for face in solid.faces:
                c = face.center
                x = c.x - (shift if kind == 'ext' else 0.0)
                if sector and abs(c.y) < 1e-9:
                    face.name = 'sector_low'
                elif sector and abs(math.sin(2.0*math.pi/args.fold)*x - math.cos(2.0*math.pi/args.fold)*c.y) < 1e-9:
                    face.name = 'sector_high'
                elif half_only and abs(c.z) < 1e-9:
                    face.name = 'midplane'
        gauge = min(exterior.vertices, key=lambda vertex: abs(vertex.p[2]))
        gauge.name = 'GND'
        record['gauge_point'] = [gauge.p[0], gauge.p[1], gauge.p[2]]
        geo = Glue([iron, air, exterior])
        inner_faces = [face for solid in geo.solids for face in solid.faces if face.name == 'kelvin_int']
        outer_faces = [face for solid in geo.solids for face in solid.faces if face.name == 'kelvin_ext']
        record['kelvin_faces'] = (len(inner_faces), len(outer_faces))
        for face in outer_faces:
            c = face.center
            partner = min(inner_faces, key=lambda other: (other.center - Pnt(c.x - shift, c.y, c.z)).Norm())
            partner.Identify(face, 'kelvin', IdentificationType.PERIODIC, gp_Trsf.Translation(Vec(shift, 0, 0)))
        if sector:
            angle = 360.0/args.fold
            faces = {name: [face for solid in geo.solids for face in solid.faces if face.name == name] for name in ('sector_low', 'sector_high')}
            record['sector_faces'] = {name: len(value) for name, value in faces.items()}
            for face in faces['sector_low']:
                x = face.center.x
                axis_origin = shift if x > 0.5*shift + 0.5*(args.flat + args.edge) else 0.0
                rotation = gp_Trsf.Rotation(Axis(Pnt(axis_origin, 0, 0), Z), math.radians(angle))
                rotated = Pnt(axis_origin + (x - axis_origin)*math.cos(math.radians(angle)) - face.center.y*math.sin(math.radians(angle)),
                              (x - axis_origin)*math.sin(math.radians(angle)) + face.center.y*math.cos(math.radians(angle)), face.center.z)
                partner = min(faces['sector_high'], key=lambda other: (other.center - rotated).Norm())
                face.Identify(partner, 'sector', IdentificationType.PERIODIC, rotation)
        with ng.TaskManager():
            parameters = MeshingParameters(maxh=args.maxh, grading=.3)
            parameters.RestrictH(shift, 0.0, 0.0, args.center_maxh)
            occ = OCCGeometry(geo)
            mesh = ng.Mesh(occ.GenerateMesh(mp=parameters))
            mesh.Curve(2)
            space = ng.H1(mesh, order=args.order)
            source = ng.GridFunction(space)
            source.Set(-applied*ng.z)
            zero = ng.GridFunction(space)
            mu = star_kelvin_mu_cf(mesh, boundary, (shift, 0.0, 0.0), mu_r_by_material={'iron': mu_r})
            options = dict(reduced_dirichlet_boundary='midplane', total_dirichlet_boundary='midplane') if half_only else {}
            result = solve(mesh, ng.CF((0, 0, applied)), source, args.flat + args.edge, (shift, 0.0, 0.0), mu_cf=mu, reduced_materials=('air',),
                           total_materials=('iron', 'kelvin'), interface_boundary='interface', kelvin_interface_boundary='kelvin_int',
                           reduced_normal_flux=MU_0*applied*ng.specialcf.normal(3)[2], reduced_flux_boundary='kelvin_int',
                           kelvin_source_potential=zero, order=args.order, **options)
            volume = ng.Integrate(ng.CF(1.0), mesh, definedon=mesh.Materials('iron'))
            mean_h = [ng.Integrate(result['H_cf'][k], mesh, definedon=mesh.Materials('iron'), order=2*args.order + 2)/volume for k in range(3)]
        record.update(elements=mesh.ne, iron_volume=volume, mean_h=mean_h, mean_mz=(mu_r - 1.0)*mean_h[2], boundaries=sorted(set(mesh.GetBoundaries())),
                      residual=result.get('linear_residual'), elapsed_seconds=time.perf_counter() - start)
        record['passed'] = True
    except Exception as error:
        record['failure'] = repr(error)
        raise
    finally:
        args.out.write_text(json.dumps(record, indent=2, default=str) + '\n')
    print(json.dumps({key: record.get(key) for key in ('reduction', 'order', 'maxh', 'elements', 'iron_volume', 'mean_mz', 'elapsed_seconds')}))


if __name__ == '__main__':
    main()
