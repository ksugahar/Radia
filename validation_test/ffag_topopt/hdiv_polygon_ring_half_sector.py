"""The polygonal ring yoke of `mixed_omega_star_kelvin_ring.py` on HDiv-MMM with the midplane mirror and the rotation group.

Refinement splits each cell of the 16-gon along its chord, its radius and its height, so the geometry stays the same polygon
at every level. Reports the mean axial magnetization of the full domain. No limit is set here: acceptance belongs to the
user.
"""
import argparse
import json
import math
from pathlib import Path
import time


def half_sector(path, stations, fold, chord_split, radial_split, height_split, inner=0.03, outer=0.05, half_height=0.004):
    """One sector of z >= 0: `stations` per ring, each chord split `chord_split` times, with identified seam faces."""
    import numpy as np
    from netgen.meshing import Element2D, Element3D, FaceDescriptor, IdentificationType, Mesh, MeshPoint, Pnt
    if stations % fold:
        raise ValueError(f"stations={stations} must be a multiple of fold={fold}")
    mesh = Mesh(dim=3)
    mesh.SetMaterial(1, "yoke")
    skin = mesh.Add(FaceDescriptor(surfnr=1, domin=1, domout=0, bc=1))
    low = mesh.Add(FaceDescriptor(surfnr=2, domin=1, domout=0, bc=2))
    high = mesh.Add(FaceDescriptor(surfnr=3, domin=1, domout=0, bc=3))
    mesh.SetBCName(0, "skin")
    mesh.SetBCName(1, "sector_low")
    mesh.SetBCName(2, "sector_high")
    chords = stations//fold
    total = chords*chord_split
    point = {}
    for it, ir, iz in np.ndindex(total + 1, radial_split + 1, height_split + 1):
        chord, fraction = divmod(it, chord_split)
        a, b = 2.0*math.pi*chord/stations, 2.0*math.pi*(chord + 1)/stations
        u, v = 1.0 - fraction/chord_split, fraction/chord_split
        radius = inner + (outer - inner)*ir/radial_split
        point[it, ir, iz] = mesh.Add(MeshPoint(Pnt(radius*(u*math.cos(a) + v*math.cos(b)), radius*(u*math.sin(a) + v*math.sin(b)), half_height*iz/height_split)))
    for it, ir, iz in np.ndindex(total, radial_split, height_split):
        corners = [point[it, ir, iz], point[it, ir + 1, iz], point[it + 1, ir + 1, iz], point[it + 1, ir, iz],
                   point[it, ir, iz + 1], point[it, ir + 1, iz + 1], point[it + 1, ir + 1, iz + 1], point[it + 1, ir, iz + 1]]
        mesh.Add(Element3D(1, corners))
        if iz == 0:
            mesh.Add(Element2D(skin, [corners[k] for k in (0, 3, 2, 1)]))
        if iz == height_split - 1:
            mesh.Add(Element2D(skin, [corners[k] for k in (4, 5, 6, 7)]))
        if ir == 0:
            mesh.Add(Element2D(skin, [corners[k] for k in (0, 4, 7, 3)]))
        if ir == radial_split - 1:
            mesh.Add(Element2D(skin, [corners[k] for k in (1, 2, 6, 5)]))
        if it == 0:
            mesh.Add(Element2D(low, [corners[k] for k in (0, 1, 5, 4)]))
        if it == total - 1:
            mesh.Add(Element2D(high, [corners[k] for k in (2, 3, 7, 6)]))
    for ir, iz in np.ndindex(radial_split + 1, height_split + 1):
        mesh.AddPointIdentification(point[0, ir, iz], point[total, ir, iz], identnr=1, type=IdentificationType.PERIODIC)
    mesh.Save(str(path))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--stations', type=int, default=16)
    parser.add_argument('--fold', type=int, default=8)
    parser.add_argument('--chord', type=int, default=1)
    parser.add_argument('--radial', type=int, default=1)
    parser.add_argument('--height', type=int, default=1)
    parser.add_argument('--order', type=int, choices=(1, 2), default=2)
    args = parser.parse_args()
    import ngsolve as ng
    from radia import vim
    mu_r, applied = 10.0, 1e5
    start = time.perf_counter()
    record = dict(schema='hdiv.polygon-ring-half-sector.v1', passed=False, stations=args.stations, fold=args.fold, chord=args.chord, radial=args.radial, height=args.height,
                  order=args.order, mu_r=mu_r, applied=applied, ngsolve=ng.__version__)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.work.mkdir(parents=True, exist_ok=True)
    try:
        path = args.work/f'half-sector-{args.stations}-{args.chord}-{args.radial}-{args.height}.vol'
        half_sector(path, args.stations, args.fold, args.chord, args.radial, args.height)
        mesh = ng.Mesh(str(path))
        with ng.TaskManager():
            result = vim.Solve(mesh, mu_r=mu_r, H_ext=ng.CF((0, 0, applied)), order=args.order, tol=1e-10, image='-z', image_cyclic=args.fold,
                               cyclic_periodic_boundaries=('sector_low', 'sector_high'))
        record.update(elements=mesh.ne, volume=ng.Integrate(ng.CF(1.0), mesh), ndof=int(result['ndof']), mean_m=[float(v) for v in result['M_avg']],
                      mean_mz=float(result['M_avg'][2]), elapsed_seconds=time.perf_counter() - start)
        record['passed'] = True
    except Exception as error:
        record['failure'] = repr(error)
        raise
    finally:
        args.out.write_text(json.dumps(record, indent=2, default=str) + '\n')
    print(json.dumps({key: record.get(key) for key in ('chord', 'radial', 'height', 'order', 'elements', 'ndof', 'volume', 'mean_mz', 'elapsed_seconds')}))


if __name__ == '__main__':
    main()
