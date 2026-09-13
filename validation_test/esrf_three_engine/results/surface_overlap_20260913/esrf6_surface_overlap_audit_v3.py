"""Read-only boundary overlap audit; mapped tessellation is not a proof for curves."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import time

import ngsolve as ng
import numpy as np
from scipy.spatial import cKDTree
from shapely.geometry import Polygon
import shapely


def tessellation(kind, n):
    if kind == ng.ET.TRIG:
        refs = [(i / n, j / n) for i in range(n + 1) for j in range(n + 1 - i)]
    elif kind == ng.ET.QUAD:
        refs = [(i / n, j / n) for i in range(n + 1) for j in range(n + 1)]
    else:
        raise ValueError(f'Unsupported boundary type {kind}')
    lookup = {(round(x*n), round(y*n)): k for k, (x, y) in enumerate(refs)}
    tris = []
    for (i, j), a in lookup.items():
        if (i+1,j) in lookup and (i,j+1) in lookup:
            b, c = lookup[i+1,j], lookup[i,j+1]
            tris.append((a,b,c))
            if (i+1,j+1) in lookup:
                tris.append((b,lookup[i+1,j+1],c))
    return np.array(refs), np.array(tris)


def triangle_contact(a, b, length_tol, area_tol):
    na = np.cross(a[1]-a[0], a[2]-a[0])
    nb = np.cross(b[1]-b[0], b[2]-b[0])
    aa, ab = np.linalg.norm(na), np.linalg.norm(nb)
    if min(aa,ab) <= area_tol:
        return None
    na, nb = na / aa, nb / ab
    da, db = (b-a[0]) @ na, (a-b[0]) @ nb
    if np.all(da > length_tol) or np.all(da < -length_tol):
        return None
    if np.all(db > length_tol) or np.all(db < -length_tol):
        return None
    if np.max(np.abs(da)) <= length_tol and np.max(np.abs(db)) <= length_tol:
        # Avoid a rotated almost-coincident vertex producing an invalid GEOS overlay.
        # A separating-axis check independently rejects edge/vertex-only contacts.
        axis=int(np.argmax(np.abs(na)))
        aa2,bb2=np.delete(a,axis,axis=1),np.delete(b,axis,axis=1)
        for poly in (aa2,bb2):
            for k in range(3):
                edge=poly[(k+1)%3]-poly[k]
                normal=np.array([-edge[1],edge[0]])/np.linalg.norm(edge)
                pa2,pb2=aa2@normal,bb2@normal
                if min(pa2.max(),pb2.max())-max(pa2.min(),pb2.min()) <= length_tol:
                    return None
        pa,pb=Polygon(aa2),Polygon(bb2)
        area=shapely.intersection(pa,pb,grid_size=length_tol*.01).area/abs(na[axis])
        if area > area_tol:
            return {'kind':'coplanar_positive_area', 'overlap_area_m2':float(area)}
        return None
    # Strict interior piercings exclude normal common-edge/common-vertex contacts.
    hits = []
    for source, target, normal in ((a,b,nb),(b,a,na)):
        distances = (source-target[0]) @ normal
        basis = np.column_stack((target[1]-target[0],target[2]-target[0]))
        for k in range(3):
            l = (k+1)%3
            if distances[k]*distances[l] >= -length_tol**2:
                continue
            t = distances[k]/(distances[k]-distances[l])
            p = source[k] + t*(source[l]-source[k])
            uv = np.linalg.lstsq(basis,p-target[0],rcond=None)[0]
            if min(uv[0],uv[1],1-uv.sum()) > 1e-8:
                hits.append(p.tolist())
    if hits:
        return {'kind':'transverse_interior_intersection','points_m':hits}
    return None


def controls():
    a=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
    shifted=a+np.array([.2,.1,0.])
    split=np.array([[0.,0.,0.],[.5,0.,0.],[0.,.5,0.]])
    adjacent=np.array([[1.,0.,0.],[1.,1.,0.],[0.,1.,0.]])
    cross=np.array([[.2,.2,-1.],[.2,.2,1.],[.8,.2,0.]])
    separated=a+np.array([0.,0.,.01])
    values={name:triangle_contact(a,b,1e-10,1e-16) for name,b in
            [('partial',shifted),('different_subdivision',split),('adjacent',adjacent),
             ('transverse',cross),('separated',separated)]}
    assert values['partial']['overlap_area_m2']>0
    assert values['different_subdivision']['overlap_area_m2']>0
    assert values['adjacent'] is None and values['separated'] is None
    assert values['transverse']['kind']=='transverse_interior_intersection'
    square=np.array([[0.,0.,0.],[1.,0.,0.],[1.,1.,0.],[0.,1.,0.]])
    diag_a=[square[[0,1,2]],square[[0,2,3]]]
    diag_b=[square[[0,1,3]],square[[1,2,3]]]
    overlap=sum(triangle_contact(x,y,1e-10,1e-16)['overlap_area_m2']
                for x in diag_a for y in diag_b)
    assert abs(overlap-1.) < 1e-14
    values['opposite_diagonal_tessellations_area']=overlap
    grazing=np.array([
        [[.0007499999999999998,.0979132938463439,.0010816555197760997],
         [.0014999999999999497,.09899494936612,0.],
         [.00224999999999995,.0979132938463439,.0010816555197760997]],
        [[-.0007499999999999998,.09791329384634387,.0010816555197760997],
         [0.,.09899494936612,0.],
         [.0007500000000000002,.0979132938463439,.0010816555197760997]]])
    values['reported_grazing_contact']=triangle_contact(*grazing,1e-10,1e-16)
    assert values['reported_grazing_contact'] is None
    return values


def audit(path, subdivisions):
    started=time.perf_counter()
    m=ng.Mesh(str(path))
    owners=Counter(f.nr for e in m.Elements(ng.VOL) for f in e.faces)
    boundary=list(m.Elements(ng.BND))
    boundary_faces=Counter(f.nr for e in boundary for f in e.faces)
    topology={
        'volume_face_owner_histogram':dict(Counter(owners.values())),
        'boundary_face_owner_histogram':dict(Counter(owners[f] for f in boundary_faces)),
        'repeated_boundary_face_ids':[f for f,c in boundary_faces.items() if c>1],
        'unregistered_one_owner_faces':[f for f,c in owners.items() if c==1 and f not in boundary_faces],
        'nonmanifold_volume_faces':[f for f,c in owners.items() if c>2],
    }
    points=[]; parents=[]; planarity=[]; edge_deviations=[]; ranges=[]
    xyz=ng.CF((ng.x,ng.y,ng.z))
    for e in boundary:
        refs,indices=tessellation(e.type,subdivisions)
        ir=ng.IntegrationRule(refs.tolist(),[1.]*len(refs))
        mapped=np.asarray(xyz(m.GetTrafo(e)(ir)),dtype=float)
        _,_,vh=np.linalg.svd(mapped-mapped.mean(axis=0),full_matrices=False)
        planarity.append(float(np.max(np.abs((mapped-mapped.mean(axis=0))@vh[-1]))))
        corner_refs=([(0,0),(1,0),(0,1)] if e.type==ng.ET.TRIG
                     else [(0,0),(1,0),(1,1),(0,1)])
        edge_error=0.
        for k,left in enumerate(corner_refs):
            right=corner_refs[(k+1)%len(corner_refs)]
            left,right=np.asarray(left),np.asarray(right)
            direction=right-left
            t=((refs-left)@direction)/(direction@direction)
            edge=np.linalg.norm(refs-left-t[:,None]*direction,axis=1)<1e-12
            a=mapped[np.argmin(np.linalg.norm(refs-left,axis=1))]
            b=mapped[np.argmin(np.linalg.norm(refs-right,axis=1))]
            physical=mapped[edge]-a
            edge_error=max(edge_error,float(np.max(np.linalg.norm(np.cross(physical,b-a),axis=1)/np.linalg.norm(b-a))))
        edge_deviations.append(edge_error)
        begin=len(points)
        points.extend(mapped[indices])
        parents.extend([e.nr]*len(indices))
        ranges.append((begin,len(points)))
    triangles=np.asarray(points)
    parents=np.array(parents)
    centers=triangles.mean(axis=1)
    radii=np.linalg.norm(triangles-centers[:,None,:],axis=2).max(axis=1)
    mins,maxs=triangles.min(axis=1),triangles.max(axis=1)
    tree=cKDTree(centers)
    tol=1e-10; area_tol=1e-16
    found={}; tested=0; candidates=0; degenerate=0
    # Radius bins avoid a single large Kelvin facet inflating every search.
    edges=np.unique(np.quantile(radii,np.linspace(0,1,9)))
    groups=[]
    for low,high in zip(edges[:-1],edges[1:]):
        ids=np.where((radii>=low)&(radii<=high if high==edges[-1] else radii<high))[0]
        if len(ids): groups.append((ids,cKDTree(centers[ids]),float(radii[ids].max())))
    for i,a in enumerate(triangles):
        if i%10000==0:
            print(f'{path.name} n={subdivisions} triangle {i}/{len(triangles)} contacts={len(found)}',flush=True)
        for ids,subtree,maxradius in groups:
            js=ids[subtree.query_ball_point(centers[i],radii[i]+maxradius+tol)]
            js=js[(js>i)&(parents[js]!=parents[i])]
            if not len(js): continue
            js=js[np.all(maxs[js]>=mins[i]-tol,axis=1)&np.all(mins[js]<=maxs[i]+tol,axis=1)]
            candidates+=len(js)
            for j in js:
                pair=tuple(sorted((int(parents[i]),int(parents[j]))))
                if pair in found: continue
                tested+=1
                result=triangle_contact(a,triangles[j],tol,area_tol)
                if result:
                    result.update(boundary_element_ids=pair,
                                  boundary_names=[boundary[k].mat for k in pair],
                                  parent_planarity_m=[planarity[k] for k in pair],
                                  mapped_subtriangle_ids=[int(i),int(j)],
                                  triangle_coordinates_m=[a.tolist(),triangles[j].tolist()])
                    found[pair]=result
    return {
        'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'volume_elements':m.ne,'boundary_elements':len(boundary),'curve_order':m.GetCurveOrder(),
        'subdivisions':subdivisions,'mapped_triangles':len(triangles),'topology':topology,
        'length_tolerance_m':tol,'area_tolerance_m2':area_tol,
        'coplanar_overlay_grid_m':tol*.01,
        'nonplanar_boundary_elements':sum(p>tol for p in planarity),
        'maximum_sampled_parent_planarity_m':max(planarity),
        'maximum_sampled_edge_chord_deviation_m':max(edge_deviations),
        'nonstraight_boundary_edge_parent_count':sum(v>tol for v in edge_deviations),
        'aabb_candidate_pairs':candidates,'tested_triangle_pairs':tested,
        'detected_parent_pairs':len(found),'contacts':list(found.values()),
        'elapsed_s':time.perf_counter()-started,
        'scope':'All registered boundary faces, including once-registered two-owner interfaces. '
                'Coordinates use NGSolve curved mapping. Triangle intersections are on a '
                'finite tessellation, not a certified curved-patch intersection proof. '
                'Positive-area coplanar overlaps and strict transverse piercings are reported; '
                'tangencies and very small features can remain unresolved.'}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--mesh',required=True)
    p.add_argument('--subdivisions',type=int,default=2)
    p.add_argument('--output',required=True)
    args=p.parse_args()
    ng.SetNumThreads(2)
    report={'controls':controls(),'ngsolve_version':ng.__version__,'shapely_version':shapely.__version__,
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    with ng.TaskManager():
        report['audit']=audit(Path(args.mesh),args.subdivisions)
    with Path(args.output).open('x',encoding='utf-8') as f:
        json.dump(report,f,indent=2,allow_nan=False)
    print(json.dumps({k:v for k,v in report['audit'].items() if k!='contacts'}),flush=True)


if __name__=='__main__': main()
