"""Rename labels on the public frozen P2 unit sphere without changing its maps."""
from pathlib import Path
import hashlib
import importlib.metadata
import json


def generate(source,destination,receipt):
    import numpy as np
    import ngsolve as ng
    from netgen.meshing import Mesh
    source,destination,receipt=map(Path,(source,destination,receipt))
    def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
    expected='a339966a361391ad2b64e7f57645336b3ac70fa549c7f08b2007875523d02236'
    if sha(source)!=expected or destination.exists() or source.resolve()==destination.resolve():
        raise ValueError('Frozen public source and a new output are required')
    original=Mesh(dim=3);original.Load(str(source))
    def topology(mesh):
        return dict(points=[list(p.p) for p in mesh.Points()],
                    volume=[(e.index,[v.nr for v in e.vertices]) for e in mesh.Elements3D()],
                    surface=[(e.index,[v.nr for v in e.vertices]) for e in mesh.Elements2D()])
    before=topology(original)
    if len(before['volume'])!=13 or len(before['surface'])!=18 or {e[0] for e in before['volume']}!={1} or {e[0] for e in before['surface']}!={1}:
        raise ValueError('Expected the single-region public sphere with one exterior group')
    def mapped(mesh):
        space=ng.Mesh(mesh)
        if space.GetCurveOrder()!=2:raise ValueError('Expected frozen P2 geometry')
        position=ng.CoefficientFunction((ng.x,ng.y,ng.z));jacobian=ng.specialcf.JacobianMatrix(3)
        arrays=[]
        for element in space.Elements(ng.VOL):
            trafo=space.GetTrafo(element)
            points=trafo(ng.IntegrationRule(element.type,8))
            arrays.extend((np.asarray(position(points)).copy(),np.asarray(jacobian(points)).copy()))
        for element in space.Elements(ng.BND):
            trafo=space.GetTrafo(element)
            points=trafo(ng.IntegrationRule(element.type,8))
            arrays.append(np.asarray(position(points)).copy())
        if any(not np.all(np.isfinite(a)) for a in arrays):raise ValueError('Nonfinite mapped geometry')
        return arrays
    before_maps=mapped(original)
    old_material,old_boundary=original.GetMaterial(1),original.GetBCName(0)
    original.SetMaterial(1,'magnetic_body');original.SetBCName(0,'exterior')
    original.Save(str(destination))
    loaded=Mesh(dim=3);loaded.Load(str(destination))
    after_maps=mapped(loaded)
    if topology(loaded)!=before or len(before_maps)!=len(after_maps) or any(not np.array_equal(a,b) for a,b in zip(before_maps,after_maps)):
        raise ValueError('Topology, coordinates or sampled P2 maps changed during label-only save')
    if loaded.GetMaterial(1)!='magnetic_body' or loaded.GetBCName(0)!='exterior' or sha(source)!=expected:
        raise ValueError('Semantic labels or original source changed')
    result=dict(accepted=True,scope='Label-only input preparation; no solver acceptance',source_sha256=expected,output_sha256=sha(destination),
                original_material=old_material,original_boundary=old_boundary,material='magnetic_body',boundary='exterior',length_unit='m',
                unit_source='Caller declaration; original independent CSG unit sphere',netgen_version=importlib.metadata.version('netgen-mesher'),
                curve_order=2,cells=13,surface_elements=18,geometry_rule_order=8,mapped_array_count=len(before_maps),
                topology_coordinates_unchanged=True,sampled_mapped_coordinates_and_jacobians_equal=True)
    receipt.write_text(json.dumps(result,indent=2),encoding='utf-8')
    return result


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path);parser.add_argument('destination',type=Path);parser.add_argument('receipt',type=Path)
    args=parser.parse_args();generate(args.source,args.destination,args.receipt)
