"""Pure NGSolve point/rule volume probe; no Radia native code or CUDA.

Diagnostic only. Pinned 6.2.2607/2608 rule-wise PYR assembly loses ~1e-10
relative volume versus analytic geometry; matching it is parity, not increased
geometry accuracy. No epsilon correction is applied here. No external send.
"""
import argparse
import json
import math
from pathlib import Path
import ngsolve as ng
from netgen.csg import Pnt
from netgen.meshing import Mesh, MeshPoint, Element3D


def element_mesh(family):
    coordinates = {
        "tet": [(0,0,0),(1,0,0),(0,1,0),(0,0,1)],
        "hex": [(0,0,0),(1,0,0),(1,1,0),(0,1,0),
                (0,0,1),(1,0,1),(1,1,1),(0,1,1)],
        "pyramid": [(0,0,0),(1,0,0),(1,1,0),(0,1,0),(.5,.5,1)],
    }
    mesh=Mesh(dim=3)
    mesh.SetMaterial(1,"body")
    ids=[mesh.Add(MeshPoint(Pnt(*p))) for p in coordinates[family]]
    mesh.Add(Element3D(1,ids))
    return ng.Mesh(mesh)


def element_probe(family, order=12):
    mesh=element_mesh(family)
    element_type={"tet":ng.ET.TET,"hex":ng.ET.HEX,"pyramid":ng.ET.PYRAMID}[family]
    exact={"tet":1/6,"hex":1.,"pyramid":1/3}[family]
    rule=ng.IntegrationRule(element_type,order)
    trafo=mesh.GetTrafo(ng.ElementId(ng.VOL,0))
    point_volume=math.fsum(float(ip.weight)*float(trafo(ip).measure) for ip in rule)
    # Python trafo(rule) exposes coordinates, not mapped measures. Use the
    # assembler's scalar and SIMD whole-rule mapping with this SAME rule object,
    # avoiding Integrate's rule-selection ambiguity. L2 P0 has one constant DOF;
    # verify its value is one, rather than assuming a normalized basis.
    space=ng.L2(mesh,order=0)
    one=ng.GridFunction(space,novisual=True)
    one.vec[:]=1
    assert space.ndof==1 and all(abs(one(trafo(ip))-1)<1e-14 for ip in rule)
    whole={}
    for simd in (False,True):
        form=ng.LinearForm(space)
        form+=ng.SymbolicLFI(space.TestFunction(),intrule=rule,simd_evaluate=simd)
        form.Assemble()
        whole["simd" if simd else "scalar"]=float(form.vec[0])
    explicit_dx=float(ng.Integrate(ng.CF(1)*ng.dx(intrules={element_type:rule}),mesh))
    integrate=float(ng.Integrate(1,mesh,order=order))
    weights=math.fsum(float(ip.weight) for ip in rule)
    values=dict(rule_weight_sum=weights,point_mapped_volume=point_volume,
                whole_rule_scalar_volume=whole["scalar"],whole_rule_simd_volume=whole["simd"],
                explicit_rule_dx_volume=explicit_dx,integrate_volume=integrate)
    return dict(family=family,analytic_volume=exact,order=order,rule_points=len(rule),
        **values,signed_relative_errors_vs_analytic={key:(value-exact)/exact for key,value in values.items()},
        signed_relative_errors_vs_rule_weights={key:(value-weights)/weights for key,value in values.items()},
        whole_rule_method="SymbolicLFI P0 constant-one, explicit same intrule, scalar and SIMD",
        point_vs_whole_scalar_relative=(point_volume-whole["scalar"])/exact)


def probe():
    with ng.TaskManager():
        rows=[element_probe(family) for family in ("tet","hex","pyramid")]
    return dict(schema="ngsolve.mapping-probe.v2",ngsolve=ng.__version__,rows=rows,
        scope="Unit straight tet/hex/PYR, explicit order12 shared rule. Diagnostic only; parity is not geometric exactness.")


def main():
    parser=argparse.ArgumentParser(__doc__)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    ng.SetNumThreads(1)
    result=probe()
    with args.output.open("x",encoding="utf-8") as stream:
        json.dump(result,stream,indent=2)
        stream.write("\n")


if __name__=="__main__":
    main()
