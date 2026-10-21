"""NGSolve-only mixed HDiv constant projection diagnostic, not MMM acceptance.

Uses a checked exported mesh and sparsecholesky, with no Radia native or GPU.
A small assembled residual need not imply accurate pointwise reproduction.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import ngsolve as ng


def probe(mesh_path, contract_path):
    checked=subprocess.run([sys.executable,"-m","cubit_mesh_export.check",str(mesh_path),
        "--strict-labels","--contract",str(contract_path),"--format","json"],capture_output=True,text=True)
    if checked.returncode:
        raise RuntimeError("Mesh check failed: "+checked.stderr+checked.stdout)
    assert json.loads(checked.stdout)["passed"]
    mesh=ng.Mesh(str(mesh_path))
    assert {len(el.vertices) for el in mesh.Elements(ng.VOL)}=={4,5,8}
    rows=[]
    with ng.TaskManager():
        for order in (1,2):
            space=ng.HDiv(mesh,order=order)
            assert space.ndof<20000
            u,v=space.TnT()
            mass=ng.BilinearForm(space,symmetric=True)
            mass+=u*v*ng.dx(bonus_intorder=12)
            mass.Assemble()
            constant=ng.CF((.3,-.7,1.1))
            rhs=ng.LinearForm(space)
            rhs+=constant*v*ng.dx(bonus_intorder=12)
            rhs.Assemble()
            field=ng.GridFunction(space,novisual=True)
            field.vec.data=mass.mat.Inverse(space.FreeDofs(),inverse="sparsecholesky")*rhs.vec
            residual=rhs.vec.CreateVector()
            residual.data=mass.mat*field.vec-rhs.vec
            elements=[]
            for el in mesh.Elements(ng.VOL):
                vertex_count=len(el.vertices)
                et={4:ng.ET.TET,5:ng.ET.PYRAMID,8:ng.ET.HEX}[vertex_count]
                rule=ng.IntegrationRule(et,12)
                trafo=mesh.GetTrafo(ng.ElementId(ng.VOL,el.nr))
                point=np.asarray([field(trafo(ip)) for ip in rule])
                bulk=np.asarray(field(trafo(rule)))
                target=np.array([.3,-.7,1.1])
                elements.append(dict(vertex_count=vertex_count,
                    point_max_error=float(np.max(abs(point-target))),
                    bulk_max_error=float(np.max(abs(bulk-target))),
                    point_bulk_max_difference=float(np.max(abs(point-bulk)))))
            maximum=max(row["point_max_error"] for row in elements)
            rows.append(dict(order=order,dofs=space.ndof,
                assembled_mass_relative_residual=float(ng.Norm(residual)/ng.Norm(rhs.vec)),
                integrated_l2_error_squared=float(ng.Integrate(ng.InnerProduct(field-constant,field-constant),mesh,order=20)),
                elements=elements,pointwise_limit=1e-11,pointwise_accepted=maximum<=1e-11))
    return dict(schema="ngsolve.mixed-projection-probe.v1",ngsolve=ng.__version__,
        mesh_sha256=hashlib.sha256(mesh_path.read_bytes()).hexdigest(),
        contract_sha256=hashlib.sha256(contract_path.read_bytes()).hexdigest(),
        rows=rows,scope="Independent NG mass projection only; no charge-Gram, Radia native or CUDA. Diagnostic, not relaxed acceptance.")


def main():
    parser=argparse.ArgumentParser(__doc__)
    parser.add_argument("--mesh",type=Path,required=True)
    parser.add_argument("--contract",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    ng.SetNumThreads(1)
    result=probe(args.mesh,args.contract)
    with args.output.open("x",encoding="utf-8") as stream:
        json.dump(result,stream,indent=2)
        stream.write("\n")


if __name__=="__main__":
    main()
