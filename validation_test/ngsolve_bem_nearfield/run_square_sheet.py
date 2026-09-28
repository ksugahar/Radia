"""Analytic square-sheet near-field gate for direct and local-expansion Laplace SL."""
import argparse,json,math,platform
from pathlib import Path
import ngsolve as ng
from ngsolve.bem import LaplaceSL
from netgen.occ import WorkPlane,OCCGeometry,Axes,X,Z
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
ng.SetNumThreads(1)
rows=[]
with ng.TaskManager():
    source=ng.Mesh(OCCGeometry(WorkPlane().RectangleC(2,2).Face()).GenerateMesh(maxh=.4))
    fes=ng.SurfaceL2(source,order=0);g=ng.GridFunction(fes);g.Set(1,definedon=source.Boundaries('.*'))
    op=LaplaceSL(fes.TrialFunction()*ng.ds(bonus_intorder=10))
    for z in (.1,.001,1e-6):
        face=WorkPlane(Axes((0,0,z),Z,X)).RectangleC(.1,.1).Face(); face.faces.name='target'
        target=ng.Mesh(OCCGeometry(face).GenerateMesh(maxh=.04)); region=target.Boundaries('target')
        exact=ng.CF(0)
        for xx,sx in ((1-ng.x,1),(-1-ng.x,-1)):
            for yy,sy in ((1-ng.y,1),(-1-ng.y,-1)):
                rr=ng.sqrt(xx*xx+yy*yy+z*z)
                exact+=sx*sy*(xx*ng.log(yy+rr)+yy*ng.log(xx+rr)-z*ng.atan(xx*yy/(z*rr)))/(4*math.pi)
        ref=ng.Integrate(exact*exact,target,definedon=region)
        row={'distance':z}
        for name,potential in [('direct',op(g)),('local_expansion',op(g,region))]:
            row[name+'_relative_l2']=math.sqrt(ng.Integrate((potential-exact)**2,target,definedon=region)/ref)
        rows.append(row)
limits={'direct_relative_l2':1e-10,'local_expansion_relative_l2':1e-8}
passed=all(math.isfinite(row[key]) and row[key]<limit for row in rows for key,limit in limits.items())
out={'ngsolve':ng.__version__,'host':platform.node(),'bonus_intorder':10,'thresholds':limits,'passed':passed,'results':rows}
a.output.parent.mkdir(parents=True,exist_ok=True)
a.output.write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
raise SystemExit(0 if passed else 1)