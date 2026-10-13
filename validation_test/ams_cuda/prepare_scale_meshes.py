"""Generate reproducible public HCurl meshes and checked versioned labels."""
from pathlib import Path
import argparse,hashlib,json,subprocess,sys

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main(a):
    import ngsolve as ng
    from netgen.csg import unit_cube
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False);ng.SetNumThreads(1)
    labels=dict(schema='radia.vol-label-contract.v1',application='native-hcurl-cuda-unit-cube',strict_labels=True)
    labels['required']=dict(materials=['domain'],boundaries=['back','front','left','right','top','bottom']);labels['allowed']=labels['required']
    contract=out/'cube-labels.json';contract.write_text(json.dumps(labels,indent=2)+'\n');series=[]
    for target,h in [(10000,.12),(30000,.085),(100000,.056)]:
        trials=[]
        for attempt in range(10):
            mesh=ng.Mesh(unit_cube.GenerateMesh(maxh=h));mesh.ngmesh.SetMaterial(1,'domain')
            fes=ng.HCurl(mesh,order=1,nograds=True,dirichlet='.*');n=fes.ndof
            trials.append(dict(maxh=h,dofs=n,elements=mesh.ne))
            if abs(n-target)/target<=.08:break
            h*= (n/target)**(1/3)
        else:raise RuntimeError('Mesh does not meet declared DoF interval')
        path=out/(str(target)+'.vol');report=out/(str(target)+'.check.json');mesh.ngmesh.Save(str(path))
        subprocess.run([sys.executable,str(a.checker.resolve()),str(path),'--strict-labels','--contract',str(contract),'--report-json',str(report)],check=True)
        receipt=json.loads(report.read_text());assert receipt['passed']
        receipt.update(input_sha256=sha(path),contract_sha256=sha(contract),checker_sha256=sha(a.checker));report.write_text(json.dumps(receipt,indent=2)+'\n')
        series.append(dict(target_dofs=target,actual_dofs=n,free_dofs=sum(fes.FreeDofs()),elements=mesh.ne,vertices=mesh.nv,mesh=path.name,check=report.name,sha256=sha(path),maxh=h,trials=trials))
    (out/'meshes.json').write_text(json.dumps(series,indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);p.add_argument('--checker',type=Path,required=True);main(p.parse_args())
