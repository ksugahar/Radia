"""Independent NGSolve FE staging parity; not HDiv charge/solve acceptance."""
from pathlib import Path
import argparse, ctypes, hashlib, json, os, sys

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--projection',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    a=parser.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    result=dict(accepted=False,scope='NGSolve high-order mapped geometry and mass staging only',
                general_hdiv_validation_complete=False,true_residual=None,
                true_residual_reason='No charge operator or material solve in this subgate',cases=[],
                binary_sha256=sha(a.projection/'radia_cuda_projection.mexw64'),source_sha256=sha(__file__))
    eng=None
    try:
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
        import numpy as np
        import ngsolve as ng
        from netgen.csg import CSGeometry,Sphere,Pnt
        from scipy.sparse import coo_matrix
        from scipy.io import savemat,loadmat
        from radia.vim._vim import _charge_basis,_charge_basis_curved
        root=Path(__file__).resolve().parents[2]
        mesh=ng.Mesh(str(root/'tests/fixtures/beam/affine_field_tetra.vol'))
        cases=[('affine',mesh,[1,2,3])]
        geometry=CSGeometry();geometry.Add(Sphere(Pnt(0,0,0),1))
        curved=ng.Mesh(geometry.GenerateMesh(maxh=3,curvaturesafety=.5,segmentsperedge=1));curved.Curve(2)
        cases.append(('curved',curved,[1,2]))
        prefix=Path(sys.base_prefix)
        os.environ['PATH']=os.pathsep.join(str(prefix/x) for x in ('','bin','Lib/site-packages/netgen','Lib/site-packages/netgen_mesher.libs','Lib/site-packages/ngsolve.libs'))+os.pathsep+os.environ['PATH']
        import matlab.engine
        eng=matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread')
        eng.addpath(str(a.projection),nargout=0);result['matlab']=eng.version();result['ngsolve']=ng.__version__
        refs=[(0,0,0),(1,0,0),(0,1,0),(0,0,1),(.5,0,0),(.5,.5,0),(0,.5,0),(0,0,.5),(.5,0,.5),(0,.5,.5)]
        for name,mesh,orders in cases:
            path=a.out/(name+'.vol');mesh.ngmesh.Save(str(path));mesh=ng.Mesh(str(path))
            for order in orders:
                with ng.TaskManager():
                    space=ng.HDiv(mesh,order=order);u,v=space.TnT()
                    mass=ng.BilinearForm(space);mass+=ng.InnerProduct(u,v)*ng.dx;mass.Assemble()
                    rows,cols,values=mass.mat.COO();matrix=coo_matrix((values,(rows,cols)),shape=(space.ndof,space.ndof)).toarray()
                    node_rule=ng.IntegrationRule(refs,[1.]*len(refs))
                    node_rows=[]
                    for e in mesh.Elements(ng.VOL):
                        trafo=mesh.GetTrafo(ng.ElementId(ng.VOL,e.nr))
                        node_rows.append([list(trafo(r).point) for r in node_rule])
                    nodes=np.array(node_rows).reshape(-1,30)
                    charges=(_charge_basis_curved(space,4) if name=='curved' else _charge_basis(space,4))
                    B=charges['B'].toarray()
                cpu_path=a.out/f'{name}-{order}-cpu.mat'
                savemat(cpu_path,dict(mass=matrix,cell_nodes=nodes,charge_map=B,host=charges['host'],kind=charges['kind'],expo=charges['expo']))
                eng.workspace['mesh_path']=str(path);eng.workspace['degree']=float(order)
                output=a.out/f'{name}-{order}.mat';eng.workspace['output_path']=str(output)
                eng.eval("staging=radia_cuda_projection('inspect_hdiv_elements',mesh_path,degree);save(output_path,'-struct','staging','-v7');",nargout=0)
                actual=loadmat(output)
                error=float(np.linalg.norm(actual['mass']-matrix)/np.linalg.norm(matrix))
                geometry_error=float(np.max(np.abs(actual['cell_nodes']-nodes)))
                charge_error=float(np.linalg.norm(actual['charge_map']-B)/np.linalg.norm(B))
                for key in ('host','kind','expo'):
                    if not np.array_equal(actual[key].reshape(-1),np.asarray(charges[key]).reshape(-1)):raise RuntimeError(f'Charge ordering mismatch: {key}')
                case=dict(geometry=name,order=order,cells=mesh.ne,dofs=space.ndof,mass_relative_error=error,mapped_node_max_absolute_error=geometry_error,curved_cells=int(actual['curved_cells'][0,0]),mesh_sha256=sha(path),cpu_reference_sha256=sha(cpu_path),native_output_sha256=sha(output))
                result['cases'].append(case)
                case['charge_map_relative_error']=charge_error
                if error>1e-12 or geometry_error>1e-13:raise RuntimeError(f'FE staging parity failed: {case}')
                if charge_error>1e-11:raise RuntimeError(f'Charge map parity failed: {case}')
                if name=='curved' and not case['curved_cells']:raise RuntimeError('Curved fixture lost its mapped geometry')
        result['accepted']=True
    except BaseException as exc:
        result['error']=f'{type(exc).__name__}: {exc}';raise
    finally:
        if eng is not None:eng.quit()
        (a.out/'result.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':main()
