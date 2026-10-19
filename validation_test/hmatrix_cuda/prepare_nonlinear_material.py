"""Independent CPU NGSolve energy/gradient/tangent fixture; no GPU execution."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--mesh',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False,parents=True)
    import ngsolve as ng
    import radia
    from radia.vim._nonlinear import _EnergyMaterialQuadrature
    from scipy.sparse import coo_matrix
    from hdiv_hdf5 import write_hdiv_h5
    ng.SetNumThreads(1);mu0=4e-7*np.pi;h=np.array([0.,1.,10.,100.,1000.,1e4,1e5]);ms=1e6;chi=1000.
    table=np.column_stack([h,mu0*(h+ms*chi*h/(ms+chi*h))]);result=dict(accepted=False,scope='Independent CPU material quadrature, not complete nonlinear solve',units='SI')
    try:
        with ng.TaskManager():
            mesh=ng.Mesh(str(a.mesh));fes=ng.HDiv(mesh,order=2);material=_EnergyMaterialQuadrature(fes,table)
            rng=np.random.default_rng(20261018);m=rng.standard_normal(fes.ndof);material.m_field.vec.FV().NumPy()[:]=m;material.samples.Interpolate(material.m_field)
            sample=material.samples.vec.FV().NumPy().reshape(3,material.npoints).T.copy();m*=.3*ms/np.max(np.linalg.norm(sample,axis=1))
            energy=material.update(m);values=material.samples.vec.FV().NumPy().reshape(3,material.npoints).T.copy()
            u,v=fes.TnT();lf=ng.LinearForm(fes);lf+=material.field*v*material.measure;lf.Assemble()
            af=ng.BilinearForm(fes);af+=material.bilinear_integrator(ng.InnerProduct(material.tangent*u,v));af.Assemble();i,j,x=af.mat.COO();W=coo_matrix((x,(i,j)),shape=(fes.ndof,fes.ndof)).tocsr();W.sort_indices()
            xyz=ng.GridFunction(ng.VectorValued(material.space,dim=3));xyz.Interpolate(ng.CoefficientFunction((ng.x,ng.y,ng.z)));points=xyz.vec.FV().NumPy().reshape(3,material.npoints).T.copy()
            direction=rng.standard_normal(fes.ndof);direction/=np.linalg.norm(direction)
            payload=dict(table=table,m=m,values=values,weights=material.weights,points=points,gradient=lf.vec.FV().NumPy().copy(),tangent=W,energy=np.array([energy]),direction=direction)
            write_hdiv_h5(a.out/'cpu.h5',payload,dict(units='SI',labels=list(mesh.GetBoundaries())))
            result.update(accepted=True,dofs=fes.ndof,points=material.npoints,integration_order=material.order,energy=energy,mesh_sha256=hashlib.sha256(a.mesh.read_bytes()).hexdigest(),radia_source=str(Path(radia.__file__).resolve()),ngsolve=ng.__version__)
    except BaseException as ex:result['error']=repr(ex);raise
    finally:(a.out/'result.json').write_text(json.dumps(result,indent=2))
if __name__=='__main__':main()
