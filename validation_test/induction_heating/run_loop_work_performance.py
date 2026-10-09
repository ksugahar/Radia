"""Phase-resolved P0 loop-product study on a self-authored flat ring.

old reproduces the production COO/dense pipeline; exact uses native matrix
products with the identical assembly. fmm is a separately labelled candidate.
Use separate processes on the same host/thread controls for memory comparisons.
"""
import argparse,ctypes,json,time,platform
from pathlib import Path
import numpy as np
import ngsolve as ng
from scipy.sparse import coo_matrix
from ngsolve.bem import LaplaceSL
from run_loop_work_ring import ring_mesh

FMM_OPTIONS = dict(fmm_minorder=20, fmm_maxdirect=100, fmm_separation=2,
                   fmm_eval_separation=3, fmm_maxlevel=20,
                   fmm_split_kr=5, fmm_order_factor=2)


def memory():
    if platform.system()!='Windows':
        import resource
        return dict(peak_working_set_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
    class Counters(ctypes.Structure):
        _fields_=[('cb',ctypes.c_ulong),('faults',ctypes.c_ulong)]+[(name,ctypes.c_size_t) for name in ('peak','current','paged_peak','paged','nonpaged_peak','nonpaged','pagefile','pagefile_peak')]
    c=Counters();c.cb=ctypes.sizeof(c)
    ctypes.windll.kernel32.GetCurrentProcess.restype=ctypes.c_void_p
    ctypes.windll.psapi.GetProcessMemoryInfo.argtypes=[ctypes.c_void_p,ctypes.POINTER(Counters),ctypes.c_ulong]
    handle=ctypes.windll.kernel32.GetCurrentProcess()
    if not ctypes.windll.psapi.GetProcessMemoryInfo(handle,ctypes.byref(c),c.cb):raise RuntimeError('Memory counter failed')
    return dict(peak_working_set_bytes=c.peak,working_set_bytes=c.current,peak_commit_bytes=c.pagefile_peak)


def run(level,mode):
    phases=[]
    def phase(name,t):phases.append(dict(name=name,seconds=time.perf_counter()-t,**memory()))
    with ng.TaskManager():
        t=time.perf_counter();mesh,points,tri=ring_mesh(level)
        space=ng.SurfaceL2(mesh,order=0,dual_mapping=False)
        dofs=np.array([space.GetDofNrs(el)[0] for el in mesh.Elements(ng.BND)])
        assert len(dofs)==space.ndof and np.array_equal(np.sort(dofs),np.arange(space.ndof))
        u,v=space.TnT();measure=ng.ds(bonus_intorder=4)
        phase('mesh_space',t)
        t=time.perf_counter();operator=LaplaceSL(u*measure,use_fmm=mode=='fmm',**FMM_OPTIONS)*v*measure
        matrix=operator.mat;phase('native_single_layer',t)
        if mode=='old':
            t=time.perf_counter();rows,cols,values=matrix.COO();phase('coo_export',t)
            t=time.perf_counter();dense=coo_matrix((values,(rows,cols)),shape=(space.ndof,space.ndof)).toarray();phase('coo_dense_conversion',t)
            t=time.perf_counter();dense=np.ascontiguousarray(dense[np.ix_(dofs,dofs)]);phase('face_order_copy',t)
        angle=np.arange(len(dofs))*.317
        source=np.c_[np.sin(angle),np.cos(.71*angle),np.sin(.37*angle+.2)]
        forward=[];reverse=[];t=time.perf_counter()
        for component in range(3):
            if mode=='old':
                forward.append(dense@source[:,component]);reverse.append(dense.T@source[:,component])
            else:
                x=matrix.CreateColVector();y=matrix.CreateRowVector()
                x.FV().NumPy()[dofs]=source[:,component]
                matrix.Mult(x,y);forward.append(y.FV().NumPy()[dofs].copy())
                matrix.T.Mult(x,y);reverse.append(y.FV().NumPy()[dofs].copy())
        phase('six_products',t)
        forward=np.array(forward).T;reverse=np.array(reverse).T
        assert np.all(np.isfinite(forward)) and np.all(np.isfinite(reverse))
        diagonal=np.sum(source*forward)
        assert diagonal>0
        return dict(schema='radia.loop-work-product-profile.v1',mode=mode,level=level,faces=len(tri),vertices=len(points),
            quadrature_bonus=4,fmm_order=20 if mode=='fmm' else None,
            fmm_parameters=FMM_OPTIONS.copy() if mode=='fmm' else None,ngsolve_version=ng.__version__,
            platform_class=platform.system(),threads=1,dense_matrix_bytes=8*len(tri)**2,
            total_seconds=sum(p['seconds'] for p in phases),phases=phases,**memory()),forward,reverse


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--level',type=int,required=True);p.add_argument('--mode',choices=['old','exact','fmm'],required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--reference',type=Path)
    a=p.parse_args();ng.SetNumThreads(1)
    result,f,r=run(a.level,a.mode)
    if a.reference:
        reference=np.load(a.reference)
        errors={name:float(np.linalg.norm(actual-reference[name])/np.linalg.norm(reference[name])) for name,actual in [('forward',f),('reverse',r)]}
        result['reference_relative_errors']=errors
        assert max(errors.values()) <= (1e-12 if a.mode=='exact' else 1e-5),errors
    np.savez(a.output.with_suffix('.npz'),forward=f,reverse=r)
    a.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result),flush=True)
