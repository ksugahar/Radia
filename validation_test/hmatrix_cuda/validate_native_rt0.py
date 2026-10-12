"""Bounded MATLAB-native RT0 setup + CUDA material-solve physical subgate."""
from pathlib import Path
import argparse,hashlib,json,os,sys,time,ctypes

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser()
    for name in ('projection','hmatrix','mesh','reference','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--dependency-root',type=Path);a=p.parse_args();a.out.mkdir(exist_ok=False,parents=True);eng=None
    result=dict(accepted=False,scope='Public affine single-tetra RT0 MATLAB-native setup, exact dense CUDA Gram and linear material solve; not general adaptive HACApK acceptance',source_sha256=sha(Path(__file__)),hashes={name:sha(path) for name,path in [('projection',a.projection/'radia_cuda_projection.mexw64'),('hmatrix',a.hmatrix/'radia_cuda_hmatrix.mexw64'),('mesh',a.mesh),('cpu_reference',a.reference)]})
    try:
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
        prefix=Path(sys.base_prefix);paths=[prefix,prefix/'bin',prefix/'Lib/site-packages/netgen',prefix/'Lib/site-packages/netgen_mesher.libs',prefix/'Lib/site-packages/ngsolve.libs',prefix/'Lib/site-packages/ngsolve_openblas']
        if a.dependency_root:paths=[a.dependency_root.parents[1]/'bin',a.dependency_root.parents[1]]+[a.dependency_root/x for x in ('netgen','netgen_mesher.libs','ngsolve.libs','ngsolve_openblas')]+paths
        os.environ['PATH']=os.pathsep.join(map(str,paths))+os.pathsep+os.environ['PATH']
        import matlab.engine
        eng=matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread');result['matlab']=eng.version()
        eng.addpath(str(a.projection),str(a.hmatrix),str(Path(__file__).resolve().parents[2]/'packages/radia-gpu/matlab'),nargout=0)
        eng.workspace['mesh_path']=str(a.mesh);eng.workspace['reference_path']=str(a.reference);eng.workspace['output_path']=str(a.out/'matlab.json')
        eng.eval("inputs=load(reference_path);test=inputs.solver_cases(end);setup=radia_cuda_projection('assemble_hdiv_rt0',mesh_path,[1 .2 .1]);mass_error=norm(setup.mass-full(test.mass),'fro')/norm(full(test.mass),'fro');rhs_error=norm(setup.rhs-test.rhs)/norm(test.rhs);mean_error=norm(setup.mean_operator-test.mean_operator,'fro')/norm(test.mean_operator,'fro');assert(max([mass_error rhs_error mean_error])<1e-12);",nargout=0)
        eng.eval("solution=radia.gpu.solveHDivRT0(mesh_path,[1 .2 .1],test.inv_chi,1e-11,100);coefficient_error=norm(solution.coefficients-test.expected)/norm(test.expected);physical_error=norm(solution.mean_magnetization-test.mean_operator*test.expected)/norm(test.mean_operator*test.expected);cpu_true_residual=norm(test.rhs-test.operator*solution.coefficients)/norm(test.rhs);assert(max([coefficient_error physical_error cpu_true_residual])<1e-10);assert(solution.true_relative_residual<=1e-11);assert(radia_cuda_hmatrix('live_count')==0);",nargout=0)
        eng.eval("r=struct('accepted',true,'mass_relative_error',mass_error,'rhs_relative_error',rhs_error,'mean_operator_relative_error',mean_error,'coefficient_relative_error',coefficient_error,'mean_magnetization_relative_error',physical_error,'independent_true_relative_residual',cpu_true_residual,'solution',solution);file=fopen(output_path,'w');assert(file>=0);fprintf(file,'%s',jsonencode(r,PrettyPrint=true));fclose(file);",nargout=0)
        result['result']=json.loads((a.out/'matlab.json').read_text(encoding='utf-8'));result['accepted']=True
    except BaseException as e:result['error']=f'{type(e).__name__}: {e}';raise
    finally:
        if eng is not None:eng.quit()
        (a.out/'result.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
if __name__=='__main__':main()
