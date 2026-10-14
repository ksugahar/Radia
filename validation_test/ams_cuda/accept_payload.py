"""Stage isolated native DLL payload; verify same-artifact MATLAB acceptance.

Third-party DLLs are a private acceptance bundle, not a published wheel payload.
No installation or source/runtime replacement occurs.
"""
from pathlib import Path
import argparse,hashlib,json,re,shutil,subprocess,sys
import numpy as np
from scipy.io import loadmat,savemat

ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def stage(out):
    out.mkdir(parents=True,exist_ok=False)
    def copy(p,name):
        dst=out/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dst)
    ams=ROOT/'build-ams-cuda';fmm=ROOT/'build-fmm-cuda/job-2rqwcq2f'
    a=json.loads((ams/'build.json').read_text());f=json.loads((fmm/'build.json').read_text())
    if sha(ams/'radia_cuda_ams.mexw64')!=a['sha256']['build-ams-cuda/radia_cuda_ams.mexw64']:raise RuntimeError('AMS build identity mismatch')
    if sha(fmm/'radia_cuda_fmm.mexw64')!=f['binary_sha256']:raise RuntimeError('FMM build identity mismatch')
    copy(ams/'radia_cuda_ams.mexw64','radia_cuda_ams.mexw64');copy(fmm/'radia_cuda_fmm.mexw64','radia_cuda_fmm.mexw64')
    copy(ams/'build.json','ams-build.json');copy(fmm/'build.json','fmm-build.json')
    for p in (ROOT/'packages/radia-gpu/matlab/+radia/+gpu').glob('*.m'):copy(p,Path('matlab/+radia/+gpu')/p.name)
    for name in ['accept_payload.m','run_payload.ps1','native_hcurl.m','fmm_numeric.m','fmm_physics.m','fmm_boundary.m']:copy(Path(__file__).with_name(name),name)
    cases=Path('C:/temp')
    copy(cases/'ams-native-hcurl-warm-20261013/model.vol','model.vol')
    copy(cases/'ams-native-hcurl-warm-20261013/vol_check.json','vol-check.json')
    copy(cases/'ams-boundary-unitcube-20261013.vol','boundary.vol')
    copy(cases/'ams-boundary-fmm3-20261013/vol_check.json','boundary-vol-check.json')
    copy(cases/'ams-fmm-physics3-20261013/input.mat','physics-input.mat')
    q=loadmat(cases/'ams-fmm-numeric-20261013/input-8.mat');q.update(separation=1.5,leaf=16.)
    savemat(out/'numeric-input.mat',{k:v for k,v in q.items() if not k.startswith('__')})
    # Resolve recursively from supported installed native libraries. MATLAB,
    # CUDA, VC runtime and Windows DLLs remain external platform dependencies.
    import netgen
    net=Path(netgen.__file__).parent
    roots=[net,Path(sys.prefix)/'bin',net.parent/'ngsolve_openblas',Path(sys.base_prefix)]
    available={p.name.lower():p for r in roots for p in r.glob('*.dll') if not p.name.lower().startswith(('vcruntime','msvcp','vcomp'))}
    dump=Path('C:/BuildTools/VS2022-CUDA/VC/Tools/MSVC/14.44.35207/bin/Hostx64/x64/dumpbin.exe')
    queue=list(out.glob('*.mexw64'));imports={};external=set();done=set()
    while queue:
        p=queue.pop();name=p.name.lower()
        if name in done:continue
        done.add(name)
        raw=subprocess.check_output([str(dump),'/dependents',str(p)],text=True,errors='replace')
        dependencies=re.findall(r'^\s+([\w.\-]+\.dll)\s*$',raw,re.M|re.I);imports[name]=dependencies
        for dependency in dependencies:
            key=dependency.lower()
            if key in available:
                target=out/'runtime'/available[key].name
                if not target.exists():copy(available[key],Path('runtime')/available[key].name)
                queue.append(target)
            else:external.add(key)
    (out/'modules.ps1').write_text("param([int]$ProcessId,[string]$OutputPath)\n$ErrorActionPreference='Stop'\n@(Get-Process -Id $ProcessId -Module | ForEach-Object { [PSCustomObject]@{name=$_.ModuleName;path=$_.FileName;sha256=([BitConverter]::ToString([Security.Cryptography.SHA256]::Create().ComputeHash([IO.File]::ReadAllBytes($_.FileName)))).Replace('-','').ToLower()} }) | ConvertTo-Json -Depth 4 | Set-Content -LiteralPath $OutputPath -Encoding UTF8\n")
    (out/'dependency-imports.json').write_text(json.dumps(dict(imports=imports,external_platform_dlls=sorted(external)),indent=2))
    (out/'owner.json').write_text(json.dumps(dict(owner='ams-fem-fmm',job='ams-payload-20261013',release_accepted=False)))
    manifest={p.relative_to(out).as_posix():sha(p) for p in out.rglob('*') if p.is_file()}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Staged',len(manifest),'files; native DLLs',len(list((out/'runtime').glob('*.dll'))))

def verify(job,out,baseline_root):
    r=dict(schema='radia.ams-fmm-payload-acceptance.v1',accepted=False,release_accepted=False,timing_accepted=False)
    try:
        manifest=json.loads((job/'manifest.json').read_text())
        for n,h in manifest.items():
            if sha(job/n)!=h:raise RuntimeError('Payload hash mismatch: '+n)
        modules=json.loads((job/'modules.json').read_text(encoding='utf-8-sig'));loaded={x['name'].lower():x for x in modules}
        gateway=Path(loaded['radia_cuda_ams.mexw64']['path']).parent
        for name in ['radia_cuda_ams.mexw64','radia_cuda_fmm.mexw64']:
            if loaded[name]['sha256']!=manifest[name]:raise RuntimeError('Loaded MEX identity mismatch: '+name)
        for p in (job/'runtime').glob('*.dll'):
            m=loaded.get(p.name.lower())
            if m is None:raise RuntimeError('Required native dependency not loaded: '+p.name)
            if m['sha256']!=sha(p):raise RuntimeError('Loaded native dependency hash mismatch: '+p.name)
            if Path(m['path']).parent!=gateway/'runtime':raise RuntimeError('Native dependency loaded outside payload: '+p.name)
        r['matlab']=json.loads((job/'matlab-result.json').read_text());r['payload_sha256']=manifest
        # Cross-host checks against the existing independently certified cases.
        pairs=[('hcurl.mat','ams-native-hcurl-warm-20261013/native.mat',['solution','rhs','mass_values']),('physics.mat','ams-fmm-physics3-20261013/output.mat',['slPotential','dlPotential','B','A','cloudField']),('numeric.mat','ams-fmm-numeric-20261013/output-8.mat',['potential','field']),('boundary.mat','ams-boundary-fmm3-20261013/native.mat',['slResult','dlResult'])]
        differences={}
        for result,baseline,keys in pairs:
            x=loadmat(job/result,simplify_cells=True);y=loadmat(baseline_root/baseline,simplify_cells=True)
            for key in keys:
                differences[result+':'+key]=float(np.linalg.norm(np.asarray(x[key])-np.asarray(y[key]))/max(np.linalg.norm(y[key]),np.finfo(float).tiny))
        if max(differences.values())>1e-8:raise RuntimeError('Same-artifact cross-host gate failed')
        # Re-evaluate the original nonlinear operator independently on CPU.
        import ngsolve as ng
        mesh=ng.Mesh(str(job/'model.vol'));fes=ng.HCurl(mesh,order=1,nograds=True,dirichlet='.*');u,v=fes.TnT()
        b2=ng.curl(u)*ng.curl(u);a=ng.BilinearForm(fes);a+=ng.Variation((.5*b2+4*.25*b2*b2+.5*u*u)*ng.dx)
        f=ng.LinearForm(fes);f+=ng.CoefficientFunction((20,0,0))*v*ng.dx
        x=ng.GridFunction(fes);x.vec.FV().NumPy()[:]=loadmat(job/'hcurl.mat',simplify_cells=True)['solution'];res=x.vec.CreateVector()
        with ng.TaskManager():f.Assemble();a.Apply(x.vec,res)
        free=np.asarray(list(fes.FreeDofs()),dtype=bool);residual=float(np.linalg.norm((res.FV().NumPy()-f.vec.FV().NumPy())[free])/np.linalg.norm(f.vec.FV().NumPy()[free]))
        if residual>1e-9:raise RuntimeError('Independent nonlinear true residual failed')
        r.update(accepted=True,cross_host_relative=differences,independent_cpu_true_residual=residual,native_dll_sha256={p.name:sha(p) for p in (job/'runtime').glob('*.dll')},scope='Exact native candidate payload; not complete installed radia-gpu wheel')
    except BaseException as e:r['error']=repr(e);raise
    finally:out.write_text(json.dumps(r,indent=2)+'\n')
    print('Accepted; independent CPU true residual',residual)

if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='mode',required=True);a=s.add_parser('stage');a.add_argument('--out',type=Path,required=True);a=s.add_parser('verify');a.add_argument('--job',type=Path,required=True);a.add_argument('--out',type=Path,required=True);a.add_argument('--baseline-root',type=Path,default=Path('C:/temp'));a=p.parse_args()
    stage(a.out) if a.mode=='stage' else verify(a.job,a.out,a.baseline_root)
