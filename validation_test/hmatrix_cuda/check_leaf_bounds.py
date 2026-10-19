"""Compile the symbolic leaf coefficient contract; no CUDA or numerical operator."""
from pathlib import Path
import argparse,hashlib,json,subprocess,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'packages/radia-simulink/tools'))
from vs2022_toolchain import discover_vs,compiler_bin

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=False);repo=Path(__file__).resolve().parents[2]
    source=Path(__file__).with_suffix('.cpp');header=repo/'packages/radia-simulink/native/hmatrix_leaf_bounds.h'
    result=dict(accepted=False,scope='Symbolic coefficient/ACA payload lifetime bounds; production callback integration remains unbuilt and excludes host metadata/runtime',native_production_built=False,performance_accepted=False,source_commit=subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],text=True).strip(),source_sha256={p.relative_to(repo).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,header,Path(__file__),repo/'packages/radia-simulink/tools/vs2022_toolchain.py']})
    result['source_dirty']=bool(subprocess.check_output(['git','-C',str(repo),'status','--porcelain'],text=True).strip())
    try:
        vs=discover_vs();compiler=compiler_bin(vs);exe=args.out/'leaf-bounds.exe';batch=args.out/'compile.cmd'
        batch.write_text('@echo off\ncall "'+str(vs/'VC/Auxiliary/Build/vcvars64.bat')+'" >nul\n"'+str(compiler/'cl.exe')+'" /nologo /std:c++17 /EHsc /MD /W4 /WX /I"'+str(header.parent)+'" "'+str(source)+'" /Fo:"'+str(args.out/'leaf-bounds.obj')+'" /Fe:"'+str(exe)+'"\n',encoding='utf-8')
        with (args.out/'compile.log').open('w') as log:subprocess.run(['cmd.exe','/d','/c',str(batch)],cwd=args.out,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=45)
        run=subprocess.run([str(exe)],capture_output=True,text=True,check=True,timeout=5);(args.out/'run.log').write_text(run.stdout+run.stderr)
        result['contract']=json.loads(run.stdout);assert result['contract']['accepted'];result['accepted']=True
        result['executable_sha256']=hashlib.sha256(exe.read_bytes()).hexdigest()
    except BaseException as ex:result['error']=repr(ex);raise
    finally:(args.out/'result.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
