"""Execute solver-heavy validation or verify an unchanged-source trusted receipt.

This reuses source validation, NEVER certification of a newly built artifact.
Ordinary tests, native builds and release gates remain separate obligations.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import io
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import urllib.request
import urllib.parse
import zipfile
import xml.etree.ElementTree as ET
from datetime import datetime, timezone, timedelta

PROFILE = 'solver-heavy'
WORKFLOW = '.github/workflows/sparsesolv.yml'
ARTIFACT = 'solver-heavy-validation-v1'
SCHEMA = 'radia.validated-source.v1'
MAX_AGE = timedelta(days=30)
JOB_CONTEXT_KEYS = frozenset({'RADIA_CI_VENV', 'RADIA_CI_OUTPUT_DIR'})
NATIVE_NAMES = ('_radia_pybind.pyd', 'sparsesolv_ngsolve.pyd', 'peec_matrices.pyd', 'axifem.pyd')
BUILD_METADATA = frozenset({'CMakeCache.txt', '.ninja_deps', '.ninja_log',
    'build.ninja', 'build_temp_msvc.bat', 'build.log', 'cmake_install.cmake',
    'CTestTestfile.cmake', 'DartConfiguration.tcl'})
ROOT_INPUTS = {'CMakeLists.txt','Build.ps1','setup.py','pyproject.toml','MANIFEST.in',
               'tools/run_test_tier.py','tools/native_build_provenance.ps1',
               'tools/validation_source_gate.py','tests/conftest.py','conftest.py','.gitignore'}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def canonical(data):
    return json.dumps(data, sort_keys=True, separators=(',',':')).encode()


def numerical_environment():
    """Normalize only reviewed per-run paths inside the runner scratch root."""
    values={key.upper():value for key,value in sorted(os.environ.items())
            if key.upper().startswith(('RADIA_','NGS','MKL_','KMP_','OMP_'))}
    temporary=os.environ.get('RUNNER_TEMP')
    if temporary and Path(temporary).is_absolute():
        scratch=Path(temporary).resolve()
        for key in JOB_CONTEXT_KEYS:
            value=values.get(key)
            if value and Path(value).is_absolute():
                resolved=Path(value).resolve()
                if resolved!=scratch and resolved.is_relative_to(scratch):
                    values[key]='<runner-job-context>'
    return values


def git(root, *args):
    return subprocess.check_output(['git','-c',f'safe.directory={root}',*args],cwd=root).decode('utf-8').strip()



def validate_native_sidecar(root, name='_radia_pybind.pyd'):
    """Verify current native provenance; this is not prior source certification."""
    if name not in NATIVE_NAMES: raise ValueError('unreviewed native module')
    binary=root/'src/radia'/name
    sidecar=Path(str(binary)+'.build.json')
    if binary.is_symlink() or sidecar.is_symlink():
        raise ValueError('linked native provenance')
    record=json.loads(sidecar.read_text(encoding='utf-8'))
    data=binary.read_bytes()
    if (record.get('schema')!='radia.native-build-provenance.v1'
            or record.get('source_commit')!=git(root,'rev-parse','HEAD')
            or record.get('source_dirty') is not False
            or record.get('source_change_fingerprint_sha256')!=digest(b'\0')
            or record.get('binary_name')!=binary.name
            or type(record.get('binary_bytes')) is not int
            or record['binary_bytes']!=len(data)
            or record.get('binary_sha256')!=digest(data)):
        raise ValueError('current native provenance mismatch')
    return record


def validated_native_inputs(root):
    """Bind the exact required modules, build directory and copied outputs."""
    build=root/'build-msvc'
    if build.is_symlink() or build.is_junction() or build.resolve()!=root.resolve()/'build-msvc':
        raise ValueError('linked build directory')
    allowed=set()
    for name in NATIVE_NAMES:
        record=validate_native_sidecar(root,name)
        if record.get('build_directory')!='build-msvc':
            raise ValueError('native build directory missing/mismatched')
        binary=root/'src/radia'/name
        stem=Path(name).stem
        choices=[parent/file for parent in (build,build/'Release')
                 for file in (name,stem+'.cp312-win_amd64.pyd')]
        outputs=[p for p in dict.fromkeys(choices) if p.exists()]
        if not outputs: raise ValueError('native build output missing')
        for output in outputs:
            if (output.is_symlink() or output.parent.is_symlink() or output.parent.is_junction()
                    or not output.is_file() or digest(output.read_bytes())!=record['binary_sha256']):
                raise ValueError('native build/copy mismatch')
            allowed.add(output.relative_to(root).as_posix())
        allowed.update({binary.relative_to(root).as_posix(),binary.relative_to(root).as_posix()+'.build.json'})
    return allowed


def mkl_build_copy(root):
    """Accept one CMake runtime copy only against its installed venv owner."""
    package=importlib.metadata.distribution('mkl')
    files=[file for file in (package.files or ()) if Path(str(file)).name.lower()=='mkl_rt.3.dll']
    if len(files)!=1 or not package.version: raise ValueError('MKL runtime owner unknown')
    source=Path(package.locate_file(files[0]))
    prefix=Path(sys.prefix).resolve()
    if (source.is_symlink() or not source.is_file() or not source.resolve().is_relative_to(prefix)):
        raise ValueError('MKL runtime is not owned by the active venv')
    copied=root/'build-msvc/mkl_rt.3.dll'
    if copied.is_symlink() or not copied.is_file(): raise ValueError('MKL build copy missing/linked')
    sha=digest(source.read_bytes())
    if digest(copied.read_bytes())!=sha: raise ValueError('MKL build copy mismatch')
    return {'package':'mkl','version':package.version,'sha256':sha}


def generated_build_metadata(root, path):
    parts=Path(path).parts
    if not parts or parts[0]!='build-msvc': return False
    if Path(path).suffix.lower() in {'.pyd','.dll','.pth','.py'}: return False
    cmake_metadata=re.fullmatch(
        r'CMakeFiles/(?:(?:[0-9][^/]+/)?(?:CMake(?:C|CXX|RC)Compiler\.cmake|CMakeSystem\.cmake|'
        r'CMakeDetermineCompilerABI_(?:C|CXX)\.bin)|'
        r'[0-9][^/]+/CompilerId(?:C|CXX)/CMake(?:C|CXX)CompilerId\.(?:c|cpp|exe)|'
        r'ShowIncludes/(?:foo\.h|main\.c)|CMakeConfigureLog\.yaml|'
        r'CMake(?:Output|Error)\.log|TargetDirectories\.txt|cmake\.check_cache|rules\.ninja)',
        '/'.join(parts[1:]))
    if not (len(parts)==2 and parts[1] in BUILD_METADATA or cmake_metadata):
        return False
    file=root/path
    if not file.is_file() or any(p.is_symlink() or p.is_junction() for p in (file,*file.parents) if p!=root.parent):
        raise ValueError('linked build metadata')
    return True


def environment_contract(root):
    validated_native_inputs(root)
    mkl_runtime=mkl_build_copy(root)
    compiler = []
    for file in (root/'build-msvc/CMakeFiles').glob('*/CMakeCXXCompiler.cmake'):
        text=file.read_text(encoding='utf-8')
        compiler.append({key:re.findall(r'set\('+key+r' "([^"\n]+)"\)',text)
                         for key in ('CMAKE_CXX_COMPILER_ID','CMAKE_CXX_COMPILER_VERSION')})
    if len(compiler)!=1 or any(len(v)!=1 for v in compiler[0].values()):
        raise ValueError('compiler identity unknown')
    cache=(root/'build-msvc/CMakeCache.txt').read_text(encoding='utf-8')
    # Keep value contents (including quoted flag spacing) intact. CMake cache
    # comments and entry ordering do not define the compiler contract.
    build_options={}
    for line in cache.splitlines():
        match=re.match(r'(CMAKE_CXX_FLAGS(?:_[A-Z0-9_]+)?|CMAKE_BUILD_TYPE):[^=]+=(.*)$',line)
        if match:
            key,value=match.groups()
            if key in build_options: raise ValueError('duplicate compiler option')
            build_options[key]=value.strip()
    if 'CMAKE_CXX_FLAGS' not in build_options or 'CMAKE_BUILD_TYPE' not in build_options:
        raise ValueError('compiler options unknown')
    sdk_versions=sorted(set(re.findall(r'10\.0\.\d+\.\d+',cache)))
    if not sdk_versions: raise ValueError('Windows SDK version unknown')
    if os.environ.get('RADIA_TESTS_ALLOW_PARTIAL')!='0': raise ValueError('partial validation context')
    threads={key:os.environ.get(key) for key in ('RADIA_TEST_NGSOLVE_THREADS','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')}
    if set(threads.values())!={'1'}: raise ValueError('unbounded or changed thread contract')
    build_tools={}
    for name,key in [('cmake','CMAKE_COMMAND'),('ninja','CMAKE_MAKE_PROGRAM')]:
        matches=re.findall(r'^'+key+r':[^=]+=(.+)$',cache,re.MULTILINE)
        if len(matches)!=1: raise ValueError('build tool identity unknown: '+name)
        executable=Path(matches[0].strip())
        if not executable.is_absolute() or not executable.is_file() or executable.is_symlink():
            raise ValueError('missing/linked build tool: '+name)
        version=subprocess.run([str(executable),'--version'],check=True,capture_output=True,
                               text=True,encoding='utf-8',errors='strict',timeout=10).stdout.strip()
        if not version: raise ValueError('empty build tool version: '+name)
        build_tools[name]={'sha256':digest(executable.read_bytes()),'version':version}
    return {'python':sys.version,'platform':platform.platform(),'machine':platform.machine(),
            'processor':platform.processor(),'compiler':compiler,'windows_sdk_versions':sdk_versions,'threads':threads,
            'compiler_options':build_options,'build_tools':build_tools,'numerical_environment':numerical_environment(),
            'mkl_runtime':mkl_runtime,
            'packages':sorted((d.metadata['Name'].lower(),d.version) for d in importlib.metadata.distributions())}


def source_fingerprint(root, environment, numerical_plan=None):
    manifest=json.loads((root/'tests/test_tier_manifest.json').read_text(encoding='utf-8'))
    profiles={}
    name=PROFILE
    while name:
        if name in profiles: raise ValueError('cyclic validation profile')
        profiles[name]=manifest['profiles'][name]
        name=profiles[name].get('extends')
    tracked=git(root,'ls-files','-z').split('\0')
    # Inspect ignored as well as ordinary untracked input. Only reconstructible
    # compiler/Python cache outputs are excluded; unknown inputs force execution.
    untracked=git(root,'ls-files','--others','-z').split('\0')
    generated={'.lib','.obj','.exp','.pdb','.ilk','.pyc'}
    native_inputs=None
    rejected=[]
    for path in untracked:
        parts=Path(path).parts
        suffix=Path(path).suffix.lower()
        if (not path or '__pycache__' in parts and suffix=='.pyc'
                or '.pytest_cache' in parts and suffix not in {'.py','.pth','.pyd','.dll'}):
            continue
        if path.startswith('build-msvc/') or path in {f'src/radia/{name}{suffix}' for name in NATIVE_NAMES for suffix in ('','.build.json')}:
            if native_inputs is None: native_inputs=validated_native_inputs(root)
            if path in native_inputs or generated_build_metadata(root,path): continue
            if path=='build-msvc/mkl_rt.3.dll':
                mkl_build_copy(root)
                continue
        if not (parts[0]=='build-msvc' and suffix in generated):
            rejected.append(path)
    if rejected:
        error=ValueError('untracked validation inputs: '+', '.join(rejected))
        error.paths=rejected
        raise error
    selected=[]
    def test_input(path):
        return path.startswith('tests/') and not (path.startswith(('tests/policy/','tests/ci_routing/','tests/validation_gate/','tests/test_ci_'))
            or path in {'tests/test_local_push_checks.py','tests/test_pre_push_check.py','tests/test_sparsesolv_ci_contract.py','tests/test_tier_manifest.json'})
    for path in tracked:
        if (path in ROOT_INPUTS or test_input(path) or path.startswith(('src/','validation_test/'))
                or (path.startswith('packages/') and 'docs' not in path.split('/'))):
            selected.append(path)
    # All production, validation inputs/helpers/fixtures and relevant build files.
    # CI workflow routing and unrelated compact-test changes are not numerical inputs.
    if not selected or not any(p.startswith('src/') for p in selected): raise ValueError('source closure unavailable')
    for profile in profiles.values():
        for path in profile['paths']:
            if path not in selected: selected.append(path)
    if numerical_plan is not None:
        heavy = [row for row in numerical_plan['cases'] if row['profile'] == PROFILE and row['selected']]
        if not heavy: raise ValueError('empty planned heavy validation')
        tracked_set = set(tracked)
        scoped = {path for row in heavy for path in row['dependencies'] if path in tracked_set}
        # Native solver builds remain indivisible, while Python consumers/fixtures are per-case.
        scoped.update(path for path in tracked if path.startswith(('src/core/', 'src/lib/', 'src/ext/', 'cmake/')))
        scoped.update(path for path in ROOT_INPUTS if path in tracked_set)
        scoped.update(path for path in ('tools/ci_numerical_plan.py', 'tools/ci_numerical_inputs.json',
                                      'tools/ci_numerical_runtime.py') if path in tracked_set)
        conservative = (bool(numerical_plan.get('unclassified_inputs'))
                        or bool(numerical_plan.get('global_selection'))
                        or bool(numerical_plan.get('global_reasons'))
                        or any(row.get('unknowns') for row in heavy))
        # Fallback selection is not a proved closure. Its causing input (including
        # workflow, runtime JSON or docs) must invalidate a prior main receipt.
        # Hash ALL tracked files and the complete environment in that case.
        selected = sorted(tracked_set - {''}) if conservative else sorted(scoped)
        profiles = {'selection_scope': 'full-tree' if conservative else 'proved-case-closure',
                    'selected_paths': sorted(row['path'] for row in heavy), 'profile_contracts': profiles}
        environment = dict(environment)
        # Do not hash the entire venv. Retain numerical/toolchain identities and
        # every declared solver dependency; Gmsh participates only for its consumers.
        relevant = {'numpy','scipy','ngsolve','netgen-mesher','mkl','mkl-devel','mkl-include',
                    'tbb','tbb-devel','intel-openmp','build123d','cadquery-ocp','h5py','pybind11','pytest','pluggy',
                    'cmake','ninja'}
        if any('gmsh' in row.get('components', []) for row in heavy): relevant.add('gmsh')
        aliases = {'OCP':'cadquery-ocp','skimage':'scikit-image','sklearn':'scikit-learn',
                   'kan':'pykan','yaml':'pyyaml','netgen':'netgen-mesher'}
        for row in heavy:
            for component in row.get('components', []):
                if component.startswith('package:'):
                    module = component.split(':', 1)[1]
                    relevant.add(aliases.get(module, module).lower().replace('_','-'))
        if not conservative:
            environment['packages'] = [item for item in environment.get('packages', []) if item[0].replace('_','-') in relevant]
    rows={}
    for path in sorted(selected):
        file=root/path
        if not file.is_file() or file.is_symlink(): raise ValueError('missing/linked validation input: '+path)
        rows[path]=digest(file.read_bytes())
    changes=git(root,'diff','--name-only','HEAD','--').splitlines()
    if set(changes).intersection(rows): raise ValueError('validation inputs dirty')
    payload={'schema':SCHEMA,'profile':PROFILE,'sources':rows,'profiles':profiles,'environment':environment}
    return digest(canonical(payload)),len(rows)


def check_junit(payload):
    cases=list(ET.fromstring(payload).iter('testcase'))
    if not cases or any(c.find(tag) is not None for c in cases for tag in ('failure','error','skipped')):
        raise ValueError('validation must execute every required case successfully')
    return len(cases)


def verified_receipt(root, artifact, run, archive, fingerprint, now):
    """Authoritative GitHub API metadata is separate from artifact contents."""
    if artifact.get('expired') is not False or artifact.get('name')!=ARTIFACT: raise ValueError('expired/wrong artifact')
    if (run.get('path')!=WORKFLOW or run.get('head_branch')!='main' or run.get('status')!='completed'
            or run.get('conclusion')!='success' or run.get('event') not in ('push','workflow_dispatch')):
        raise ValueError('untrusted run provenance')
    if artifact.get('workflow_run',{}).get('id')!=run['id']: raise ValueError('artifact run mismatch')
    if artifact.get('digest')!='sha256:'+digest(archive): raise ValueError('archive digest absent/mismatch')
    expires=datetime.fromisoformat(artifact['expires_at'].replace('Z','+00:00'))
    if expires<=now: raise ValueError('expired artifact')
    updated=datetime.fromisoformat(run['updated_at'].replace('Z','+00:00'))
    if updated>now or now-updated>MAX_AGE: raise ValueError('receipt too old/future')
    # A successful PR/fork artifact can never become a reusable main receipt.
    subprocess.run(['git','-c',f'safe.directory={root}','merge-base','--is-ancestor',run['head_sha'],'origin/main'],
                   cwd=root,check=True,capture_output=True)
    with zipfile.ZipFile(io.BytesIO(archive)) as z:
        if set(z.namelist())!={'receipt.json','solver-heavy.xml'}: raise ValueError('unexpected evidence archive')
        if any(i.file_size>16*1024*1024 for i in z.infolist()): raise ValueError('oversized evidence')
        original_receipt=z.read('receipt.json')
        receipt=json.loads(original_receipt)
        xml=z.read('solver-heavy.xml')
    if (receipt.get('schema')!=SCHEMA or receipt.get('profile')!=PROFILE or receipt.get('accepted') is not True
            or receipt.get('reusable') is not True or receipt.get('fingerprint')!=fingerprint
            or receipt.get('source_sha')!=run['head_sha'] or receipt.get('run_id')!=run['id']
            or receipt.get('run_attempt')!=run['run_attempt'] or receipt.get('junit_sha256')!=digest(xml)
            or type(receipt.get('case_count')) is not int or receipt.get('case_count')!=check_junit(xml)):
        raise ValueError('receipt/input/evidence mismatch')
    return receipt,xml,original_receipt


class _SafeRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, request, fp, code, message, headers, newurl):
        if urllib.parse.urlsplit(newurl).scheme != 'https':
            raise ValueError('non-HTTPS artifact redirect')
        redirected=super().redirect_request(request,fp,code,message,headers,newurl)
        if urllib.parse.urlsplit(newurl).hostname != 'api.github.com':
            redirected.remove_header('Authorization')
        return redirected


class GitHub:
    def __init__(self, repository, token):
        self.base='https://api.github.com/repos/'+repository
        self.token=token
    def get(self, path, binary=False):
        request=urllib.request.Request(self.base+path,headers={'Authorization':'Bearer '+self.token,
            'Accept':'application/vnd.github+json','X-GitHub-Api-Version':'2022-11-28'})
        with urllib.request.build_opener(_SafeRedirect()).open(request,timeout=30) as response:
            data=response.read(32*1024*1024+1)
        if len(data)>32*1024*1024: raise ValueError('oversized API response')
        return data if binary else json.loads(data)
    def find(self, root, fingerprint):
        # Bounded search; pagination/cache miss safely executes validation.
        artifacts=self.get('/actions/artifacts?name='+ARTIFACT+'&per_page=30')['artifacts']
        for artifact in artifacts[:30]:
            try:
                run=self.get('/actions/runs/'+str(artifact['workflow_run']['id']))
                if run.get('head_branch')!='main' or run.get('conclusion')!='success': continue
                archive=self.get('/actions/artifacts/'+str(artifact['id'])+'/zip',binary=True)
                receipt,xml,original_receipt=verified_receipt(root,artifact,run,archive,fingerprint,datetime.now(timezone.utc))
                return artifact,receipt,xml,original_receipt
            except Exception:
                continue  # Untrusted/unknown/expired is a miss, never PASS.
        return None


def run(root, out, force=False, client=None, numerical_plan=None):
    planned = None
    if numerical_plan is not None:
        from ci_numerical_runtime import verified_plan, selected_paths
        planned = verified_plan(root, numerical_plan)
        selected_paths(planned, PROFILE)
    out.mkdir(parents=True,exist_ok=True)
    identity={'source_sha':git(root,'rev-parse','HEAD'),'run_id':int(os.environ.get('GITHUB_RUN_ID','0')),
              'run_attempt':int(os.environ.get('GITHUB_RUN_ATTEMPT','0'))}
    result={'schema':SCHEMA,'profile':PROFILE,'current':identity,'accepted':False,
            'current_artifact_certified':False,'mode':'execute','reuse_reason':None,
            'verified_at':datetime.now(timezone.utc).isoformat()}
    # Operational paths remain observable without making each run a new
    # numerical environment. No other RADIA/NGS/library variable is normalized.
    result['job_context']={key:os.environ[key] for key in JOB_CONTEXT_KEYS if key in os.environ}
    fingerprint=None
    try:
        environment=environment_contract(root)
        result['mkl_runtime']=environment.get('mkl_runtime')
        fingerprint,count=(source_fingerprint(root,environment,planned) if planned is not None
                           else source_fingerprint(root,environment))
        result.update(fingerprint=fingerprint,source_count=count)
    except Exception as error:
        result['reuse_reason']='unknown fingerprint: '+type(error).__name__
        result['reuse_detail']=str(error)
        result['rejected_untracked_paths']=getattr(error,'paths',[])
        print(result['reuse_reason']+': '+result['reuse_detail'],flush=True)
    found=None
    if not force and fingerprint and client:
        try: found=client.find(root,fingerprint)
        except Exception as error: result['reuse_reason']='receipt lookup unavailable: '+type(error).__name__
    if found:
        artifact,receipt,xml,original_receipt=found
        # Preserve original proof verbatim. Verification has its own current identity.
        (out/'receipt.json').write_bytes(original_receipt)
        (out/'solver-heavy.xml').write_bytes(xml)
        result.update(accepted=True,mode='verified-reuse',original=receipt,artifact_id=artifact['id'])
    else:
        evidence=out/'solver-heavy.xml'
        if evidence.is_symlink(): raise ValueError('linked output evidence')
        evidence.unlink(missing_ok=True)  # Never accept a previous attempt's XML.
        command=[sys.executable,'tools/run_test_tier.py','--profile',PROFILE,'--junitxml',str(out/'solver-heavy.xml')]
        if numerical_plan is not None: command += ['--numerical-plan', str(numerical_plan)]
        process=subprocess.run(command,cwd=root,check=False)
        if process.returncode==0:
            try:
                xml=(out/'solver-heavy.xml').read_bytes()
                count=check_junit(xml)
                receipt={'schema':SCHEMA,'profile':PROFILE,**identity,'fingerprint':fingerprint,
                         'accepted':True,'reusable':bool(fingerprint),'case_count':count,'junit_sha256':digest(xml)}
                (out/'receipt.json').write_bytes(canonical(receipt))
                result.update(accepted=True,case_count=count)
            except Exception as error:
                result['failure']='invalid validation evidence: '+type(error).__name__
        else: result['exit_code']=process.returncode
    (out/'verification.json').write_bytes(canonical(result))
    return 0 if result['accepted'] else 1


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--force',action='store_true')
    parser.add_argument('--numerical-plan',type=Path)
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    token=os.environ.get('GH_TOKEN')
    repo=os.environ.get('GITHUB_REPOSITORY')
    client=GitHub(repo,token) if token and repo else None
    sys.exit(run(root,args.output,args.force,client,args.numerical_plan))
