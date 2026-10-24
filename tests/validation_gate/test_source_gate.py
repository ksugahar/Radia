"""Guards: untrusted receipts, changed numerical inputs and incorrect scope reuse fail.

Assumes: source-contract fixtures and explicit synthetic environments; no numerical
accuracy, native compilation, heavy validation or network is certified here.
"""
import copy
from datetime import datetime,timezone,timedelta
import importlib.util
import io
import json
from pathlib import Path
import zipfile
import pytest

ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('validation_gate',ROOT/'tools/validation_source_gate.py')
gate=importlib.util.module_from_spec(spec);spec.loader.exec_module(gate)
NOW=datetime(2026,10,17,tzinfo=timezone.utc)
XML=b'<testsuite><testcase name="physics"/></testsuite>'


def proof():
    receipt={'schema':gate.SCHEMA,'profile':gate.PROFILE,'source_sha':'a'*40,'run_id':12,'run_attempt':1,
             'fingerprint':'fingerprint','accepted':True,'reusable':True,'case_count':1,'junit_sha256':gate.digest(XML)}
    run={'id':12,'run_attempt':1,'path':gate.WORKFLOW,'head_branch':'main','status':'completed','conclusion':'success',
         'event':'push','head_sha':'a'*40,'updated_at':NOW.isoformat()}
    return receipt,run


def archive(receipt,xml=XML):
    data=io.BytesIO()
    with zipfile.ZipFile(data,'w') as z:
        z.writestr('receipt.json',json.dumps(receipt,indent=2));z.writestr('solver-heavy.xml',xml)
    return data.getvalue()


def artifact(data):
    return {'name':gate.ARTIFACT,'expired':False,'digest':'sha256:'+gate.digest(data),'workflow_run':{'id':12},
            'expires_at':(NOW+timedelta(days=1)).isoformat()}


def test_trusted_receipt_and_verbatim_evidence(monkeypatch,tmp_path):
    receipt,run=proof();data=archive(receipt)
    calls=[];monkeypatch.setattr(gate.subprocess,'run',lambda command,**kw:calls.append(command))
    actual,xml,raw=gate.verified_receipt(tmp_path,artifact(data),run,data,'fingerprint',NOW)
    assert actual==receipt and xml==XML and raw==json.dumps(receipt,indent=2).encode()
    assert '--is-ancestor' in calls[0] and calls[0][-1]=='origin/main'


@pytest.mark.parametrize('field,value',[('head_branch','feature'),('conclusion','failure'),('event','pull_request'),
                                         ('path','.github/workflows/untrusted.yml'),('run_attempt',2)])
def test_untrusted_run_rejected(monkeypatch,tmp_path,field,value):
    receipt,run=proof();run[field]=value;data=archive(receipt)
    monkeypatch.setattr(gate.subprocess,'run',lambda *a,**k:None)
    with pytest.raises(ValueError):gate.verified_receipt(tmp_path,artifact(data),run,data,'fingerprint',NOW)


@pytest.mark.parametrize('change',['archive','expired','expirytime','old','wrong-input','wrong-sha','failed-evidence','missing-digest','boolean-count','unknown-expiry','nonboolean-expiry'])
def test_bad_artifact_rejected(monkeypatch,tmp_path,change):
    receipt,run=proof();xml=XML
    if change=='wrong-sha':receipt['source_sha']='b'*40
    if change=='boolean-count':receipt['case_count']=True
    if change=='failed-evidence':xml=b'<testsuite><testcase><skipped/></testcase></testsuite>';receipt['junit_sha256']=gate.digest(xml)
    data=archive(receipt,xml);meta=artifact(data)
    if change=='archive':data+=b'corrupt'
    if change=='expired':meta['expired']=True
    if change=='unknown-expiry':meta.pop('expired')
    if change=='nonboolean-expiry':meta['expired']=0
    if change=='expirytime':meta['expires_at']=NOW.isoformat()
    if change=='old':run['updated_at']=(NOW-timedelta(days=31)).isoformat()
    if change=='missing-digest':meta.pop('digest')
    fingerprint='different' if change=='wrong-input' else 'fingerprint'
    monkeypatch.setattr(gate.subprocess,'run',lambda *a,**k:None)
    with pytest.raises(ValueError):gate.verified_receipt(tmp_path,meta,run,data,fingerprint,NOW)


def test_main_ancestry_required(monkeypatch,tmp_path):
    receipt,run=proof();data=archive(receipt)
    def fail(*a,**k):raise gate.subprocess.CalledProcessError(1,['git'])
    monkeypatch.setattr(gate.subprocess,'run',fail)
    with pytest.raises(gate.subprocess.CalledProcessError):gate.verified_receipt(tmp_path,artifact(data),run,data,'fingerprint',NOW)


def test_source_scope_includes_physics_inputs_but_not_workflow_routing(tmp_path,monkeypatch):
    files={'src/solver.py':'physics','validation_test/case/input.vol':'mesh','tests/conftest.py':'fixture',
           'packages/demo/tests/fixture.h5':'fixture-bytes','Build.ps1':'build','.github/workflows/sparsesolv.yml':'routing','tools/validation_source_gate.py':'gate'}
    for path,text in files.items():
        p=tmp_path/path;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(text)
    manifest=tmp_path/'tests/test_tier_manifest.json'
    manifest.write_text(json.dumps({'profiles':{'solver-heavy':{'paths':['validation_test/case/test.py']}}}))
    test=tmp_path/'validation_test/case/test.py';test.write_text('assert True');files[str(test.relative_to(tmp_path))]='assert True'
    monkeypatch.setattr(gate,'git',lambda root,*args:'\0'.join(files) if args[0]=='ls-files' and '--others' not in args else '')
    initial,_=gate.source_fingerprint(tmp_path,{'deps':'1'})
    (tmp_path/'.github/workflows/sparsesolv.yml').write_text('new routing')
    assert gate.source_fingerprint(tmp_path,{'deps':'1'})[0]==initial
    for path in ['src/solver.py','validation_test/case/input.vol','Build.ps1','tests/conftest.py','packages/demo/tests/fixture.h5']:
        p=tmp_path/path;p.write_text(files[path]+'changed')
        assert gate.source_fingerprint(tmp_path,{'deps':'1'})[0]!=initial
        p.write_text(files[path])
    assert gate.source_fingerprint(tmp_path,{'deps':'2'})[0]!=initial


@pytest.mark.parametrize('lookup', ['miss','error','force','unknown'])
def test_miss_unknown_force_execute_validation(tmp_path,monkeypatch,lookup):
    monkeypatch.setattr(gate,'git',lambda *a:'a'*40)
    monkeypatch.setattr(gate,'environment_contract',lambda root:{})
    def fingerprint(*a):
        if lookup=='unknown':raise ValueError('unknown')
        return 'fingerprint',1
    monkeypatch.setattr(gate,'source_fingerprint',fingerprint)
    class Client:
        def find(self,*a):
            if lookup=='force':pytest.fail('force must not look up cache')
            if lookup=='error':raise OSError('offline')
            return None
    calls=[]
    def execute(command,**kw):
        calls.append(command);(tmp_path/'out/solver-heavy.xml').write_bytes(XML)
        return type('Result',(),{'returncode':0})()
    monkeypatch.setattr(gate.subprocess,'run',execute)
    assert gate.run(tmp_path,tmp_path/'out',lookup=='force',Client())==0
    report=json.loads((tmp_path/'out/verification.json').read_text())
    assert calls and report['mode']=='execute' and not report['current_artifact_certified']


def test_failed_execution_never_passes(tmp_path,monkeypatch):
    monkeypatch.setattr(gate,'git',lambda *a:'a'*40)
    monkeypatch.setattr(gate,'environment_contract',lambda root:{})
    monkeypatch.setattr(gate,'source_fingerprint',lambda *a:('fingerprint',1))
    monkeypatch.setattr(gate.subprocess,'run',lambda *a,**k:type('Result',(),{'returncode':1})())
    assert gate.run(tmp_path,tmp_path/'out')==1
    assert not json.loads((tmp_path/'out/verification.json').read_text())['accepted']


def test_workflow_keeps_ordinary_and_release_gates_separate():
    workflow=(ROOT/'.github/workflows/sparsesolv.yml').read_text()
    assert '--profile solver-numerics' in workflow and '-RadiaOnly -BuildJobs 2' in workflow
    assert 'tools/validation_source_gate.py' in workflow and 'force_validation:' in workflow
    assert "steps.heavy-validation.outputs.mode == 'execute'" in workflow
    assert 'actions: read' in workflow


def test_validation_job_fetches_main_ancestry():
    import yaml
    workflow=yaml.safe_load((ROOT/'.github/workflows/sparsesolv.yml').read_text())
    checkout=next(step for step in workflow['jobs']['ams-regression']['steps']
                  if step.get('uses','').startswith('actions/checkout@'))
    assert checkout['with']['fetch-depth']==0

def test_redirect_does_not_forward_token_to_artifact_storage():
    request=gate.urllib.request.Request('https://api.github.com/repos/example/repo/actions/artifacts/1/zip',headers={'Authorization':'Bearer secret'})
    redirect=gate._SafeRedirect().redirect_request(request,None,302,'redirect',{},'https://storage.example/signed-archive')
    assert not redirect.has_header('Authorization')
    with pytest.raises(ValueError):gate._SafeRedirect().redirect_request(request,None,302,'redirect',{},'http://storage.example/archive')

def test_verified_reuse_keeps_original_identity_and_bytes(tmp_path,monkeypatch):
    receipt,prior=proof();data=archive(receipt);raw=json.dumps(receipt,indent=2).encode()
    monkeypatch.setattr(gate,'git',lambda *a:'b'*40)
    monkeypatch.setattr(gate,'environment_contract',lambda root:{})
    monkeypatch.setattr(gate,'source_fingerprint',lambda *a:('fingerprint',1))
    class Client:
        def find(self,*a):return {'id':99},receipt,XML,raw
    monkeypatch.setattr(gate.subprocess,'run',lambda *a,**k:pytest.fail('reused validation executed'))
    assert gate.run(tmp_path,tmp_path/'out',client=Client())==0
    report=json.loads((tmp_path/'out/verification.json').read_text())
    assert report['current']['source_sha']=='b'*40 and report['original']['source_sha']=='a'*40
    assert not report['current_artifact_certified'] and report['mode']=='verified-reuse'
    assert (tmp_path/'out/receipt.json').read_bytes()==raw

def test_old_xml_cannot_turn_empty_execution_into_pass(tmp_path,monkeypatch):
    out=tmp_path/'out';out.mkdir();(out/'solver-heavy.xml').write_bytes(XML)
    monkeypatch.setattr(gate,'git',lambda *a:'a'*40)
    monkeypatch.setattr(gate,'environment_contract',lambda root:{})
    monkeypatch.setattr(gate,'source_fingerprint',lambda *a:('fingerprint',1))
    monkeypatch.setattr(gate.subprocess,'run',lambda *a,**k:type('Result',(),{'returncode':0})())
    assert gate.run(tmp_path,out)==1
    assert not json.loads((out/'verification.json').read_text())['accepted']


@pytest.mark.parametrize('input_path', ['src/extra.py', 'validation_test/case/generated.h5',
                                       'packages/demo/tests/new_fixture.vol',
                                       'tests/local/conftest.py', 'local.pth',
                                       'src/radia/radia_motor_rom.dll', 'src/radia/extra.pyd'])
def test_unknown_untracked_inputs_prevent_reuse(tmp_path, monkeypatch, input_path):
    manifest=tmp_path/'tests/test_tier_manifest.json'
    manifest.parent.mkdir()
    manifest.write_text(json.dumps({'profiles':{'solver-heavy':{'paths':[]}}}))
    def files(root, *args):
        return input_path if '--others' in args else 'src/solver.py'
    monkeypatch.setattr(gate, 'git', files)
    with pytest.raises(ValueError, match='untracked validation input'):
        gate.source_fingerprint(tmp_path, {})


@pytest.mark.parametrize('failure', [None, 'schema', 'source_commit', 'source_dirty',
                                     'source_change_fingerprint_sha256', 'binary_name',
                                     'binary_bytes', 'binary_sha256'])
def test_current_native_sidecar_is_hash_and_source_bound(tmp_path, monkeypatch, failure):
    binary=tmp_path/'src/radia/_radia_pybind.pyd'
    binary.parent.mkdir(parents=True)
    data=b'current test binary';binary.write_bytes(data)
    record=dict(schema='radia.native-build-provenance.v1',source_commit='a'*40,
                source_dirty=False,source_change_fingerprint_sha256=gate.digest(b'\0'),
                binary_name=binary.name,binary_bytes=len(data),binary_sha256=gate.digest(data))
    if failure:record[failure]='wrong'
    Path(str(binary)+'.build.json').write_text(json.dumps(record))
    monkeypatch.setattr(gate, 'git', lambda *a:'a'*40)
    if failure:
        with pytest.raises(ValueError, match='native provenance mismatch'):
            gate.validate_native_sidecar(tmp_path)
    else:
        gate.validate_native_sidecar(tmp_path)


@pytest.mark.parametrize('ignored', [False, True])
def test_git_scan_covers_whole_repository_and_ignored_inputs(tmp_path, ignored):
    import subprocess
    subprocess.run(['git','init','-q',str(tmp_path)],check=True)
    if ignored:
        (tmp_path/'.gitignore').write_text('tests/local/\n')
    path=tmp_path/'tests/local/conftest.py'
    path.parent.mkdir(parents=True);path.write_text('fixture code')
    # Real Git, including an ignored file, must return an input outside the old
    # src/validation_test/packages pathspec. No native build is required.
    manifest=tmp_path/'tests/test_tier_manifest.json'
    manifest.write_text(json.dumps({'profiles':{'solver-heavy':{'paths':[]}}}))
    with pytest.raises(ValueError,match='untracked validation input'):
        gate.source_fingerprint(tmp_path,{})
    assert 'tests/local/conftest.py' in gate.git(tmp_path,'ls-files','--others','-z')


def build_environment(tmp_path, monkeypatch, flags=''):
    monkeypatch.setattr(gate,'validated_native_inputs',lambda root:set())
    monkeypatch.setattr(gate,'mkl_build_copy',lambda root:{'package':'mkl','version':'test','sha256':'a'*64})
    monkeypatch.setattr(gate.importlib.metadata,'distributions',lambda:[])
    for key in ('RADIA_TEST_NGSOLVE_THREADS','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
        monkeypatch.setenv(key,'1')
    monkeypatch.setenv('RADIA_TESTS_ALLOW_PARTIAL','0')
    compiler=tmp_path/'build-msvc/CMakeFiles/4.0/CMakeCXXCompiler.cmake'
    compiler.parent.mkdir(parents=True,exist_ok=True)
    compiler.write_text('set(CMAKE_CXX_COMPILER_ID "MSVC")\nset(CMAKE_CXX_COMPILER_VERSION "19.44")\n')
    cache=tmp_path/'build-msvc/CMakeCache.txt'
    cache.write_text('SDK:STRING=10.0.26100.0\nCMAKE_BUILD_TYPE:STRING=Release\n'
                     'CMAKE_CXX_FLAGS:STRING='+flags+'\nCMAKE_CXX_FLAGS_RELEASE:STRING=/O2\n')
    return gate.environment_contract(tmp_path)


@pytest.mark.parametrize('setting', ['/fp:fast','/arch:AVX2'])
def test_compiler_flags_change_environment_fingerprint(tmp_path,monkeypatch,setting):
    before=build_environment(tmp_path,monkeypatch)
    after=build_environment(tmp_path,monkeypatch,setting)
    assert gate.digest(gate.canonical(before)) != gate.digest(gate.canonical(after))


def test_build_type_and_config_flags_are_bound(tmp_path,monkeypatch):
    before=build_environment(tmp_path,monkeypatch)
    cache=tmp_path/'build-msvc/CMakeCache.txt'
    for old,new in [('=Release','=Debug'),('/O2','/O1 /fp:fast')]:
        original=cache.read_text();cache.write_text(original.replace(old,new))
        assert gate.environment_contract(tmp_path)['compiler_options'] != before['compiler_options']
        cache.write_text(original)


@pytest.mark.parametrize('key', ['RADIA_OPTION','NGSOLVE_NUM_THREADS','MKL_CBWR','KMP_AFFINITY','OMP_PROC_BIND'])
def test_numerical_environment_changes_are_bound(tmp_path,monkeypatch,key):
    before=build_environment(tmp_path,monkeypatch)
    monkeypatch.setenv(key,'changed')
    after=gate.environment_contract(tmp_path)
    assert before['numerical_environment'] != after['numerical_environment']
    assert after['numerical_environment'][key]=='changed'


@pytest.mark.parametrize('key', sorted(gate.JOB_CONTEXT_KEYS))
def test_reviewed_job_paths_normalize_only_inside_runner_temp(tmp_path,monkeypatch,key):
    monkeypatch.setenv('RUNNER_TEMP',str(tmp_path))
    monkeypatch.setenv(key,str(tmp_path/'run-one'))
    before=gate.numerical_environment()
    monkeypatch.setenv(key,str(tmp_path/'run-two'))
    assert gate.numerical_environment()==before
    monkeypatch.setenv(key,str(tmp_path.parent/'outside-runner-temp'))
    assert gate.numerical_environment()!=before


def test_other_numerical_paths_are_not_normalized(tmp_path,monkeypatch):
    monkeypatch.setenv('RUNNER_TEMP',str(tmp_path))
    monkeypatch.setenv('RADIA_OTHER_PATH',str(tmp_path/'one'))
    before=gate.numerical_environment()
    monkeypatch.setenv('RADIA_OTHER_PATH',str(tmp_path/'two'))
    assert gate.numerical_environment()!=before


@pytest.mark.parametrize('mode', ['hit','miss','bad-then-hit'])
def test_find_executes_real_filter_and_verification(tmp_path,monkeypatch,mode):
    receipt,run=proof();run['updated_at']=datetime.now(timezone.utc).isoformat()
    data=archive(receipt);meta=artifact(data)
    meta['id']=99
    meta['expires_at']=(datetime.now(timezone.utc)+timedelta(days=1)).isoformat()
    bad=copy.deepcopy(meta);bad['id']=98;bad['workflow_run']={'id':11}
    calls=[]
    class Client(gate.GitHub):
        def get(self,path,binary=False):
            calls.append(path)
            if path.startswith('/actions/artifacts?'):
                return {'artifacts':[bad,meta] if mode=='bad-then-hit' else [meta]}
            if path=='/actions/runs/11':raise OSError('unavailable run')
            if path=='/actions/runs/12':
                return {**run,'head_branch':'feature'} if mode=='miss' else run
            if path=='/actions/artifacts/99/zip':return data
            raise AssertionError(path)
    monkeypatch.setattr(gate.subprocess,'run',lambda *a,**kw:None)
    found=Client('example/repo','unused').find(tmp_path,'fingerprint')
    if mode=='miss':
        assert found is None
        assert not any(p.endswith('/zip') for p in calls)
    else:
        assert found[1]==receipt and found[2]==XML
        if mode=='bad-then-hit':assert '/actions/runs/11' in calls


def test_find_bounds_search_to_thirty_candidates(tmp_path):
    calls=[]
    class Client(gate.GitHub):
        def get(self,path,binary=False):
            if path.startswith('/actions/artifacts?'):
                assert 'per_page=30' in path
                return {'artifacts':[{'workflow_run':{'id':i}} for i in range(35)]}
            calls.append(path)
            raise OSError('untrusted or unavailable')
    assert Client('example/repo','unused').find(tmp_path,'fingerprint') is None
    assert len(calls)==30


def native_tree(tmp_path,monkeypatch):
    monkeypatch.setattr(gate,'git',lambda *a:'a'*40)
    build=tmp_path/'build-msvc';build.mkdir()
    for name in gate.NATIVE_NAMES:
        binary=tmp_path/'src/radia'/name;binary.parent.mkdir(parents=True,exist_ok=True)
        binary.write_bytes(('current '+name).encode())
        (build/name).write_bytes(binary.read_bytes())
        record=dict(schema='radia.native-build-provenance.v1',source_commit='a'*40,
            source_dirty=False,source_change_fingerprint_sha256=gate.digest(b'\0'),
            binary_name=name,binary_bytes=binary.stat().st_size,
            binary_sha256=gate.digest(binary.read_bytes()),build_directory='build-msvc')
        Path(str(binary)+'.build.json').write_text(json.dumps(record))
    return build


def test_validated_build_copies_and_generated_metadata(tmp_path,monkeypatch):
    build=native_tree(tmp_path,monkeypatch)
    allowed=gate.validated_native_inputs(tmp_path)
    assert all('build-msvc/'+name in allowed for name in gate.NATIVE_NAMES)
    metadata=build/'CMakeFiles/4.0/CMakeCXXCompiler.cmake'
    metadata.parent.mkdir(parents=True);metadata.write_text('generated')
    assert gate.generated_build_metadata(tmp_path,metadata.relative_to(tmp_path).as_posix())
    for name in ('extra.pyd','extra.dll','inject.pth','conftest.py'):
        assert not gate.generated_build_metadata(tmp_path,'build-msvc/CMakeFiles/'+name)
    assert not gate.generated_build_metadata(tmp_path,'build-msvc/CMakeFiles/arbitrary-input.h')
    assert not gate.generated_build_metadata(tmp_path,'other/CMakeFiles/generated.txt')


@pytest.mark.parametrize('change',['old-sidecar','different-dir','copied-bytes','missing-copy','different-native'])
def test_unqualified_build_context_cannot_exempt_scratch(tmp_path,monkeypatch,change):
    build=native_tree(tmp_path,monkeypatch)
    name=gate.NATIVE_NAMES[0];sidecar=tmp_path/'src/radia'/f'{name}.build.json'
    record=json.loads(sidecar.read_text())
    if change=='old-sidecar':record.pop('build_directory')
    if change=='different-dir':record['build_directory']='other-build'
    if change=='copied-bytes':(build/name).write_bytes(b'different copied native')
    if change=='missing-copy':(build/name).unlink()
    if change=='different-native':(tmp_path/'src/radia'/gate.NATIVE_NAMES[1]).write_bytes(b'stale')
    sidecar.write_text(json.dumps(record))
    with pytest.raises(ValueError):gate.validated_native_inputs(tmp_path)


def test_unknown_native_inside_declared_build_forces_execution(tmp_path,monkeypatch):
    build=native_tree(tmp_path,monkeypatch)
    (build/'unexpected.pyd').write_bytes(b'unknown module')
    tests=tmp_path/'tests';tests.mkdir()
    (tests/'test_tier_manifest.json').write_text(json.dumps({'profiles':{'solver-heavy':{'paths':[]}}}))
    monkeypatch.setattr(gate,'git',lambda root,*a:'build-msvc/unexpected.pyd' if '--others' in a else 'a'*40)
    with pytest.raises(ValueError,match='untracked validation input'):
        gate.source_fingerprint(tmp_path,{})


@pytest.mark.parametrize('change',[None,'copy','outside-venv','unknown-owner'])
def test_mkl_copy_is_bound_to_installed_venv_package(tmp_path,monkeypatch,change):
    prefix=tmp_path/'venv';source=prefix/'Library/bin/mkl_rt.3.dll'
    source.parent.mkdir(parents=True);source.write_bytes(b'installed MKL runtime')
    copied=tmp_path/'build-msvc/mkl_rt.3.dll'
    copied.parent.mkdir();copied.write_bytes(source.read_bytes())
    class Package:
        version='2026.test'
        files=['Library/bin/mkl_rt.3.dll']
        def locate_file(self,file):return prefix/file
    package=Package()
    monkeypatch.setattr(gate.sys,'prefix',str(prefix))
    monkeypatch.setattr(gate.importlib.metadata,'distribution',lambda name:package)
    if change=='copy':copied.write_bytes(b'unverified runtime')
    if change=='outside-venv':monkeypatch.setattr(gate.sys,'prefix',str(tmp_path/'other-venv'))
    if change=='unknown-owner':package.files=[]
    if change:
        with pytest.raises(ValueError):gate.mkl_build_copy(tmp_path)
    else:
        assert gate.mkl_build_copy(tmp_path)==dict(package='mkl',version='2026.test',sha256=gate.digest(source.read_bytes()))


def test_untracked_rejection_reports_all_paths_including_root_compiler_objects(tmp_path,monkeypatch):
    manifest=tmp_path/'tests/test_tier_manifest.json';manifest.parent.mkdir()
    manifest.write_text(json.dumps({'profiles':{'solver-heavy':{'paths':[]}}}))
    paths=['tests/inject.pth','root.lib','root.pyc']
    monkeypatch.setattr(gate,'git',lambda root,*args:'\0'.join(paths) if '--others' in args else 'src/solver.py')
    with pytest.raises(ValueError) as caught:gate.source_fingerprint(tmp_path,{})
    assert caught.value.paths==paths



def test_planned_fingerprint_keeps_solver_fixture_native_and_environment_inputs(tmp_path, monkeypatch):
    files = {'src/solver.py': 'solver', 'src/viewer.py': 'viewer', 'src/core/native.cpp': 'native',
             'tests/conftest.py': 'fixtures', 'validation_test/case/test.py': 'case',
             'validation_test/case/input.vol': 'mesh'}
    for name, text in files.items():
        path = tmp_path / name; path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding='utf-8')
    (tmp_path / 'tests/test_tier_manifest.json').write_text(json.dumps({'profiles': {
        'solver-heavy': {'paths': ['validation_test/case/test.py']}}}), encoding='utf-8')
    monkeypatch.setattr(gate, 'git', lambda root, *args: '\0'.join(files)
                        if args[0] == 'ls-files' and '--others' not in args else '')
    value = {'cases': [dict(profile='solver-heavy', selected=True, path='validation_test/case/test.py',
                           dependencies=['src/solver.py', 'validation_test/case/test.py',
                                         'validation_test/case/input.vol', 'tests/conftest.py'],
                           components=['solver'], unknowns=[])]}
    env = {'packages': [('ngsolve', '6.2.2608'), ('gmsh', '5.dev1'), ('mcp', '1')],
           'compiler': 'MSVC14.44', 'threads': 2}
    initial = gate.source_fingerprint(tmp_path, env, value)[0]
    (tmp_path / 'src/viewer.py').write_text('viewer changed', encoding='utf-8')
    assert gate.source_fingerprint(tmp_path, env, value)[0] == initial
    other = copy.deepcopy(env); other['packages'][1] = ('gmsh', '5.dev2')
    assert gate.source_fingerprint(tmp_path, other, value)[0] == initial
    for name in ('src/solver.py', 'src/core/native.cpp', 'validation_test/case/input.vol'):
        (tmp_path / name).write_text(files[name] + 'changed', encoding='utf-8')
        assert gate.source_fingerprint(tmp_path, env, value)[0] != initial
        (tmp_path / name).write_text(files[name], encoding='utf-8')
    other['packages'][0] = ('ngsolve', '6.2.next')
    assert gate.source_fingerprint(tmp_path, other, value)[0] != initial
    gmsh_case = copy.deepcopy(value); gmsh_case['cases'][0]['components'].append('gmsh')
    assert gate.source_fingerprint(tmp_path, env, gmsh_case)[0] != gate.source_fingerprint(
        tmp_path, dict(env, packages=[('ngsolve', '6.2.2608'), ('gmsh', '5.dev2'), ('mcp', '1')]), gmsh_case)[0]
    with pytest.raises(ValueError, match='empty planned heavy'):
        gate.source_fingerprint(tmp_path, env, {'cases': []})
