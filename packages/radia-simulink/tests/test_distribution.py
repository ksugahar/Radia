"""Guard the installed resource contract without requiring a MATLAB license."""
import hashlib
import importlib.util
import json
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]


def load_api():
    spec=importlib.util.spec_from_file_location("gpu_contract_api",ROOT/"src/radia_simulink/__init__.py")
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize('mutation', [
    'none', 'notice-missing', 'ngsolve-license-missing', 'netgen-license-missing',
    'license-truncated', 'license-modified', 'source-link-missing',
    'dependency-missing', 'dependency-unpinned', 'dependency-duplicate',
    'dependency-not-gpu', 'solver-abi-mismatch', 'crlf-license',
])
def test_dependency_only_distribution_audit(tmp_path, mutation):
    """Audit bytes in the archive, not a receipt claiming licensing acceptance."""
    import tomllib
    import zipfile
    from email.message import Message
    spec = importlib.util.spec_from_file_location('dependency_audit', ROOT/'verify_wheel.py')
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    metadata = Message()
    gpu = tomllib.loads((ROOT/'pyproject.toml').read_text(encoding='utf-8'))[
        'project']['optional-dependencies']['gpu']
    for requirement in gpu:
        if requirement.startswith('nvidia-cusparse-cu12'):
            if mutation == 'dependency-missing':
                continue
            if mutation == 'dependency-unpinned':
                requirement = requirement.replace('==', '>=')
            if mutation == 'dependency-duplicate':
                metadata['Requires-Dist'] = requirement + '; extra == "gpu"'
            if mutation == 'dependency-not-gpu':
                metadata['Requires-Dist'] = requirement + '; extra == "test"'
                continue
        if requirement.startswith('netgen-mesher') and mutation == 'solver-abi-mismatch':
            requirement = 'netgen-mesher==0.0.0'
        metadata['Requires-Dist'] = requirement + (' and ' if ';' in requirement else '; ') + 'extra == "gpu"'
    files = {name: (ROOT/name).read_bytes() for name in (
        'THIRD_PARTY_NOTICES.md', 'licenses/NGSOLVE-LICENSE.txt', 'licenses/NETGEN-LICENSE.txt')}
    if mutation == 'notice-missing':
        files.pop('THIRD_PARTY_NOTICES.md')
    if mutation == 'ngsolve-license-missing':
        files.pop('licenses/NGSOLVE-LICENSE.txt')
    if mutation == 'netgen-license-missing':
        files.pop('licenses/NETGEN-LICENSE.txt')
    if mutation == 'license-truncated':
        files['licenses/NGSOLVE-LICENSE.txt'] = files['licenses/NGSOLVE-LICENSE.txt'][:100]
    if mutation == 'license-modified':
        files['licenses/NETGEN-LICENSE.txt'] += b'changed'
    if mutation == 'source-link-missing':
        files['THIRD_PARTY_NOTICES.md'] = files['THIRD_PARTY_NOTICES.md'].replace(
            b'https://github.com/NGSolve/netgen', b'https://example.invalid/source')
    if mutation == 'crlf-license':
        for name in ('licenses/NGSOLVE-LICENSE.txt', 'licenses/NETGEN-LICENSE.txt'):
            files[name] = files[name].replace(b'\r\n', b'\n').replace(b'\n', b'\r\n')
    wheel = tmp_path/'dependency-audit.whl'
    with zipfile.ZipFile(wheel, 'w') as archive:
        for name, value in files.items():
            archive.writestr('radia_simulink/' + name, value)
    with zipfile.ZipFile(wheel) as archive:
        if mutation in ('none', 'crlf-license'):
            report = verifier.audit_dependency_distribution(archive, metadata)
            assert report['passed'] and not report['dependency_dlls_bundled']
            assert report['full_lgpl21_texts_checked']
        else:
            with pytest.raises(ValueError):
                verifier.audit_dependency_distribution(archive, metadata)


def test_installed_manifest_rejects_tampering_and_path_escape(tmp_path,monkeypatch):
    api=load_api()
    monkeypatch.setattr(api,"files",lambda package:tmp_path)
    (tmp_path/"matlab").mkdir()
    payload=tmp_path/"matlab/kernel.m"
    payload.write_text("checked payload")
    manifest=dict(schema="radia.gpu.distribution.v1",version=api.__version__,source_commit="a"*40,
                  capabilities=["test payload"],runtime_requirements=[],
                  files={"matlab/kernel.m":hashlib.sha256(payload.read_bytes()).hexdigest()})
    (tmp_path/"manifest.json").write_text(json.dumps(manifest))
    checked = api.verify()
    assert checked["passed"]
    assert checked['manifest_sha256'] == hashlib.sha256((tmp_path/'manifest.json').read_bytes()).hexdigest()
    payload.write_text("changed payload")
    with pytest.raises(RuntimeError,match="hash mismatch"):
        api.verify()
    manifest["files"]={"../escape":"0"*64}
    (tmp_path/"manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(RuntimeError,match="invalid payload path"):
        api.verify()


def test_package_keeps_solver_runtime_dependencies_optional_and_lazy():
    import tomllib
    metadata=tomllib.loads((ROOT/"pyproject.toml").read_text())
    assert metadata["project"]["dependencies"]==[]
    extras=metadata["project"]["optional-dependencies"]
    assert set(extras["exchange"]).issubset(extras["test"])
    api=load_api()
    assert "ngsolve" not in api.__dict__ and "matlab" not in api.__dict__
    assert "checkout" not in api.matlab_path.__doc__.lower()


@pytest.mark.parametrize('member',['../escape','radia_simulink/matlab/radia_mex.mexw64',
                                  'radia_simulink/matlab/external.dll',
                                  'radia_simulink/matlab/external.DLL','radia/__init__.py'])
def test_wheel_rejects_escape_cpu_solver_and_external_payload(tmp_path,member):
    import zipfile
    spec=importlib.util.spec_from_file_location('gpu_wheel_contract',ROOT/'verify_wheel.py')
    verifier=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    wheel=tmp_path/'radia_simulink-0.1.0a1-py3-none-win_amd64.whl'
    with zipfile.ZipFile(wheel,'w') as archive:
        archive.writestr(member,b'forbidden')
    with pytest.raises(ValueError):
        verifier.verify(wheel)


def test_installed_acceptance_rejects_a_different_build_of_the_same_source():
    spec=importlib.util.spec_from_file_location('gpu_installed_identity',ROOT/'verify_wheel.py')
    verifier=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    wheel = dict(source_commit='a'*40, version='0.1.0', manifest_sha256='b'*64)
    verifier.verify_installed_identity(dict(wheel), wheel)
    installed = dict(wheel, manifest_sha256='c'*64)
    with pytest.raises(ValueError, match='differs from release wheel'):
        verifier.verify_installed_identity(installed, wheel)
def test_runtime_and_package_versions_agree():
    import ast
    import tomllib
    metadata = tomllib.loads((ROOT/'pyproject.toml').read_text())
    module = ast.parse((ROOT/'src/radia_simulink/__init__.py').read_text())
    assignment = next(node for node in module.body if isinstance(node, ast.Assign)
                      and any(isinstance(target, ast.Name) and target.id=='__version__' for target in node.targets))
    assert ast.literal_eval(assignment.value) == metadata['project']['version']


@pytest.mark.parametrize('component,field', [('hmatrix','fault_injection'),('projection','fault_injection'),('projection','lifecycle_diagnostic')])
@pytest.mark.parametrize('mode', [False, True, None, 0, 'False'])
def test_hdiv_mex_wheel_requires_explicit_false_diagnostic_mode(component, field, mode):
    spec=importlib.util.spec_from_file_location('hdiv_mex_receipt_verifier',ROOT/'verify_wheel.py')
    verifier=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    manifest=dict(source_commit='a'*40,files={
        'matlab/radia_cuda_hmatrix.mexw64':'b'*64,
        'matlab/radia_cuda_projection.mexw64':'c'*64})
    artifacts={name:dict(source_dirty=False,source_commit='a'*40,
        taskmanager_started=False,fault_injection=False,
        binary_sha256=manifest['files']['matlab/radia_cuda_'+name+'.mexw64'])
        for name in ('hmatrix','projection')}
    artifacts['projection']['lifecycle_diagnostic']=False
    if mode is None:
        artifacts[component].pop(field)
    else:
        artifacts[component][field]=mode
    receipt=dict(schema='radia.gpu.hdiv-payload.v1',artifacts=artifacts)
    if mode is False:
        verifier.verify_hdiv_receipt(manifest,receipt)
    else:
        with pytest.raises(ValueError,match='HDiv CUDA binary receipt mismatch'):
            verifier.verify_hdiv_receipt(manifest,receipt)





@pytest.mark.parametrize('wheel_tag,native_tag,accepted',[
    ('cp312-cp312','cp312',True),('py3-none','cp312',False),
    ('cp313-cp313','cp312',False),('cp312-cp311','cp312',False)])
def test_python_cuda_wheel_requires_exact_native_abi(tmp_path,wheel_tag,native_tag,accepted):
    import zipfile
    spec=importlib.util.spec_from_file_location('wheel_native_contract',ROOT/'verify_wheel.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    name=f'_csr_cuda.{native_tag}-win_amd64.pyd'
    wheel=tmp_path/f'radia_simulink-0.1.5-{wheel_tag}-win_amd64.whl'
    digest=hashlib.sha256(b'packaging-only synthetic binary').hexdigest()
    manifest=dict(source_commit='a'*40,files={name:digest})
    receipt=dict(schema='radia.csr-python-build.v1',source_commit='a'*40,binary=name,
                 binary_sha256=digest,matlab_engine_user_dependency=False)
    with zipfile.ZipFile(wheel,'w') as z:
        z.writestr('radia_simulink/'+name,b'packaging-only synthetic binary')
        z.writestr('radia_simulink/python-build.json',json.dumps(receipt))
    with zipfile.ZipFile(wheel) as z:
        if accepted:module.verify_python_binding(wheel,z.namelist(),z,manifest)
        else:
            with pytest.raises(ValueError,match='ABI'):
                module.verify_python_binding(wheel,z.namelist(),z,manifest)


def test_native_payload_uses_platform_library_root():
    # Load the actual native-distribution policy without importing the build backend.
    import ast
    spec = importlib.util.spec_from_file_location("native_payload_policy", ROOT / "native_payload.py")
    policy = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(policy)
    assert policy.NativePayloadDistribution().has_ext_modules() is True
    assert "include native_payload.py" in (ROOT / "MANIFEST.in").read_text().splitlines()
    # setup.py must apply that exact policy to the backend Distribution.
    tree = ast.parse((ROOT / "setup.py").read_text())
    distribution = next(node for node in tree.body
                        if isinstance(node, ast.ClassDef) and node.name == "BinaryDistribution")
    assert [ast.unparse(base) for base in distribution.bases] == [
        "payload_module_native().NativePayloadDistribution", "Distribution"]
    assert not any(isinstance(node, ast.FunctionDef) and node.name == "has_ext_modules"
                   for node in distribution.body)


@pytest.mark.parametrize('mutation', ['none', 'abi', 'wheel-abi', 'hash', 'source', 'engine', 'receipt-missing'])
def test_csr_wheel_requires_matching_abi_and_build_receipt(tmp_path, mutation):
    """Guards: CSR adapter identity. Assumes: synthetic packaging bytes, no runtime qualification."""
    import zipfile
    spec=importlib.util.spec_from_file_location('wheel_csr_contract', ROOT/'verify_wheel.py')
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    wheel_abi='cp313' if mutation=='wheel-abi' else 'cp312'
    wheel=tmp_path/f'radia_simulink-0.1.5-cp312-{wheel_abi}-win_amd64.whl'
    abi='cp313' if mutation=='abi' else 'cp312'
    name=f'_csr_cuda.{abi}-win_amd64.pyd'
    payload=b'packaging-only synthetic csr'
    digest=hashlib.sha256(payload).hexdigest()
    manifest=dict(source_commit='a'*40, files={name:digest})
    receipt=dict(schema='radia.csr-python-build.v1', source_commit='a'*40,
                 binary=name, binary_sha256=digest, matlab_engine_user_dependency=False)
    if mutation=='hash': receipt['binary_sha256']='b'*64
    if mutation=='source': receipt['source_commit']='c'*40
    if mutation=='engine': receipt['matlab_engine_user_dependency']=True
    with zipfile.ZipFile(wheel, 'w') as archive:
        archive.writestr('radia_simulink/'+name,payload)
        if mutation!='receipt-missing':
            archive.writestr('radia_simulink/python-build.json',json.dumps(receipt))
    with zipfile.ZipFile(wheel) as archive:
        if mutation=='none': module.verify_python_binding(wheel,archive.namelist(),archive,manifest)
        else:
            with pytest.raises(ValueError):
                module.verify_python_binding(wheel,archive.namelist(),archive,manifest)


@pytest.mark.parametrize('companion', ['none', 'csr'])
def test_simulink_wheel_rejects_radia_owned_hdiv_python_binding(tmp_path, companion):
    """Guards: one-way package ownership. Assumes: synthetic native bytes, not solver evidence."""
    import zipfile
    spec=importlib.util.spec_from_file_location('wheel_hdiv_boundary', ROOT/'verify_wheel.py')
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    wheel=tmp_path/'radia_simulink-0.1.5-cp312-cp312-win_amd64.whl'
    manifest=dict(source_commit='a'*40,files={})
    with zipfile.ZipFile(wheel,'w') as archive:
        for component in (('csr','hdiv') if companion=='csr' else ('hdiv',)):
            name=f'_{component}_cuda.cp312-win_amd64.pyd';payload=b'synthetic packaging bytes'
            digest=hashlib.sha256(payload).hexdigest();manifest['files'][name]=digest
            receipt=dict(schema=f'radia.{component}-python-build.v1',source_commit='a'*40,
                         binary=name,binary_sha256=digest,matlab_engine_user_dependency=False)
            archive.writestr('radia_simulink/'+name,payload)
            record='python-build.json' if component=='csr' else 'hdiv-python-build.json'
            archive.writestr('radia_simulink/'+record,json.dumps(receipt))
    with zipfile.ZipFile(wheel) as archive:
        with pytest.raises(ValueError,match='Python CUDA binding wheel ABI mismatch'):
            module.verify_python_binding(wheel,archive.namelist(),archive,manifest)
    with pytest.raises(ValueError,match='CPU solver or external executable'):
        module.verify(wheel)


@pytest.mark.parametrize('binary,receipt', [(True, False), (False, True)])
def test_resident_wheel_rejects_unpaired_payload(tmp_path, binary, receipt):
    import zipfile
    spec = importlib.util.spec_from_file_location('resident_pair_contract', ROOT/'verify_wheel.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    wheel = tmp_path/'radia_simulink-0.1.5-cp312-cp312-win_amd64.whl'
    name = '_resident_amg.cp312-win_amd64.pyd'
    with zipfile.ZipFile(wheel, 'w') as archive:
        if binary:
            archive.writestr('radia_simulink/'+name, b'synthetic packaging fixture')
        if receipt:
            archive.writestr('radia_simulink/resident-amg-build.json', '{}')
    with zipfile.ZipFile(wheel) as archive:
        with pytest.raises(ValueError, match='Resident AMG'):
            module.verify_python_binding(wheel, archive.namelist(), archive,
                                         dict(source_commit='a'*40, files={}))


@pytest.mark.parametrize('companion', ['pure', 'csr'])
def test_resident_orphan_receipt_cannot_hide_in_other_wheel(tmp_path, companion):
    import zipfile
    spec = importlib.util.spec_from_file_location('resident_orphan_contract', ROOT/'verify_wheel.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    tag = 'py3-none' if companion == 'pure' else 'cp312-cp312'
    wheel = tmp_path/f'radia_simulink-0.1.5-{tag}-win_amd64.whl'
    name = '_csr_cuda.cp312-win_amd64.pyd'
    manifest = dict(source_commit='a'*40, files={name: 'b'*64})
    with zipfile.ZipFile(wheel, 'w') as archive:
        archive.writestr('radia_simulink/resident-amg-build.json', '{}')
        if companion == 'csr':
            archive.writestr('radia_simulink/'+name, b'synthetic CSR packaging fixture')
            archive.writestr('radia_simulink/python-build.json', json.dumps(dict(
                schema='radia.csr-python-build.v1', source_commit='a'*40,
                binary=name, binary_sha256='b'*64, matlab_engine_user_dependency=False)))
    with zipfile.ZipFile(wheel) as archive:
        with pytest.raises(ValueError, match='Resident AMG receipt requires'):
            module.verify_python_binding(wheel, archive.namelist(), archive, manifest)
