"""Exercise checked-input failure boundaries without FE or a GPU runtime."""
import importlib.util
import hashlib
import json
import sys
import types
import os
from pathlib import Path

import numpy as np
import pytest


@pytest.fixture
def inputs():
    path=Path(__file__).resolve().parents[1]/"validation_test/hmatrix_cuda/hdiv_validation_inputs.py"
    spec=importlib.util.spec_from_file_location("hdiv_validation_inputs_test",path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def test_oracle_columns_and_shape_mismatch(inputs):
    data=dict(expected=np.array([[1.],[2.]]),rhs=np.array([1.,2.]),system=np.eye(2),mean_operator=np.ones((3,2)))
    expected,_,rhs=inputs.oracle_vectors(data)
    assert expected.shape==rhs.shape==(2,)
    data['system']=np.eye(3)
    with pytest.raises(ValueError,match='shapes mismatch'):inputs.oracle_vectors(data)


def _stub_receipt_verifier(monkeypatch,verify):
    package=types.ModuleType('radia');package.__path__=[]
    module=types.ModuleType('radia.hdiv_gpu');module._verify_mesh_receipt=verify
    package.hdiv_gpu=module
    monkeypatch.setitem(sys.modules,'radia',package)
    monkeypatch.setitem(sys.modules,'radia.hdiv_gpu',module)


def _checked_inputs(tmp_path):
    mesh=tmp_path/'mesh.vol';mesh.write_bytes(b'input')
    contract=tmp_path/'contract.json';contract.write_bytes(b'{}')
    receipt=tmp_path/'receipt.json';receipt.write_bytes(b'{}')
    provenance=tmp_path/'provenance.json';provenance.write_text('{}')
    digest=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
    return mesh,contract,digest(mesh),digest(contract),receipt,digest(receipt),provenance


def test_checked_hash_failure_never_invokes_receipt_verifier(inputs,tmp_path,monkeypatch):
    def forbidden(*args,**kwargs):raise AssertionError('Hash failure must precede receipt verification')
    _stub_receipt_verifier(monkeypatch,forbidden)
    mesh,contract,_,contract_sha,receipt,receipt_sha,provenance=_checked_inputs(tmp_path)
    report=tmp_path/'report.json'
    with pytest.raises(ValueError,match='hash mismatch: mesh'):
        inputs.checked_mesh(mesh,contract,'0'*64,contract_sha,receipt,receipt_sha,provenance,report)
    assert not report.exists()


def test_checked_mesh_rejects_non_tetrahedral_report(inputs,tmp_path,monkeypatch):
    quality=dict(volume_element_count=8,tetrahedron_count=0)
    report_bytes=json.dumps(dict(quality=quality)).encode()
    import base64
    _stub_receipt_verifier(monkeypatch,lambda *args:dict(
        report_base64=base64.b64encode(report_bytes).decode(),checker='stub'))
    args=_checked_inputs(tmp_path);report=tmp_path/'report.json'
    with pytest.raises(ValueError,match='exclusively tetrahedral'):inputs.checked_mesh(*args,report)
    assert not report.exists()
    quality['tetrahedron_count']=8;report_bytes=json.dumps(dict(quality=quality)).encode()
    assert inputs.checked_mesh(*args,report)['report']['quality']==quality
    assert report.read_bytes()==report_bytes


def test_explicit_dependency_directory_required(inputs,tmp_path):
    with pytest.raises(ValueError,match='Explicit'):inputs.dependency_directories([])
    missing=tmp_path/'missing'
    with pytest.raises(FileNotFoundError):inputs.dependency_directories([missing])


def test_cpu_receipt_must_match_loaded_binary(inputs,tmp_path):
    binary=tmp_path/'reference.pyd';binary.write_bytes(b'actual native')
    receipt=tmp_path/'receipt.json'
    receipt.write_text(json.dumps(dict(source_dirty=False,binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest())))
    assert inputs.cpu_receipt(receipt,binary)['source_dirty'] is False
    binary.write_bytes(b'wrong native')
    with pytest.raises(ValueError,match='build receipt'):inputs.cpu_receipt(receipt,binary)


def test_runtime_blas_change_is_rejected(inputs,monkeypatch):
    pools=[dict(user_api='blas',num_threads=1)]
    monkeypatch.setitem(sys.modules,'threadpoolctl',types.SimpleNamespace(threadpool_info=lambda:pools))
    assert inputs.observed_blas()==pools
    pools[0]['num_threads']=18
    with pytest.raises(RuntimeError,match='single-thread'):inputs.observed_blas()


def test_oracle_csr_sorts_but_never_coalesces(inputs):
    from scipy.sparse import csr_matrix
    unordered=csr_matrix((np.array([2.,1.]),np.array([1,0]),np.array([0,2])),shape=(1,2))
    ordered=inputs.sorted_oracle_csr(unordered)
    assert ordered.has_canonical_format
    np.testing.assert_array_equal(ordered.toarray(),unordered.toarray())
    assert not unordered.has_sorted_indices
    duplicates=csr_matrix((np.array([1.,2.]),np.array([0,0]),np.array([0,2])),shape=(1,2))
    with pytest.raises(ValueError,match='duplicate'):inputs.sorted_oracle_csr(duplicates)


@pytest.mark.skipif(os.name!='nt',reason='Windows process memory counters')
def test_process_memory_reports_owned_process_samples():
    path=Path(__file__).resolve().parents[1]/'validation_test/hmatrix_cuda/hdiv_process_memory.py'
    spec=importlib.util.spec_from_file_location('hdiv_memory_test',path)
    memory=importlib.util.module_from_spec(spec);spec.loader.exec_module(memory)
    with memory.ProcessMemory(os.getpid()) as sample:
        assert sample.samples[0]['rss_bytes']>0
    result=sample.result()
    assert result['sample_count']>=2 and result['sampled_peak_private_bytes']>0

@pytest.mark.skipif(os.name!='nt',reason='Windows process affinity')
def test_observed_process_affinity(inputs):
    assert inputs.process_affinity(os.getpid()) > 0


def test_memory_checkpoint_survives_failed_gate(tmp_path):
    import importlib.util
    import json
    path=Path(__file__).resolve().parents[1]/'validation_test/hmatrix_cuda/hdiv_process_memory.py'
    spec=importlib.util.spec_from_file_location('memory_checkpoint',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    output=tmp_path/'samples.jsonl'
    def reject(sample):
        assert json.loads(output.read_text()) == sample
        raise AssertionError('memory gate failed')
    with output.open('x') as stream:
        with pytest.raises(AssertionError,match='memory gate failed'):
            module.checkpoint_sample(stream,dict(iteration=99,host_bytes=12054528),reject)
    assert json.loads(output.read_text())['iteration']==99
