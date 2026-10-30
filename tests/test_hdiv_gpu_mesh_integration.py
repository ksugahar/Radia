"""Guards: receipt content/options/provenance fail before FE/GPU entry and relocation.

Assumes: reports are synthetic transport fixtures, not checker or GPU acceptance.
The build-host oracle separately executes the authoritative installed checker.
"""
import base64
import copy
import gzip
import hashlib
import importlib.util
import json
import sys
import types
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]

def api():
    spec=importlib.util.spec_from_file_location("hdiv_receipt_contract",ROOT/"src/radia/hdiv_gpu.py")
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def fixture(tmp_path,compressed=False):
    data=(ROOT/"tests/fixtures/hdiv/hex_pair.vol").read_bytes()
    if compressed:data=gzip.compress(data,mtime=0)
    mesh=tmp_path/("relocated.vol.gz" if compressed else "relocated.vol");mesh.write_bytes(data)
    contract=tmp_path/"caller-labels.json"
    labels=dict(schema="radia.vol-label-contract.v1",required=dict(materials=["body"],boundaries=["exterior"]),allowed=dict(materials=["body"],boundaries=["exterior"]))
    contract.write_text(json.dumps(labels))
    producer="C:/producer/input/mesh.vol"+(".gz" if compressed else "")
    report=dict(schema="cubit-mesh-export.vol-check.v1",passed=True,warnings=[],vol_file=producer,json_file=None,
        cad_reference=dict(available=False,auto_discovered=False),
        labels=dict(schema=labels['schema'],passed=True,strict_labels=True,warnings=[],required=labels['required'],allowed=labels['allowed']),
        quality=dict(passed=True,warnings=[]),boundary_domain_ownership=dict(issues=[]),
        versions=dict(cubit_mesh_export="2.2.6",ngsolve="6.2.2608"))
    checker=dict(distribution="cubit-mesh-export",version="2.2.6",record_sha256="b"*64)
    value=dict(schema="radia.checked-vol-receipt.v1",checker=checker,ngsolve_version="6.2.2608",python_version="recorded-only",
        argv=["C:/producer/python.exe","-m","cubit_mesh_export.check",producer,"--strict-labels","--contract","C:/producer/input/labels.json","--format","json"],
        returncode=0,mesh_sha256=hashlib.sha256(data).hexdigest(),contract_sha256=hashlib.sha256(contract.read_bytes()).hexdigest())
    receipt=tmp_path/"check.receipt.json"
    return mesh,contract,receipt,value,report


def pins(receipt,value,report):
    raw=json.dumps(report).encode();value['report_base64']=base64.b64encode(raw).decode();value['report_sha256']=hashlib.sha256(raw).hexdigest()
    receipt.write_text(json.dumps(value))
    return dict(mesh_check_receipt=receipt,mesh_check_receipt_sha256=hashlib.sha256(receipt.read_bytes()).hexdigest(),
        mesh_check_provenance=copy.deepcopy({k:value[k] for k in ('checker','argv','ngsolve_version')}))


@pytest.mark.parametrize('compressed',[False,True])
def test_valid_receipt_relocates_and_preserves_compressed_snapshot(tmp_path,monkeypatch,compressed):
    module=api();mesh,contract,receipt,value,report=fixture(tmp_path,compressed);options=pins(receipt,value,report)
    monkeypatch.setattr(module.metadata,'version',lambda name:'6.2.2608')
    expected=mesh.read_bytes()
    def initialize(self,snapshot,law,order):
        assert snapshot.name==('mesh.vol.gz' if compressed else 'mesh.vol')
        assert snapshot.read_bytes()==expected
    monkeypatch.setattr(module.PreparedHDivNonlinear,'_initialize',initialize)
    module.PreparedHDivNonlinear(mesh,{},backend='gpu',label_contract=contract,**options)


@pytest.mark.parametrize('mutation',['mesh','contract','report','options','checker','quality','strict','warnings','digest','cad','abi','provenance','duplicate','nonfinite'])
def test_bad_receipt_cannot_reach_workflow_native_or_fe(tmp_path,monkeypatch,mutation):
    module=api();mesh,contract,receipt,value,report=fixture(tmp_path)
    if mutation=='quality':report['quality']['passed']=False
    if mutation=='strict':report['labels']['strict_labels']=False
    if mutation=='warnings':report['warnings']=['finding']
    if mutation=='cad':report['cad_reference']['available']=True
    if mutation=='options':value['argv']+=['--no-quality']
    if mutation=='abi':value['ngsolve_version']='0.0.0'
    options=pins(receipt,value,report)
    if mutation=='mesh':mesh.write_bytes(b'modified')
    if mutation=='contract':contract.write_bytes(b'{}')
    if mutation=='report':
        value['report_base64']=base64.b64encode(b'{}').decode();receipt.write_text(json.dumps(value));options['mesh_check_receipt_sha256']=hashlib.sha256(receipt.read_bytes()).hexdigest()
    if mutation=='checker':options['mesh_check_provenance']['checker']['record_sha256']='c'*64
    if mutation=='digest':options['mesh_check_receipt_sha256']='0'*64
    if mutation=='provenance':options['mesh_check_provenance']={}
    if mutation in ('duplicate','nonfinite'):
        receipt.write_text('{"schema":"radia.checked-vol-receipt.v1","schema":"again"}' if mutation=='duplicate' else '{"value":NaN}')
        options['mesh_check_receipt_sha256']=hashlib.sha256(receipt.read_bytes()).hexdigest()
    monkeypatch.setattr(module.metadata,'version',lambda name:'6.2.2608')
    touched=[]
    monkeypatch.setattr(module,'_prepare_host_components',lambda *a,**kw:touched.append('FE-staging'))
    def native_constructor(*args, **kwargs):
        touched.append('native-constructor')
        raise AssertionError('Native constructor before receipt')
    native = types.SimpleNamespace(_ResidentHDiv=native_constructor,
                                   _ResidentBHMaterial=native_constructor)
    def native_import(name):
        touched.append(name)
        raise AssertionError('Native import before receipt')
    gpu = types.ModuleType('radia.gpu')
    gpu._native = native_import
    monkeypatch.setitem(sys.modules, 'radia.gpu', gpu)
    monkeypatch.setitem(sys.modules, '_hdiv_cuda', native)
    import builtins
    original=builtins.__import__
    def guarded(name,*args,**kwargs):
        if name in ('ngsolve','radia.gpu','_hdiv_cuda') or (name == 'gpu' and '_native' in (args[2] if len(args)>2 else ())):
            touched.append(name);raise AssertionError('Workflow import before receipt')
        return original(name,*args,**kwargs)
    monkeypatch.setattr(builtins,'__import__',guarded)
    with pytest.raises(ValueError):module.PreparedHDivNonlinear(mesh,{},backend='gpu',label_contract=contract,**options)
    with pytest.raises(ValueError):module.PreparedHDivNonlinear.from_bh_table(mesh,[],backend='gpu',label_contract=contract,**options)
    with pytest.raises(ValueError):module.prepare_hdiv_host_setup(mesh,backend='gpu',label_contract=contract,**options)
    assert touched==[],touched
