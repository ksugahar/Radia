"""Build headers stay byte-exact and independent of the installed SDK."""
import hashlib
import importlib.util
from pathlib import Path


def test_header_snapshot_is_exact_and_does_not_share_storage(tmp_path):
    spec=importlib.util.spec_from_file_location('native_builder',Path(__file__).parents[1]/'tools/build_csr_python.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    source=tmp_path/'installed-sdk';header=source/'pybind11/detail/common.h';header.parent.mkdir(parents=True)
    original=b'// public synthetic header\n#define CONTRACT_VERSION 3\n';header.write_bytes(original)
    destination=tmp_path/'job/dependency-headers/netgen'
    hashes=module.snapshot_headers(source,destination)
    assert hashes=={'pybind11/detail/common.h':hashlib.sha256(original).hexdigest()}
    copied=destination/'pybind11/detail/common.h';assert copied.read_bytes()==original
    copied.write_bytes(b'changed owned snapshot')
    assert header.read_bytes()==original


def test_logical_mapping_preserves_code_and_original_line_numbers(tmp_path):
    spec=importlib.util.spec_from_file_location('source_paths',Path(__file__).parents[1]/'tools/native_source_paths.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    original=tmp_path/'original';original.mkdir();header=original/'nested.h'
    body=b'#pragma once\nstatic int contract(){return __LINE__;}\n';header.write_bytes(body)
    destination=tmp_path/'compiled'
    manifest=module.stage_logical_sources([header],original,destination,'ngsolve/include')
    assert header.read_bytes()==body
    assert (destination/'nested.h').read_bytes()==b'#line 1 "ngsolve/include/nested.h"\n'+body
    assert manifest['nested.h']['original_sha256']==hashlib.sha256(body).hexdigest()


def test_resident_binding_stages_immutable_sdk_and_project_diagnostics(tmp_path, monkeypatch):
    import json
    tools = Path(__file__).parents[1]/'tools'
    monkeypatch.syspath_prepend(str(tools))
    spec = importlib.util.spec_from_file_location('resident_source_builder', tools/'build_resident_amg_python.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    sdk = tmp_path/'private-installed-sdk'
    header = sdk/'pybind11/detail/common.h'; header.parent.mkdir(parents=True)
    body = b'#pragma once\nstatic const char* diagnostic_file = __FILE__;\n'
    header.write_bytes(body)
    package = tmp_path/'private-checkout/packages/radia-simulink'; package.mkdir(parents=True)
    native = package.parents[1]/'src/core/gpu'; native.mkdir(parents=True)
    project = {'resident_amg_python.cpp': b'\xef\xbb\xbf#include "resident_amg_host.hpp"\n',
               'resident_amg_host.hpp': b'#pragma once\nstatic_assert(sizeof(int)>=2);\n'}
    for name, data in project.items(): (native/name).write_bytes(data)
    out = tmp_path/'private-build'; out.mkdir()
    headers, source, provenance = module.stage_binding_sources(package, out, sdk)
    assert headers.is_relative_to(out) and source.is_relative_to(out)
    assert header.read_bytes() == body
    assert (out/'dependency-headers-original/netgen/pybind11/detail/common.h').read_bytes() == body
    assert (headers/'pybind11/detail/common.h').read_bytes() == b'#line 1 "ngsolve/include/pybind11/detail/common.h"\n'+body
    for name, data in project.items():
        assert (native/name).read_bytes() == data
        bom = b'\xef\xbb\xbf' if data.startswith(b'\xef\xbb\xbf') else b''
        assert (out/'logical-native'/name).read_bytes() == bom+f'#line 1 "radia/native/{name}"\n'.encode()+data[len(bom):]
    for key, digest in provenance.items():
        record = out/(key.removesuffix('_sha256').replace('_', '-')+'.json')
        assert hashlib.sha256(record.read_bytes()).hexdigest() == digest
    manifest = json.loads((out/'logical-header-paths.json').read_bytes())
    assert manifest['pybind11/detail/common.h']['original_sha256'] == hashlib.sha256(body).hexdigest()
