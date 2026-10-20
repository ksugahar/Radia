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
