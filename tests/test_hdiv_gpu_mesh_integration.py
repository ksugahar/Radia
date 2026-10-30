"""Guards: real checked mesh snapshots preserve content and compressed filenames.

Assumes: the checker package requires its pinned Netgen runtime; native lane only.
"""
import importlib.util
import json
from pathlib import Path
import pytest
ROOT = Path(__file__).resolve().parents[1]
API = "src/radia/hdiv_gpu.py"

def load(relative):
    spec = importlib.util.spec_from_file_location("hdiv_mesh_integration", ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("compressed", [False, True])
def test_real_label_checker_is_filename_neutral_and_preserves_compression(tmp_path,monkeypatch,compressed):
    import gzip,hashlib,os
    module=load(API)
    contents=(ROOT/"tests/fixtures/hdiv/hex_pair.vol").read_bytes()
    if compressed:contents=gzip.compress(contents,mtime=0)
    path=tmp_path/("caller-name.vol.gz" if compressed else "caller-name.vol");path.write_bytes(contents)
    contract=tmp_path/"caller-labels.json"
    contract.write_text(json.dumps(dict(schema="radia.vol-label-contract.v1",
        required=dict(materials=["body"],boundaries=["exterior"]),
        allowed=dict(materials=["body"],boundaries=["exterior"]))))
    monkeypatch.setenv("PYTHONPATH",str(ROOT/"packages/cubit-mesh-export/src")+os.pathsep+os.environ.get("PYTHONPATH",""))
    def initialize(self,snapshot,law,order):
        assert snapshot.name == ("mesh.vol.gz" if compressed else "mesh.vol")
        assert snapshot.read_bytes() == contents
        assert self.mesh_sha256 == hashlib.sha256(contents).hexdigest()
    monkeypatch.setattr(module.PreparedHDivNonlinear,"_initialize",initialize)
    module.PreparedHDivNonlinear(path,{},backend="gpu",label_contract=contract)
