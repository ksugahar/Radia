import importlib.util,json
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[3]
spec=importlib.util.spec_from_file_location('payload_provenance', ROOT/'packages/cubit-mesh-export/src/cubit_mesh_export/_native_provenance.py')
provenance=importlib.util.module_from_spec(spec)
spec.loader.exec_module(provenance)

@pytest.mark.parametrize('sdk', ['2025.12','2026.8'])
def test_record_replaces_previous_sdk_with_configured_sdk(tmp_path,monkeypatch,sdk):
    package=tmp_path/'package'
    package.mkdir()
    (package/'native_payloads.json').write_text(json.dumps({'payloads':{name:{'cubit_version':'2025.12'} for name in provenance.REQUIRED_PAYLOADS}}))
    for name in provenance.REQUIRED_PAYLOADS:
        (package/name).write_bytes(b'candidate')
    for build in ['build-pyd','build-ccm']:
        folder=tmp_path/'src/cubit_plugin'/build
        folder.mkdir(parents=True)
        (folder/'CMakeCache.txt').write_text(f'Cubit_DIR:PATH=C:/Program Files/Coreform Cubit {sdk}/cmake\n')
    monkeypatch.setattr(provenance,'_source_commit',lambda root:'test')
    monkeypatch.setattr(provenance.metadata,'version',lambda name:'6.2.2607')
    manifest=provenance.record_manifest(tmp_path,package)
    assert {payload['cubit_version'] for payload in manifest['payloads'].values()} == {sdk}
