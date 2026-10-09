"""Privacy source bindings reject numerical and literal changes."""
import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location("privacy_source_contract", Path(__file__).resolve().parents[1] / "validation_test/radia_mcp/privacy_source_contract.py")
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)

def test_cpp_redaction_preserves_literals_and_directives():
    original = '#define K 2\nconst char *s="//literal"; auto x=R"tag(/*literal*/)tag"; double y=1.0; // private runtime\n'
    redacted = original.replace('private runtime', 'execution runtime')
    assert c.cpp_tokens(original) == c.cpp_tokens(redacted)
    assert c.cpp_tokens(original) != c.cpp_tokens(original.replace('1.0', '1.1'))
    assert c.cpp_tokens(original) != c.cpp_tokens(original.replace('//literal', '//changed'))
    assert c.cpp_tokens('#define K 2\nK') != c.cpp_tokens('#define K 2 K')
    assert c.cpp_tokens('x + + y') != c.cpp_tokens('x ++ y')

def test_updated_byte_hash_cannot_authorize_numerical_drift():
    path = 'src/core/rad_hacapk_hdiv.cpp'
    historical = 'double value=1.0; // original comment\n'
    redacted = 'double value=1.0; // neutral comment\n'
    record = {'method': c.AUDITED_PATHS[path], 'historical_text_sha256': c.text_sha256(historical), 'published_text_sha256': c.text_sha256(redacted), 'semantic_sha256': c.semantic_sha256(path, historical)}
    assert c.verifies_redaction(path, c.text_sha256(historical), redacted, record)
    drift = redacted.replace('1.0', '1.1')
    record['published_text_sha256'] = c.text_sha256(drift)
    assert not c.verifies_redaction(path, c.text_sha256(historical), drift, record)
    assert not c.verifies_redaction('unknown.cpp', c.text_sha256(historical), drift, record)

def test_matlab_exact_runtime_pair_only():
    path = 'validation_test/radia_mcp/generate_motor_angle_family_mex_artifact.m'
    old = 'metadata = struct("hostname", hostName); power = 2;\n'
    new = 'metadata = struct("platform_class", string(computer(\'arch\'))); power = 2;\n'
    assert c.semantic_sha256(path, old) == c.semantic_sha256(path, new)
    assert c.semantic_sha256(path, old) != c.semantic_sha256(path, new.replace('power = 2', 'power = 3'))


def test_private_same_host_proof_and_public_payload():
    import pytest
    spec = importlib.util.spec_from_file_location("privacy_assembler", Path(__file__).resolve().parents[1] / "validation_test/optimization/assemble_optuna_release_evidence.py")
    assembler = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(assembler)
    for left, right in (({}, {}), ({"host": "worker-one"}, {}), ({"host": "worker-one"}, {"host": "worker-two"})):
        with pytest.raises(RuntimeError):
            assembler._require_private_same_host(left, right)
    assembler._require_private_same_host({"host": "worker-one"}, {"host": "WORKER-ONE"})
    private = {"host": "worker-one", "runtime": {"hostname": "worker-one", "platform_class": "Windows", "threads": 8}, "power": [1.0, 2.0]}
    public = assembler._public_runtime_payload(private)
    assert public == {"runtime": {"platform_class": "Windows", "threads": 8}, "power": [1.0, 2.0]}
    assert private["host"] == "worker-one"

    assembler._validate_public_runtime_payload(public)
    for raw in ({"note": "192.168.1.2"}, {"nested": [r"W:\private\file"]}, {"details": "mdx1"}, {"nested": [r"\\storage\share"]}):
        with pytest.raises(RuntimeError):
            assembler._validate_public_runtime_payload(raw)
