"""Provenance checks shared by legacy and current evidence contracts."""
import importlib

import pytest

from radia.ltspice._artifact_identity import (
    digest_is_sha256,
    legacy_digest_is_sha256,
    same_generation,
)


@pytest.mark.parametrize('value,expected', [
    ('a' * 64, True), ('ABCDEF01' * 8, True), ('0' * 63, False),
    ('0' * 65, False), ('g' * 64, False), (' ' + '0' * 63, False),
    ('', False), (None, False), (False, False), (b'a' * 64, False),
    (int('1' * 64), False),
])
def test_digest_requires_exact_hex_string(value, expected):
    assert digest_is_sha256(value) is expected


def test_old_coercion_does_not_leak_into_new_contract():
    value = int('1' * 64)
    assert legacy_digest_is_sha256(value)
    assert not digest_is_sha256(value)


@pytest.mark.parametrize('contract,expected', [
    ({'generation_id': 'run', 'source': 'run', 'result': 'run'}, True),
    ({'generation_id': 'run', 'source': 'other', 'result': 'run'}, False),
    ({'generation_id': 'run', 'source': 'run'}, False),
    ({'generation_id': '', 'source': '', 'result': ''}, False),
    ({'generation_id': 12, 'source': '12', 'result': '12'}, True),
    ({'generation_id': 12, 'source': 12, 'result': 12}, False),
])
def test_generation_requires_nonempty_matching_fields(contract, expected):
    assert same_generation(contract, 'source', 'result') is expected


@pytest.mark.parametrize('version', range(43, 56))
def test_numbered_gates_use_shared_predicates(version):
    module = importlib.import_module(f'radia.ltspice.ltspice_v{version}_gates')
    if version <= 44:
        assert module._sha256 is legacy_digest_is_sha256
    else:
        assert module._digest is digest_is_sha256
    if version >= 48:
        assert module._generation is same_generation
