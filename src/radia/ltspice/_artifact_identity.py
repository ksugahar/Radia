"""Shared artifact provenance predicates, independent of numbered gate APIs."""
from collections.abc import Mapping


def digest_is_sha256(value: object) -> bool:
    """Require a hexadecimal string; accept either letter case."""
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def legacy_digest_is_sha256(value: object) -> bool:
    """Preserve v43/v44 string coercion without weakening later contracts."""
    return digest_is_sha256(str(value or ""))


def same_generation(contract: Mapping[str, object], *names: str) -> bool:
    """Require every named field to equal the nonempty generation identifier."""
    generation = str(contract.get("generation_id") or "")
    return bool(generation) and all(contract.get(name) == generation for name in names)
