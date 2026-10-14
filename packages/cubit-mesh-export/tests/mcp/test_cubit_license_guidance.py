"""Retired authentication advice must not return through old MCP topic names."""

import pytest

from cubit_mesh_export.mcp.knowledge.license import (
    LICENSE_PER_USER_RULE, LICENSE_2025_12_TOKEN_AUTH, get_license_documentation,
)


@pytest.mark.parametrize("topic", ["", "per_user"])
def test_per_user_topic_preserves_current_guidance(topic):
    text = get_license_documentation(topic)
    assert text == LICENSE_PER_USER_RULE
    assert "official Coreform Cubit shortcut" in text
    assert "is retired" in text
    assert "does not establish" in text


def test_token_cache_topic_does_not_restore_login_automation():
    text = get_license_documentation("token_cache")
    assert text == LICENSE_2025_12_TOKEN_AUTH
    assert "not a verified" in text
    assert "official Coreform UI" in text


def test_all_topics_exclude_credential_commands_and_preserve_routing():
    text = get_license_documentation("all")
    assert LICENSE_PER_USER_RULE in text
    assert LICENSE_2025_12_TOKEN_AUTH in text
    for obsolete in ("--login", "--logout", "-ForceLogin", "takeown", "icacls"):
        assert obsolete not in text
    assert get_license_documentation(None) == LICENSE_PER_USER_RULE
    assert get_license_documentation(" TOKEN_CACHE ") == LICENSE_2025_12_TOKEN_AUTH
    assert "Unknown license topic" in get_license_documentation("missing")
    for retired in ("force_login", "retry_backoff", "portal_retirement",
                    "admin_overwrite", "token_auth_2025_12"):
        assert "Unknown license topic" in get_license_documentation(retired)
