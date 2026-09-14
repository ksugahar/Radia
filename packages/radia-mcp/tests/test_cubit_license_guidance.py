"""Retired authentication advice must not return through old MCP topic names."""

import pytest

from radia_mcp.cubit.knowledge.license import (
    LICENSE_PER_USER_RULE, LICENSE_2025_12_TOKEN_AUTH, get_license_documentation,
)


@pytest.mark.parametrize("topic", [
    "", "per_user_rule", "per_user", "user_rule", "admin_overwrite",
    "admin", "renewals", "warning",
])
def test_per_user_aliases_preserve_current_guidance(topic):
    text = get_license_documentation(topic)
    assert text == LICENSE_PER_USER_RULE
    assert "official Coreform Cubit shortcut" in text
    assert "is retired" in text
    assert "does not establish" in text


@pytest.mark.parametrize("topic", [
    "token_auth_2025_12", "token_auth", "token", "2025_12", "2025.12",
    "twin_issue", "retry_backoff", "force_login", "portal_retirement",
    "login_tokens", "login_tokens_json",
])
def test_historical_token_aliases_do_not_restore_login_automation(topic):
    text = get_license_documentation(topic)
    assert text == LICENSE_2025_12_TOKEN_AUTH
    assert "not a verified" in text
    assert "rather\nthan prescribe automatic login retries or logout" in text


def test_all_topics_exclude_credential_commands_and_preserve_routing():
    text = get_license_documentation("all")
    assert LICENSE_PER_USER_RULE in text
    assert LICENSE_2025_12_TOKEN_AUTH in text
    for obsolete in ("--login", "--logout", "-ForceLogin", "takeown", "icacls"):
        assert obsolete not in text
    assert get_license_documentation(None) == LICENSE_PER_USER_RULE
    assert get_license_documentation(" TOKEN ") == LICENSE_2025_12_TOKEN_AUTH
    assert "Unknown license topic" in get_license_documentation("missing")
