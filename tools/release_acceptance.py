"""Shared solver acceptance roles; physical identities are private configuration.

Exact-wheel release acceptance runs on the validation role only. The CPU build
host prepares the release candidate and is not an acceptance target.
"""
import os

RELEASE_ACCEPTANCE_HOSTS = ("validation",)
PUBLIC_RELEASE_ACCEPTANCE_ROLES = {"validation": "validation"}
CURRENT_RELEASE_ACCEPTANCE_HOSTNAMES = {
    role: os.environ.get("RADIA_ACCEPTANCE_HOSTNAME_" + role.upper(), "").strip().lower()
    for role in RELEASE_ACCEPTANCE_HOSTS
}
RELEASE_ACCEPTANCE_HOSTNAMES = {
    role: frozenset({name}) if name else frozenset()
    for role, name in CURRENT_RELEASE_ACCEPTANCE_HOSTNAMES.items()
}
