"""Shared solver acceptance roles; physical identities are private configuration."""
import os

RELEASE_ACCEPTANCE_HOSTS = ("validation", "100")
PUBLIC_RELEASE_ACCEPTANCE_ROLES = {"validation": "validation", "100": "build"}
CURRENT_RELEASE_ACCEPTANCE_HOSTNAMES = {
    role: os.environ.get("RADIA_ACCEPTANCE_HOSTNAME_" + role.upper(), "").strip().lower()
    for role in RELEASE_ACCEPTANCE_HOSTS
}
RELEASE_ACCEPTANCE_HOSTNAMES = {
    role: frozenset({name}) if name else frozenset()
    for role, name in CURRENT_RELEASE_ACCEPTANCE_HOSTNAMES.items()
}
RELEASE_ACCEPTANCE_HOST_LABELS = {"validation": "validation", "100": "build"}
