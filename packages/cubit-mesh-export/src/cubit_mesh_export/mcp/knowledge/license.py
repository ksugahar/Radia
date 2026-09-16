"""Current laboratory Cubit activation guidance, not an authentication client."""

LICENSE_PER_USER_RULE = """
# Cubit activation belongs to the intended Windows user

## Current operating policy

Use the official Coreform Cubit shortcut and its login/license UI in the
intended user's own Windows session. The former CoreformCubit.ps1 launcher
is retired; do not recreate it or deploy a replacement credential launcher.
An administrator's successful SSH probe does not establish that another
user's interactive Cubit session is licensed.

Keep credentials and tokens out of chat, command arguments, source files,
logs, and shared folders. Do not copy another user's credential/cache files,
rewrite their ownership, or delete them as a routine repair. If authentication
fails, ask the intended user to complete the official login UI; retain the
non-secret error and consult vendor support if it persists.

## Read-only diagnosis

Separate these observations:
- Installation and plugin files are present.
- The intended account can authenticate through the official UI.
- An authorized command actually completes and its output passes its checks.

A file's presence, age, owner, or a successful version probe alone proves
neither an active seat nor successful mesh export. Inspect only necessary
metadata; do not print token/credential contents. Report which Windows
account and which execution route was actually checked.

GUI startup and GUI tests belong only to cubit-mesh-export on LAB and
100号機. The bundled Cubit MCP supports AI-driven headless operation.
Radia solver workflows normally consume checked .vol files. Batch generation
uses APREPRO or Cubit's Python API through cubit-mesh-export; Radia CI
consumes fixtures and does not launch Cubit's GUI.

Preserve human-owned sessions. A license diagnosis does not authorize
logging users out, releasing a seat, terminating Cubit, or reinstalling
plugins. Perform state-changing recovery only as an explicitly coordinated
operation, then verify it in the affected user's context.
"""

LICENSE_2025_12_TOKEN_AUTH = """
# Cubit 2025.12 cache observations and current recovery boundary

Historical LAB observations found login_tokens.json and cubit_creds under
the current user's local Coreform application-data directory. Older
investigations referenced a nested licenses/renewals cache. These are
version-specific observations, not a stable API or an activation recipe.

Do not infer token expiry from file modification time. The previously
suggested 7-day refresh rule and 10-14-day token lifetime were not a verified
vendor contract and must not drive automation. Likewise, an HTTP 500 is not
proof of duplicate-token prevention or a reason to force logout. Record the
actual non-secret error, runtime version, and affected account.

Authentication is handled through the official Coreform UI by the intended
user. Do not revive the retired custom launcher, embed passwords in scripts,
or copy cache files.
Do not promise that a cached seat remains valid after an authentication error.

Consult current vendor guidance for the installed Cubit version when a
UI login does not resolve the issue. The historical portal-retirement date
does not establish the status of a service today. Confirm actual command
completion separately from license metadata and preserve active sessions.
"""


def get_license_documentation(topic: str = "per_user") -> str:
    """Return current license guidance by its canonical topic name."""
    topic = (topic or "").lower().strip()
    if topic == "all":
        return (
            f"# ===== Per-user activation rule =====\n\n"
            f"{LICENSE_PER_USER_RULE}\n\n"
            f"# ===== 2025.12 token-based auth =====\n\n"
            f"{LICENSE_2025_12_TOKEN_AUTH}"
        )
    if topic == "token_cache":
        return LICENSE_2025_12_TOKEN_AUTH
    if topic in ("", "per_user"):
        return LICENSE_PER_USER_RULE
    return (
        f"Unknown license topic '{topic}'.  Available: per_user, "
        f"token_cache, all."
    )
