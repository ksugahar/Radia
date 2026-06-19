"""Token-aware GitHub REST helper (gh-free).

The unauthenticated GitHub API allows only 60 requests/hour per IP, which
release-day / CI-watch polling exhausts in minutes (then every check 403s).
This helper reads a Personal Access Token and uses the authenticated
5000 req/hour limit instead.  With no token it falls back to
unauthenticated (still works for the public ksugahar/Radia repo).

Token source (first hit wins) -- NEVER committed:
  1. $GH_TOKEN
  2. $GITHUB_TOKEN
  3. ~/.radia/gh_token          (gitignored; create it yourself: a file
                                 containing a single classic/fine-grained
                                 PAT with public_repo read scope)
  4. `gh auth token`            (GitHub CLI browser login, OS keyring;
                                 set RADIA_GH to a gh executable if needed)

So a one-time `setx GH_TOKEN ghp_xxx` (or writing ~/.radia/gh_token) lifts
the whole repo's CI tooling -- tools/check_ci.py, release_triple ci-verify,
ad-hoc polling -- to 5000 req/hr.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import urllib.request

_API = "https://api.github.com/"


def token() -> str | None:
    for k in ("GH_TOKEN", "GITHUB_TOKEN"):
        v = os.environ.get(k)
        if v and v.strip():
            return v.strip()
    for p in (os.path.expanduser("~/.radia/gh_token"),
              os.path.expanduser("~/.config/radia/gh_token")):
        try:
            with open(p, encoding="utf-8") as f:
                t = f.read().strip()
            if t:
                return t
        except OSError:
            pass
    gh_token = _gh_cli_token()
    if gh_token:
        return gh_token
    return None


def _gh_cli_token() -> str | None:
    candidates = []
    configured = os.environ.get("RADIA_GH")
    if configured:
        candidates.append(configured)
    found = shutil.which("gh")
    if found:
        candidates.append(found)
    if not candidates:
        return None
    for gh in candidates:
        for args in (
            [gh, "auth", "token", "--hostname", "github.com"],
            [gh, "auth", "token"],
        ):
            try:
                proc = subprocess.run(
                    args,
                    check=False,
                    capture_output=True,
                    text=True,
                    timeout=5,
                )
            except (OSError, subprocess.TimeoutExpired):
                continue
            if proc.returncode == 0 and proc.stdout.strip():
                return proc.stdout.strip()
    return None


def authenticated() -> bool:
    return token() is not None


def gh_get(path: str, timeout: int = 30):
    """GET the GitHub API.  `path` may be a full URL or an api-relative path
    (e.g. "repos/ksugahar/Radia/actions/runs?per_page=8").  Returns
    (parsed_json, headers_dict).  Adds the Authorization header iff a token
    is available."""
    url = path if path.startswith("http") else _API + path.lstrip("/")
    headers = {"Accept": "application/vnd.github+json", "User-Agent": "radia-ci"}
    tok = token()
    if tok:
        headers["Authorization"] = f"Bearer {tok}"
    req = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r), dict(r.headers)


def rate_limit() -> dict:
    """Return the core rate-limit dict {limit, remaining, reset}.  The
    /rate_limit endpoint is itself exempt from the limit."""
    data, _ = gh_get("rate_limit")
    return data["resources"]["core"]
