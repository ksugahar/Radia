"""License pre-warm for Coreform Cubit (Learn-edition flow).

Mirrors the user's `public-safe curated corpus` flow:

  1. Read the cached renewals file at
     %LOCALAPPDATA%\\Coreform\\Cubit\\Coreform\\licenses\\renewals
  2. Parse the JSON, find the best remaining-days license.
  3. Skip login if (cache age <= 3 days) AND (best remaining > 7 days).
  4. Otherwise run `rlm_activate.exe --login <email> <pass>` before
     Cubit launch.

**Why**: on a cold license checkout the Cubit daemon's "ready" marker
takes > 60 s to appear (Cubit blocks on RLM server contact).  After
this warmup, the first command usually comes back in 2 – 3 s.

**Credentials**: read EXCLUSIVELY from environment variables
``RADIA_CUBIT_LEARN_EMAIL`` and ``RADIA_CUBIT_LEARN_PASSWORD`` at
call time.  No hardcoded defaults -- if either is unset, warmup is
skipped and we rely on whatever state is already in the Learn
cache.  This package is published to PyPI so it must not embed
real credentials; lab-side setup is responsible for setting the
env vars (e.g. via the launcher script that wraps this code).

**Fail-open**: any failure here MUST NOT block Cubit launch.  The
warmup is an optimisation; if it can't run (rlm_activate missing,
RLM server unreachable, credentials wrong), fall through silently
and let Cubit do its own license checkout.
"""
from __future__ import annotations

import datetime as _dt
import json as _json
import os
import re
import subprocess
from pathlib import Path
from typing import Optional


# Credentials are read from env vars only.  Do NOT hard-code defaults
# in this file -- this package is published to PyPI as radia-mcp and
# any value here would be world-public the moment a release is cut.
# Lab-side launchers are responsible for setting RADIA_CUBIT_LEARN_*
# before invoking the warmup (see the Cubit launcher .ps1 / SKILL).
_CACHE_AGE_DAYS = 3      # re-login if cache older than this
_MIN_REMAINING_DAYS = 7  # re-login if best license has fewer days left


def _renewals_path() -> Path:
    base = os.environ.get("LOCALAPPDATA")
    if not base:
        return Path()  # no LOCALAPPDATA on this platform; warmup a no-op
    return Path(base) / "Coreform" / "Cubit" / "Coreform" / "licenses" / "renewals"


def _credentials() -> tuple[Optional[str], Optional[str]]:
    email = os.environ.get("RADIA_CUBIT_LEARN_EMAIL")
    password = os.environ.get("RADIA_CUBIT_LEARN_PASSWORD")
    if not email or not password:
        return None, None
    return email, password


def _needs_login() -> tuple[bool, dict]:
    """Return (needs_login, diagnostic_info).

    Follows the ps1 logic: parse renewals cache, check cache age vs
    `_CACHE_AGE_DAYS`, and best-license expiry vs `_MIN_REMAINING_DAYS`.
    """
    info: dict = {"renewals_path": None, "cache_age_days": None,
                  "best_remaining_days": None}
    p = _renewals_path()
    if not p or not p.is_file():
        info["reason"] = "no renewals cache"
        return True, info
    info["renewals_path"] = str(p)

    try:
        raw = p.read_bytes().decode("utf-8", errors="replace")
    except OSError:
        info["reason"] = "cannot read renewals cache"
        return True, info

    m = re.search(r'\{"Valid".*', raw, re.DOTALL)
    if not m:
        info["reason"] = "renewals cache has no Valid section"
        return True, info
    try:
        data = _json.loads(m.group(0))
    except _json.JSONDecodeError:
        info["reason"] = "renewals Valid JSON unparseable"
        return True, info

    valid = data.get("Valid", {})
    today = _dt.date.today()

    age_days = 999
    cached_today = valid.get("today")
    if cached_today:
        try:
            cd = _dt.datetime.strptime(cached_today, "%Y-%m-%d").date()
            age_days = (today - cd).days
        except ValueError:
            pass
    info["cache_age_days"] = age_days

    best_remaining = 0
    for lic in valid.get("licenses", []) or []:
        exp = lic.get("expiration")
        if not exp:
            continue
        try:
            ed = _dt.datetime.strptime(exp, "%Y-%m-%d").date()
            days = (ed - today).days
            if days > best_remaining:
                best_remaining = days
        except ValueError:
            pass
    info["best_remaining_days"] = best_remaining

    if age_days <= _CACHE_AGE_DAYS and best_remaining > _MIN_REMAINING_DAYS:
        info["reason"] = (f"cache fresh ({age_days}d old, "
                          f"{best_remaining}d remaining)")
        return False, info
    info["reason"] = (f"cache stale or expiring "
                      f"(age={age_days}d, remain={best_remaining}d)")
    return True, info


def warmup_license(bin_dir: Path, timeout_s: float = 30.0,
                    verbose: bool = False, force: bool = True) -> dict:
    """Pre-warm Cubit Learn license before launching Cubit.

    `force=True` (default, 2026-04-25): always run `--logout` then `--login`
    to clear any stale signed-renewals state that Cubit's runtime would
    reject. Empirical observation: skipping rlm_activate when the cache
    looks fresh is intermittently followed by a license-warning modal
    dialog in GUI mode (which blocks `-input` playback). Forcing logout
    + login eliminates the dialog.

    `force=False`: legacy behavior — skip rlm_activate when cache is
    fresh (`age <= _CACHE_AGE_DAYS` and `best_remaining > _MIN_REMAINING_DAYS`).

    Returns a dict with keys:
      status : "ok" | "skipped" | "error"
      reason : human-readable detail
      info   : diagnostic from _needs_login() (cache age etc.)
      action : "login_ran" | "login_skipped" | "no_rlm_activate" | ...

    `bin_dir`: Cubit `bin/` directory (contains rlm_activate.exe).

    Never raises.  On any error, returns status="error" and the
    caller (CubitSession launch) continues without warmup.
    """
    result = {"status": "skipped", "reason": "", "info": {},
              "action": "login_skipped"}

    needs_login, info = _needs_login()
    result["info"] = info
    if not force and not needs_login:
        result["reason"] = info.get("reason", "cache fresh")
        return result

    email, password = _credentials()
    if not email or not password:
        result["reason"] = "no credentials (RADIA_CUBIT_LEARN_EMAIL not set)"
        result["action"] = "login_skipped"
        return result

    rlm_exe = Path(bin_dir) / "rlm_activate.exe"
    if not rlm_exe.is_file():
        result["status"] = "error"
        result["reason"] = f"rlm_activate.exe not found in {bin_dir}"
        result["action"] = "no_rlm_activate"
        return result

    # Force logout first to release any stuck server-side seat / clear
    # stale renewals signature. Logout failure is non-fatal — proceed.
    if force:
        try:
            subprocess.run(
                [str(rlm_exe), "--logout"],
                capture_output=True, text=True, timeout=timeout_s,
            )
        except (subprocess.TimeoutExpired, OSError):
            pass

    try:
        proc = subprocess.run(
            [str(rlm_exe), "--login", email, password],
            capture_output=True, text=True, timeout=timeout_s,
            # Detached so any stdin doesn't hang the parent
        )
    except subprocess.TimeoutExpired:
        result["status"] = "error"
        result["reason"] = f"rlm_activate timeout after {timeout_s}s"
        result["action"] = "login_timeout"
        return result
    except OSError as exc:
        result["status"] = "error"
        result["reason"] = f"rlm_activate launch failed: {exc}"
        result["action"] = "login_launch_failed"
        return result

    if proc.returncode == 0:
        result["status"] = "ok"
        result["reason"] = ("forced logout+login OK" if force
                            else info.get("reason",
                                          "login succeeded (cache was stale)"))
        result["action"] = "login_ran"
    else:
        result["status"] = "error"
        result["reason"] = (f"rlm_activate exit {proc.returncode}: "
                            f"{proc.stdout[-200:] or proc.stderr[-200:]}")
        result["action"] = "login_failed"
    return result
