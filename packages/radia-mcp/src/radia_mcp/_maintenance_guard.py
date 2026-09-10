"""Host-wide cooperative exclusion and fail-closed maintenance receipts.

This does not intercept direct pip calls or non-cooperating file writers.
The lock file is permanent; unresolved operations survive process death.
"""

from __future__ import annotations

import json
import os
import re
import socket
import sys
import tempfile
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

STATE_ROOT = (
    Path(r"C:\temp") if os.name == "nt" else Path(tempfile.gettempdir())
) / "radia-mcp-maintenance"


def identity(path: Path) -> str:
    return os.path.normcase(str(path.resolve(strict=True)))


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _atomic(path: Path, record: dict) -> None:
    temporary = path.with_name(path.name + "." + uuid.uuid4().hex + ".tmp")
    try:
        with temporary.open("x", encoding="utf-8") as stream:
            json.dump(record, stream, ensure_ascii=True, indent=2, allow_nan=False)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


@contextmanager
def _exclusive():
    # A fixed host-wide lock deliberately serializes even different interpreter
    # generations: the same client config can refer to either interpreter.
    STATE_ROOT.mkdir(parents=True, exist_ok=True)
    with (STATE_ROOT / "deployment.lock").open("a+b") as lock:
        lock.seek(0, os.SEEK_END)
        if lock.tell() == 0:
            lock.write(b"\0")
            lock.flush()
        lock.seek(0)
        try:
            if os.name == "nt":
                import msvcrt

                msvcrt.locking(lock.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl

                fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            raise ValueError("Another maintenance owner holds the deployment lock") from exc
        try:
            yield
        finally:
            lock.seek(0)
            if os.name == "nt":
                msvcrt.locking(lock.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


def _state() -> dict:
    path = STATE_ROOT / "state.json"
    if not path.exists():
        return {}
    value = json.loads(path.read_text(encoding="utf-8"))
    if (
        not isinstance(value, dict)
        or value.get("schema") != "radia-mcp.activation.v1"
        or not re.fullmatch(r"[0-9a-f]{32}", str(value.get("change_id", "")))
        or not isinstance(value.get("scope"), dict)
    ):
        raise ValueError("Invalid maintenance state; reconcile manually")
    return value


def _save(record: dict) -> None:
    # Publish state first: a crash between the two writes stays fail-closed.
    _atomic(STATE_ROOT / "state.json", record)
    receipts = STATE_ROOT / "receipts"
    receipts.mkdir(exist_ok=True)
    _atomic(receipts / (record["change_id"] + ".json"), record)


def _authorization(owner: str, reason: str, clients_idle: bool) -> None:
    if not owner.strip() or not reason.strip():
        raise ValueError("Explicit owner and reason are required")
    if clients_idle is not True:
        raise ValueError("All affected clients must be confirmed idle before activation")


@contextmanager
def change(
    *, owner: str, reason: str, clients_idle: bool, scope: dict, expected: dict, proposed: dict
):
    _authorization(owner, reason, clients_idle)
    with _exclusive():
        previous = _state()
        if previous and previous.get("status") not in {"completed", "reconciled"}:
            raise ValueError("Unresolved maintenance change; explicit reconciliation required")
        record = {
            "schema": "radia-mcp.activation.v1",
            "change_id": uuid.uuid4().hex,
            "previous_change_id": previous.get("change_id"),
            "owner": owner,
            "reason": reason,
            "host": socket.gethostname(),
            "interpreter": identity(Path(sys.executable)),
            "pid": os.getpid(),
            "scope": scope,
            "expected": expected,
            "proposed": proposed,
            "started_at": _now(),
            "status": "running",
            "clients_idle": "operator-attested; not automatically discovered",
            "live_status": "unverified",
            "client_reconnect": "not-performed",
            "client_results": [
                {"target": target.strip(), "status": "unverified", "observed_at": None}
                for target in scope.get("targets", "").split(",")
                if target.strip()
            ],
        }
        _save(record)
        try:
            yield record
        except BaseException as exc:
            record.update(status="failed", finished_at=_now(), error_type=type(exc).__name__)
            _save(record)  # Never persist exception text / configuration secrets.
            raise
        else:
            record.update(status="completed", finished_at=_now())
            _save(record)


def reconcile(change_id: str, *, owner: str, reason: str, clients_idle: bool, observe) -> dict:
    """Acknowledge an audited interrupted operation; never retry or roll back."""
    _authorization(owner, reason, clients_idle)
    with _exclusive():
        previous = _state()
        if previous.get("change_id") != change_id or previous.get("status") not in {
            "running",
            "failed",
        }:
            raise ValueError("No matching unresolved change")
        if previous.get("interpreter") != identity(Path(sys.executable)):
            raise ValueError("Reconcile using the original target interpreter")
        observed = observe(previous["scope"])
        # Retain the failed receipt, including a process-death 'running' record.
        receipts = STATE_ROOT / "receipts"
        receipts.mkdir(exist_ok=True)
        _atomic(receipts / (previous["change_id"] + ".json"), previous)
        result = {
            **previous,
            "change_id": uuid.uuid4().hex,
            "previous_change_id": change_id,
            "owner": owner,
            "reason": reason,
            "status": "reconciled",
            "observed": observed,
            "finished_at": _now(),
            "live_status": "unverified",
            "client_reconnect": "not-performed",
        }
        _save(result)
        return result
