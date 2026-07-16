from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import urllib.error

import pytest


_SCRIPT = Path(__file__).resolve().parents[1] / "tools" / "upload_release_asset.py"
if str(_SCRIPT.parent) not in sys.path:
    sys.path.insert(0, str(_SCRIPT.parent))
_SPEC = importlib.util.spec_from_file_location("upload_release_asset", _SCRIPT)
assert _SPEC is not None and _SPEC.loader is not None
upload_release_asset = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(upload_release_asset)


def _http_error(code: int) -> urllib.error.HTTPError:
    return urllib.error.HTTPError(
        "https://api.github.com/test", code, "temporary", {}, None
    )


def test_urlopen_retries_transient_http_error(monkeypatch):
    response = object()
    outcomes = [_http_error(503), _http_error(502), response]
    sleeps = []

    def fake_urlopen(_req, *, timeout):
        assert timeout == 30
        outcome = outcomes.pop(0)
        if isinstance(outcome, Exception):
            raise outcome
        return outcome

    monkeypatch.setattr(upload_release_asset.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(upload_release_asset.time, "sleep", sleeps.append)

    assert upload_release_asset._urlopen_with_retry(object(), timeout=30) is response
    assert sleeps == [5, 10]


def test_urlopen_does_not_retry_nontransient_http_error(monkeypatch):
    calls = 0

    def fake_urlopen(_req, *, timeout):
        nonlocal calls
        calls += 1
        raise _http_error(404)

    monkeypatch.setattr(upload_release_asset.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(
        upload_release_asset.time,
        "sleep",
        lambda _delay: pytest.fail("404 must not be retried"),
    )

    with pytest.raises(urllib.error.HTTPError) as exc_info:
        upload_release_asset._urlopen_with_retry(object(), timeout=30)
    assert exc_info.value.code == 404
    assert calls == 1
