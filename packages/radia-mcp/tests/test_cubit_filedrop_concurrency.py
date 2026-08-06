"""Multi-client concurrency stress for the GUI file-drop transport.

No Cubit needed: a fake responder thread plays the bootstrap.py role
(scan *.req, echo a response into out/<stem>.resp).  What is stressed
is the CLIENT side under concurrency: two attached CubitSession clients
(distinct `_request_stem` client ids -- the collision codex fixed) fire
interleaved calls from multiple threads each; every response must reach
the caller that sent the request with its own payload (no cross-talk),
and the drop dir must end clean (no request/response debris).
"""

import json
import os
import threading
import time

from radia_mcp.cubit import session as cubit_session


class FakeBootstrapResponder(threading.Thread):
    """Minimal stand-in for bootstrap.py's QTimer poll loop."""

    def __init__(self, drop, outbox):
        super().__init__(daemon=True)
        self.drop = drop
        self.outbox = outbox
        self.stop_flag = threading.Event()
        self.handled = 0

    def run(self):
        while not self.stop_flag.is_set():
            for req_path in sorted(self.drop.glob("*.req")):
                stem = req_path.stem
                try:
                    req = json.loads(req_path.read_text(encoding="utf-8"))
                except (OSError, json.JSONDecodeError):
                    continue
                payload = {"id": req.get("id"), "ok": True,
                           "result": {"echo": req.get("args"),
                                      "op": req.get("op")}}
                tmp = self.outbox / f"{stem}.resp.tmp"
                dst = self.outbox / f"{stem}.resp"
                tmp.write_text(json.dumps(payload), encoding="utf-8")
                os.replace(tmp, dst)
                req_path.unlink(missing_ok=True)
                self.handled += 1
            time.sleep(0.002)


def _attached_client(drop, outbox):
    sess = cubit_session.CubitSession.__new__(cubit_session.CubitSession)
    sess._bin_dir = drop
    sess._mode = "gui"
    sess._proc = None
    sess._next_id = 1
    sess._lock = threading.Lock()
    sess._ready_info = {"ready": True, "protocol_version": 2,
                        "pid": os.getpid()}
    sess._drop_dir = drop
    sess._outbox = outbox
    sess._owned = False
    sess._job_handle = None
    sess._last_license_warmup = {}
    sess._command_history = []
    sess._command_history_max = 100
    return sess


def test_two_clients_concurrent_no_crosstalk(tmp_path):
    drop = tmp_path / "cubit-session"
    outbox = drop / "out"
    outbox.mkdir(parents=True)
    (drop / "pid.lock").write_text(str(os.getpid()), encoding="utf-8")
    (drop / "ready").write_text(json.dumps(
        {"ready": True, "protocol_version": 2, "pid": os.getpid()}),
        encoding="utf-8")

    responder = FakeBootstrapResponder(drop, outbox)
    responder.start()
    try:
        clients = [_attached_client(drop, outbox) for _ in range(2)]
        assert clients[0]._request_stem(1) != clients[1]._request_stem(1)

        errors = []
        n_threads, n_calls = 4, 25

        def worker(client_idx, thread_idx):
            sess = clients[client_idx]
            for k in range(n_calls):
                token = f"c{client_idx}-t{thread_idx}-{k}"
                try:
                    resp = sess.call("probe", [token], timeout_s=20.0,
                                     _recover=False)
                    assert resp["ok"], resp
                    assert resp["result"]["echo"] == [token], resp
                except Exception as exc:  # noqa: BLE001
                    errors.append(f"{token}: {type(exc).__name__}: {exc}")

        threads = [threading.Thread(target=worker, args=(ci, ti))
                   for ci in range(2) for ti in range(n_threads)]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=120)
        assert not errors, errors[:5]
        total = 2 * n_threads * n_calls
        assert responder.handled >= total
    finally:
        responder.stop_flag.set()
        responder.join(timeout=10)

    # transport debris check: no unconsumed requests or responses
    assert list(drop.glob("*.req")) == []
    assert list(drop.glob("*.req.tmp")) == []
    assert list(outbox.glob("*.resp")) == []


def test_response_id_mismatch_fails_loud(tmp_path):
    """A response carrying the WRONG id must raise, not be delivered."""
    import pytest

    drop = tmp_path / "cubit-session"
    outbox = drop / "out"
    outbox.mkdir(parents=True)
    sess = _attached_client(drop, outbox)

    def rogue_responder():
        for _ in range(400):
            reqs = list(drop.glob("*.req"))
            if reqs:
                stem = reqs[0].stem
                (outbox / f"{stem}.resp").write_text(
                    json.dumps({"id": 999999, "ok": True, "result": "x"}),
                    encoding="utf-8")
                return
            time.sleep(0.005)

    t = threading.Thread(target=rogue_responder, daemon=True)
    t.start()
    with pytest.raises(cubit_session.CubitSessionError) as ei:
        sess._call_via_filedrop({"id": 1, "op": "ping", "args": [],
                                 "protocol_version": 2}, timeout_s=10.0)
    assert "id mismatch" in str(ei.value)
