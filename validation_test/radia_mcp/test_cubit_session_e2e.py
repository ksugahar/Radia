"""Live-Cubit session E2E validation for the radia_mcp.cubit MCP stack.

Promoted from C:/temp/cubit_probe_selftest (2026-08-05) per the
promotion ladder: these are the executable truths behind the MathWorks
pattern port (waves 1-6) -- the GUI file-drop path, multi-client
concurrency, the batch stdio path, and the cubit<->build123d probe
contract on one geometry.

Driving-policy note (Sugahara, 2026-08-05): APREPRO/Python on the
HEADLESS route is the primary agent path; the GUI session is the USER's
visual-debugging aid.  The GUI tests here exist to protect that
user-facing debug surface (including cubit_snapshot, which needs a
rendering window), not to promote GUI driving for automation.

Requirements: a local Coreform Cubit install (a license seat is consumed
briefly per test; sessions use RADIA_CUBIT_SESSION_MODE=new so the
shared per-user daemon is never touched, and every session is shut down
in a finally block).  Skipped wholesale when Cubit is absent.
"""

import json
import os
import threading
import time

import pytest
from radia_mcp.cubit import session as cs

pytestmark = pytest.mark.skipif(
    cs.find_cubit_install() is None,
    reason="Coreform Cubit not installed")


@pytest.fixture()
def gui_session(monkeypatch):
    """A fresh private-daemon GUI session, always shut down."""
    monkeypatch.setenv("RADIA_CUBIT_SESSION_MODE", "new")
    sess = cs.CubitSession(mode="gui")
    try:
        yield sess
    finally:
        try:
            sess.shutdown()
        except Exception:
            pass


def test_gui_full_stack(gui_session):
    """Spawn -> diagnostics -> cmd -> shared probes -> verified snapshot
    -> journal -> ownership shutdown -> private-dir cleanup."""
    sess = gui_session
    info = sess.ensure_started()
    drop = sess._drop_dir
    assert info.get("protocol_version") == 2
    assert sess._owned is True
    assert drop.name != "cubit-session"          # private dir (mode=new)
    assert not (drop / "startup_error.txt").exists()

    r = sess.call("cmd", ["reset", "create brick x 10",
                          "create cylinder radius 3 height 20",
                          "move volume 2 x 12",
                          "volume all size 2", "mesh volume all"])
    assert r["ok"], r

    bad = sess.call("cmd", ["draw volume 99"])
    per = bad.get("result") or []
    assert per and per[-1]["ok"] is False
    assert per[-1]["error_count_delta"] >= 1     # silent-error detection

    smy = sess.call("probe", ["summary"])["result"]
    ent = sess.call("probe", ["entities"])["result"]
    pv = sess.call("probe", ["per_volume"])["result"]
    assert smy["volumes"] == 2 and smy["hexes"] > 0
    assert len(ent["volumes"]) == 2
    assert all("bbox_min" in v and "extent" in v for v in ent["volumes"])
    assert len(pv) == 2 and all(row["meshed"] for row in pv)

    png = str(drop / "e2e_view.png")
    snap = sess.call("snapshot", [png])["result"]
    assert snap.get("ok") and snap.get("bytes", 0) > 1000, snap
    assert os.path.getsize(png) == snap["bytes"]

    hist = sess._command_history
    assert any(h["line"] == "create brick x 10" and h["ok"] for h in hist)
    assert any(h["line"] == "draw volume 99" and not h["ok"] for h in hist)

    report = sess.shutdown()
    assert report["stopped"] == "owned-child"
    assert not drop.exists()                     # private dir removed


def test_gui_two_client_concurrency(gui_session):
    """One private daemon, two attached clients, concurrent probes with
    zero cross-talk; state changes propagate between clients."""
    a = gui_session
    a.ensure_started()
    a.call("cmd", ["reset", "create brick x 10", "mesh volume all"])

    b = cs.CubitSession.__new__(cs.CubitSession)
    b._bin_dir = a._bin_dir
    b._mode = "gui"
    b._proc = None
    b._next_id = 1
    b._lock = threading.Lock()
    b._ready_info = dict(a._ready_info or {"ready": True})
    b._drop_dir = a._drop_dir
    b._outbox = a._outbox
    b._owned = False
    b._job_handle = None
    b._last_license_warmup = {}
    b._command_history = []
    b._command_history_max = 100

    errors, results = [], []

    def worker(sess, tag, n=10):
        for k in range(n):
            try:
                r = sess.call("probe", ["summary"], timeout_s=30.0,
                              _recover=False)
                assert r["ok"] and r["result"]["volumes"] == 1, r
                results.append(tag)
            except Exception as exc:  # noqa: BLE001
                errors.append(f"{tag}#{k}: {type(exc).__name__}: {exc}")

    threads = [threading.Thread(target=worker, args=(s, f"{nm}{i}"))
               for nm, s in (("A", a), ("B", b)) for i in range(3)]
    for t in threads:
        t.start()
    deadline = time.monotonic() + 120.0
    for t in threads:
        t.join(timeout=max(0.0, deadline - time.monotonic()))
    assert not [t.name for t in threads if t.is_alive()]
    assert not errors, errors[:5]
    assert len(results) == 60

    b.call("cmd", ["create sphere radius 2"])
    smy = a.call("probe", ["summary"])["result"]
    assert smy["volumes"] == 2                   # B's change visible to A


def test_batch_probe_stack():
    """Batch stdio daemon: entities/labels probes + label audit +
    silent block no-op detection (the Cubit 2025.12 phantom-block
    hazard)."""
    sess = cs.CubitSession(mode="batch")
    try:
        ready = sess.ensure_started()
        assert ready.get("ready"), ready

        r = sess.call("cmd", ["reset", "create brick x 10",
                              "volume all size 2", "mesh volume all",
                              'block 1 add volume 1', 'block 1 name "iron"'])
        assert r["ok"], r

        ent = sess.call("probe", ["entities"])["result"]
        assert ent["volume_count"] == 1
        v = ent["volumes"][0]
        assert v["extent"] == [10.0, 10.0, 10.0]

        lab = sess.call("probe", ["labels"])["result"]
        assert lab["blocks"][0]["name"] == "iron"
        assert lab["blocks"][0]["volume_elems"] > 0
        assert lab["audit"]["passed"] is True

        # Phantom-block hazard: wrong-kind add silently no-ops
        sess.call("cmd", ['block 2 add tri in surface 1'])
        lab2 = sess.call("probe", ["labels"])["result"]
        ids = [b["id"] for b in lab2["blocks"]]
        if 2 in ids:
            b2 = next(b for b in lab2["blocks"] if b["id"] == 2)
            assert b2["surface_elems"] == 0      # membership absent
    finally:
        sess.shutdown()


def test_cad_mesh_probe_contract(tmp_path):
    """The cross-server probe contract on one geometry: a labeled
    build123d assembly, probed on the CAD side and meshed+probed in the
    Cubit batch daemon, must agree per body in the SAME vocabulary."""
    pytest.importorskip("build123d")
    from build123d import Box, Compound, Cylinder, Pos, export_step
    from radia_mcp.build123d.server import build123d_probe
    from radia_mcp.common.server_hardening import PROBE_SOLID_CORE_KEYS

    step = tmp_path / "compat_asm.step"
    box = Box(1, 1, 2)
    box.label = "yoke"
    cyl = Pos(3, 0, 0) * Cylinder(0.4, 1.0)
    cyl.label = "core"
    export_step(Compound(children=[box, cyl]), str(step))

    cad = json.loads(build123d_probe(str(step), query="entities"))
    assert cad["status"] == "ok" and cad["solid_count"] == 2

    sess = cs.CubitSession(mode="batch")
    try:
        sess.ensure_started()
        fwd = str(step).replace("\\", "/")
        r = sess.call("cmd", ["reset", f'import step "{fwd}"',
                              "volume all size 0.2", "mesh volume all"])
        assert r["ok"], r
        mesh = sess.call("probe", ["entities"])["result"]
    finally:
        sess.shutdown()

    for s in cad["solids"]:
        assert PROBE_SOLID_CORE_KEYS <= set(s)
    for v in mesh["volumes"]:
        assert PROBE_SOLID_CORE_KEYS <= set(v)

    def closest(mesh_vols, centroid):
        return min(mesh_vols, key=lambda v: sum(
            (x - y) ** 2 for x, y in zip(v["centroid"], centroid)))

    for s in cad["solids"]:
        m = closest(mesh["volumes"], s["centroid"])
        rel = abs(m["volume"] - s["volume"]) / s["volume"]
        assert rel < 5e-3, (s["label"], s["volume"], m["volume"])
