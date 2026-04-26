"""cubit_exec_safely end-to-end: build user state, try AI recipe safely, verify."""
from __future__ import annotations
import subprocess, sys, json
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
subprocess.run(['taskkill', '/F', '/IM', 'coreform_cubit.exe'],
               capture_output=True, check=False)

from radia_mcp.cubit import session as cs
from radia_mcp.cubit import server as srv

cs.CubitSession.reset()
sess = cs.CubitSession.get()
sess.ensure_started()
print("=== 1. Build 'user' state in live GUI (box + sphere) ===")
r = sess.call('cmd', [
    'brick x 10',
    'sphere radius 3',
], timeout_s=60)
for line in r['result']:
    print(f"  {'ok ' if line['ok'] else 'FAIL'}: {line['line']}")
s0 = sess.call('probe', ['summary'])['result']
print(f"  pre-AI: {s0}")

print("\n=== 2. AI proposes a clean recipe (mesh all volumes) ===")
good = ['volume all size 1.0', 'volume all scheme tetmesh', 'mesh volume all']
r = json.loads(srv.cubit_exec_safely(good, timeout_s=600))
print(f"  status: {r['status']}")
print(f"  checkpoint_label: {r['checkpoint_label']}")
print(f"  checkpoint_path: {r['checkpoint_path']}")
print(f"  applied_to_gui: {r['applied_to_gui']}")
if 'live' in r:
    print(f"  live summary: {r['live'].get('summary')}")

print("\n=== 3. AI proposes a broken recipe (non-existent volume) ===")
bad = ['volume 99 size 0.5', 'mesh volume 99']  # volume 99 does not exist
r = json.loads(srv.cubit_exec_safely(bad, timeout_s=300))
print(f"  status: {r['status']}")
print(f"  stage: {r.get('stage')}")
print(f"  checkpoint_label: {r['checkpoint_label']}")
print(f"  applied_to_gui: {r['applied_to_gui']}")
print(f"  note: {r.get('note','')[:200]}")

print("\n=== 4. GUI state is preserved (recipe #2 applied, #3 blocked) ===")
s1 = sess.call('probe', ['summary'])['result']
print(f"  post-safe: {s1}")
print(f"  hexes+tets change vs pre-AI: "
      f"{s1['tets']-s0['tets']} tets, {s1['hexes']-s0['hexes']} hexes")

print("\n=== 5. Checkpoints left on disk (rollback available) ===")
r = json.loads(srv.cubit_list_checkpoints())
labels = [e.get('label','?') for e in r.get('checkpoints', [])]
autos = [l for l in labels if l.startswith('_autosafe_')]
print(f"  total: {r.get('count',0)}")
print(f"  _autosafe_* (produced by exec_safely): {len(autos)}")
for l in autos[:3]:
    print(f"    {l}")
