"""Push C++ build artifacts (.pyd) from LAB to mdx via base64-over-ssh.

mdx has no MSVC + Intel MKL + NGSolve build environment, so it cannot
build the C++ extensions itself. The standard `download_binaries.sh`
GitHub-Releases route doesn't work either because mdx has no `gh` CLI
installed. This tool pushes the LAB-built `.pyd` files directly via
base64-over-ssh.

Run AFTER every `Build.ps1` on LAB whose output mdx should see.
Idempotent and safe to re-run.

The deploy skill (`.claude/skills/deploy/SKILL.md` Stage 2 mdx) calls
this from "Step B".
"""

from __future__ import annotations

import base64
import subprocess
import sys
from pathlib import Path

LAB_REPO = Path(r"S:\Radia\01_GitHub")
MDX_USER_HOST = "mdx"

# (LAB source path, list of remote destinations on mdx)
PYDS = [
    (LAB_REPO / "src" / "radia" / "_radia_pybind.pyd",
     [r"C:\Radia\01_GitHub\src\radia\_radia_pybind.pyd"]),
    (LAB_REPO / "src" / "radia" / "cln_core.pyd",
     [r"C:\Radia\01_GitHub\src\radia\cln_core.pyd"]),
    (LAB_REPO / "src" / "radia" / "peec_matrices.pyd",
     [r"C:\Radia\01_GitHub\src\radia\peec_matrices.pyd"]),
    # CRITICAL: cubit_mesh_curver.pyd lives in TWO paths on the local clone.
    # LAB Build.ps1 writes both; on mdx (no build) we mirror both.
    (LAB_REPO / "packages" / "cubit-mesh-export" / "src" / "cubit_mesh_export"
        / "cubit_mesh_curver.pyd",
     [r"C:\Radia\01_GitHub\src\radia\cubit_mesh_curver.pyd",
      r"C:\Radia\01_GitHub\packages\cubit-mesh-export\src\cubit_mesh_export"
      r"\cubit_mesh_curver.pyd"]),
]


def push_one(src: Path, dst: str) -> bool:
    if not src.exists():
        print(f"  [SKIP] {src} (LAB has not built this .pyd)")
        return False
    b64 = base64.b64encode(src.read_bytes()).decode()
    cmd = [
        "ssh",
        MDX_USER_HOST,
        f"pwsh -Command \"[IO.File]::WriteAllBytes('{dst}', "
        f"[Convert]::FromBase64String([Console]::In.ReadToEnd()))\"",
    ]
    p = subprocess.run(cmd, input=b64, text=True, capture_output=True)
    if p.returncode != 0:
        print(f"  [FAIL rc={p.returncode}] {dst}")
        if p.stderr:
            print(f"    stderr: {p.stderr[:200]}")
        return False
    print(f"  [OK] {dst}  ({src.stat().st_size} bytes)")
    return True


def stop_locks_on_mdx() -> None:
    """Stop processes that hold .pyd handles open on mdx.

    Any python.exe / pythonw.exe / jupyter* process that has imported
    radia will hold _radia_pybind.pyd open and cause WriteAllBytes to
    fail with WinError 32. Same for mcp-server-* console scripts.

    This tool is destructive (kill -9 equivalent on Windows). The mdx
    box is a verification target only, so dropping in-flight Python
    sessions is acceptable.
    """
    cmd = [
        "ssh",
        MDX_USER_HOST,
        "pwsh -Command \""
        "Get-Process -ErrorAction SilentlyContinue | Where-Object { "
        "$_.Name -in 'python','pythonw','jupyter','jupyter-lab','ipykernel' -or "
        "$_.Name -like 'mcp-server*' "
        "} | ForEach-Object { Stop-Process -Id $_.Id -Force -ErrorAction SilentlyContinue }; "
        "Start-Sleep -Seconds 3\"",
    ]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        print(f"  [WARN] stop-locks rc={p.returncode}: {p.stderr[:200]}")


def main() -> int:
    if not LAB_REPO.exists():
        print(f"ERROR: LAB repo not found at {LAB_REPO}")
        return 2
    print(f"Pushing C++ binaries from {LAB_REPO} to {MDX_USER_HOST}:")
    print(f"  (stopping mdx mcp-server* + python.exe processes first)")
    stop_locks_on_mdx()
    n_total = 0
    n_ok = 0
    for src, dsts in PYDS:
        for dst in dsts:
            n_total += 1
            if push_one(src, dst):
                n_ok += 1
    print(f"\nDone: {n_ok} / {n_total} succeeded.")
    return 0 if n_ok == n_total else 1


if __name__ == "__main__":
    sys.exit(main())
