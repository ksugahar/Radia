"""Three-way update of the vendored compact Netgen snapshot.

    python update_snapshot.py --netgen-repo <netgen clone> --old <commit> --new <tag>

For every file under netgen_src/ and every *_patched.cpp (a trimmed copy of
an upstream file), merge
    ours   = the vendored file with its local patches,
    base   = the upstream file at --old (the current snapshot),
    theirs = the upstream file at --new,
with `git merge-file`, writing LF line endings.  Conflicts are left in the
files with standard markers and listed at the end; resolve them, then build.
A trimmed *_patched.cpp does not receive upstream functions added in regions
it omits: the link step reports them, and they must be copied by hand.
Update netgen_version.hpp to the new tag afterwards.
"""
import argparse
import subprocess
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
SNAP = HERE / "netgen_src"
PATCHED = {
    "meshclass_patched.cpp": "libsrc/meshing/meshclass.cpp",
    "meshtype_patched.cpp": "libsrc/meshing/meshtype.cpp",
    "basegeom_patched.cpp": "libsrc/meshing/basegeom.cpp",
    "flags_patched.cpp": "libsrc/core/flags.cpp",
    "utils_patched.cpp": "libsrc/core/utils.cpp",
    "version_patched.cpp": "libsrc/core/version.cpp",
    "statushandler_patched.cpp": "libsrc/core/statushandler.cpp",
}


def show(repo, rev, path):
    r = subprocess.run(["git", "-C", str(repo), "show", f"{rev}:{path}"], capture_output=True)
    return None if r.returncode else r.stdout.replace(b"\r\n", b"\n")


def merge(repo, old, new, target, upath):
    base, theirs = show(repo, old, upath), show(repo, new, upath)
    ours = target.read_bytes().replace(b"\r\n", b"\n")
    if theirs is None:
        return "removed-upstream"
    if base == theirs:
        return "unchanged"
    if ours == base:
        target.write_bytes(theirs)
        return "fast-forward"
    with tempfile.TemporaryDirectory() as td:
        o, b, t = (Path(td) / n for n in ("ours", "base", "theirs"))
        o.write_bytes(ours); b.write_bytes(base); t.write_bytes(theirs)
        r = subprocess.run(["git", "merge-file", "-L", "ours", "-L", "base", "-L", "theirs",
                            str(o), str(b), str(t)], capture_output=True)
        target.write_bytes(o.read_bytes())
        return "merged" if r.returncode == 0 else f"CONFLICT({r.returncode})"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--netgen-repo", type=Path, required=True)
    ap.add_argument("--old", required=True, help="upstream commit of the current snapshot")
    ap.add_argument("--new", required=True, help="upstream tag or commit to update to")
    a = ap.parse_args()
    results = {}
    for f in sorted(SNAP.rglob("*")):
        if f.is_file():
            rel = f.relative_to(SNAP).as_posix()
            results[rel] = merge(a.netgen_repo, a.old, a.new, f, "libsrc/" + rel)
    for name, upath in PATCHED.items():
        results[name] = merge(a.netgen_repo, a.old, a.new, HERE / name, upath)
    for k, v in results.items():
        if v != "unchanged":
            print(f"{v:18s} {k}")
    conflicts = [k for k, v in results.items() if v.startswith("CONFLICT")]
    print(f"unchanged: {sum(v == 'unchanged' for v in results.values())}, conflicts: {len(conflicts)}")
    return 1 if conflicts else 0


if __name__ == "__main__":
    raise SystemExit(main())
