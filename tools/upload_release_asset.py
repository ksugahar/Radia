"""Upload one or more files to a GitHub Release — gh-CLI-free.

Mirror of tools/download_release_asset.py for the WRITE side: replaces
`gh release upload <tag> <file>... --clobber`.

Requires `GITHUB_TOKEN` (or `GH_TOKEN`) in env with `contents: write`
scope. For the pre-push hook on LAB this is a personal access token
that the developer sets in their shell rc or .git/safe.directory
config; the workflow uses ${{ secrets.GITHUB_TOKEN }} as usual.

If no token is set, this exits 0 with a warning (so it never blocks
push) — same graceful-fail semantics as the original bash hook.

Usage:
    python tools/upload_release_asset.py \
        --repo ksugahar/Radia \
        --tag  binaries \
        src/radia/foo.pyd src/radia/bar.pyd

Exit codes:
    0  upload succeeded, OR no token in env (graceful skip)
    1  asset push failed after retries (rare; release missing etc.)
"""
from __future__ import annotations
import argparse
import json
import os
import sys
import urllib.error
import urllib.request
from pathlib import Path


API_BASE = "https://api.github.com"
UPLOAD_BASE = "https://uploads.github.com"


def _headers() -> dict[str, str]:
    token = os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN")
    if not token:
        return {}
    return {
        "Authorization": f"Bearer {token}",
        "Accept": "application/vnd.github+json",
        "X-GitHub-Api-Version": "2022-11-28",
        "User-Agent": "radia-asset-upload",
    }


def fetch_release(repo: str, tag: str) -> dict | None:
    """Returns release dict or None on 404 / network error."""
    url = f"{API_BASE}/repos/{repo}/releases/tags/{tag}"
    try:
        req = urllib.request.Request(url, headers=_headers())
        with urllib.request.urlopen(req, timeout=30) as r:
            return json.load(r)
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None
        raise
    except (urllib.error.URLError, TimeoutError):
        return None


def delete_asset(repo: str, asset_id: int) -> None:
    """DELETE /repos/{repo}/releases/assets/{id}."""
    url = f"{API_BASE}/repos/{repo}/releases/assets/{asset_id}"
    req = urllib.request.Request(url, method="DELETE", headers=_headers())
    try:
        with urllib.request.urlopen(req, timeout=30):
            pass
    except urllib.error.HTTPError as e:
        if e.code != 404:
            raise


def upload_asset(release: dict, file_path: Path) -> None:
    """POST file to release.upload_url, --clobber semantics."""
    name = file_path.name
    # First, delete any existing asset with the same name (the --clobber
    # part of the original gh command).
    for existing in release.get("assets", []):
        if existing["name"] == name:
            print(f"  Removing existing {name} (id={existing['id']})")
            delete_asset(release["url"].split("/repos/")[1].split("/")[0] + "/"
                            + release["url"].split("/repos/")[1].split("/")[1],
                            existing["id"])

    # GitHub upload_url has a template suffix like "{?name,label}" — strip it.
    upload_template = release["upload_url"]
    base = upload_template.split("{")[0]
    upload_url = f"{base}?name={urllib_quote(name)}"

    body = file_path.read_bytes()
    headers = _headers()
    headers["Content-Type"] = "application/octet-stream"
    headers["Content-Length"] = str(len(body))
    req = urllib.request.Request(upload_url, data=body, method="POST",
                                    headers=headers)
    with urllib.request.urlopen(req, timeout=300) as r:
        result = json.load(r)
    print(f"  Uploaded {name}: {result.get('size')} bytes "
            f"(state={result.get('state')})")


def urllib_quote(s: str) -> str:
    """Tiny shim to avoid importing urllib.parse just for one quote."""
    import urllib.parse
    return urllib.parse.quote(s, safe="")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--repo", required=True,
                    help="owner/repo (e.g. ksugahar/Radia)")
    p.add_argument("--tag",  required=True,
                    help="release tag name (e.g. 'binaries')")
    p.add_argument("files", nargs="+", help="one or more files to upload")
    args = p.parse_args()

    token = os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN")
    if not token:
        print("Warning: GITHUB_TOKEN/GH_TOKEN not set in env. "
                "Skipping binary upload (push will continue).")
        return 0

    release = fetch_release(args.repo, args.tag)
    if release is None:
        print(f"::warning::Release '{args.tag}' on {args.repo} not found. "
                f"Create it with: python tools/create_release.py "
                f"--tag {args.tag} --title 'Binary Files'  "
                f"(then re-run push)")
        return 0  # graceful, don't block push

    failed = []
    for f in args.files:
        p = Path(f)
        if not p.exists():
            print(f"  Skip missing: {f}")
            continue
        try:
            upload_asset(release, p)
        except urllib.error.HTTPError as e:
            print(f"  HTTP {e.code} uploading {p.name}: {e.reason}")
            failed.append(f)
        except (urllib.error.URLError, TimeoutError) as e:
            print(f"  Network error uploading {p.name}: {e}")
            failed.append(f)

    if failed:
        print(f"::error::Failed to upload {len(failed)} asset(s).")
        return 1
    print("Binary upload complete.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
