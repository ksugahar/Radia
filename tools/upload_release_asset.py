"""Upload one or more files to a GitHub Release -- gh-CLI-free.

Mirror of tools/download_release_asset.py for the WRITE side: replaces
`gh release upload <tag> <file>... --clobber`.

Requires a token with `contents: write` scope. It uses the same local
discovery order as tools/gh_api.py: `GH_TOKEN`, `GITHUB_TOKEN`,
`~/.radia/gh_token`, then `gh auth token` (GitHub CLI browser login).
The workflow uses ${{ secrets.GITHUB_TOKEN }} as usual.

If no token is set, this exits 0 with a warning (so it never blocks
push) -- same graceful-fail semantics as the original bash hook.

Usage:
    python tools/upload_release_asset.py \
        --repo ksugahar/Radia \
        --tag  binaries \
        src/radia/foo.pyd src/radia/bar.pyd

Exit codes:
    0  upload succeeded, OR no token found locally (graceful skip)
    1  asset push failed after retries (rare; release missing etc.)
"""
from __future__ import annotations
import argparse
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

from gh_api import token as github_token


API_BASE = "https://api.github.com"
UPLOAD_BASE = "https://uploads.github.com"


def _headers() -> dict[str, str]:
    token = github_token()
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


def create_release(repo: str, tag: str, title: str = "Binary Files",
                   notes: str = "Pre-built modules. Updated by CI.") -> dict:
    """POST /repos/{repo}/releases -- gh-free `gh release create`."""
    url = f"{API_BASE}/repos/{repo}/releases"
    body = json.dumps({"tag_name": tag, "name": title, "body": notes,
                       "make_latest": "false"}).encode()
    headers = _headers()
    headers["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=body, method="POST", headers=headers)
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.load(r)


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

    # GitHub upload_url has a template suffix like "{?name,label}" -- strip it.
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
    p.add_argument("--create-if-missing", action="store_true",
                    help="create the release (POST /releases) if the tag has "
                         "no release yet -- gh-free `gh release create`")
    p.add_argument("files", nargs="+", help="one or more files to upload")
    args = p.parse_args()

    token = github_token()
    if not token:
        print("Warning: no GitHub token found "
                "(GH_TOKEN/GITHUB_TOKEN, ~/.radia/gh_token, or gh auth). "
                "Skipping binary upload (push will continue).")
        return 0

    release = fetch_release(args.repo, args.tag)
    if release is None:
        if args.create_if_missing:
            print(f"Release '{args.tag}' not found -- creating it")
            try:
                release = create_release(args.repo, args.tag)
            except urllib.error.HTTPError as e:
                print(f"::error::failed to create release '{args.tag}': "
                        f"{e.code} {e.reason}")
                return 1
        else:
            print(f"::warning::Release '{args.tag}' on {args.repo} not found. "
                    f"Re-run with --create-if-missing to create it.")
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
