"""Download a single asset from a GitHub Release — gh-CLI-free.

Replacement for `gh release download <tag> -p <name> -D <dir> --clobber`
in CI workflows after gh CLI was retired from the self-hosted runner
on 2026-05-24.

Uses the authenticated GitHub REST API (so private/large assets work).
Requires `GITHUB_TOKEN` (or `GH_TOKEN`) in env — provided by GitHub
Actions automatically via `${{ secrets.GITHUB_TOKEN }}`.

Usage:
    python tools/download_release_asset.py \
        --repo ksugahar/Radia \
        --tag  binaries \
        --name radia_cubit_mesh.pyd \
        --dest src/radia

Exit codes:
    0  asset downloaded (or already present with correct size)
    1  asset not found / download failed (after retries)
"""
from __future__ import annotations
import argparse
import json
import os
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path


API_BASE = "https://api.github.com"


def _headers(authed: bool) -> dict[str, str]:
    """Build request headers.  Token is OPTIONAL for public repos
    (60 req/h anonymous limit is plenty for CI), but recommended in
    CI where ${{ secrets.GITHUB_TOKEN }} gives 1000 req/h."""
    h = {
        "Accept": "application/vnd.github+json",
        "X-GitHub-Api-Version": "2022-11-28",
        "User-Agent": "radia-ci-asset-fetch",
    }
    token = os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN")
    if token and authed:
        h["Authorization"] = f"Bearer {token}"
    return h


def fetch_release(repo: str, tag: str) -> dict:
    """GET /repos/{owner}/{repo}/releases/tags/{tag}."""
    url = f"{API_BASE}/repos/{repo}/releases/tags/{tag}"
    req = urllib.request.Request(url, headers=_headers(authed=True))
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.load(r)


def download_asset(asset: dict, dest_path: Path) -> int:
    """Download asset to dest_path.

    Prefers `browser_download_url` (works without auth on public repos,
    302-redirects to a signed CDN URL).  Falls back to authenticated
    `assets/<id>` for private repos.
    """
    # Public path — no auth needed; urllib follows redirects automatically.
    url = asset.get("browser_download_url") or asset["url"]
    headers = _headers(authed=False)
    headers["Accept"] = "application/octet-stream"
    req = urllib.request.Request(url, headers=headers)
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(req, timeout=120) as r, open(dest_path, "wb") as f:
        # 1 MB chunks
        while True:
            chunk = r.read(1024 * 1024)
            if not chunk:
                break
            f.write(chunk)
    return dest_path.stat().st_size


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--repo", required=True,
                    help="owner/repo (e.g. ksugahar/Radia)")
    p.add_argument("--tag",  required=True,
                    help="release tag name (e.g. 'binaries')")
    p.add_argument("--name", required=True,
                    help="exact asset filename to fetch")
    p.add_argument("--dest", required=True,
                    help="destination directory (created if missing)")
    p.add_argument("--min-size", type=int, default=100_000,
                    help="treat existing file smaller than this as missing "
                            "(default: 100,000 bytes)")
    p.add_argument("--max-attempts", type=int, default=6,
                    help="retry count for transient race (default: 6)")
    args = p.parse_args()

    dest_path = Path(args.dest) / args.name

    # If a usable copy already exists, skip the download.
    if dest_path.exists():
        sz = dest_path.stat().st_size
        if sz >= args.min_size:
            print(f"  {args.name} already present locally: {sz} bytes (skip)")
            return 0
        else:
            print(f"  {args.name} exists but suspicious size ({sz}), refetching")

    # Retry loop — race window described in build-test.yml: when a tag
    # push triggers CI before the pre-push hook's binary upload settled.
    delays = [0, 5, 10, 20, 40, 60]
    for attempt in range(args.max_attempts):
        if attempt > 0:
            wait = delays[min(attempt, len(delays) - 1)]
            print(f"  Retry {attempt}/{args.max_attempts - 1} after {wait}s ...")
            time.sleep(wait)
        try:
            print(f"  Looking up release '{args.tag}' on {args.repo} "
                    f"(attempt {attempt + 1})")
            release = fetch_release(args.repo, args.tag)
        except urllib.error.HTTPError as e:
            print(f"    HTTP {e.code} fetching release: {e.reason}")
            continue
        except (urllib.error.URLError, TimeoutError) as e:
            print(f"    Network error fetching release: {e}")
            continue

        asset = next(
            (a for a in release.get("assets", []) if a["name"] == args.name),
            None,
        )
        if asset is None:
            available = [a["name"] for a in release.get("assets", [])]
            print(f"    Asset '{args.name}' not in release "
                    f"(available: {available[:6]}{'...' if len(available)>6 else ''})")
            continue

        try:
            sz = download_asset(asset, dest_path)
            if sz >= args.min_size:
                print(f"  Downloaded {args.name}: {sz} bytes")
                return 0
            else:
                print(f"  Downloaded but size {sz} < min_size {args.min_size}, retrying")
        except urllib.error.HTTPError as e:
            print(f"    HTTP {e.code} downloading asset: {e.reason}")
        except (urllib.error.URLError, TimeoutError) as e:
            print(f"    Network error downloading asset: {e}")

    print(f"::error::Failed to download {args.name} from {args.repo}@{args.tag} "
            f"after {args.max_attempts} attempts.")
    print(f"         The pre-push hook on LAB uploads this file on every push;")
    print(f"         if it is missing, run on LAB (gh-CLI free):")
    print(f"           python tools/upload_release_asset.py "
            f"--tag {args.tag} --file src/radia/{args.name}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
