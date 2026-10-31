"""Read-only publication proof for one tagged source and immutable CI wheel.

No build, install, native import, tag mutation or publication happens here.
"""
from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import tomllib
from urllib import error, parse, request
import zipfile

PACKAGE = Path(__file__).resolve().parents[1]
MAX_ARCHIVE = 256 * 1024 * 1024
MAX_JSON = 8 * 1024 * 1024
REQUIRED_JOBS = ("build-test", "radia-optuna installed-wheel MATLAB/Simulink E2E")
IDENTITY = ("schema", "repository", "ref", "tag", "tag_object_sha", "source_sha",
            "ci_run_id", "run_attempt", "job_ids", "artifact_id", "wheel",
            "wheel_sha256", "version", "archive_sha256")


class GateError(ValueError):
    pass


def require(value, message):
    if not value:
        raise GateError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def identifier(value):
    require(type(value) is int and value > 0, "Expected a positive service ID")
    return value


def checked_url(url, *, signed=False, public=False):
    u = parse.urlsplit(url)
    require(u.scheme == "https" and u.username is None and u.password is None
            and u.port in (None, 443) and not u.fragment, "Unsafe service URL")
    if signed:
        allowed = (u.hostname in {"release-assets.githubusercontent.com", "objects.githubusercontent.com"}
                   or bool(u.hostname and u.hostname.endswith(".blob.core.windows.net")))
    else:
        allowed = u.hostname == ("pypi.org" if public else "api.github.com")
    require(allowed, "Unapproved API/download host")
    return url


class NoRedirect(request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


class Remote:
    def __init__(self, token, opener=None):
        self.token = token
        self.opener = opener or request.build_opener(NoRedirect())

    def read(self, url, *, binary=False, public=False, absent=False):
        checked_url(url, public=public)
        limit = MAX_ARCHIVE if binary else MAX_JSON
        authorization = not public
        for hop in range(4):
            headers = {"Accept": "application/octet-stream" if binary else
                       ("application/json" if public else "application/vnd.github+json")}
            if authorization:
                headers["Authorization"] = "Bearer " + self.token
                headers["X-GitHub-Api-Version"] = "2022-11-28"
            try:
                with self.opener.open(request.Request(url, headers=headers), timeout=30) as response:
                    length = response.headers.get("Content-Length")
                    require(length is None or 0 <= int(length) <= limit, "Response too large")
                    body = response.read(limit + 1)
                    require(len(body) <= limit, "Response too large")
                    return body
            except error.HTTPError as failure:
                if failure.code == 404 and absent and hop == 0:
                    return None
                if binary and failure.code in (301, 302, 303, 307, 308):
                    url = checked_url(failure.headers.get("Location", ""), signed=True)
                    authorization = False  # Never forward the GitHub bearer token to signed hosts.
                    continue
                raise GateError("Service request failed; absence is not established") from None
            except (error.URLError, TimeoutError, OSError):
                raise GateError("Service unavailable; absence is not established") from None
        raise GateError("Too many download redirects")

    def json(self, path, *, absent=False, public=False):
        url = path if public else "https://api.github.com" + path
        body = self.read(url, absent=absent, public=public)
        if body is None:
            return None
        try:
            value = json.loads(body)
        except (ValueError, UnicodeError):
            raise GateError("Invalid service JSON") from None
        require(isinstance(value, dict), "Expected a service object")
        return value

    def pages(self, path, key):
        result = []
        for page in range(1, 101):
            separator = "&" if "?" in path else "?"
            value = self.json(f"{path}{separator}per_page=100&page={page}").get(key)
            require(isinstance(value, list), "Invalid paginated response")
            result.extend(value)
            if len(value) < 100:
                return result
        raise GateError("Pagination limit exceeded")


def ref_identity(context, version):
    repo = context["repository"]
    require(re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", repo)
            and not any(x in {".", ".."} for x in repo.split("/")), "Invalid repository")
    require(re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+", version), "Invalid release version")
    tag = "radia-optuna-v" + version
    require(context["ref_type"] == "tag" and context["ref"] == "refs/tags/" + tag,
            "Dispatch must use the exact package-version tag ref")
    require(re.fullmatch(r"[0-9a-f]{40}", context["sha"]), "Invalid triggering SHA")
    return repo, tag


def peel(reference, fetch):
    node = reference["object"]
    seen = set()
    for _ in range(8):
        sha = node.get("sha", "")
        require(re.fullmatch(r"[0-9a-f]{40}", sha) and sha not in seen, "Invalid/cyclic tag object")
        seen.add(sha)
        if node.get("type") == "commit":
            return sha
        require(node.get("type") == "tag", "Tag does not peel to a commit")
        value = fetch(sha)
        require(value.get("sha") == sha, "Annotated tag object changed")
        node = value["object"]
    raise GateError("Tag nesting limit exceeded")


def run_identity(run, jobs, repository, source_sha, run_id):
    for key, expected in {"id": run_id, "status": "completed", "conclusion": "success",
                          "head_sha": source_sha, "head_branch": "main",
                          "path": ".github/workflows/radia-optuna.yml"}.items():
        require(run.get(key) == expected, f"CI {key}: got {run.get(key)!r}, expected {expected!r}")
    require(run.get("event") in {"push", "workflow_dispatch"}, "Unapproved CI event")
    require(run.get("head_repository", {}).get("full_name") == repository, "Foreign CI repository")
    selected = {}
    for name in REQUIRED_JOBS:
        matches = [job for job in jobs if job.get("name") == name]
        require(len(matches) == 1, f"Required job count for {name!r}: {len(matches)}, expected 1")
        job = matches[0]
        require(job.get("status") == "completed" and job.get("conclusion") == "success",
                "Required CI job did not pass: " + name)
        selected[name] = identifier(job.get("id"))
    return {"run_attempt": identifier(run.get("run_attempt")), "job_ids": selected}


def choose_artifact(artifacts):
    matches = [item for item in artifacts if item.get("name") == "radia-optuna-wheel"]
    require(len(matches) == 1, "Expected exactly one CI wheel artifact")
    item = matches[0]
    require(item.get("expired") is False and 0 < item.get("size_in_bytes", 0) <= MAX_ARCHIVE,
            "Expired or oversized CI artifact")
    identifier(item.get("id"))
    return item


def unpack_wheel(data, version):
    require(len(data) <= MAX_ARCHIVE, "Oversized artifact archive")
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        items = archive.infolist()
        require(len(items) == 1, "Ambiguous/duplicate artifact ZIP members")
        item = items[0]
        require(not item.is_dir() and PurePosixPath(item.filename).name == item.filename
                and "\\" not in item.filename and ":" not in item.filename
                and stat.S_IFMT(item.external_attr >> 16) != stat.S_IFLNK,
                "Unsafe artifact ZIP path")
        require(re.fullmatch(r"radia_optuna-" + re.escape(version) + r"-(cp\d+)-\1-win_amd64\.whl",
                             item.filename), "Unexpected wheel name/ABI")
        require(not item.flag_bits & 1 and item.file_size <= MAX_ARCHIVE
                and item.file_size <= max(1, item.compress_size) * 100, "Unsafe ZIP size/encryption")
        return item.filename, archive.read(item)


def wheel_source(data, source_sha):
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        items = archive.infolist()
        names = [item.filename for item in items]
        require(len(names) <= 10000 and len(names) == len(set(names)), "Duplicate/oversized wheel member list")
        require(sum(item.file_size for item in items) <= MAX_ARCHIVE, "Oversized wheel expansion")
        for item in items:
            path = PurePosixPath(item.filename)
            require(not path.is_absolute() and ".." not in path.parts
                    and "\\" not in item.filename and ":" not in item.filename
                    and stat.S_IFMT(item.external_attr >> 16) != stat.S_IFLNK,
                    "Unsafe wheel member path")
        build = json.loads(archive.read("radia_optuna/native_build.json"))
        require(build.get("source_sha") == source_sha, "Native receipt/source SHA mismatch")
        require(build.get("build_route") == "nvcc-shared-library" and build.get("wrapper") == "plain-msvc-mex",
                "Unapproved native build route")
        records = [name for name in names if name.endswith(".dist-info/RECORD")]
        require(len(records) == 1, "Missing/ambiguous wheel RECORD")
        seen = set()
        for name, value, size in csv.reader(io.StringIO(archive.read(records[0]).decode("utf8"))):
            require(name not in seen and name in names, "Invalid wheel RECORD member")
            seen.add(name)
            if name == records[0]:
                require(value == size == "", "Invalid RECORD self entry")
            else:
                body = archive.read(name)
                expected = "sha256=" + base64.urlsafe_b64encode(hashlib.sha256(body).digest()).decode().rstrip("=")
                require(value == expected and size == str(len(body)), "Wheel RECORD hash mismatch")
        require(seen == set(names), "Unrecorded wheel members")


def pypi_resume(value, filename, sha256, version):
    if value is None:
        return False
    require(value.get("info", {}).get("version") == version, "PyPI version identity mismatch")
    matches = [item for item in value.get("urls", []) if item.get("filename") == filename]
    require(len(matches) == 1, "Existing PyPI version does not have the exact candidate filename")
    item = matches[0]
    require(item.get("digests", {}).get("sha256") == sha256 and item.get("yanked", False) is False,
            "Existing PyPI bytes differ or are yanked; overwrite prohibited")
    return True


def compare_proofs(actual, frozen):
    for key in IDENTITY:
        require(key in frozen and actual.get(key) == frozen[key], "Approved proof changed: " + key)


def verify(remote, context, run_id, expected_hash, destination, frozen=None):
    version = tomllib.loads((PACKAGE / "pyproject.toml").read_text(encoding="utf8"))["project"]["version"]
    repo, tag = ref_identity(context, version)
    require(re.fullmatch(r"[0-9a-fA-F]{64}", expected_hash), "Invalid candidate SHA256")
    expected_hash = expected_hash.lower()
    prefix = "/repos/" + repo
    reference = remote.json(prefix + "/git/ref/tags/" + tag)
    require(reference.get("ref") == context["ref"], "API tag ref differs from dispatch")
    source_sha = peel(reference, lambda sha: remote.json(prefix + "/git/tags/" + sha))
    require(source_sha == context["sha"], "Triggering SHA differs from peeled tag commit")
    run = remote.json(prefix + f"/actions/runs/{identifier(run_id)}")
    jobs = remote.pages(prefix + f"/actions/runs/{run_id}/jobs?filter=latest", "jobs")
    identity = run_identity(run, jobs, repo, source_sha, run_id)
    artifact = choose_artifact(remote.pages(prefix + f"/actions/runs/{run_id}/artifacts", "artifacts"))
    raw = remote.read("https://api.github.com" + prefix + f"/actions/artifacts/{artifact['id']}/zip", binary=True)
    archive_digest = artifact.get("digest")
    require(isinstance(archive_digest, str) and re.fullmatch(r"sha256:[0-9a-f]{64}", archive_digest),
            "Missing or invalid artifact archive digest")
    require(archive_digest == "sha256:" + digest(raw), "Artifact archive digest mismatch")
    filename, data = unpack_wheel(raw, version)
    require(digest(data) == expected_hash, "Wheel differs from accepted candidate SHA256")
    wheel_source(data, source_sha)
    proof = {"schema": "radia-optuna.publication-proof.v1", "repository": repo,
             "ref": context["ref"], "tag": tag, "tag_object_sha": reference["object"]["sha"],
             "source_sha": source_sha, "ci_run_id": run_id, **identity,
             "artifact_id": artifact["id"], "wheel": filename,
             "wheel_sha256": expected_hash, "archive_sha256": digest(raw), "version": version}
    if frozen is not None:
        compare_proofs(proof, frozen)
    require(not destination.is_symlink() and (not destination.exists() or not any(destination.iterdir())),
            "Artifact output must be a fresh empty directory")
    destination.mkdir(parents=True, exist_ok=True)
    wheel = destination / filename
    with wheel.open("xb") as stream:
        stream.write(data)
    # Reuse the package's existing cross-checkout verifier; no extraction/build/native import.
    spec = importlib.util.spec_from_file_location("optuna_publication_wheel_verifier", PACKAGE / "verify_wheel.py")
    require(spec is not None and spec.loader is not None, "Missing package verifier")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    verified = module.verify(wheel, release_candidate=True)
    require(verified.get("ok") is True, "Package wheel verification did not pass")
    proof["pypi_exists"] = pypi_resume(remote.json(
        f"https://pypi.org/pypi/radia-optuna/{version}/json", public=True, absent=True), filename, expected_hash, version)
    release = remote.json(prefix + "/releases/tags/" + tag, absent=True)
    proof["release_exists"] = release is not None
    proof["github_asset_exists"] = False
    if release is not None:
        require(release.get("tag_name") == tag and release.get("draft") is False
                and release.get("prerelease") is False, "Unexpected GitHub release state")
        assets = [item for item in release.get("assets", []) if item.get("name") == filename]
        require(len(assets) <= 1, "Ambiguous existing release assets")
        if assets:
            asset = assets[0]
            require(asset.get("state") == "uploaded" and asset.get("size") == len(data), "Existing asset identity mismatch")
            actual = remote.read("https://api.github.com" + prefix + f"/releases/assets/{identifier(asset.get('id'))}", binary=True)
            require(digest(actual) == expected_hash, "Existing release asset bytes differ; overwrite prohibited")
            proof["github_asset_exists"] = True
    return proof


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-id", type=int, required=True)
    p.add_argument("--candidate-sha256", required=True)
    p.add_argument("--destination", type=Path, required=True)
    p.add_argument("--receipt", type=Path, required=True)
    p.add_argument("--approved-receipt", type=Path)
    args = p.parse_args()
    context = {"repository": os.environ["GITHUB_REPOSITORY"], "ref": os.environ["GITHUB_REF"],
               "ref_type": os.environ["GITHUB_REF_TYPE"], "sha": os.environ["GITHUB_SHA"]}
    frozen = json.loads(args.approved_receipt.read_text(encoding="utf8")) if args.approved_receipt else None
    require(not args.receipt.exists(), "Proof output collision")
    proof = verify(Remote(os.environ["GH_TOKEN"]), context, args.run_id,
                   args.candidate_sha256, args.destination, frozen)
    args.receipt.write_text(json.dumps(proof, indent=2)+"\n", encoding="utf8", newline="\n")
    if os.environ.get("GITHUB_OUTPUT"):
        with open(os.environ["GITHUB_OUTPUT"], "a", encoding="utf8", newline="\n") as output:
            for key in ("artifact_id", "source_sha", "tag", "version", "wheel", "wheel_sha256",
                        "pypi_exists", "release_exists", "github_asset_exists"):
                value = str(proof[key]).lower() if isinstance(proof[key], bool) else str(proof[key])
                output.write(f"{key}={value}\n")
    print(json.dumps({key: proof[key] for key in ("source_sha", "artifact_id", "wheel_sha256", "version")}))


if __name__ == "__main__":
    main()
