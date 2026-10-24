"""Scan public source and explicit artifacts with approved external criteria."""
import argparse
import hashlib
import zipfile
import json
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
APPROVED_RULES_SHA256 = "08019894d2a840682ff78b6a7015a29c14a82549b3337a638ed0f131b31fbbdc"
FORBIDDEN = None


class PolicyInputError(ValueError):
    pass


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise PolicyInputError("Duplicate private-policy field")
        result[key] = value
    return result


def load_rules(explicit_file=None):
    """Require exact reviewed bytes; no embedded or missing-input fallback."""
    filename = explicit_file or os.environ.get("RADIA_PUBLIC_NAME_RULES_FILE")
    payload = os.environ.get("RADIA_PUBLIC_NAME_RULES_JSON")
    if explicit_file:
        payload = None
    elif bool(filename) == bool(payload):
        raise PolicyInputError("Exactly one private-policy source is required")
    try:
        raw = Path(filename).read_bytes() if filename else payload.encode("utf-8")
    except (OSError, AttributeError):
        raise PolicyInputError("Private-policy source is unavailable") from None
    if hashlib.sha256(raw).hexdigest() != APPROVED_RULES_SHA256:
        raise PolicyInputError("Private-policy fingerprint does not match the reviewed policy")
    try:
        text = raw.decode("utf-8")
        if os.environ.get("GITHUB_ACTIONS") == "true":
            # Workflow command escaping prevents multiline command injection.
            masked = text.replace("%", "%25").replace("\r", "%0D").replace("\n", "%0A")
            print("::add-mask::" + masked, flush=True)
        data = json.loads(text, object_pairs_hook=_unique_object)
        if not isinstance(data, dict) or set(data) != {"schema", "flags", "patterns"}:
            raise PolicyInputError("Invalid private-policy fields")
        if data["schema"] != "radia.private-public-name-rules.v1" or data["flags"] != ["IGNORECASE"]:
            raise PolicyInputError("Unsupported private-policy schema or flags")
        if not isinstance(data["patterns"], list) or len(data["patterns"]) != 1 or not isinstance(data["patterns"][0], str) or not data["patterns"][0]:
            raise PolicyInputError("Invalid private-policy pattern set")
        return re.compile(data["patterns"][0], re.IGNORECASE)
    except (UnicodeError, json.JSONDecodeError, re.error):
        raise PolicyInputError("Private-policy data cannot be decoded or compiled") from None


def notebook_text(document):
    for cell in document.get("cells", []):
        yield "".join(cell.get("source", []))
        for output in cell.get("outputs", []):
            yield "".join(output.get("text", []))
            for mime, value in output.get("data", {}).items():
                if mime.startswith("text/"):
                    yield "".join(value) if isinstance(value, list) else str(value)


def binary_hits(data):
    # Scan ASCII/UTF-8 plus both UTF-16 byte alignments and endian orders.
    decoded = [data.decode("latin1"), data.decode("utf-8", errors="ignore")]
    for encoding in ("utf-16-le", "utf-16-be"):
        decoded.extend(data[offset:].decode(encoding, errors="ignore") for offset in (0, 1))
    return any(FORBIDDEN.search(text) for text in decoded)


def artifact_hits(path):
    bad = []
    if FORBIDDEN.search(path.name):
        bad.append(str(path) + ": name")
    if path.suffix == ".whl":
        with zipfile.ZipFile(path) as wheel:
            for member in wheel.infolist():
                if FORBIDDEN.search(member.filename) or binary_hits(wheel.read(member)):
                    bad.append(str(path) + ": member " + member.filename)
    elif binary_hits(path.read_bytes()):
        bad.append(str(path) + ": bytes")
    return bad


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, action="append", default=[])
    parser.add_argument("--report", type=Path)
    parser.add_argument("--rules-file", type=Path)
    args = parser.parse_args()
    global FORBIDDEN
    try:
        FORBIDDEN = load_rules(args.rules_file)
    except PolicyInputError as error:
        raise SystemExit("FAIL: " + str(error)) from None
    paths = subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT).decode().split("\0")
    bad = []
    semantic_count = 0
    semantic_hash = hashlib.sha256()
    for name in filter(None, paths):
        path = ROOT / name
        if FORBIDDEN.search(name):
            bad.append(name + ": path")
        if not path.is_file():
            continue
        data = path.read_bytes()
        if b"\0" in data:
            continue
        try:
            text = data.decode("utf-8")
        except UnicodeDecodeError:
            continue
        texts = notebook_text(json.loads(text)) if path.suffix == ".ipynb" else [text]
        semantic_count += 1
        semantic_hash.update(name.encode("utf-8") + b"\0" + data)
        if any(FORBIDDEN.search(part) for part in texts):
            bad.append(name + ": text")
    artifacts = []
    for path in args.artifact:
        data = path.read_bytes()
        artifacts.append({"path": str(path), "sha256": hashlib.sha256(data).hexdigest()})
        bad.extend(artifact_hits(path))
    report = {"criteria_sha256": APPROVED_RULES_SHA256, "status": "FAIL" if bad else "PASS", "source_semantic_files": semantic_count,
              "source_semantic_sha256": semantic_hash.hexdigest(),
              "source_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
              "source_scope": "tracked paths and UTF-8 semantic text; notebook text excludes binary image payloads",
              "artifact_scope": "only explicitly supplied artifacts; no release-artifact result inferred",
              "artifacts": artifacts, "findings": bad}
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if bad:
        raise SystemExit("\n".join(bad))
    print("PASS: public product names; explicit artifacts=" + str(len(artifacts)))


if __name__ == "__main__":
    main()
