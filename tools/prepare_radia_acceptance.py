"""Export public numeric acceptance after verifying private current-machine records.

Keep the raw input directories and one shared random batch key in private
release storage. Public exports contain opaque distinct-machine tokens and
hashes of the verified originals; they contain no hostnames or runtime paths.
"""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile

from release_acceptance import RELEASE_ACCEPTANCE_HOSTS
from verify_radia_promotion import prepare_public_host, strict_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--target", choices=RELEASE_ACCEPTANCE_HOSTS, required=True)
    parser.add_argument("--private-batch-key", type=Path, required=True)
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--ci-run", type=int, required=True)
    args = parser.parse_args()
    with zipfile.ZipFile(args.wheel) as archive:
        identity = dict(
            wheel_sha256=hashlib.sha256(args.wheel.read_bytes()).hexdigest(),
            native_sha256=hashlib.sha256(archive.read("radia/_radia_pybind.pyd")).hexdigest(),
            source_commit=args.source_commit, ci_run=args.ci_run,
            sources_verified=sum(n.startswith("radia/") and n.endswith(".py") for n in archive.namelist()),
            version=args.wheel.name.split("-")[1])
    a, f, xml = prepare_public_host(
        strict_json((args.input_dir / "acceptance.json").read_bytes()),
        strict_json((args.input_dir / "full6.json").read_bytes()),
        (args.input_dir / "focused.xml").read_bytes(), identity, args.target,
        args.private_batch_key.read_bytes())
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, record in (("acceptance.json", a), ("full6.json", f)):
        (args.output_dir / name).write_text(json.dumps(record, indent=2, allow_nan=False) + "\n",
                                         encoding="utf-8")
    (args.output_dir / "focused.xml").write_bytes(xml)
    print("Verified public acceptance export complete")


if __name__ == "__main__":
    main()
