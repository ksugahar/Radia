"""Fail closed when a radia-optuna wheel is partial or crosses its boundary."""

from __future__ import annotations

import argparse
import json
import re
import sys
import tomllib
import zipfile
from email.parser import BytesParser
from pathlib import Path, PurePosixPath

PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parents[1]
PACKAGE_PREFIX = PurePosixPath("radia_optuna")
MATLAB_PREFIX = PACKAGE_PREFIX / "matlab"
OPTUNA_PREFIX = MATLAB_PREFIX / "+radia" / "+optuna"
FORBIDDEN_NAME_PARTS = (
    "radia_mex",
    "ngsolve",
    "netgen",
    "mkl_",
    "mkl.",
    "libmkl",
    "radia_pybind",
)
TEXT_PAYLOAD_SUFFIXES = {".json", ".m", ".md", ".py"}


def _fail(messages: list[str]) -> None:
    for message in messages:
        print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def _metadata(archive: zipfile.ZipFile, names: set[str]):
    candidates = sorted(name for name in names if name.endswith(".dist-info/METADATA"))
    if len(candidates) != 1:
        _fail([f"expected one METADATA entry, found {len(candidates)}"])
    return BytesParser().parsebytes(archive.read(candidates[0]))


def _payload_matches(
    member: str,
    wheel_payload: bytes,
    source_payload: bytes,
    *,
    release_candidate: bool,
) -> bool:
    """Compare canonical sources without confusing checkout/build variance."""
    if not release_candidate:
        return wheel_payload == source_payload
    if member.endswith(".mexw64"):
        # The exact main-CI wheel has already passed strict same-workspace
        # fidelity. A separately rebuilt MEX is not byte reproducible.
        return True
    if PurePosixPath(member).suffix.lower() in TEXT_PAYLOAD_SUFFIXES:
        return wheel_payload.replace(b"\r\n", b"\n") == source_payload.replace(
            b"\r\n", b"\n"
        )
    return wheel_payload == source_payload


def _is_windows_x64_pe(payload: bytes) -> bool:
    if len(payload) < 1024 or payload[:2] != b"MZ" or len(payload) < 64:
        return False
    pe_offset = int.from_bytes(payload[60:64], "little")
    if pe_offset < 64 or pe_offset + 6 > len(payload):
        return False
    return (
        payload[pe_offset : pe_offset + 4] == b"PE\0\0"
        and int.from_bytes(payload[pe_offset + 4 : pe_offset + 6], "little")
        == 0x8664
    )


def _source_payloads(source_manifest: dict[str, object]) -> dict[str, Path]:
    matlab_root = REPO_ROOT / "matlab"
    payloads = {
        str(PACKAGE_PREFIX / path.name): path
        for path in (PACKAGE_ROOT / "src" / "radia_optuna").iterdir()
        if path.is_file() and path.suffix in {".py", ".json"}
    }
    for source in (matlab_root / "+radia" / "+optuna").rglob("*"):
        if not source.is_file() or source.suffix not in {".m", ".bin"}:
            continue
        member = MATLAB_PREFIX / PurePosixPath(source.relative_to(matlab_root).as_posix())
        payloads[str(member)] = source
    for name in ("optuna_upstream_compatibility.json", "optuna49_api_coverage.json"):
        payloads[str(MATLAB_PREFIX / name)] = matlab_root / name
    payloads[str(MATLAB_PREFIX / "optuna_mex.mexw64")] = (
        matlab_root / "optuna_mex.mexw64"
    )
    for relative in source_manifest["simulink_entry_points"]:
        member = MATLAB_PREFIX / PurePosixPath(str(relative))
        payloads[str(member)] = matlab_root / str(relative)
    payloads[str(MATLAB_PREFIX / "README.md")] = PACKAGE_ROOT / "MATLAB_README.md"
    payloads[str(MATLAB_PREFIX / "LICENSE")] = REPO_ROOT / "LICENSE"
    payloads[str(MATLAB_PREFIX / "THIRD_PARTY_NOTICES.md")] = (
        PACKAGE_ROOT / "THIRD_PARTY_NOTICES.md"
    )
    return payloads


def verify(wheel: Path, *, release_candidate: bool = False) -> dict[str, object]:
    errors: list[str] = []
    if not re.fullmatch(r"radia_optuna-[^-]+-py3-none-win_amd64\.whl", wheel.name):
        errors.append(f"wheel must be tagged py3-none-win_amd64: {wheel.name}")

    source_project = tomllib.loads(
        (PACKAGE_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    )["project"]
    source_version = source_project["version"]
    root_project = tomllib.loads(
        (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    )["project"]
    source_manifest = json.loads(
        (PACKAGE_ROOT / "src" / "radia_optuna" / "manifest.json").read_text(
            encoding="utf-8"
        )
    )
    if source_manifest["version"] != source_version:
        errors.append(
            "package source versions are not synchronized: "
            f"package={source_version}, manifest={source_manifest['version']}"
        )
    expected_base = [f"radia-optuna=={source_version}"]
    expected_upstream = [f"radia-optuna[upstream]=={source_version}"]
    extras = root_project["optional-dependencies"]
    if extras.get("optuna") != expected_base:
        errors.append(f"radia[optuna] must pin {expected_base[0]}")
    if extras.get("optuna-upstream") != expected_upstream:
        errors.append(
            f"radia[optuna-upstream] must pin {expected_upstream[0]}"
        )

    with zipfile.ZipFile(wheel) as archive:
        names = set(archive.namelist())
        source_payloads = _source_payloads(source_manifest)
        metadata = _metadata(archive, names)
        wheel_version = metadata.get("Version")
        if wheel_version != source_version:
            errors.append(
                f"wheel version {wheel_version!r} does not match {source_version!r}"
            )

        init_source = archive.read(str(PACKAGE_PREFIX / "__init__.py")).decode(
            "utf-8"
        )
        init_match = re.search(
            r'^__version__\s*=\s*["\']([^"\']+)["\']', init_source, re.MULTILINE
        )
        init_version = init_match.group(1) if init_match else None
        if init_version != source_version:
            errors.append(
                f"wheel __version__ {init_version!r} does not match {source_version!r}"
            )

        wheel_metadata_entries = sorted(
            name for name in names if name.endswith(".dist-info/WHEEL")
        )
        if len(wheel_metadata_entries) != 1:
            errors.append(
                "expected one WHEEL metadata entry, found "
                f"{len(wheel_metadata_entries)}"
            )
        else:
            wheel_tags = [
                line.removeprefix("Tag: ")
                for line in archive.read(wheel_metadata_entries[0])
                .decode("utf-8")
                .splitlines()
                if line.startswith("Tag: ")
            ]
            if wheel_tags != ["py3-none-win_amd64"]:
                errors.append(f"unexpected internal wheel tags: {wheel_tags}")

        expected_fixed = {
            str(MATLAB_PREFIX / "optuna_mex.mexw64"),
            str(MATLAB_PREFIX / "optuna_upstream_compatibility.json"),
            str(MATLAB_PREFIX / "optuna49_api_coverage.json"),
            str(MATLAB_PREFIX / "README.md"),
            str(MATLAB_PREFIX / "LICENSE"),
            str(MATLAB_PREFIX / "THIRD_PARTY_NOTICES.md"),
            str(
                OPTUNA_PREFIX
                / "+internal"
                / "sobol_direction_numbers.bin"
            ),
            str(PACKAGE_PREFIX / "manifest.json"),
        }
        missing = sorted(expected_fixed - names)
        if missing:
            errors.append("missing required entries: " + ", ".join(missing))

        stale_payloads = sorted(
            member
            for member, source in source_payloads.items()
            if member in names
            and not _payload_matches(
                member,
                archive.read(member),
                source.read_bytes(),
                release_candidate=release_candidate,
            )
        )
        missing_source_payloads = sorted(set(source_payloads).difference(names))
        if missing_source_payloads:
            errors.append(
                "missing checked source payloads: "
                + ", ".join(missing_source_payloads)
            )
        if stale_payloads:
            errors.append(
                "wheel payload differs from the checked monorepo source: "
                + ", ".join(stale_payloads)
            )

        if str(MATLAB_PREFIX / "THIRD_PARTY_NOTICES.md") in names:
            notices = archive.read(
                str(MATLAB_PREFIX / "THIRD_PARTY_NOTICES.md")
            ).decode("utf-8")
            required_notices = (
                "Copyright (c) 2018 Preferred Networks, Inc.",
                "Copyright (c) 2025 Preferred Networks, Inc.",
                "Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.",
                "Joe--Kuo criterion-6 Sobol direction-number table",
                (
                    "Optuna, the Optuna logo and any related marks are trademarks "
                    "of Preferred Networks, Inc."
                ),
                "independent, unofficial project",
            )
            missing_notices = [
                notice for notice in required_notices if notice not in notices
            ]
            if missing_notices:
                errors.append(
                    "third-party notice is incomplete: "
                    + ", ".join(missing_notices)
                )

        matlab_files = sorted(
            name
            for name in names
            if name.startswith(f"{OPTUNA_PREFIX}/") and name.endswith(".m")
        )
        expected_count = int(source_manifest["matlab_file_count"])
        if len(matlab_files) != expected_count:
            errors.append(
                f"expected {expected_count} MATLAB files, found {len(matlab_files)}"
            )

        mex_files = sorted(name for name in names if name.endswith(".mexw64"))
        if mex_files != [str(MATLAB_PREFIX / "optuna_mex.mexw64")]:
            errors.append(f"unexpected MEX inventory: {mex_files}")
        elif not _is_windows_x64_pe(archive.read(mex_files[0])):
            errors.append("optuna_mex is not a valid Windows x64 PE binary")

        forbidden = sorted(
            name
            for name in names
            if any(part in name.lower() for part in FORBIDDEN_NAME_PARTS)
        )
        if forbidden:
            errors.append("forbidden solver/runtime entries: " + ", ".join(forbidden))

        wheel_manifest = json.loads(
            archive.read(str(PACKAGE_PREFIX / "manifest.json")).decode("utf-8")
        )
        if wheel_manifest != source_manifest:
            errors.append("wheel manifest differs from the checked source manifest")
        if wheel_manifest.get("radia_runtime_required") is not False:
            errors.append("standalone manifest must set radia_runtime_required=false")
        if wheel_manifest.get("simulink_standalone") is not True:
            errors.append("standalone manifest must set simulink_standalone=true")
        if wheel_manifest.get("native_sobol_max_dimension") != 21201:
            errors.append("standalone manifest must declare native Sobol dimension 21201")

        simulink_entries = {
            str(MATLAB_PREFIX / PurePosixPath(relative))
            for relative in source_manifest["simulink_entry_points"]
        }
        missing_simulink = sorted(simulink_entries - names)
        if missing_simulink:
            errors.append(
                "missing standalone Simulink entries: " + ", ".join(missing_simulink)
            )
        packaged_simulink_block_entries = {
            name
            for name in names
            if name.startswith(f"{MATLAB_PREFIX}/+radia/+simulink/")
            or name == str(MATLAB_PREFIX / "radia_optuna_sfun.m")
        }
        declared_simulink_block_entries = {
            name
            for name in simulink_entries
            if "/+radia/+simulink/" in name
            or name == str(MATLAB_PREFIX / "radia_optuna_sfun.m")
        }
        if packaged_simulink_block_entries != declared_simulink_block_entries:
            errors.append(
                "standalone Simulink block inventory differs from the manifest"
            )

        adapters = set(source_manifest["radia_integration_adapters"])
        adapter_entries = {
            str(MATLAB_PREFIX / PurePosixPath(relative)) for relative in adapters
        }
        if not adapter_entries.issubset(names):
            errors.append("one or more declared Radia integration adapters are absent")

        requirements = metadata.get_all("Requires-Dist", [])
        if not any(
            re.search(r"^optuna\s*==\s*4\.9\.0\s*;.*extra\s*==\s*['\"]upstream['\"]", req)
            for req in requirements
        ):
            errors.append("the upstream extra does not pin optuna==4.9.0")

    if errors:
        _fail(errors)

    result = {
        "schema": "radia-optuna.wheel-verification.v1",
        "ok": True,
        "wheel": str(wheel.resolve()),
        "version": source_version,
        "platform_tag": "py3-none-win_amd64",
        "matlab_file_count": len(matlab_files),
        "native_gateway": "optuna_mex",
        "native_command_count": source_manifest["native_command_count"],
        "simulink_standalone": source_manifest["simulink_standalone"],
        "simulink_entry_count": len(simulink_entries),
        "radia_integration_adapter_count": len(adapters),
        "source_fidelity_verified": True,
        "source_fidelity_mode": (
            "release-candidate" if release_candidate else "same-workspace-strict"
        ),
        "native_gateway_verification": (
            "exact-main-ci+pe-x64"
            if release_candidate
            else "same-workspace-byte-exact+pe-x64"
        ),
        "source_fidelity_file_count": len(source_payloads),
    }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wheel", type=Path)
    parser.add_argument("--json", action="store_true", dest="as_json")
    parser.add_argument(
        "--release-candidate",
        action="store_true",
        help=(
            "verify an exact successful main-CI artifact across checkouts: "
            "normalize text line endings and validate the MEX as PE/x64"
        ),
    )
    args = parser.parse_args()
    result = verify(args.wheel, release_candidate=args.release_candidate)
    if args.as_json:
        print(json.dumps(result, indent=2))
    else:
        print(
            "radia-optuna wheel PASS: "
            f"{result['matlab_file_count']} MATLAB files, "
            f"{result['native_command_count']} native commands"
        )


if __name__ == "__main__":
    main()
