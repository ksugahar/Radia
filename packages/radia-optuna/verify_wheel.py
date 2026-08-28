"""Fail closed when a radia-optuna wheel is partial or crosses its boundary."""

from __future__ import annotations

import argparse
import json
import re
import struct
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
TEXT_PAYLOAD_SUFFIXES = {
    ".json", ".m", ".md", ".ps1", ".py", ".toml", ".txt",
}
TEXT_PAYLOAD_NAMES = {"LICENSE"}


def _fail(messages: list[str]) -> None:
    for message in messages:
        print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def _metadata(archive: zipfile.ZipFile, names: set[str]):
    candidates = sorted(name for name in names if name.endswith(".dist-info/METADATA"))
    if len(candidates) != 1:
        _fail([f"expected one METADATA entry, found {len(candidates)}"])
    return BytesParser().parsebytes(archive.read(candidates[0]))


def _normalize_pe_timestamps(payload: bytes) -> bytes:
    """Remove reproducibility-neutral PE timestamps from a comparison copy."""
    if len(payload) < 0x40 or payload[:2] != b"MZ":
        return payload
    pe_offset = struct.unpack_from("<I", payload, 0x3C)[0]
    if pe_offset + 24 > len(payload) or payload[pe_offset:pe_offset + 4] != b"PE\0\0":
        return payload

    normalized = bytearray(payload)
    coff_offset = pe_offset + 4
    section_count = struct.unpack_from("<H", payload, coff_offset + 2)[0]
    optional_size = struct.unpack_from("<H", payload, coff_offset + 16)[0]
    normalized[coff_offset + 4:coff_offset + 8] = b"\0" * 4

    optional_offset = coff_offset + 20
    if optional_offset + optional_size > len(payload) or optional_size < 2:
        return bytes(normalized)
    optional_magic = struct.unpack_from("<H", payload, optional_offset)[0]
    directory_layout = {
        0x10B: (optional_offset + 92, optional_offset + 96),
        0x20B: (optional_offset + 108, optional_offset + 112),
    }.get(optional_magic)
    if directory_layout is None:
        return bytes(normalized)
    directory_count_offset, data_directory_offset = directory_layout
    if directory_count_offset + 4 > optional_offset + optional_size:
        return bytes(normalized)
    if struct.unpack_from("<I", payload, directory_count_offset)[0] <= 6:
        return bytes(normalized)

    debug_directory_entry = data_directory_offset + 6 * 8
    if debug_directory_entry + 8 > optional_offset + optional_size:
        return bytes(normalized)
    debug_rva, debug_size = struct.unpack_from("<II", payload, debug_directory_entry)
    if not debug_rva or debug_size < 28:
        return bytes(normalized)

    section_offset = optional_offset + optional_size
    debug_file_offset = None
    for section_index in range(section_count):
        entry = section_offset + section_index * 40
        if entry + 40 > len(payload):
            break
        virtual_size, virtual_address, raw_size, raw_offset = struct.unpack_from(
            "<IIII", payload, entry + 8
        )
        section_span = max(virtual_size, raw_size)
        if virtual_address <= debug_rva < virtual_address + section_span:
            candidate = raw_offset + (debug_rva - virtual_address)
            if candidate + debug_size <= len(payload):
                debug_file_offset = candidate
            break
    if debug_file_offset is None:
        return bytes(normalized)

    for entry in range(debug_file_offset, debug_file_offset + debug_size, 28):
        if entry + 28 > len(payload):
            break
        normalized[entry + 4:entry + 8] = b"\0" * 4
    return bytes(normalized)


def _normalized_payload(member: str, payload: bytes) -> bytes:
    path = PurePosixPath(member)
    if path.suffix.lower() == ".mexw64":
        return _normalize_pe_timestamps(payload)
    if (
        path.suffix.lower() in TEXT_PAYLOAD_SUFFIXES
        or path.name in TEXT_PAYLOAD_NAMES
    ):
        return payload.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return payload


def _source_payloads(source_manifest: dict[str, object]) -> dict[str, Path]:
    matlab_root = REPO_ROOT / "matlab"
    payloads = {
        str(PACKAGE_PREFIX / path.name): path
        for path in (PACKAGE_ROOT / "src" / "radia_optuna").iterdir()
        if path.is_file() and path.suffix in {".py", ".json"}
    }
    for source in (matlab_root / "+radia" / "+optuna").rglob("*.m"):
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


def verify(wheel: Path, *, source_fidelity: bool = True) -> dict[str, object]:
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
        source_payloads = _source_payloads(source_manifest) if source_fidelity else {}
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
            str(PACKAGE_PREFIX / "manifest.json"),
        }
        missing = sorted(expected_fixed - names)
        if missing:
            errors.append("missing required entries: " + ", ".join(missing))

        if source_fidelity:
            stale_payloads = sorted(
                member
                for member, source in source_payloads.items()
                if member in names
                and _normalized_payload(member, archive.read(member))
                != _normalized_payload(member, source.read_bytes())
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
        "artifact_integrity_verified": True,
        "source_fidelity_verified": source_fidelity,
        "source_fidelity_file_count": len(source_payloads),
    }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wheel", type=Path)
    parser.add_argument("--json", action="store_true", dest="as_json")
    parser.add_argument(
        "--artifact-only",
        action="store_true",
        help="verify the wheel's internal release contract without local-source bytes",
    )
    args = parser.parse_args()
    result = verify(args.wheel, source_fidelity=not args.artifact_only)
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
