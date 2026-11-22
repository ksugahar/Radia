"""Version-pinned Onshape CAD export with local provenance."""

import base64
import hashlib
import json
import os
import re
import tempfile
import time
from http.client import HTTPException
from pathlib import Path
from typing import Annotated, Literal
from urllib.error import HTTPError, URLError
from urllib.parse import urljoin, urlsplit
from urllib.request import HTTPRedirectHandler, Request, build_opener

from mcp.server.fastmcp import FastMCP

from radia_mcp.common import register_status_tool
from pydantic import Field


mcp = FastMCP("onshape-research-mcp")
_ORIGIN = "https://cad.onshape.com"
_API = "/api/v11"
_ELEMENTS = {"part_studio": "partstudios", "assembly": "assemblies"}
_TRANSLATION_FORMATS = {"x_t": "PARASOLID", "sat": "ACIS"}
_ID = re.compile(r"[0-9a-f]{24}")
_STORAGE = re.compile(r"(?:[a-z0-9.-]+\.)?s3(?:[.-][a-z0-9-]+)*\.amazonaws\.com")
_DEFAULT_MAX_BYTES = 256 * 1024 * 1024
_STEP_UNITS = {"METER", "CENTIMETER", "MILLIMETER", "INCH", "FOOT", "YARD"}


class _NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


def _open(opener, url, authorization=None, payload=None, *, accept="application/json", timeout=30):
    headers = {"Accept": accept}
    if authorization is not None:
        headers["Authorization"] = authorization
    if payload is not None:
        headers["Content-Type"] = "application/json"
    request = Request(url, data=payload, headers=headers)
    return opener.open(request, timeout=timeout)


def _validate_download(first, last, file_format):
    head, tail = first.lstrip(), last.rstrip()
    if file_format == "step":
        valid = head.startswith(b"ISO-10303-21;") and tail.endswith(b"END-ISO-10303-21;")
    elif file_format == "sat":
        valid = (re.match(rb"\d+[ \t]+\d+[ \t]+\d+[ \t]+\d+[ \t]*\r?\n", head)
                 and tail.endswith(b"End-of-ACIS-data"))
    else:
        valid = not head.startswith(b"PK") and re.search(rb"(?m)^\*\*PARASOLID\b", head[:4096])
    if not valid:
        raise ValueError(f"Download does not match the expected .{file_format} file.")


def _download(opener, url, authorization, target, file_format="step", max_bytes=_DEFAULT_MAX_BYTES):
    for _ in range(5):
        location = urlsplit(url)
        same_origin = location.netloc == "cad.onshape.com"
        storage = bool(_STORAGE.fullmatch(location.netloc))
        if location.scheme != "https" or not (same_origin or storage):
            raise ValueError("Unrecognized Onshape download redirect destination")
        try:
            response = _open(opener, url, authorization if same_origin else None,
                             accept="application/octet-stream")
        except HTTPError as exc:
            if exc.code not in {301, 302, 303, 307, 308} or not exc.headers.get("Location"):
                raise
            with exc:
                url = urljoin(url, exc.headers["Location"])
            continue
        with response, target.open("xb") as output:
            length = getattr(response, "headers", {}).get("Content-Length")
            if length is not None and not re.fullmatch(r"[0-9]{1,20}", length):
                raise ValueError("Invalid download Content-Length.")
            if length is not None and int(length) > max_bytes:
                raise ValueError(f"Download exceeds the {max_bytes}-byte limit.")
            digest = hashlib.sha256()
            first, last, received = b"", b"", 0
            while chunk := response.read(min(1024 * 1024, max_bytes - received + 1)):
                received += len(chunk)
                if received > max_bytes:
                    raise ValueError(f"Download exceeds the {max_bytes}-byte limit.")
                if not first:
                    first = chunk[:4096]
                last = (last + chunk)[-4096:]
                digest.update(chunk)
                output.write(chunk)
            _validate_download(first, last, file_format)
            return digest.hexdigest()
    raise ValueError("Too many Onshape download redirects")


def _error(stage, message):
    return {"status": "error", "stage": stage, "error": message}


@mcp.tool()
async def onshape_export_cad(
    document_id: str,
    element_id: str,
    version_id: str,
    element_type: str,
    export_dir: str,
    file_format: Literal["step", "x_t", "sat"],
    configuration: str = "",
    step_unit: str | None = None,
    max_download_bytes: Annotated[int, Field(strict=True, gt=0)] = _DEFAULT_MAX_BYTES,
    poll_timeout_s: Annotated[int, Field(strict=True, gt=0)] = 120,
) -> dict:
    """Export one local CAD file from a research Onshape version.

    Use part_studio or assembly, and STEP, Parasolid text or ACIS text via
    file_format. Supply the Onshape configuration encodedId, or omit it to
    use the fixed version's default. Only cad.onshape.com is supported.
    Set ONSHAPE_ACCESS_KEY and ONSHAPE_SECRET_KEY in the server environment.
    STEP declaration readback requires the build123d extra. Optional STEP
    units are METER, CENTIMETER, MILLIMETER, INCH, FOOT and YARD. Omitting
    step_unit uses the server default. Download and polling default to
    256 MiB and 120 seconds, respectively. Both accept positive integers.
    """
    import anyio

    return await anyio.to_thread.run_sync(
        _export_cad, document_id, element_id, version_id, element_type,
        export_dir, file_format, configuration, step_unit, max_download_bytes, poll_timeout_s,
    )


def _export_cad(document_id, element_id, version_id, element_type, export_dir,
                file_format, configuration, step_unit, max_download_bytes, poll_timeout_s):
    if not all(_ID.fullmatch(value) for value in (document_id, element_id, version_id)):
        return _error("input", "Provide document, element and immutable version IDs (24 lowercase hex digits).")
    if element_type not in _ELEMENTS or not export_dir.strip():
        return _error("input", "Provide a part_studio or assembly and a local export directory.")
    if file_format not in {"step", *_TRANSLATION_FORMATS}:
        return _error("input", "Use file_format step, x_t or sat.")
    if step_unit is not None and (file_format != "step" or step_unit not in _STEP_UNITS):
        return _error("input", "step_unit must be a supported Onshape STEP unit.")
    if any(type(value) is not int or value <= 0
           for value in (max_download_bytes, poll_timeout_s)):
        return _error("input", "Download and polling limits must be positive integers.")
    try:
        poll_timeout_s = float(poll_timeout_s)
    except OverflowError:
        return _error("input", "Polling deadline is outside the supported range.")
    access, secret = os.getenv("ONSHAPE_ACCESS_KEY"), os.getenv("ONSHAPE_SECRET_KEY")
    if not access or not secret:
        return _error("auth", "Set ONSHAPE_ACCESS_KEY and ONSHAPE_SECRET_KEY in the server environment.")
    authorization = "Basic " + base64.b64encode(f"{access}:{secret}".encode()).decode("ascii")
    opener = build_opener(_NoRedirect())
    source = f"/{_ELEMENTS[element_type]}/d/{document_id}/v/{version_id}/e/{element_id}"
    options = {"storeInDocument": False, "grouping": True,
               "notifyUser": False, "triggerAutoDownload": False}
    if file_format == "step":
        endpoint = source + "/export/step"
        options["isYAxisUp"] = False
        if configuration:
            options["advancedParams"] = {"configuration": configuration}
        if step_unit is not None:
            options["stepUnit"] = step_unit
    else:
        endpoint = source + "/translations"
        options["formatName"] = _TRANSLATION_FORMATS[file_format]
        if configuration:
            options["configuration"] = configuration
        if file_format == "x_t":
            options["parasolidMode"] = "TEXT"
    stage = "export"
    translation_id = None
    try:
        with _open(opener, _ORIGIN + _API + endpoint, authorization,
                   json.dumps(options).encode()) as response:
            state = json.load(response)
        translation_id = state["id"]
        if not isinstance(translation_id, str) or not _ID.fullmatch(translation_id):
            return _error(stage, "Onshape returned an invalid translation ID.")
        stage = "poll"
        deadline = time.monotonic() + poll_timeout_s
        poll_url = f"{_ORIGIN}{_API}/translations/{translation_id}"
        while state.get("requestState") == "ACTIVE":
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                return {**_error(stage, "Onshape polling deadline exceeded."),
                        "translation_id": translation_id}
            with _open(opener, poll_url, authorization, timeout=min(30, remaining)) as response:
                state = json.load(response)
            if time.monotonic() >= deadline:
                return {**_error(stage, "Onshape polling deadline exceeded."),
                        "translation_id": translation_id}
            if state.get("requestState") == "ACTIVE":
                time.sleep(min(2, max(0, deadline - time.monotonic())))
        if state.get("requestState") != "DONE":
            return {**_error(stage, "Onshape translation failed or did not finish within the polling window."),
                    "translation_id": translation_id}
        files = state.get("resultExternalDataIds")
        if not isinstance(files, list) or len(files) != 1 or not isinstance(files[0], str) or not _ID.fullmatch(files[0]):
            return _error("download", f"Expected one external .{file_format} result from Onshape.")
        stage = "download"
        root = Path(export_dir).resolve()
        root.mkdir(parents=True, exist_ok=True)
        url = f"{_ORIGIN}{_API}/documents/d/{document_id}/externaldata/{files[0]}"
        record = {
            "schema_version": 1, "source": "onshape", "document_id": document_id,
            "element_id": element_id, "version_id": version_id, "element_type": element_type,
            "configuration": configuration or None, "translation_id": translation_id,
            "configuration_mode": "explicit" if configuration else "version_default",
            "api_version": _API.rsplit("/", 1)[-1], "export_options": options,
            "format": file_format,
        }
        with tempfile.TemporaryDirectory(prefix=".onshape-", dir=root) as staging:
            folder = Path(staging) / "bundle"
            folder.mkdir()
            artifact = folder / f"model.{file_format}"
            digest = _download(opener, url, authorization, artifact, file_format, max_download_bytes)
            stage = "provenance"
            record.update(file=artifact.name, file_sha256=digest)
            if file_format == "step":
                from radia_mcp.build123d.step_units import read_step_length_units

                record.update(requested_step_unit=step_unit,
                              declared_length_units=read_step_length_units(artifact),
                              measurement_length_unit="mm")
                record.update(step=artifact.name, step_sha256=digest)
            manifest = folder / "provenance.json"
            manifest.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8", newline="\n")
            published = root / f"onshape-{file_format}-{Path(staging).name.removeprefix('.onshape-')}"
            folder.rename(published)
        result = {"status": "ok", **record, "file": str(published / artifact.name),
                  "provenance": str(published / manifest.name)}
        if file_format == "step":
            result["step"] = result["file"]
        return result
    except HTTPError as exc:
        with exc:
            diagnostic = {
                400: "Check version, configuration and export options.",
                401: "Check the API credentials.",
                403: "Check export permissions for the fixed version.",
                404: "Check the supplied IDs or whether the export is still available.",
                429: "Onshape rate limit reached.",
            }.get(exc.code, "Onshape rejected the request.")
            if stage == "download" and exc.code == 403:
                diagnostic = "Download access denied. The signed download link may have expired."
            result = {**_error(stage, diagnostic), "http_status": exc.code}
            if stage == "poll":
                result["translation_id"] = translation_id
            return result
    except (URLError, OSError, HTTPException):
        if stage == "poll":
            return {**_error(stage, "Onshape polling connection failed or timed out."),
                    "translation_id": translation_id}
        return _error(stage, "Onshape connection or local file operation failed.")
    except ImportError:
        return _error(stage, "STEP unit inspection requires the build123d extra.")
    except ValueError as exc:
        return _error(stage, str(exc))
    except (KeyError, AttributeError, TypeError):
        return _error(stage, "Unexpected Onshape response.")


register_status_tool(
    mcp,
    server_name="onshape-research-mcp", tool_name="onshape_status",
    description="Research-only, version-pinned Onshape STEP, Parasolid and ACIS export",
    subpackage="radia_mcp.onshape",
    related_servers=["build123d-mcp"],
    optional_deps=["build123d"],
)


def main():
    import argparse
    import asyncio

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selftest", action="store_true", help="Check local tool registration without network.")
    if parser.parse_args().selftest:
        tools = asyncio.run(mcp.list_tools())
        if {tool.name for tool in tools} != {"onshape_export_cad", "onshape_status"}:
            raise SystemExit("Unexpected Onshape research tool catalog")
        print("Onshape research tool registration passed")
        return
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
