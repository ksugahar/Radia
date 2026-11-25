"""Version-pinned Onshape CAD export with local provenance."""

import base64
import hashlib
import json
import os
import re
import shutil
import tempfile
import time
from http.client import HTTPException
from pathlib import Path
from typing import Annotated, Literal, get_args
from urllib.error import HTTPError, URLError
from urllib.parse import urljoin, urlsplit
from urllib.request import HTTPRedirectHandler, Request, build_opener

from mcp.server.fastmcp import FastMCP
from mcp.types import ToolAnnotations

from pydantic import Field

from radia_mcp.common import register_status_tool
from radia_mcp.common.server_hardening import ANN_READONLY_WEB

from .assembly import inspect_assembly
from .discovery import ID as _ID, ORIGIN as _ORIGIN, discover, parse_document_url
from .parameters import UnknownConfigurationParameters, read_parameters
from .shaded_views import ViewName, render_shaded_views


mcp = FastMCP("onshape-mcp")
# Every tool calls cad.onshape.com; export and shaded views also write local files.
_ANN_WRITES_WEB = ToolAnnotations(
    readOnlyHint=False, destructiveHint=False, idempotentHint=False, openWorldHint=True,
)
_API = "/api/v11"
_ELEMENTS = {"part_studio": "partstudios", "assembly": "assemblies"}
_TRANSLATION_FORMATS = {"x_t": "PARASOLID", "sat": "ACIS"}
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


def _authorization():
    access, secret = os.getenv("ONSHAPE_ACCESS_KEY"), os.getenv("ONSHAPE_SECRET_KEY")
    if not access or not secret:
        return None
    return "Basic " + base64.b64encode(f"{access}:{secret}".encode()).decode("ascii")


def _auth_error():
    return _error("auth", "Set ONSHAPE_ACCESS_KEY and ONSHAPE_SECRET_KEY in the server environment.")


@mcp.tool(annotations=ANN_READONLY_WEB)
async def onshape_discover(
    operation: Literal["parse_url", "documents", "elements", "versions"],
    url: str = "",
    document_id: str = "",
    search: str = "",
    offset: Annotated[int, Field(strict=True, ge=0)] = 0,
    limit: Annotated[int, Field(strict=True, gt=0, le=100)] = 20,
) -> dict:
    """Explore authorized cad.onshape.com documents, without writing to Onshape.

    parse_url is local. documents lists one page of My Documents (optional search).
    versions accepts a document_id or URL and returns version URLs.
    elements accepts a workspace, version or microversion URL and lists CAD tabs.
    For a version URL each tab supplies export_args for onshape_export_cad. Add
    file_format and export_dir when exporting. A tab URL selects only
    that tab. Only version URLs provide export_args.
    """
    import anyio

    return await anyio.to_thread.run_sync(
        _discover, operation, url, document_id, search, offset, limit
    )


def _discover(operation, url, document_id, search, offset, limit):
    try:
        reference = parse_document_url(url) if url else None
        if operation == "parse_url":
            if reference is None:
                return _error("input", "Provide a document URL to parse.")
            return {"status": "ok", "reference": reference}
        if operation == "elements" and (reference is None or not reference["reference_id"]):
            return _error("input", "Elements require a workspace, version or microversion URL.")
        if operation == "versions" and reference is None and not _ID.fullmatch(document_id):
            return _error("input", "Versions require a valid document ID or URL.")
        if reference is not None and document_id and document_id != reference["document_id"]:
            return _error("input", "Document ID does not match the URL.")
        if operation == "documents" and (url or document_id):
            return _error("input", "Document search accepts a search string rather than a document reference.")
    except ValueError:
        return _error("input", "Provide a valid cad.onshape.com document URL with one configuration at most.")

    return _read_onshape("discovery", lambda request: discover(
        request, operation, reference, document_id, search, offset, limit
    ))


def _read_onshape(stage, action):
    authorization = _authorization()
    if authorization is None:
        return _auth_error()
    opener = build_opener(_NoRedirect())

    def request(path, options=None):
        payload = None if options is None else json.dumps(options).encode("utf-8")
        with _open(opener, _ORIGIN + _API + path, authorization, payload) as response:
            ceiling = 4 * 1024 * 1024
            payload = response.read(ceiling + 1)
            if len(payload) > ceiling:
                raise ValueError("Onshape read response too large.")
            return json.loads(payload)

    try:
        return {"status": "ok", **action(request)}
    except HTTPError as exc:
        exc.close()
        return {**_error(stage, "Onshape rejected the read request."), "http_status": exc.code}
    except (URLError, OSError, HTTPException):
        return _error(stage, "Onshape read or local file operation failed.")
    except ImportError:
        return _error(stage, "This operation requires the build123d extra.")
    except UnknownConfigurationParameters as exc:
        return {**_error("input", str(exc)), "unknown_parameter_ids": exc.ids}
    except (ValueError, KeyError, TypeError, AttributeError):
        return _error(stage, "Unexpected or oversized Onshape read response.")


@mcp.tool(annotations=ANN_READONLY_WEB)
async def onshape_read_parameters(
    document_id: str,
    element_id: str,
    version_id: str,
    element_type: Literal["part_studio", "assembly"],
    configuration: str = "",
    parameter_values: dict[str, str] | None = None,
) -> dict:
    """Read parameters from a fixed Onshape version, using discovery export_args.

    Part Studios return formatted variable values, parsed numeric values and
    units where available. Assemblies return configuration metadata.
    Pass an existing configuration encodedId, or supply parameter IDs and
    option values to compute one with Onshape's version-pinned encoder. The
    encoder uses POST and does not change the document. Returned export_args
    can be passed to onshape_export_cad with file_format and export_dir.
    """
    import anyio

    return await anyio.to_thread.run_sync(
        _read_parameters, document_id, element_id, version_id, element_type,
        configuration, parameter_values,
    )


def _read_parameters(document_id, element_id, version_id, element_type,
                     configuration, parameter_values):
    if not all(_ID.fullmatch(value) for value in (document_id, element_id, version_id)):
        return {**_error("input", "Provide document, element and immutable version IDs (24 lowercase hex digits)."),
                "cache": "not_eligible"}
    if element_type not in _ELEMENTS or (configuration and parameter_values):
        return _error("input", "Provide a CAD element type and only one configuration representation.")
    source = {"document_id": document_id, "element_id": element_id, "version_id": version_id,
              "element_type": element_type, "configuration": configuration}
    return _read_onshape("parameters", lambda request: read_parameters(request, source, parameter_values))


@mcp.tool(annotations=_ANN_WRITES_WEB)
async def onshape_render_shaded_views(
    document_id: str,
    element_id: str,
    version_id: str,
    element_type: Literal["part_studio", "assembly"],
    export_dir: str,
    configuration: str = "",
    views: list[ViewName] | None = None,
) -> dict:
    """Save shaded PNG previews for a fixed version on cad.onshape.com.

    Use discovery or parameter export_args. Default views are isometric and
    front; other presets are back, top, bottom, left and right. Each view
    costs one read request. Requires Pillow from the build123d extra.
    Use STEP readback for geometry measurements and comparisons.
    """
    import anyio

    selected = ["isometric", "front"] if views is None else views
    if not all(_ID.fullmatch(value) for value in (document_id, element_id, version_id)):
        return {**_error("input", "Provide document, element and immutable version IDs (24 lowercase hex digits)."),
                "cache": "not_eligible"}
    if (element_type not in _ELEMENTS or not export_dir.strip() or not selected
            or len(set(selected)) != len(selected)
            or any(view not in get_args(ViewName) for view in selected)):
        return _error("input", "Provide a CAD element type, output directory and unique named views.")
    source = {"document_id": document_id, "element_id": element_id, "version_id": version_id,
              "element_type": element_type, "configuration": configuration}

    def action(request):
        images = render_shaded_views(
            request, _ELEMENTS[element_type], document_id, element_id, version_id,
            configuration, selected, export_dir,
        )
        return {**source, "views": images}

    return await anyio.to_thread.run_sync(_read_onshape, "shaded_views", action)


@mcp.tool(annotations=ANN_READONLY_WEB)
async def onshape_inspect_assembly(
    document_id: str,
    element_id: str,
    version_id: str,
    element_type: Literal["assembly"],
    configuration: str = "",
    max_part_requests: Annotated[int, Field(strict=True, gt=0, le=100)] = 32,
) -> dict:
    """Read fixed-version assembly structure, mates and mass properties.

    Pass discovery export_args for an assembly. Part occurrences carry the
    top-level row-major transform (translation in m), volume in mm³ and
    transformed centroid in mm. Part properties use each source microversion
    and its effective configuration. Repeated sources reuse one GET per call.
    Missing or budget-limited part properties are explicitly unresolved.
    Match occurrences to STEP labels explicitly when comparing measurements.
    """
    import anyio

    if not all(_ID.fullmatch(value) for value in (document_id, element_id, version_id)):
        return {**_error("input", "Provide document, element and immutable version IDs (24 lowercase hex digits)."),
                "cache": "not_eligible"}
    if element_type != "assembly" or type(max_part_requests) is not int or not 1 <= max_part_requests <= 100:
        return _error("input", "Provide an assembly and a part request limit from 1 through 100.")
    source = {"document_id": document_id, "element_id": element_id, "version_id": version_id,
              "element_type": element_type, "configuration": configuration}

    def action(request):
        result = inspect_assembly(request, document_id, element_id, version_id,
                                  configuration, max_part_requests)
        return {**source, **result}

    return await anyio.to_thread.run_sync(_read_onshape, "assembly", action)


def _step_cache_key(endpoint, options):
    request = [_ORIGIN, _API, endpoint, "step", options]
    return hashlib.sha256(json.dumps(request, sort_keys=True).encode()).hexdigest()


def _publish_cached_step(entry, key, root, max_bytes):
    """Verify the bytes copied to the caller's destination before publication."""
    manifest = entry / "provenance.json"
    record = json.loads(manifest.read_text(encoding="utf-8"))
    if (record.get("cache_key") != key or record.get("format") != "step"
            or record.get("file") != "model.step"):
        raise ValueError("Cache provenance mismatch")
    root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".onshape-", dir=root) as staging:
        folder = Path(staging) / "bundle"
        folder.mkdir()
        artifact = folder / "model.step"
        digest = hashlib.sha256()
        size = 0
        with (entry / artifact.name).open("rb") as reader, artifact.open("wb") as writer:
            for chunk in iter(lambda: reader.read(1024 * 1024), b""):
                size += len(chunk)
                if size > max_bytes:
                    raise ValueError("Cached STEP exceeds the requested size limit")
                digest.update(chunk)
                writer.write(chunk)
        if digest.hexdigest() != record.get("file_sha256") or record.get("step_sha256") != record["file_sha256"]:
            raise ValueError("Cached STEP hash mismatch")
        (folder / manifest.name).write_text(json.dumps(record, indent=2) + "\n",
                                            encoding="utf-8", newline="\n")
        published = root / f"onshape-step-{Path(staging).name.removeprefix('.onshape-')}"
        folder.rename(published)
    return _export_result(record, published, "hit")


def _store_cached_step(entry, folder):
    entry.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".onshape-", dir=entry.parent) as staging:
        bundle = Path(staging) / "bundle"
        shutil.copytree(folder, bundle)
        try:
            bundle.rename(entry)
        except FileExistsError:
            # A concurrent export already published this immutable request.
            pass


def _export_result(record, published, cache):
    result = {**record, "status": "ok", "cache": cache,
              "file": str(published / record["file"]),
              "provenance": str(published / "provenance.json")}
    if record["format"] == "step":
        result["step"] = result["file"]
    return result


@mcp.tool(annotations=_ANN_WRITES_WEB)
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
    Set ONSHAPE_STEP_CACHE_DIR to enable a private, user-local fixed-version
    STEP cache. A verified hit can be served without Onshape credentials.
    """
    import anyio

    return await anyio.to_thread.run_sync(
        _export_cad, document_id, element_id, version_id, element_type,
        export_dir, file_format, configuration, step_unit, max_download_bytes, poll_timeout_s,
    )


def _export_cad(document_id, element_id, version_id, element_type, export_dir,
                file_format, configuration, step_unit, max_download_bytes, poll_timeout_s):
    if not all(_ID.fullmatch(value) for value in (document_id, element_id, version_id)):
        return {**_error("input", "Provide document, element and immutable version IDs (24 lowercase hex digits)."),
                "cache": "not_eligible"}
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
    cache_root = os.environ.get("ONSHAPE_STEP_CACHE_DIR") if file_format == "step" else None
    entry = None
    cache = "disabled" if file_format == "step" else "not_eligible"
    if cache_root:
        key = _step_cache_key(endpoint, options)
        entry = Path(cache_root) / key
        cache = "miss"
        try:
            if entry.exists():
                return _publish_cached_step(entry, key, Path(export_dir).resolve(), max_download_bytes)
        except (OSError, ValueError, TypeError, AttributeError, KeyError):
            return {**_error("cache", "Cached STEP failed verification or could not be copied."),
                    "cache": "error"}
    authorization = _authorization()
    if authorization is None:
        return {**_auth_error(), "cache": cache}
    opener = build_opener(_NoRedirect())
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
            if entry is not None:
                record["cache_key"] = key
            manifest = folder / "provenance.json"
            manifest.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8", newline="\n")
            published = root / f"onshape-{file_format}-{Path(staging).name.removeprefix('.onshape-')}"
            if entry is not None:
                stage = "cache"
                _store_cached_step(entry, folder)
            folder.rename(published)
        return _export_result(record, published, cache)
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
    server_name="onshape-mcp", tool_name="onshape_status",
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
        # onshape_reload_code exists only in editable development installs.
        if {tool.name for tool in tools} - {"onshape_reload_code"} != {
            "onshape_export_cad", "onshape_discover", "onshape_read_parameters",
            "onshape_render_shaded_views", "onshape_inspect_assembly", "onshape_status",
        }:
            raise SystemExit("Unexpected Onshape research tool catalog")
        print("Onshape research tool registration passed")
        return
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
