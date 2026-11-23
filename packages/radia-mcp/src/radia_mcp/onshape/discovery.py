"""Bounded, read-only Onshape document discovery."""

import re
from urllib.parse import parse_qs, urlencode, urlsplit


ORIGIN = "https://cad.onshape.com"
ID = re.compile(r"[0-9a-f]{24}")
_DOCUMENT_URL = re.compile(
    rf"/documents/(?P<document_id>{ID.pattern})"
    rf"(?:/(?P<reference_kind>[wvm])/(?P<reference_id>{ID.pattern})"
    rf"(?:/e/(?P<element_id>{ID.pattern}))?)?/?"
)
_ELEMENT_TYPES = {"PARTSTUDIO": "part_studio", "ASSEMBLY": "assembly"}


def parse_document_url(url):
    parsed = urlsplit(url)
    match = _DOCUMENT_URL.fullmatch(parsed.path)
    if parsed.scheme != "https" or parsed.netloc != "cad.onshape.com" or not match:
        raise ValueError("Provide a valid https://cad.onshape.com/documents URL.")
    configurations = parse_qs(parsed.query, keep_blank_values=True).get("configuration", [])
    if len(configurations) > 1:
        raise ValueError("Provide at most one Onshape configuration in the URL.")
    fields = match.groupdict()
    return {**fields, "immutable": fields["reference_kind"] in {"v", "m"},
            "configuration": configurations[0] if configurations else ""}


def _id(value):
    if not isinstance(value, str) or not ID.fullmatch(value):
        raise ValueError("Invalid Onshape ID.")
    return value


def _items(value):
    if not isinstance(value, list) or not all(isinstance(item, dict) for item in value):
        raise ValueError("Unexpected Onshape list response.")
    return value


def _document_url(document_id):
    return f"{ORIGIN}/documents/{document_id}"


def discover(get, operation, reference, document_id, search, offset, limit):
    """Use one GET for a typed operation. The transport never follows next URLs."""
    if operation == "documents":
        query = urlencode({"filter": 0, "q": search, "offset": offset, "limit": limit})
        page = get(f"/documents?{query}")
        items = _items(page["items"])
        documents = []
        for item in items:
            if item.get("resourceType", "document").casefold() != "document":
                continue
            did = _id(item["id"])
            documents.append({"document_id": did, "name": item.get("name"), "url": _document_url(did)})
        return {"documents": documents, "next_offset": offset + len(items) if items and page.get("next") else None}

    did = _id(reference["document_id"] if reference is not None else document_id)
    if operation == "versions":
        query = urlencode({"offset": offset, "limit": limit})
        items = _items(get(f"/documents/d/{did}/versions?{query}"))
        versions = []
        for item in items:
            vid = _id(item["id"])
            versions.append({
                "version_id": vid, "name": item.get("name"), "microversion_id": item.get("microversion"),
                "url": f"{_document_url(did)}/v/{vid}",
            })
        return {"document_id": did, "versions": versions,
                "next_offset": offset + len(items) if len(items) >= limit else None}

    kind, rid = reference["reference_kind"], reference["reference_id"]
    items = _items(get(f"/documents/d/{did}/{kind}/{rid}/elements"))
    elements = []
    selected = reference["element_id"]
    for item in items:
        if selected is not None and item.get("id") != selected:
            continue
        element_type = _ELEMENT_TYPES.get(item.get("elementType"))
        if element_type is None or item.get("deleted"):
            continue
        eid = _id(item["id"])
        row = {"element_id": eid, "element_type": element_type, "name": item.get("name"),
               "microversion_id": item.get("microversionId")}
        if kind == "v":
            row["export_args"] = {"document_id": did, "element_id": eid, "version_id": rid,
                                  "element_type": element_type, "configuration": reference["configuration"]}
        elements.append(row)
    return {"reference": reference, "elements": elements}
