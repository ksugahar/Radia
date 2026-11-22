"""Read the structure and part properties of an immutable Onshape assembly."""

from math import isfinite
from urllib.parse import quote, urlencode

from .discovery import ID


def _source_key(row):
    return (row.get("documentId"), row.get("elementId"), row.get("documentMicroversion"),
            row.get("fullConfiguration") or row.get("configuration"))


def _definitions(root, nested):
    entries = {}
    for definition in [root, *nested]:
        key = _source_key(definition)
        if key in entries:
            raise ValueError("Duplicate Onshape assembly definition.")
        instances = {}
        for row in definition["instances"]:
            instance_id = row["id"]
            if instance_id in instances:
                raise ValueError("Duplicate instance ID in an Onshape assembly definition.")
            instances[instance_id] = row
        entries[key] = instances
    return entries


def _resolve(root, definitions, path):
    current = root
    for index, instance_id in enumerate(path):
        instance = definitions.get(current, {}).get(instance_id)
        if instance is None or instance.get("suppressed") or index == len(path) - 1:
            return instance, current
        if instance["type"] != "Assembly":
            return None, current
        current = _source_key(instance)
    return None, current


def _numbers(values, count):
    if not isinstance(values, list) or len(values) < count:
        return None
    result = values[:count]
    return result if all(type(v) in (int, float) and isfinite(v) for v in result) else None


def _properties(body, transform=None):
    volumes = _numbers(body.get("volume"), 1)
    center = _numbers(body.get("centroid"), 3)
    if volumes is None or center is None or volumes[0] <= 0:
        return None
    if transform is not None:
        center = [sum(transform[i + j] * center[j] for j in range(3)) + transform[i + 3]
                  for i in (0, 4, 8)]
    mass = _numbers(body.get("mass"), 1) if body.get("hasMass") else None
    return {"volume_mm3": volumes[0] * 1e9, "centroid_mm": [v * 1e3 for v in center],
            "mass_kg": mass[0] if mass else None}


def _part_properties(request, instance, linking_document):
    did, eid, mid, _ = _source_key(instance)
    part_id, configuration = instance.get("partId"), instance.get("configuration")
    if (not all(isinstance(v, str) and ID.fullmatch(v) for v in (did, eid, mid))
            or not isinstance(part_id, str) or part_id in {"", ".", ".."}
            or not isinstance(configuration, str)):
        return None
    query = {"configuration": configuration, "useMassPropertyOverrides": "false"}
    if linking_document != did:
        if not isinstance(linking_document, str) or not ID.fullmatch(linking_document):
            return None
        query["linkDocumentId"] = linking_document
    endpoint = f"/parts/d/{did}/m/{mid}/e/{eid}/partid/{quote(part_id, safe='')}/massproperties"
    response = request(f"{endpoint}?{urlencode(query)}")
    if response.get("microversionId") not in (None, mid):
        raise ValueError("Part mass properties came from a different microversion.")
    return response["bodies"].get(part_id)


def inspect_assembly(request, document_id, element_id, version_id, configuration, max_part_requests):
    endpoint = f"/assemblies/d/{document_id}/v/{version_id}/e/{element_id}"
    configured = {"configuration": configuration} if configuration else {}
    query = urlencode({**configured, "includeMateFeatures": "true",
                       "includeNonSolids": "false", "excludeSuppressed": "false"})
    result = request(f"{endpoint}?{query}")
    root, nested = result["rootAssembly"], result.get("subAssemblies", [])
    root_key = _source_key(root)
    definitions = _definitions(root, nested)
    parts, unresolved, cache = [], [], {}
    for occurrence in root["occurrences"]:
        path = occurrence["path"]
        instance, owner = _resolve(root_key, definitions, path)
        if instance is None:
            unresolved.append({"path": path, "reason": "unknown_occurrence"})
            continue
        if instance.get("suppressed") or instance["type"] != "Part":
            continue
        transform = _numbers(occurrence["transform"], 16)
        if transform is None or len(occurrence["transform"]) != 16 or transform[12:] != [0, 0, 0, 1]:
            unresolved.append({"path": path, "reason": "invalid_placement"})
            continue
        row = {"path": path, "name": instance.get("name"), "part_id": instance.get("partId"),
               "document_id": instance.get("documentId"), "element_id": instance.get("elementId"),
               "source_microversion_id": instance.get("documentMicroversion"),
               "configuration": instance.get("configuration"), "hidden": occurrence.get("hidden"),
               "transform_row_major_m": transform, "properties": None}
        parts.append(row)
        key = (_source_key(instance), instance.get("configuration"), instance.get("partId"), owner[0])
        if key not in cache:
            if len(cache) >= max_part_requests:
                unresolved.append({"path": path, "reason": "part_request_limit"})
                continue
            cache[key] = _part_properties(request, instance, owner[0])
        body = cache[key]
        row["properties"] = _properties(body, transform) if body is not None else None
        if row["properties"] is None:
            unresolved.append({"path": path, "reason": "part_properties_unavailable"})

    mass_query = f"?{urlencode(configured)}" if configured else ""
    aggregate = _properties(request(f"{endpoint}/massproperties{mass_query}"))
    return {
        "source_microversion_id": root.get("documentMicroversion"),
        "coordinate_frame": "top_level_assembly",
        "assembly_properties": aggregate,
        "definitions": [
            {"document_id": key[0], "element_id": key[1], "source_microversion_id": key[2],
             "configuration": key[3], "instances": list(definitions[key].values()),
             "mates": [feature for feature in definition.get("features", [])
                       if feature.get("featureType", "").casefold() in {"mate", "materelation"}]}
            for definition in [root, *nested] for key in [_source_key(definition)]
        ],
        "part_occurrences": parts,
        "unresolved": unresolved,
    }
