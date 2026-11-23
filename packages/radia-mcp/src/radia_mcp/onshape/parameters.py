"""Read fixed-version Onshape variable tables and configuration metadata."""

import math
import re
from urllib.parse import urlencode


_MEASURE = re.compile(
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"
    r"(?:\s+([^\d\W].*|[°%µ].*))?"
)


class UnknownConfigurationParameters(ValueError):
    def __init__(self, ids):
        self.ids = sorted(ids)
        super().__init__("Unknown Onshape configuration parameter IDs.")


def _variable(row):
    formatted = row.get("value")
    match = _MEASURE.fullmatch(formatted.strip()) if isinstance(formatted, str) else None
    number = float(match[1]) if match else None
    if number is not None and not math.isfinite(number):
        number = None
    return {
        "name": row["name"], "expression": row.get("expression"), "type": row.get("type"),
        "value_text": formatted, "value": number, "unit": match[2] if number is not None else None,
    }


def read_parameters(request, export_args, parameter_values):
    did, vid, eid = (export_args[key] for key in ("document_id", "version_id", "element_id"))
    path = f"/d/{did}/v/{vid}/e/{eid}"
    metadata = request(f"/elements{path}/configuration")
    definitions = metadata["configurationParameters"]
    if not isinstance(definitions, list) or metadata.get("microversionSkew"):
        raise ValueError("Onshape did not return configuration definitions for this fixed version.")
    encoded = export_args["configuration"]
    if parameter_values:
        available = {entry["parameterId"] for entry in definitions}
        unknown = parameter_values.keys() - available
        if unknown:
            raise UnknownConfigurationParameters(unknown)
        # This operation computes an encoding using the chosen version. It
        # does not call the workspace endpoint for changing configuration.
        options = [{"parameterId": key, "parameterValue": value}
                   for key, value in sorted(parameter_values.items())]
        result = request(f"/elements/d/{did}/e/{eid}/configurationencodings?{urlencode({'versionId': vid})}",
                         {"parameters": options})
        encoded = result["encodedId"]
        if not isinstance(encoded, str) or not encoded:
            raise ValueError("Onshape did not return a configuration encodedId.")

    tables = None
    if export_args["element_type"] == "part_studio":
        query = {"includeValuesAndReferencedVariables": "true"}
        if encoded:
            query["configuration"] = encoded
        raw = request(f"/variables{path}/variables?{urlencode(query)}")
        if not isinstance(raw, list):
            raise ValueError("Onshape did not return a variable table list.")
        tables = [{"variable_studio_reference": table.get("variableStudioReference"),
                   "variables": [_variable(row) for row in table["variables"]]} for table in raw]
    return {
        "export_args": {**export_args, "configuration": encoded},
        "source_microversion_id": metadata.get("sourceMicroversion"),
        "configuration": {"parameters": definitions, "encodedId": encoded or None,
                          "current_configuration": metadata.get("currentConfiguration")},
        "variable_tables": tables,
    }
