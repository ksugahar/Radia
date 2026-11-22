"""Read pinned Onshape shaded previews as locally verified PNG artifacts."""

import base64
import hashlib
import io
import math
import tempfile
from pathlib import Path
from typing import Literal
from urllib.parse import urlencode


ViewName = Literal["isometric", "front", "back", "top", "bottom", "left", "right"]
_DIMENSION = 512
_ISOMETRIC = (
    1 / math.sqrt(2), 1 / math.sqrt(2), 0, 0,
    -1 / math.sqrt(6), 1 / math.sqrt(6), 2 / math.sqrt(6), 0,
    1 / math.sqrt(3), -1 / math.sqrt(3), 1 / math.sqrt(3), 0,
)
_VIEW_MATRICES = {
    "isometric": ",".join(f"{component:.15g}" for component in _ISOMETRIC),
    **{name: name for name in ("front", "back", "top", "bottom", "left", "right")},
}


def _png(response, image_module):
    images = response["images"]
    if not isinstance(images, list) or len(images) != 1:
        raise ValueError("Expected one shaded PNG in the Onshape response.")
    encoded = images[0]
    if isinstance(encoded, str):
        payload = base64.b64decode(encoded, validate=True)
    elif isinstance(encoded, list) and all(isinstance(chunk, str) for chunk in encoded):
        # The generated API schema describes the image as an array of byte strings.
        payload = b"".join(base64.b64decode(chunk, validate=True) for chunk in encoded)
    else:
        raise ValueError("Onshape returned an unexpected shaded image encoding.")
    if not payload.startswith(b"\x89PNG\r\n\x1a\n"):
        raise ValueError("The shaded view is not a PNG.")
    try:
        with image_module.open(io.BytesIO(payload), formats=["PNG"]) as image:
            if not all(0 < size <= _DIMENSION for size in image.size):
                raise ValueError("The shaded PNG exceeds the requested dimensions.")
            image.verify()
    except (OSError, SyntaxError, image_module.DecompressionBombError) as exc:
        raise ValueError("Onshape returned an invalid or oversized shaded PNG.") from exc
    return payload


def render_shaded_views(request, collection, document_id, element_id, version_id,
                        configuration, views, export_dir):
    from PIL import Image  # Check the optional dependency before any API calls.

    endpoint = f"/{collection}/d/{document_id}/v/{version_id}/e/{element_id}/shadedviews"
    root = Path(export_dir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    artifacts = {}
    with tempfile.TemporaryDirectory(prefix=".onshape-views-", dir=root) as staging:
        folder = Path(staging) / "bundle"
        folder.mkdir()
        for name in views:
            options = {"viewMatrix": _VIEW_MATRICES[name], "outputWidth": _DIMENSION,
                       "outputHeight": _DIMENSION, "pixelSize": 0, "useAntiAliasing": "true"}
            if configuration:
                options["configuration"] = configuration
            payload = _png(request(f"{endpoint}?{urlencode(options)}"), Image)
            image = folder / f"{name}.png"
            image.write_bytes(payload)
            artifacts[name] = {"path": image.name, "sha256": hashlib.sha256(payload).hexdigest()}
        published = root / Path(staging).name.removeprefix(".")
        folder.rename(published)
    for artifact in artifacts.values():
        artifact["path"] = str(published / artifact["path"])
    return artifacts
