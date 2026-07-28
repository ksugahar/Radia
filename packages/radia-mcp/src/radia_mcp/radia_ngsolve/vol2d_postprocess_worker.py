"""Owned worker for the dimension-2 postprocessing MCP tool."""

from __future__ import annotations

import json
import sys

from .vol2d_postprocess import analyze_vol2d_postprocess


def main() -> None:
    request = json.loads(sys.stdin.buffer.read().decode("utf-8"))
    result = analyze_vol2d_postprocess(request)
    print(json.dumps(result, sort_keys=True, separators=(",", ":")))


if __name__ == "__main__":
    main()
