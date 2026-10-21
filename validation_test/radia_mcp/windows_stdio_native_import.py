"""Exercise lazy native imports over real stdio, saving results before judging.

Run with an isolated environment containing radia-mcp, mcp, NumPy and SciPy.
Use the externally admitted CPU affinity and memory cap for this validation.
The default reader exercises the fix; --reader sdk is a diagnostic control.
"""

from __future__ import annotations

import argparse
import asyncio
import json
import os
from pathlib import Path
import sys
import time


def serve(args):
    from mcp.server.fastmcp import FastMCP

    server = FastMCP("stdio-native-import-validation")

    @server.tool()
    def status() -> dict:
        return {"ready": True, "numpy_loaded": "numpy" in sys.modules}

    @server.tool()
    def native_import() -> dict:
        import numpy
        import scipy.linalg

        return {"numpy": numpy.__version__, "scipy": scipy.__version__}

    @server.tool()
    def echo(text: str) -> str:
        return text

    @server.tool()
    def expected_error() -> str:
        raise ValueError("expected validation error")

    if args.reader == "sdk":
        server.run(transport="stdio")
    else:
        from radia_mcp.common.windows_stdio import run_stdio

        run_stdio(server)


async def probe(args, report, save):
    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client

    env = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", NGS_NUM_THREADS="1")
    params = StdioServerParameters(
        command=args.python,
        args=[str(Path(__file__).resolve()), "--server", "--reader", args.reader],
        env=env,
    )
    with (args.out / "server-stderr.log").open("w", encoding="utf-8") as log:
        async with stdio_client(params, errlog=log) as (read, write):
            async with ClientSession(read, write) as session:
                await session.initialize()
                status = await session.call_tool("status", {})
                assert not status.isError
                cold = json.loads(status.content[0].text)
                assert cold["ready"] and not cold["numpy_loaded"]
                report["stage"] = "before_native_import"
                save()
                imported = await session.call_tool("native_import", {})
                assert not imported.isError
                report["versions"] = json.loads(imported.content[0].text)
                text = "field-確認-" + "x" * 20000
                echoed = await session.call_tool("echo", {"text": text})
                assert not echoed.isError and echoed.content[0].text == text
                error = await session.call_tool("expected_error", {})
                assert error.isError
                recovered = await session.call_tool("status", {})
                assert not recovered.isError
                report.update(accepted=True, stage="complete", cold_start=True,
                              large_utf8_frame=True, error_recovery=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--server", action="store_true")
    parser.add_argument("--reader", choices=("available", "sdk"), default="available")
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--timeout", type=float, default=20)
    args = parser.parse_args()
    if args.server:
        serve(args)
        return 0
    if args.out is None or args.timeout <= 0:
        parser.error("--out and a positive --timeout are required")
    args.out.mkdir(parents=True, exist_ok=False)
    begin = time.monotonic()
    report = dict(accepted=False, reader=args.reader, stage="starting",
                  python=sys.version.split()[0], library_threads=1)

    def save():
        report["elapsed_seconds"] = time.monotonic() - begin
        with (args.out / "result.json").open("w", encoding="utf-8") as stream:
            json.dump(report, stream, indent=2)
            stream.flush()
            os.fsync(stream.fileno())

    save()
    try:
        asyncio.run(asyncio.wait_for(probe(args, report, save), args.timeout))
    except Exception as error:
        report.update(accepted=False, error_type=type(error).__name__, error=str(error))
    finally:
        save()
    return 0 if report["accepted"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
