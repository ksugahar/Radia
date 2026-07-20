"""mathematica MCP tools — wolframscript subprocess bridge."""
from __future__ import annotations

import os
import json
import subprocess
import shutil
import tempfile
import time
from pathlib import Path


# ============================================================
# wolframscript.exe 探索
# ============================================================
# 1. 環境変数 WOLFRAMSCRIPT_PATH (override)
# 2. PATH (where wolframscript)
# 3. Windows: 既知の install dir からスキャン
# 4. Mac/Linux: 一般的な場所
def _find_wolframscript() -> str | None:
    env = os.environ.get("WOLFRAMSCRIPT_PATH")
    if env and Path(env).is_file():
        return env

    found = shutil.which("wolframscript")
    if found:
        return found

    # Windows 既知 install dir を新しいバージョン優先で探索
    # 14.x からは `Wolfram Research\Wolfram\<ver>\` にリブランド、
    # WolframScript は別パッケージとして `\WolframScript\` 直下にも入る。
    candidates: list[str] = []
    standalone = Path(r"C:\Program Files\Wolfram Research\WolframScript\wolframscript.exe")
    if standalone.is_file():
        candidates.append(str(standalone))
    for base in [
        r"C:\Program Files\Wolfram Research\Wolfram",
        r"C:\Program Files (x86)\Wolfram Research\Wolfram",
        r"C:\Program Files\Wolfram Research\Mathematica",
        r"C:\Program Files (x86)\Wolfram Research\Mathematica",
    ]:
        if Path(base).is_dir():
            for sub in sorted(Path(base).iterdir(), reverse=True):
                exe = sub / "wolframscript.exe"
                if exe.is_file():
                    candidates.append(str(exe))

    # Mac
    for p in [
        "/Applications/Mathematica.app/Contents/MacOS/wolframscript",
        "/usr/local/bin/wolframscript",
    ]:
        if Path(p).is_file():
            candidates.append(p)

    # Linux
    for p in ["/usr/bin/wolframscript", "/opt/Mathematica/Executables/wolframscript"]:
        if Path(p).is_file():
            candidates.append(p)

    return candidates[0] if candidates else None


_WOLFRAMSCRIPT_PATH: str | None = None  # cached


def _wolframscript_path() -> str:
    """Cached wolframscript locator. Raises FileNotFoundError if not found."""
    global _WOLFRAMSCRIPT_PATH
    if _WOLFRAMSCRIPT_PATH is None:
        _WOLFRAMSCRIPT_PATH = _find_wolframscript()
    if not _WOLFRAMSCRIPT_PATH:
        raise FileNotFoundError(
            "wolframscript.exe not found. Install Mathematica or set "
            "WOLFRAMSCRIPT_PATH env var."
        )
    return _WOLFRAMSCRIPT_PATH


def _read_process_text(path: str) -> str:
    try:
        return Path(path).read_text(encoding="utf-8", errors="replace").strip()
    except OSError:
        return ""


def _run_wolframscript(arguments: list[str], timeout: int) -> dict:
    """Run wolframscript without inheriting the MCP stdio transport."""
    started = time.perf_counter()
    timed_out = False
    with tempfile.TemporaryDirectory(prefix="mma_ws_") as temp_dir:
        stdout_path = os.path.join(temp_dir, "out.txt")
        stderr_path = os.path.join(temp_dir, "err.txt")
        try:
            with open(stdout_path, "wb") as stdout_file, open(stderr_path, "wb") as stderr_file:
                result = subprocess.run(
                    arguments,
                    stdin=subprocess.DEVNULL,
                    stdout=stdout_file,
                    stderr=stderr_file,
                    timeout=timeout,
                )
            exit_code = result.returncode
        except subprocess.TimeoutExpired:
            timed_out = True
            exit_code = -2
        stdout = _read_process_text(stdout_path)
        stderr = _read_process_text(stderr_path)
    return {
        "result": stdout,
        "stderr": stderr,
        "exit_code": exit_code,
        "timed_out": timed_out,
        "duration_seconds": round(time.perf_counter() - started, 6),
    }


def _parse_json_output(text: str) -> tuple[object | None, str]:
    """Parse a JSON result, tolerating informational lines before it."""
    candidates = [text.strip()]
    candidates.extend(
        line.strip()
        for line in reversed(text.splitlines())
        if line.lstrip().startswith(("{", "["))
    )
    errors: list[str] = []
    for candidate in dict.fromkeys(candidates):
        if not candidate:
            continue
        try:
            return json.loads(candidate), ""
        except json.JSONDecodeError as exc:
            errors.append(str(exc))
    return None, errors[-1] if errors else "empty output"


# ============================================================
# MCP tools (mathematica_*)
# ============================================================
def mathematica_evaluate(code: str, timeout: int = 60) -> dict:
    """Evaluate Wolfram Language code via wolframscript and return result.

    The canonical low-level bridge. Claude generates Wolfram Language code
    from natural language intent, calls this tool, parses the textual
    output, and explains the result to the user.

    Parameters
    ----------
    code : str
        Wolfram Language expression(s). Multi-line accepted (use ; or
        newlines to chain). Examples:
            "1 + 1"
            "FullSimplify[Sin[x]^2 + Cos[x]^2]"
            "Quantity[1, \"Tesla\"] // UnitConvert"
            "Curl[{0, x*y, z}, {x, y, z}]"
            "D[Sin[x], x]"
    timeout : int
        Wall-clock timeout in seconds (default 60). Long symbolic
        computations may need 120-300.

    Returns
    -------
    dict with keys:
        result      : str   — stdout (textual evaluation result)
        stderr      : str   — stderr (warnings, error messages)
        exit_code   : int   — wolframscript exit code (0 = success)
        wolframscript_path : str — actual exe used
        timed_out   : bool  — True iff killed by timeout
        code        : str   — echoed input

    Notes
    -----
    - Output is Wolfram's text form (InputForm by default for headless).
      Use ToString[expr, InputForm] inside `code` if needed for clarity.
    - For LaTeX output, wrap your expression: TeXForm[ToExpression[...]]
    - Empty stdout + zero exit + populated stderr usually means a
      license-server / kernel startup issue.
    - For multiple-line scripts, use newlines or ; in code; both work.
    - First call after boot pays a 1-2 sec kernel startup cost. Cached
      after that for repeat calls within the same session is NOT
      currently implemented — each call spawns a fresh kernel.

    Wolframscript Gotchas (hard-won — read before generating long scripts!)
    ----------------------------------------------------------------------
    1. **Underscore `_` in identifiers = Pattern**. `KphiG_num` is parsed as
       `KphiG` followed by `Pattern[num, _]` (matches anything, names it
       "num"), NOT as a single variable. This silently corrupts assignments,
       breaks Solve, and crashes JSON Export with cryptic
       "Expression Pattern cannot be exported" errors.
       → Use camelCase or no separator: `KphiGnum`, `KphiGtest`, `myVarName`.

    2. **Multi-line expressions outside brackets DON'T continue**. With:
           phiExpr = a + b + c
                  + d + e
       wolframscript treats this as TWO statements: `phiExpr = a + b + c`
       then `+d + e` (a no-op). The `+` at line start is NOT continuation.
       → Keep RHS on one line, OR wrap in brackets `( ... )`/`[ ... ]`/`{ ... }`
       (continuation works inside those).

    3. **Scientific notation `*^` is fragile in wolframscript**. `4 Pi*^-7`
       (interactive: 4π·1e-7) sometimes fails parsing in -file mode with
       "Invalid syntax in or before ...".
       → Prefer explicit `4*Pi*10^-7`.

    4. **Pattern functions can interact badly with Solve**. Defining
           phi[s_, z_, a_, b_, ...] := a + b s + ...
       and calling `phi[sa, za, a, b, ...]` inside Table+Solve can return
       an empty solution list for reasons that are hard to reproduce.
       → Prefer a plain expression + ReplaceAll:
           phiExpr = a + b s + c z + d s z;
           eqns = Table[pVars[[k]] == (phiExpr /. {s -> ..., z -> ...}), ...];
           Solve[eqns, coefVars]

    5. **Heavy `Together`/`FullSimplify` on big rational expressions OOMs**.
       Substituting a symbolic Vandermonde inverse into a 9x9 quadratic
       form and calling Together can blow up to 8+ GB and never finish.
       → Output unsimplified `Expand`'d expressions; let downstream codegen
       (C++ / Python) handle the symbolic representation. Often you don't
       need beautiful simplified output, just correct InputForm strings.

    6. **Avoid symbolic Vandermonde substitution in element-matrix
       derivations**. For a high-order FE basis with N monomials, output
       only the NxN coefficient-form matrix K_phi (one polynomial integral
       per entry). Compute the Vandermonde V (N nodes -> N coefficients)
       and apply V^{-T} K_phi V^{-1} **numerically at runtime** in the
       calling code. This keeps Mathematica fast (seconds) and avoids the
       expression explosion that comes from substituting a symbolic
       Vandermonde inverse.

    7. **JSON Export silently writes 0 bytes on Pattern leakage**. If you
       see "Expression Pattern cannot be exported as JSON" in stderr, the
       output file is empty even though "Saved: ..." prints. Cause is
       usually #1 (an underscore-as-pattern leaked into an Association
       value). Hunt down stray underscores in identifier names.

    Examples
    --------
    >>> mathematica_evaluate("1 + 1")
    {"result": "2", "stderr": "", "exit_code": 0, ...}

    >>> mathematica_evaluate("FullSimplify[Sin[x]^2 + Cos[x]^2]")
    {"result": "1", "stderr": "", "exit_code": 0, ...}

    >>> mathematica_evaluate(
    ...   "Module[{e=E0 Cos[k x - omega t], b=B0 Cos[k x - omega t]},"
    ...   " FullSimplify[Curl[{0,e,0},{x,y,z}]+D[{0,0,b},t]==0,"
    ...   " k E0==omega B0]]"
    ... )
    """
    try:
        ws = _wolframscript_path()
    except FileNotFoundError as e:
        return {
            "result": "",
            "stderr": str(e),
            "exit_code": -1,
            "wolframscript_path": "",
            "timed_out": False,
            "code": code,
        }

    # IMPORTANT (2026-06-27 fix): redirect the child's stdout/stderr to TEMP FILES, not pipes
    # (capture_output=True).  `wolframscript -code` launches a grandchild WolframKernel that
    # INHERITS the stdout/stderr write handles; if that kernel lingers (it does in the
    # console-less, job-object MCP-server launch context -- though NOT from an interactive
    # shell, which is why a bare `subprocess.run(..., capture_output=True)` works when tested by
    # hand but hangs inside mcp-server-mathematica), the pipe never reaches EOF and
    # `subprocess.run`'s internal communicate() BLOCKS until `timeout` even though wolframscript
    # itself already exited and printed the answer.  Writing to files means we wait only for the
    # wolframscript PROCESS to exit (prompt) and then read the files; a lingering grandchild
    # holding a file handle does not block us.  Also pin stdin=DEVNULL so the kernel can never
    # touch the server's JSON-RPC stdin.
    execution = _run_wolframscript([ws, "-code", code], timeout)
    return {
        **execution,
        "wolframscript_path": ws,
        "code": code,
    }


def mathematica_run_script(
    script_path: str,
    timeout: int = 300,
    result_format: str = "auto",
    require_ok: bool = False,
) -> dict:
    """Run a local Wolfram Language script and return structured results.

    This is the verification-oriented companion to ``mathematica_evaluate``.
    It runs a tracked ``.wls``, ``.wl``, or ``.m`` file with
    ``wolframscript -file`` so a large batch of identities shares one kernel
    startup.  Scripts may print a final JSON object such as
    ``{"ok": true, "checks": {...}, "failures": []}``.

    Args:
        script_path: Local path to a Wolfram Language script.
        timeout: Wall-clock timeout in seconds.
        result_format: ``"text"``, ``"json"``, or ``"auto"``.  Auto parses
            JSON when possible and otherwise preserves text output.
        require_ok: If true, the parsed JSON must contain boolean ``ok: true``.

    Returns:
        Process status, stdout/stderr, elapsed time, parsed JSON (when
        available), and ``verification_ok``.  ``ok`` is false on process,
        JSON, or requested verification failure.
    """
    if timeout < 1:
        return {"ok": False, "error": "timeout must be >= 1 second", "script_path": script_path}
    format_name = result_format.lower().strip()
    if format_name not in {"auto", "json", "text"}:
        return {
            "ok": False,
            "error": "result_format must be 'auto', 'json', or 'text'",
            "script_path": script_path,
        }
    script = Path(script_path).expanduser().resolve()
    if not script.is_file():
        return {"ok": False, "error": "script file not found", "script_path": str(script)}
    if script.suffix.lower() not in {".wls", ".wl", ".m"}:
        return {
            "ok": False,
            "error": "script must use .wls, .wl, or .m",
            "script_path": str(script),
        }
    try:
        ws = _wolframscript_path()
    except FileNotFoundError as exc:
        return {
            "ok": False,
            "error": str(exc),
            "script_path": str(script),
            "wolframscript_path": "",
        }

    execution = _run_wolframscript([ws, "-file", str(script)], timeout)
    process_ok = execution["exit_code"] == 0 and not execution["timed_out"]
    parsed = None
    json_error = ""
    if format_name != "text":
        parsed, json_error = _parse_json_output(execution["result"])
    verification_ok = parsed.get("ok") if isinstance(parsed, dict) and isinstance(parsed.get("ok"), bool) else None
    ok = process_ok
    if format_name == "json" and parsed is None:
        ok = False
    if require_ok and verification_ok is not True:
        ok = False
    error = ""
    if not process_ok:
        error = execution["stderr"] or "wolframscript failed"
    elif format_name == "json" and parsed is None:
        error = f"script output is not valid JSON: {json_error}"
    elif require_ok and verification_ok is None:
        error = "parsed JSON does not contain boolean key 'ok'"
    elif require_ok and verification_ok is False:
        failures = parsed.get("failures", []) if isinstance(parsed, dict) else []
        error = f"verification failed: {failures}"
    return {
        "ok": ok,
        "process_ok": process_ok,
        "verification_ok": verification_ok,
        "parsed": parsed,
        "json_error": json_error if parsed is None and format_name != "text" else "",
        "script_path": str(script),
        "wolframscript_path": ws,
        "error": error,
        **execution,
    }


def mathematica_status() -> dict:
    """Diagnostic: check wolframscript availability + version + license.

    Returns
    -------
    dict with keys:
        available   : bool  — wolframscript.exe found AND license OK
        path        : str   — wolframscript.exe path (or "" if not found)
        version     : str   — Mathematica $Version (or "" if license fail)
        license_ok  : bool  — True iff "1+1" returns "2"
        license_msg : str   — error message if license_ok is False
    """
    try:
        ws = _wolframscript_path()
    except FileNotFoundError:
        return {
            "available": False,
            "path": "",
            "version": "",
            "license_ok": False,
            "license_msg": "wolframscript.exe not found",
        }

    # 1+1 license check.  The FIRST wolframscript call after the MCP server
    # boots is a COLD START (kernel launch + license activation) that can take
    # 20-30 s -- far over the old 15 s budget -- and a single-seat license
    # already held by an interactive Wolfram Kernel / Mathematica GUI makes
    # wolframscript WAIT for the seat.  Both surface as an EMPTY result, which
    # must NOT be misreported as a license failure (the symptom was
    # available=False, license_msg="no output" even though `wolframscript -code
    # "1+1"` returns 2 in ~5 s once warm).  Use a generous timeout and label a
    # timeout distinctly from a genuine license/output error.
    test = mathematica_evaluate("1 + 1", timeout=60)
    license_ok = test["result"].strip() == "2"
    license_msg = ""
    if not license_ok:
        if test.get("timed_out"):
            license_msg = (
                "wolframscript produced no output within 60 s -- usually a slow "
                "cold-start (kernel + license activation) or a single-seat license "
                "already held by an open Wolfram Kernel / Mathematica GUI, NOT a "
                "license failure. Close other kernels and retry, or call "
                "mathematica_evaluate directly (it has its own timeout)."
            )
        else:
            license_msg = test["stderr"] or test["result"] or "no output"

    # version (kernel is warm now if 1+1 passed; keep headroom anyway)
    ver_result = mathematica_evaluate("$Version", timeout=30)
    version = ver_result["result"].strip().strip('"') if ver_result["exit_code"] == 0 else ""

    return {
        "available": license_ok,
        "path": ws,
        "version": version,
        "license_ok": license_ok,
        "license_msg": license_msg,
    }


# ============================================================
# High-level helpers (simplify / to_tex / vector_calc / unit_convert /
# solve / integrate / differentiate / check_identity).
# Defined in helpers.py and re-exported here so that grant_writing-style
# auto-registration in __init__.py picks them up.
# ============================================================
from .helpers import (  # noqa: E402, F401
    mathematica_simplify,
    mathematica_to_tex,
    mathematica_check_identity,
    mathematica_check_identities,
    mathematica_verification_guide,
    mathematica_vector_calc,
    mathematica_unit_convert,
    mathematica_solve,
    mathematica_integrate,
    mathematica_differentiate,
)
