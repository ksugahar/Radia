"""Evaluate a Wolfram expression once and export its result as TeX."""

from __future__ import annotations

import json
import math
import re
import unicodedata
from pathlib import Path


def mathematica_export_equation(
    expression: str,
    action: str = "tex",
    output_path: str | None = None,
    timeout: int = 60,
) -> dict:
    """Calculate a Wolfram expression and export its result as TeX.

    Actions: tex returns TeX and InputForm; save also writes a UTF-8 .tex
    file that an equation editor can open.
    Matrices use MatrixForm for display only; InputForm retains the result.
    Each request evaluates once in a fresh kernel; prior notebook variables
    are unavailable. Include definitions in expression, e.g. a Module.
    Use one expression; join statements with ';' or wrap them in Module.
    Unevaluated symbolic heads are allowed; Null and failure markers are not.
    A failed calculation never writes a file.
    """
    from .tools import _parse_json_output, mathematica_evaluate

    try:
        if not isinstance(expression, str) or not expression.strip():
            raise ValueError("expression must be a non-empty Wolfram expression")
        if any(
            unicodedata.category(char) == "Cc" and char not in "\t\r\n"
            for char in expression
        ):
            raise ValueError("expression contains unsupported control characters")
        if action not in {"tex", "save"}:
            raise ValueError("action must be tex or save")
        if isinstance(timeout, bool) or not math.isfinite(timeout) or not 0 < timeout <= 300:
            raise ValueError("timeout must be in (0, 300]")
        destination = None
        if action == "save":
            if not output_path:
                raise ValueError("output_path is required for save")
            destination = Path(output_path).expanduser().resolve()
            if destination.suffix.lower() != ".tex":
                raise ValueError("output extension must be .tex")
        elif output_path is not None:
            raise ValueError("output_path is only used by save")
    except (TypeError, ValueError, OSError) as exc:
        return {"ok": False, "stage": "input", "error": str(exc)}

    literal = json.dumps(expression.strip(), ensure_ascii=False)
    code = (
        'Module[{held, result, display, reason = ""}, Block[{$MessageList = {}},'
        f"held = Check[ToExpression[{literal}, InputForm, HoldComplete], $Failed];"
        'result = If[held === $Failed, $Failed, If[Length[held] != 1, reason = "Use one expression; join statements with semicolons or wrap them in Module"; $Failed, Check[ReleaseHold[held], $Failed]]];'
        'If[result === Null || !FreeQ[result, $Failed | $Aborted], Print[ExportString[<|"ok"->False,'
        '"error"->(If[reason === "", "", reason <> "; "] <> "Invalid equation result: " <> ToString[result, InputForm] <> "; messages: " <> ToString[$MessageList, InputForm])|>, "RawJSON", "Compact"->True]],'
        "display = If[MatrixQ[result], MatrixForm[result], result];"
        'Print[ExportString[<|"ok"->True, "tex"->ToString[TeXForm[display]],'
        '"input_form"->ToString[result, InputForm]|>, "RawJSON", "Compact"->True]]]]]'
    )
    raw = mathematica_evaluate(code, timeout=timeout)
    if raw.get("exit_code") != 0 or raw.get("timed_out"):
        return {
            "ok": False,
            "stage": "evaluation",
            "error": raw.get("stderr") or "Wolfram evaluation failed",
            "raw": raw,
        }
    payload, error = _parse_json_output(raw.get("result", ""))
    if (
        not isinstance(payload, dict)
        or payload.get("ok") is not True
        or not isinstance(payload.get("tex"), str)
        or not payload["tex"].strip()
        or not isinstance(payload.get("input_form"), str)
    ):
        return {
            "ok": False,
            "stage": "evaluation",
            "error": error
            or (payload.get("error") if isinstance(payload, dict) else None)
            or "Wolfram did not return a valid equation",
            "raw": raw,
        }
    result = {
        "ok": True,
        "action": action,
        "tex": payload["tex"],
        "input_form": payload["input_form"],
        "raw": raw,
    }
    # Centered array and matrix have the same layout; the native editor
    # otherwise renders TeXForm's {cc} column specification as visible text.
    array_start = r"\begin{array}"
    if array_start in result["tex"]:
        pattern = r"\\begin\{array\}\{c+\}"
        if len(re.findall(pattern, result["tex"])) != result["tex"].count(array_start):
            result.update(
                ok=False,
                stage="conversion",
                error="Only centered TeXForm arrays can be exported as matrix TeX",
            )
            return result
        result["wolfram_tex"] = result["tex"]
        result["tex"] = re.sub(pattern, lambda match: r"\begin{matrix}", result["tex"]).replace(
            r"\end{array}", r"\end{matrix}"
        )
    if action == "tex":
        return result
    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes((result["tex"] + "\n").encode("utf-8"))
        result["output_path"] = str(destination)
    except OSError as exc:
        result.update(ok=False, stage="save", error=str(exc))
    return result
