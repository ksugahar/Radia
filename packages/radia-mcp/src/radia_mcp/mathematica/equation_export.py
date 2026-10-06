"""Evaluate once and hand the resulting equation to the native editor."""

from __future__ import annotations

import math
import json
import re
from pathlib import Path


def mathematica_export_equation(
    expression: str,
    action: str = "tex",
    output_path: str | None = None,
    target: str = "office",
    executable: str | None = None,
    timeout: int = 60,
) -> dict:
    """Calculate a Wolfram expression and export its result through EqnEdit64.

    Actions: tex returns TeX and InputForm without changing the clipboard;
    save writes a UTF-8 .tex file that EqnEdit64 can open; copy publishes an
    editable Office equation (target=office) or another supported clipboard
    target; render writes PNG/EMF for preview. Images are not editable math.
    Matrices use MatrixForm for display only; InputForm retains the result.
    Each request evaluates once in a fresh kernel; prior notebook variables
    are unavailable. Include definitions in expression, e.g. a Module.
    Use one expression; join statements with ';' or wrap them in Module.
    Unevaluated symbolic heads are allowed; Null and failure markers are not.
    A failed calculation never invokes the editor. Backend failures retain
    the generated TeX so the caller can inspect unsupported notation.
    """
    from .tools import _parse_json_output, mathematica_evaluate

    try:
        if not isinstance(expression, str) or not expression.strip():
            raise ValueError("expression must be a non-empty Wolfram expression")
        if any(ord(char) < 32 and char not in "\t\r\n" for char in expression):
            raise ValueError("expression contains unsupported control characters")
        if action not in {"tex", "save", "copy", "render"}:
            raise ValueError("action must be tex, save, copy, or render")
        if isinstance(timeout, bool) or not math.isfinite(timeout) or not 0 < timeout <= 300:
            raise ValueError("timeout must be in (0, 300]")
        destination = None
        if action in {"save", "render"}:
            if not output_path:
                raise ValueError("output_path is required for save/render")
            destination = Path(output_path).expanduser().resolve()
            allowed = {".tex"} if action == "save" else {".png", ".emf"}
            if destination.suffix.lower() not in allowed:
                raise ValueError(f"output extension must be one of {sorted(allowed)}")
        elif output_path is not None:
            raise ValueError("output_path is only used by save/render")
        if action == "copy" and target not in {"office", "powerpoint", "google-slides", "png"}:
            raise ValueError("unsupported clipboard target")
    except (TypeError, ValueError, OSError) as exc:
        return {"ok": False, "stage": "input", "error": str(exc)}

    literal = json.dumps(expression, ensure_ascii=False)
    code = (
        'Module[{held, result, display, reason = ""}, Block[{$MessageList = {}},'
        f"held = Check[ToExpression[{literal}, InputForm, HoldComplete], $Failed];"
        'result = If[held === $Failed, $Failed, If[Length[held] != 1, reason = "Use one expression; join statements with semicolons or wrap them in Module"; $Failed, Check[ReleaseHold[held], $Failed]]];'
        'If[result === Null || !FreeQ[result, $Failed | $Aborted], Print[ExportString[<|"ok"->False,'
        '"error"->(reason <> "; Invalid equation result: " <> ToString[result, InputForm] <> "; messages: " <> ToString[$MessageList, InputForm])|>, "RawJSON", "Compact"->True]],'
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
                error="Only centered TeXForm arrays can be exported to EqnEdit64",
            )
            return result
        result["wolfram_tex"] = result["tex"]
        result["tex"] = re.sub(pattern, lambda match: r"\begin{matrix}", result["tex"]).replace(
            r"\end{array}", r"\end{matrix}"
        )
    if action == "tex":
        return result
    if action == "save":
        try:
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes((result["tex"] + "\n").encode("utf-8"))
            result["output_path"] = str(destination)
        except OSError as exc:
            result.update(ok=False, stage="save", error=str(exc))
        return result
    from radia_mcp.presentation._equation_cli import (
        presentation_copy_equation,
        presentation_render_equation,
    )

    if action == "copy":
        backend = presentation_copy_equation(result["tex"], target=target, executable=executable)
    else:
        backend = presentation_render_equation(
            result["tex"],
            str(destination),
            image_format=destination.suffix.lower()[1:],
            executable=executable,
        )
    result["backend_result"] = backend
    if not backend["ok"]:
        result.update(ok=False, stage="editor", error=backend.get("error", "EqnEdit64 failed"))
    return result
