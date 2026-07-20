"""Mathematica high-level helpers.

mathematica_evaluate (低レベル subprocess bridge) の上に、数式工程
(数式 DB との往復、paper 執筆、Maxwell 検証、物理単位変換) で頻出
する操作を 1 関数 1 目的で wrap した薄いラッパー群。

LAB / kubota / 学生がよく行う操作:
  - simplify: FullSimplify  (式の整理)
  - to_tex:   TeXForm        (paper 用 LaTeX 文字列)
  - check_identity: 等式 LHS == RHS が常に成り立つか
  - vector_calc: Maxwell 系の Curl/Div/Grad/Laplacian
  - unit_convert: A/m ↔ T 等
  - solve / integrate / differentiate: 教科書的操作

設計方針:
- 各 helper は純粋に mathematica_evaluate を呼び、Wolfram 側で
  Wrap した式を組む。Python 側で Wolfram 構文を組み立てる。
- 戻りは {ok, value (str), tex (str, on demand), code, raw}。
  失敗時は ok=False + error。
"""

from __future__ import annotations

import json


def mathematica_evaluate(code: str, timeout: int = 60) -> dict:
    """Call the low-level bridge lazily to keep both import orders valid."""
    from .tools import mathematica_evaluate as _evaluate

    return _evaluate(code, timeout=timeout)


def _parse_json_output(text: str) -> tuple[object | None, str]:
    """Use the bridge JSON parser without creating an import cycle."""
    from .tools import _parse_json_output as _parse

    return _parse(text)


# -------------------------------------------------------------------
# 共通 result builder
# -------------------------------------------------------------------
def _wrap(code: str, expect_str: bool = False, timeout: int = 60) -> dict:
    r = mathematica_evaluate(code, timeout=timeout)
    out = {
        "ok": (r["exit_code"] == 0) and not r["timed_out"],
        "value": r["result"],
        "code": code,
        "raw": r,
    }
    if not out["ok"]:
        out["error"] = r["stderr"] or r["result"] or "wolframscript failed"
    return out


# -------------------------------------------------------------------
# 1. simplify
# -------------------------------------------------------------------
def mathematica_simplify(
    expression: str,
    assumptions: str | None = None,
    timeout: int = 60,
) -> dict:
    """Wolfram FullSimplify[expression, assumptions] を実行。

    Args:
        expression: Wolfram 式 (str)。例 "Sin[x]^2 + Cos[x]^2"
        assumptions: 仮定 (Wolfram Boolean expr)。例
                     "Element[k, Reals] && k > 0"
        timeout: 秒。複雑式は 120-300 推奨。

    Returns:
        {ok, value, code, raw}.  value は簡約後の式 (Wolfram 文字列)。

    Examples:
        >>> mathematica_simplify("Sin[x]^2 + Cos[x]^2")["value"]
        '1'
        >>> mathematica_simplify("Sqrt[x^2]", "x > 0")["value"]
        'x'
    """
    if assumptions:
        code = f"FullSimplify[{expression}, {assumptions}]"
    else:
        code = f"FullSimplify[{expression}]"
    return _wrap(code, timeout=timeout)


# -------------------------------------------------------------------
# 2. to_tex
# -------------------------------------------------------------------
def mathematica_to_tex(expression: str, timeout: int = 60) -> dict:
    """Wolfram の TeXForm を文字列で取得 (paper / 数式 DB 登録用)。

    `ToString[TeXForm[...]]` で改行や余分な引用符を除去。

    Args:
        expression: Wolfram 式 (str)。

    Returns:
        {ok, tex, value (=tex), code, raw}.
    """
    code = f"ToString[TeXForm[{expression}]]"
    out = _wrap(code, expect_str=True, timeout=timeout)
    if out["ok"]:
        out["tex"] = out["value"].strip()
    return out


# -------------------------------------------------------------------
# 3. check_identity
# -------------------------------------------------------------------
def mathematica_check_identity(
    lhs: str,
    rhs: str,
    assumptions: str | None = None,
    timeout: int = 120,
) -> dict:
    """式 LHS == RHS が常に成り立つかを FullSimplify で判定。

    eqnedt32 数式 DB に登録された式の Mathematica 検証用。

    Args:
        lhs, rhs: Wolfram 式 (str)。
        assumptions: 仮定。例 "Element[x, Reals] && k > 0".
        timeout: 秒。

    Returns:
        {ok, identical (bool|None), value, code, raw}.
        identical:
            True  → FullSimplify が True を返した (恒等式)
            False → FullSimplify が False or 簡約された別式を返した
            None  → 評価できなかった or 解析失敗
    """
    expr = f"FullSimplify[({lhs}) - ({rhs}) == 0"
    if assumptions:
        expr += f", {assumptions}"
    expr += "]"
    out = _wrap(expr, timeout=timeout)
    if out["ok"]:
        v = out["value"].strip()
        if v == "True":
            out["identical"] = True
        elif v == "False":
            out["identical"] = False
        else:
            out["identical"] = None  # 簡約失敗 or partial
    else:
        out["identical"] = None
    return out


def mathematica_check_identities(
    claims: list[dict],
    timeout: int = 180,
) -> dict:
    """Verify several named identities in one Mathematica kernel.

    Each claim is a mapping with ``name``, ``lhs``, ``rhs``, and optional
    ``assumptions`` keys.  Batching avoids one Wolfram kernel cold start per
    formula and returns a named pass/fail association suitable for course,
    paper, and solver-regression checks.
    """
    if not claims:
        return {"ok": False, "error": "claims must contain at least one identity", "claims": []}
    if len(claims) > 200:
        return {"ok": False, "error": "at most 200 identities may be checked at once", "claims": []}

    entries: list[str] = []
    names: list[str] = []
    for index, claim in enumerate(claims, start=1):
        if not isinstance(claim, dict):
            return {"ok": False, "error": f"claim {index} must be an object", "claims": names}
        name = str(claim.get("name", "")).strip()
        lhs = str(claim.get("lhs", "")).strip()
        rhs = str(claim.get("rhs", "")).strip()
        assumptions = str(claim.get("assumptions", "True")).strip() or "True"
        if not name or not lhs or not rhs:
            return {
                "ok": False,
                "error": f"claim {index} requires non-empty name, lhs, and rhs",
                "claims": names,
            }
        if name in names:
            return {"ok": False, "error": f"duplicate claim name: {name}", "claims": names}
        names.append(name)
        wolfram_name = json.dumps(name, ensure_ascii=False)
        entries.append(
            f"{wolfram_name} -> TrueQ[FullSimplify[({lhs}) == ({rhs}), "
            f"Assumptions -> ({assumptions})]]"
        )

    association = ",\n    ".join(entries)
    code = (
        "Module[{checks},\n"
        "  checks = Association[\n    " + association + "\n  ];\n"
        "  Print[ExportString[<|\"ok\" -> And @@ Values[checks], "
        "\"checks\" -> checks, \"failures\" -> Keys[Select[checks, Not]]|>, "
        "\"RawJSON\", \"Compact\" -> True]]\n"
        "]"
    )
    raw = mathematica_evaluate(code, timeout=timeout)
    process_ok = raw["exit_code"] == 0 and not raw["timed_out"]
    parsed, parse_error = _parse_json_output(raw["result"]) if process_ok else (None, "")
    if not process_ok:
        return {
            "ok": False,
            "error": raw["stderr"] or raw["result"] or "wolframscript failed",
            "claims": names,
            "code": code,
            "raw": raw,
        }
    if not isinstance(parsed, dict) or not isinstance(parsed.get("ok"), bool):
        return {
            "ok": False,
            "error": f"Mathematica did not return the expected JSON result: {parse_error}",
            "claims": names,
            "code": code,
            "raw": raw,
        }
    return {
        "ok": parsed["ok"],
        "checks": parsed.get("checks", {}),
        "failures": parsed.get("failures", []),
        "claims": names,
        "code": code,
        "raw": raw,
    }


def mathematica_verification_guide(topic: str = "electromagnetics") -> dict:
    """Return the recommended Mathematica verification workflow.

    The guide keeps agents from launching one Wolfram kernel per equation and
    records what symbolic checks can and cannot establish.  It is deliberately
    executable-tool oriented: use ``mathematica_check_identities`` for a small
    collection of independent claims and ``mathematica_run_script`` for a
    tracked course, paper, or solver verification suite.
    """
    normalized = topic.strip().lower().replace("_", "-")
    aliases = {
        "em": "electromagnetics",
        "electromagnetic": "electromagnetics",
        "course": "electromagnetics",
        "differential-forms": "differential-forms",
        "differentialforms": "differential-forms",
        "forms": "differential-forms",
        "paper": "paper",
        "publication": "paper",
    }
    normalized = aliases.get(normalized, normalized)
    if normalized not in {"electromagnetics", "differential-forms", "paper"}:
        return {
            "ok": False,
            "error": "topic must be electromagnetics, differential-forms, or paper",
            "topic": topic,
        }

    common = {
        "small_batch": "mathematica_check_identities",
        "tracked_suite": "mathematica_run_script",
        "single_exploration": "mathematica_evaluate",
        "rules": [
            "Batch named identities in one kernel instead of paying one cold start per formula.",
            "State assumptions explicitly, including real-valued variables, positivity, and nonzero denominators.",
            "Test governing equations, limiting cases, dimensions, and sign conventions separately.",
            "Treat a symbolic pass as equation evidence, not as proof that a Canvas, caption, or physical interpretation is correct.",
            "Keep reusable verification in a tracked .wls/.wl/.m file and require a final JSON object with ok, checks, and failures.",
        ],
        "json_contract": {
            "ok": True,
            "checks": {"namedCheck": True},
            "failures": [],
        },
    }

    topic_checks = {
        "electromagnetics": {
            "checks": [
                "Coordinate handedness and basis-vector cross products.",
                "Gradient, curl, divergence, and gauge-invariance identities.",
                "Biot-Savart or finite-wire formulas against symmetry-axis and far-field limits.",
                "Maxwell equations, constitutive substitutions, SI dimensions, and energy derivatives.",
                "Parameter-dependent formulas at several representative and boundary values.",
            ],
            "starter_claims": [
                {
                    "name": "curlGradient",
                    "lhs": "Curl[Grad[phi[x,y,z], {x,y,z}], {x,y,z}]",
                    "rhs": "{0,0,0}",
                    "assumptions": "Element[{x,y,z}, Reals]",
                },
                {
                    "name": "divergenceCurl",
                    "lhs": "With[{aa={a1[x,y,z],a2[x,y,z],a3[x,y,z]}}, Div[Curl[aa,{x,y,z}],{x,y,z}]]",
                    "rhs": "0",
                    "assumptions": "Element[{x,y,z}, Reals]",
                },
            ],
        },
        "differential-forms": {
            "checks": [
                "Exterior derivative nilpotency and pullback/exterior-derivative commutation.",
                "Orientation-sensitive signs, especially straight versus twisted forms under Kelvin inversion.",
                "Hodge-star scaling with metric, material law, form degree, and spatial dimension.",
                "Stokes pairings, Whitney degrees of freedom, and weak-form integration by parts.",
                "Energy-density top forms separately from integrated scalar energy.",
            ],
            "companion_server": "mcp-server-differential-forms",
        },
        "paper": {
            "checks": [
                "Every displayed equation used for a numerical claim.",
                "Equivalent formulations under the exact assumptions stated in the manuscript.",
                "Asymptotic limits, parameter domains, units, and sign conventions.",
                "Numerical spot checks independent of the production implementation.",
            ],
            "reporting": "Publish the named checks and failures; keep runtime and raw output in the reproducibility record.",
        },
    }
    return {"ok": True, "topic": normalized, **common, **topic_checks[normalized]}


# -------------------------------------------------------------------
# 4. vector_calc (Maxwell 系)
# -------------------------------------------------------------------
def mathematica_vector_calc(
    expression: str,
    operation: str,
    variables: str = "{x, y, z}",
    timeout: int = 60,
) -> dict:
    """ベクトル解析: Curl / Div / Grad / Laplacian / Cross / Dot。

    Maxwell 方程式・ベクトル場の確認に使う。

    Args:
        expression: ベクトル / スカラー場の Wolfram 式。
                    例 "{0, x*y, z}", "Sin[x] * Cos[y]"
        operation: "curl" / "div" / "grad" / "laplacian"
                   / "cross" / "dot" のいずれか。
                   "cross" / "dot" のときは expression に
                   2 ベクトルをコンマ区切りで:
                   '"{0, 1, 0}", "{x, y, z}"'
        variables: 微分変数。default {x, y, z}。
                   円柱座標なら {r, theta, z}, "Cylindrical".

    Returns:
        {ok, value, code, raw}.

    Examples:
        >>> mathematica_vector_calc("{0, x*y, 0}", "curl")["value"]
        '{0, 0, y}'
        >>> mathematica_vector_calc("Sin[x]", "laplacian")["value"]
        '-Sin[x]'
    """
    op = operation.lower()
    if op == "curl":
        code = f"Curl[{expression}, {variables}]"
    elif op == "div":
        code = f"Div[{expression}, {variables}]"
    elif op == "grad":
        code = f"Grad[{expression}, {variables}]"
    elif op == "laplacian":
        code = f"Laplacian[{expression}, {variables}]"
    elif op == "cross":
        code = f"Cross[{expression}]"  # caller が "{a},{b}" を入れる
    elif op == "dot":
        code = f"Dot[{expression}]"
    else:
        return {
            "ok": False,
            "error": f"unknown operation '{operation}'. Use curl/div/grad/laplacian/cross/dot.",
            "code": "",
            "raw": None,
        }
    return _wrap(code, timeout=timeout)


# -------------------------------------------------------------------
# 5. unit_convert
# -------------------------------------------------------------------
def mathematica_unit_convert(
    quantity: str,
    target_unit: str,
    timeout: int = 60,
) -> dict:
    """物理単位変換 (Wolfram Quantity / UnitConvert)。

    Args:
        quantity: 元の値. 例 "Quantity[1, \\"Tesla\\"]"
                  または "1 Tesla" の素朴 str (内部で UnitConvert 入り口に放り込む)。
        target_unit: 変換先。例 '"Gauss"', '"Amperes/Meter"'.

    Returns:
        {ok, value, code, raw}.

    Examples:
        >>> mathematica_unit_convert('Quantity[1, "Tesla"]', '"Gauss"')["value"]
        '10000 G'
        >>> mathematica_unit_convert('Quantity[1, "Tesla"] / (4 Pi 10^-7)', '"Amperes/Meter"')["value"]
        '796178.... A/m'
    """
    code = f"UnitConvert[{quantity}, {target_unit}]"
    return _wrap(code, timeout=timeout)


# -------------------------------------------------------------------
# 6. solve
# -------------------------------------------------------------------
def mathematica_solve(
    equations: str,
    unknowns: str,
    use_reduce: bool = False,
    timeout: int = 120,
) -> dict:
    """方程式 (系) を Solve で解く。

    Args:
        equations: 等式 1 つまたはリスト。例 "x^2 == 4"
                   または "{x + y == 1, x - y == 0}".
        unknowns: 未知数。例 "x" または "{x, y}".
        use_reduce: True なら Reduce[] (より一般的、不等式・領域・cases も)。
                    Solve よりも結果が冗長になりがち。

    Returns:
        {ok, value, code, raw}.  value は Wolfram 形式の解。

    Examples:
        >>> mathematica_solve("x^2 == 4", "x")["value"]
        '{{x -> -2}, {x -> 2}}'
        >>> mathematica_solve("{x + y == 3, x - y == 1}", "{x, y}")["value"]
        '{{x -> 2, y -> 1}}'
    """
    fn = "Reduce" if use_reduce else "Solve"
    code = f"{fn}[{equations}, {unknowns}]"
    return _wrap(code, timeout=timeout)


# -------------------------------------------------------------------
# 7. integrate
# -------------------------------------------------------------------
def mathematica_integrate(
    integrand: str,
    variable: str,
    lower: str | None = None,
    upper: str | None = None,
    assumptions: str | None = None,
    timeout: int = 120,
) -> dict:
    """積分 (定積分 / 不定積分)。

    Args:
        integrand: 被積分関数。例 "Sin[x]^2"
        variable: 積分変数。例 "x"
        lower, upper: 定積分のとき両方指定。None なら不定積分。
        assumptions: Integrate の仮定。例 "k > 0".

    Returns:
        {ok, value, code, raw}.

    Examples:
        >>> mathematica_integrate("Sin[x]", "x")["value"]
        '-Cos[x]'
        >>> mathematica_integrate("Exp[-x^2]", "x", "-Infinity", "Infinity")["value"]
        'Sqrt[Pi]'
    """
    if (lower is None) ^ (upper is None):
        return {
            "ok": False,
            "error": "lower and upper must both be set or both None",
            "code": "",
            "raw": None,
        }
    if lower is None:
        body = f"Integrate[{integrand}, {variable}"
    else:
        body = f"Integrate[{integrand}, {{{variable}, {lower}, {upper}}}"
    if assumptions:
        body += f", Assumptions -> ({assumptions})"
    body += "]"
    return _wrap(body, timeout=timeout)


# -------------------------------------------------------------------
# 8. differentiate
# -------------------------------------------------------------------
def mathematica_differentiate(
    expression: str,
    variable: str,
    order: int = 1,
    timeout: int = 60,
) -> dict:
    """微分。

    Args:
        expression: 関数 (Wolfram 式)。例 "Sin[x] * Cos[k * x]"
        variable: 変数。例 "x"
        order: 階数。default 1.

    Returns:
        {ok, value, code, raw}.

    Examples:
        >>> mathematica_differentiate("Sin[x]", "x")["value"]
        'Cos[x]'
        >>> mathematica_differentiate("x^4", "x", order=2)["value"]
        '12 x^2'
    """
    if order < 1:
        return {
            "ok": False,
            "error": "order must be >= 1",
            "code": "",
            "raw": None,
        }
    if order == 1:
        code = f"D[{expression}, {variable}]"
    else:
        code = f"D[{expression}, {{{variable}, {order}}}]"
    return _wrap(code, timeout=timeout)
