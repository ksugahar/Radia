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

from .tools import mathematica_evaluate


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
