"""Put EQNEDT64's boxes next to TeX's, on the same equations at the same size.

TeX will state the width, height and depth of anything it sets, so the
comparison is against numbers rather than against a picture of numbers -- and
it needs no window, which the screen-measuring did.

THE REFERENCE MUST BE SET IN THE SAME FONT, and getting that wrong cost two
wrong conclusions.  tex_reference.tex is built with XeLaTeX and unicode-math
against Latin Modern Math -- the very .otf this editor reads.  Written with
lmodern instead, TeX quietly uses lmex10, a Type1 font whose radical ladder is
finer than the OpenType one; a root over a fraction then read as 18 percent out
when it was in fact exact, and a change made to the integral on the strength of
that comparison turned out to be wrong in the other direction.

Rebuild the reference:

    xelatex -jobname=dims validation_test/equation/tex_reference.tex

and point EQ_TEX_DIMS at the dims.txt it writes.
"""
import os
import re
import sys

sys.path.insert(0, r"S:\Radia\01_GitHub\.claude\worktrees\radia-equation\src")
import radia.equation as eq   # noqa: E402

TEX = r"C:\temp\ee3\tex\dims.txt"

CASES = {
    "x": "x", "ab": "ab", "abc": "abc",
    "frac": r"\frac{a}{b}", "frac_wide": r"\frac{abc}{d}",
    "frac_nested": r"\frac{\frac{p}{q}}{c}",
    "frac_nested2": r"\frac{\frac{\frac{1}{2}}{3}}{4}",
    "sqrt2": r"\sqrt{2}", "sqrt_frac": r"\sqrt{\frac{a}{b}}",
    "sqrt_deep": r"\sqrt{\frac{\frac{a}{b}}{c}}", "root3": r"\sqrt[3]{x}",
    "sup": "x^{2}", "sub": "x_{i}", "subsup": "x_{i}^{2}",
    "sup_sup": "x^{y^{z}}",
    "int": r"\int_{0}^{T}", "int_f": r"\int_{0}^{T}f",
    "sum": r"\sum_{n=1}^{N}", "oint": r"\oint_{C}",
    "plus": "a+b", "eq": "a=b", "comma": "a,b", "paren": "a(b)",
    "times": "a*b",
}

PT = re.compile(r"(-?[0-9.]+)pt")


def tex_dims():
    out = {}
    with open(TEX, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            name = line.split()[0]
            nums = [float(v) for v in PT.findall(line)]
            if len(nums) == 3:
                out[name] = tuple(nums)          # wd, ht, dp
    return out


def ours(latex):
    """Width, height above the baseline and depth below it, in points --
    the same three numbers TeX reports, asked for directly."""
    st = eq.SvgStyle()
    st.padding = 0.0
    return eq.tex_metrics(latex, st)


def main():
    tex = tex_dims()
    print("%-13s %-23s %-23s %s" % ("", "TeX (wd ht dp)", "EQNEDT64", "ht+dp ratio"))
    rows = []
    for name, latex in CASES.items():
        if name not in tex:
            continue
        tw, th, td = tex[name]
        ow, oh, od = ours(latex)
        rows.append((name, tw, th, td, ow, oh, od))
        r = (oh + od) / (th + td) if (th + td) else 0
        print("%-13s %6.2f %6.2f %6.2f   %6.2f %6.2f %6.2f   %.3f"
              % (name, tw, th, td, ow, oh, od, r))

    print()
    print("=== shape, which is the part that is ours to get right ===")
    d = {r[0]: r for r in rows}

    def tot(rec, side):
        return rec[2] + rec[3] if side == "tex" else rec[5] + rec[6]

    def rel(a, b, side):
        if a not in d or b not in d:
            return None
        return tot(d[a], side) / tot(d[b], side)

    checks = [
        ("a fraction against a letter", "frac", "x"),
        ("a fraction inside one", "frac_nested", "frac"),
        ("two deep", "frac_nested2", "frac"),
        ("a root against a letter", "sqrt2", "x"),
        ("a root over a fraction", "sqrt_frac", "frac"),
        ("a superscript against a letter", "sup", "x"),
        ("stacked superscripts", "sup_sup", "sup"),
        ("an integral with limits", "int", "x"),
        ("a summation with limits", "sum", "x"),
    ]
    print("%-32s %8s %8s %s" % ("", "TeX", "ours", "difference"))
    for label, a, b in checks:
        t, o = rel(a, b, "tex"), rel(a, b, "ours")
        if t is None or o is None:
            continue
        print("%-32s %8.3f %8.3f %+7.1f %%"
              % (label, t, o, (o / t - 1) * 100))

    print()
    print("=== where a fraction sits about the baseline ===")
    for name in ("frac", "frac_wide", "frac_nested"):
        if name not in d:
            continue
        _n, tw, th, td, ow, oh, od = d[name]
        print("  %-13s TeX %5.1f %% above    ours %5.1f %% above"
              % (name, 100 * th / (th + td), 100 * oh / (oh + od)))

    print()
    print("=== spacing, as extra width over the bare letters ===")
    if "ab" in d:
        base = d["ab"][1], d["ab"][4]
        for name, sym in (("plus", "+"), ("eq", "="), ("comma", ","),
                          ("times", "*"), ("paren", "(")):
            if name not in d:
                continue
            _n, tw, _th, _td, ow, _oh, _od = d[name]
            print("  %-3s  TeX +%5.2f pt   ours +%5.2f pt"
                  % (sym, tw - base[0], ow - base[1]))


if __name__ == "__main__":
    main()
