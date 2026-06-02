"""poster sub-module — 学会ポスター発表の作文・組版支援.

Promoted from the COMSOL Conference 2025 Tokyo poster (Kelvin.tex,
2025-12-05) as a reusable sub-skill. presentation/ がスライド用なのに
対し、こちらは A1/A0 ポスター特化:

    - 遠目 (1m) でも読めるフォントサイズ basement
    - 図優位 (text < 図 + caption 面積)
    - tcolorbox/multicols ベースの bxjsarticle テンプレ
    - platex + dvipdfmx 経由のコンパイル (texcompile 連携)

Tools (prefix ``poster_``):
    - poster_template_kelvin       — A1 縦・bxjsarticle ポスター雛形を返す
    - poster_lint                  — ポスター .tex の遠目可読性 lint
    - poster_compile               — platex + dvipdfmx で PDF 生成
    - poster_figures_audit         — figures/ の参照整合チェック
"""
from . import tools as _tools


def register(mcp) -> int:
    """Register poster_* tools with the given FastMCP instance.

    Returns the count of tools registered.
    """
    count = 0
    for name in dir(_tools):
        if name.startswith("poster_") and callable(getattr(_tools, name)):
            mcp.tool()(getattr(_tools, name))
            count += 1
    return count
