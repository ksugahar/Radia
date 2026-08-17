from __future__ import annotations

from radia_mcp.paper_writing.tools import (
    paper_writing_check_conclusion_first_use,
)


def test_accepts_terms_numbers_and_symbols_introduced_in_body() -> None:
    tex = r"""
\section{Method}
The PEEC model is solved by GMRES.  The surface impedance $Z_s$ is used.
\section{Results}
The computation time was reduced by 10\%.
\section{Conclusion}
The PEEC model with GMRES and $Z_s$ reduced the computation time by 10\%.
This result enables practical design exploration.
"""

    result = paper_writing_check_conclusion_first_use(tex)

    assert result["passed"] is True
    assert result["issue_count"] == 0
    assert result["new_technical_terms"] == []
    assert result["new_math_symbols"] == []
    assert result["new_numeric_claims"] == []


def test_flags_new_method_symbol_number_and_citation_in_conclusion() -> None:
    tex = r"""
\section{Method}
The PEEC model is used.
\section{Results}
The computation time decreased.
\section{Conclusion}
HACApK-compressed GMRES with $Z_s$ reduced the computation time by 50\%
\cite{NewSolver2026}.
"""

    result = paper_writing_check_conclusion_first_use(tex)

    assert result["passed"] is False
    assert "HACApK-compressed" in result["new_technical_terms"]
    assert "GMRES" in result["new_technical_terms"]
    assert "Z_s" in result["new_math_symbols"]
    assert "50%" in result["new_numeric_claims"]
    assert result["new_citation_keys"] == ["NewSolver2026"]


def test_accepts_new_synthesis_without_new_technical_artifact() -> None:
    tex = r"""
\section{Method}
The PEEC model is used for induction heating analysis.
\section{Results}
The model preserved the temperature distribution.
\section{まとめ}
以上の結果は，コイル形状の設計探索を実用的な時間で行えることを示唆する．
"""

    result = paper_writing_check_conclusion_first_use(tex)

    assert result["passed"] is True
    assert result["conclusion_heading"] == "まとめ"


def test_reports_missing_conclusion_heading() -> None:
    result = paper_writing_check_conclusion_first_use(
        "Method and results are described without a conclusion heading."
    )

    assert result["score"] is None
    assert result["conclusion_found"] is False


def test_whitelist_suppresses_domain_standard_acronym() -> None:
    tex = r"""
\section{Results}
The loss was reduced.
\section{Conclusion}
The result supports CLN-based design.
"""

    result = paper_writing_check_conclusion_first_use(tex, whitelist="CLN")

    assert result["passed"] is True
    assert result["new_technical_terms"] == []


def test_excludes_bibliography_commands_after_conclusion() -> None:
    tex = r"""
\section{Results}
The PEEC model reduced the loss.
\section{Conclusion}
The PEEC model reduced the loss.
\bibliographystyle{IEEJtran}
\bibliography{references}
\end{document}
"""

    result = paper_writing_check_conclusion_first_use(tex)

    assert result["passed"] is True
    assert result["new_technical_terms"] == []
