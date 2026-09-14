"""Optional real-TeX acceptance: generated bbl only, no manuscript-local bib.

Run explicitly with --basetemp under C:/temp; retain PDF/JSON there for visual QA.
The data is synthetic and is not a source of real citation metadata.
"""
import hashlib
import json
import shutil
import subprocess

import pytest


@pytest.mark.skipif(any(shutil.which(tool) is None for tool in ("bibtex", "kpsewhich", "pdflatex")),
                    reason="real bibliography acceptance requires a TeX installation")
@pytest.mark.parametrize("style", ["plain", "IEEEtran"])
def test_generated_bbl_renders_accents_literals_and_corporate_authors(tmp_path, monkeypatch, style):
    from radia_mcp.bibliography.plans import T14_canonical as canonical
    pypdf = pytest.importorskip("pypdf")
    parent = tmp_path / "canonical-fixture.bib"
    parent.write_text(r'''
@article{personal,
 author={M{\"u}ller, Jane and {van der Vorst}, Henry},
 title={Magnetic fields and 50\% efficiency},
 journal={Synthetic Journal of Test Fixtures}, year={2024},
 volume={12}, number={3}, pages={12--19}}
@misc{corporate,
 author={{Research and Development Group}},
 title={Validation\_2026 \& source integrity}, year={2025}}
@book{nested,
 author={Doe, John}, title={{Nested {Magnetic} Models}},
 publisher={Synthetic Test Press}, year={2023}}
''', encoding="utf-8")
    monkeypatch.setattr(canonical, "CANONICAL", parent)
    manuscript = tmp_path / "manuscript"
    manuscript.mkdir()
    tex = manuscript / "paper.tex"
    tex.write_text(r"""\documentclass[10pt,twocolumn]{article}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage[margin=2cm]{geometry}
\begin{document}
\section*{Synthetic bibliography acceptance}
These are synthetic fixtures, not real publications.
Personal names and literal symbols~\cite{personal}, organization
names~\cite{corporate}, and protected title groups~\cite{nested}
exercise the generated bibliography contract.
\bibliographystyle{""" + style + r"""}
\bibliography{references}
\end{document}
""", encoding="utf-8")
    generation = canonical.bibliography_make_bbl(str(tex))
    assert generation.startswith("bibliography_make_bbl:"), generation
    bbl = tex.with_suffix(".bbl")
    original_bbl = bbl.read_bytes()
    for _ in range(2):
        result = subprocess.run(
            [shutil.which("pdflatex"), "-no-shell-escape", "-interaction=nonstopmode",
             "-halt-on-error", tex.name], cwd=manuscript, capture_output=True, timeout=60,
        )
        assert result.returncode == 0, result.stdout.decode("utf-8", errors="replace")[-3000:]
    assert not list(manuscript.glob("*.bib"))
    assert bbl.read_bytes() == original_bbl
    log = tex.with_suffix(".log").read_text(encoding="utf-8", errors="replace")
    for failure in ("undefined", "Overfull", "Missing character"):
        assert failure not in log, log[-3000:]
    pdf = tex.with_suffix(".pdf")
    reader = pypdf.PdfReader(pdf)
    assert len(reader.pages) == 1
    text = " ".join(reader.pages[0].extract_text().split())
    for phrase in ("Müller", "van der Vorst", "Research and Development Group",
                   "50%", "Validation_2026 & source integrity", "Nested Magnetic Models"):
        assert phrase in text, (phrase, text)
    assert "\ufffd" not in text
    (manuscript / "acceptance.json").write_text(json.dumps({
        "style": style, "generation": generation, "page_count": 1,
        "pdf_sha256": hashlib.sha256(pdf.read_bytes()).hexdigest(),
        "bbl_sha256": hashlib.sha256(original_bbl).hexdigest(),
        "text_checks": "passed", "visual_review": "pending",
        "fixture": "synthetic; no canonical data changed",
    }, indent=2), encoding="utf-8")
