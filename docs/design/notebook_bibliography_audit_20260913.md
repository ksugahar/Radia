# Notebook bibliography migration: initial audit

Base: `fbfecfa6f`, 2026-09-13. All 51 `docs/**/*.ipynb` files were parsed and
their cell sources screened for bibliography headings, citation syntax,
author/year references, arXiv/DOI links and source attributions. This is not a
claim to have checked every scientific statement for sufficient citation.
Initially there were no `.bib` or `.bbl` files below docs and no notebook
declarations using the canonical bibliography pipeline.

## Eight migrated notebooks

`docs/open_boundary/open_boundary_demo.ipynb` now declares seven explicit
canonical keys and ships its generated `.bbl`. The reference cell is the actual
TeX4ht-rendered bbl body, not retyped BibTeX fields. Existing calculation source,
execution counts and outputs are unchanged. The bibliography HTML is saved and
needs no live TeX installation to read. The bibliography-only generation was run
with BibTeX and TeX4ht from TeX Live 2026; the numerical notebook was not rerun.

The seven works cover the existing Grote/Keller, Hagstrom/Warburton,
Warburg/Cauer, Kameari CLN and Freeman/Lowther attributions. Two missing parent
entries were verified and added using the publisher's metadata:

- [Grote and Keller, SIAM J. Appl. Math. 55(2), 280–297 (1995)](https://epubs.siam.org/doi/10.1137/S0036139993269266).
- [Hagstrom and Warburton, SIAM J. Numer. Anal. 47(5), 3678–3704 (2009)](https://epubs.siam.org/doi/10.1137/090745477).

This establishes the citation pipeline, not numerical equivalence to those
papers. Existing parent records were reused without inventing missing metadata.

On 2026-09-14, `docs/section_optics/section_optics_design.ipynb` was migrated
using the same pipeline and the explicit key `steinberg2024beamline`. Its
handwritten full citation was replaced by a link to the generated reference
cell; the open-access preprint link and the distinction between the paper's
Halbach arrays and this notebook's iron electromagnet remain intact. The
[APS publisher record](https://journals.aps.org/prab/abstract/10.1103/PhysRevAccelBeams.27.071601)
verifies the authors, title, volume 27, article 071601 and publication year
2024; [arXiv 2402.01120](https://arxiv.org/abs/2402.01120) verifies the preprint.
Calculation cells, saved outputs and widget state are unchanged; no numerical
revalidation is claimed. The parent addition also leaves the first notebook's
selected-source fingerprint valid.

`docs/nodal_force/nodal_force.ipynb` was migrated on 2026-09-14 with three
explicit keys and the same generated-bbl/display contract. Two missing records
were verified against [IEEE-deposited Crossref metadata](https://api.crossref.org/works/10.1109/TMAG.1983.1062812)
and the [Elsevier article](https://www.sciencedirect.com/science/article/pii/S0377042703009816).
The [Emerald eggshell record](https://doi.org/10.1108/03321640410553427) identifies
the COMPEL article as 2004, not 2003, and spells the coauthor Deliége. The
notebook now explicitly cites that journal article; no separate 2003 conference
publication is claimed. Its calculation cells and saved outputs are unchanged.

## Final five migrations (2026-09-14)

| Notebook (below docs) | Explicit canonical bibliography |
| --- | --- |
| `analytical_formulas/analytical_formulas.ipynb` | 12 records: nine IEE Japan reports, PPPL-1517, Montgomery (1969), Ortner et al. (2023 volume; online 2022). |
| `cohomology/tomega_wire.ipynb` | Kotiuga (1987), explicitly selected Bossavit book (1998), Pellikka et al. (2013). |
| `mesh_fusion/mesh_fusion.ipynb` | Five records: Becker et al., Hansbo et al., Buffa et al., and Egger preprints 2005.12020 / 2112.05572. |
| `mixed_galerkin/mixed_galerkin_results.ipynb` | Six records, explicitly distinguishing Senior's 1960 article and 1962 note, plus Mitzner, Yuferev/Ida, Kameari and Quarteroni/Valli. |
| `peec_integration/peec_showcase.ipynb` | Dowell, FASTHENRY, explicitly selected Grover 2004 reprint and PRIMA. Numeric/figure provenance is separately qualified below and in the notebook. |

The other 43 notebooks were not assigned dummy bibliographies. Screening found
internal file links, implementation dates and numerical reference comparisons;
those are not automatically scholarly citations. The absence of a detected
bibliography declaration is **unclaimed**, not an all-clear citation audit.
The Kelvin notebook's laboratory webpage link remains an ordinary source link,
not a generated scholarly reference list.

## Source reconciliation and remaining evidence limits

All eight screened notebooks now use the generated-bbl/display contract. This
closes their bibliography-pipeline migration, **not** every historical attribution
or numerical validation question. Ambiguous old labels were not silently promoted
to verified citations: newly selected background works and remaining limitations
are explicit in the notebook narrative, which supersedes retained script labels.

- **Analytical formulas:** the previously unavailable nine originals were located
  in `W:/03_文献・論文/00_電磁界解析/02_解析積分公式集/00_electromagnetic_integral_formula_series/`.
  Their front-page author lists and SA/RM identifiers were read and visually
  verified. Part 4 uses its own front-page author order, not a later part's
  reordered citation; Part 8 has Matsuo alone. Years follow report identifiers,
  not received dates (Part 2 received in 2002 belongs to 2003; Part 9 received
  in 2006 belongs to 2007). The nine source PDFs remain in the literature store,
  not copied into Git. Part 1 reference 15 and the
  [PPPL report](https://www.osti.gov/servlets/purl/6159556) identify Weissenburger
  and Christensen (1979), PPPL-1517. The
  [digitized Montgomery book record](https://books.google.com/books?id=D1138LYUoXAC)
  identifies Wiley-Interscience (1969). The
  [Ortner publisher record](https://www.mdpi.com/2673-8724/3/1/2) distinguishes
  online publication on 2022-12-30 from Magnetism 3(1), 11--31 (2023).
- **Cohomology:** [Kotiuga's own publication list](https://people.bu.edu/prk/Publications.htm)
  identifies the 1987 cuts paper (J. Appl. Phys. 61(8), 3916--3918), and the
  [SIAM record](https://epubs.siam.org/doi/10.1137/130906556) identifies Pellikka
  et al. (2013). The revision explicitly selects the
  [1998 Academic Press book](https://www.sciencedirect.com/book/monograph/9780121187101/computational-electromagnetism)
  as background, not proof of the original bare surname's intent. The lab copy
  under `01_教科書/20_数学/01_微分幾何・外微分形式/03_Bossavit/01_Computational_Electromagnetism/`
  includes a corrected author version (`00Frontpages.pdf`, correction note dated
  2003-11-20). `2.pdf`, printed p. 58, solution 2.6 discusses cuts and relative
  homology; that passage was read and visually verified. The reference still
  identifies the published 1998 book, not an invented 2003 edition.
- **Mesh fusion:** [arXiv 2005.12020](https://arxiv.org/abs/2005.12020) lists
  Egger, Harutyunyan, Merkel and Schöps, unlike the saved code-cell attribution
  containing Loescher and Steinbach. [arXiv 2112.05572](https://arxiv.org/abs/2112.05572)
  has five authors, including Löscher and Merkel, and corresponds to a 2022
  journal publication (the preprint is 2021). The
  [composite-grid candidate](https://doi.org/10.1051/m2an:2003039) has **Anita
  Hansbo, Peter Hansbo and Mats G. Larson**, not just the abbreviated pair in
  the notebook. The exact
  [Becker publisher paper](https://www.esaim-m2an.org/articles/m2an/pdf/2003/02/m2an0173.pdf)
  and [Buffa archive record](https://www.numdam.org/item/M2AN_2001__35_2_191_0/)
  also identify the other works. Bibliographic identity does not establish the notebook's specific
  corner-stabilization claim. That claim needs a matching passage/theorem;
  do not certify it merely by adding a bibliography entry.
- **Legacy mixed_galerkin:** Senior's *A note on impedance boundary conditions*
  ([publisher-deposited metadata](https://api.crossref.org/works/10.1139/p62-067),
  Canadian J. Phys. 40, 663--665, 1962) is distinct from the 1960 article under
  `Senior1962`. Both are now explicitly cited as background, not as proof of
  which paper supplied the notebook's "Senior tower" coefficients. Existing
  keys are not renamed. Yuferev/Ida uses the parent publication year 2009; the
  shared `34_SIBC/00_Comprehensive/Surface Impedance Boundary Conditions a Comprehensive Approach.pdf`
  copyright page (PDF page 4, ISBN 978-1-4200-4489-8) was visually checked and
  says 2010. Both date conventions are explained rather than inventing a second
  work. The accuracy table and coefficient derivations were not revalidated.
- **PEEC:** `validation_test/peec_integration/ngsbem_peec_demo/compute_L_final.py`
  repeats the FastHenry 21.6 nH comparison as a printed constant. The sibling
  `verify_loop_peec_vs_ngsbem.py` gives a Grover-labelled formula but no edition,
  page or equation locator. FastHenry-like Radia examples are not independent
  upstream FastHenry run evidence. Locate the actual run/model and Grover
  passage before claiming those numerical references are sourced. The
  [digitized Grover record](https://books.google.com/books?id=K3KHi9lIltsC)
  identifies the 2004 reprint of the 1946 book; that edition is explicitly
  selected only as background. General FASTHENRY/PRIMA papers do not validate
  saved constants or figure provenance. The notebook now says this prominently
  and no longer calls its opening summary an independently verified showcase.
  The four script-generated figure filenames are not evidence of an external
  publication or measured-data origin.

These findings do not change solver code or invalidate/recompute saved output;
they limit what this bibliography migration can responsibly claim.

## Reproduction and acceptance

Use the [notebook finalization procedure](../../packages/radia-mcp/docs/operations/bibliography-finalization.md#notebook-references)
to generate from `metadata.radia.bibliography.keys` through the existing
`bibliography_make_bbl` tool and render the bbl with installed TeX4ht.
No separate MCP tool, parser, notebook-local `.bib`, or JSON sidecar was added.
The fast contract checks only notebooks that explicitly declare migration;
it must not be reported as resolving the scientific evidence limitations above.
All original code cells (including outputs, execution counts and cell metadata)
and pre-existing notebook metadata/widget state are preserved. No calculation
was executed to manufacture missing provenance.

Final-five acceptance: 113 focused tests passed (docs notebook contract plus
bibliography canonical/finalization/citation-source contracts). All 27 original
code cells compare equal to the base, including saved outputs and execution
metadata; all original notebook metadata, including widget state, also compares
equal after removing only the new bibliography declaration. Five independent
bbl generations matched their rendering copies byte-for-byte after LF
normalization, and all five TeX4ht renders succeeded. The earlier three notebooks'
selected-source fingerprints remain valid after the 21 parent additions.

No shared editable installation, live client or solver binary was changed.
