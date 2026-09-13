# Notebook bibliography migration: initial audit

Base: `fbfecfa6f`, 2026-09-13. All 51 `docs/**/*.ipynb` files were parsed and
their cell sources screened for bibliography headings, citation syntax,
author/year references, arXiv/DOI links and source attributions. This is not a
claim to have checked every scientific statement for sufficient citation.
Initially there were no `.bib` or `.bbl` files below docs and no notebook
declarations using the canonical bibliography pipeline.

## Three migrated notebooks

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

## Five notebooks remain unresolved

| Notebook (below docs) | Citation work still required |
| --- | --- |
| `analytical_formulas/analytical_formulas.ipynb` | Identify and register the IEEJ Parts 1–9 series, Weissenburger/Christensen report, Montgomery book and Ortner formulation from the source cross-reference table. |
| `cohomology/tomega_wire.ipynb` | Identify Kotiuga 1987, the particular Bossavit source, and Pellikka 2013; an author surname alone does not identify a publication. |
| `mesh_fusion/mesh_fusion.ipynb` | Register the exact Becker/Hansbo/Stenberg, Hansbo/Larson, Buffa/Maday/Rapetti and Egger papers, including arXiv 2005.12020 and 2112.05572. |
| `mixed_galerkin/mixed_galerkin_results.ipynb` | Resolve legacy Senior 1962 versus parent key `Senior1962` (verified record year 1960), and Yuferev–Ida 2010 versus parent `YuferevIda2009`. Mitzner, Kameari and Quarteroni/Valli have candidate parent keys. Do not silently change publication identity to match a year. |
| `peec_integration/peec_showcase.ipynb` | Resolve Grover and FastHenry reference identities; Dowell has parent `dowell1966eddy`. Named PRIMA/source figure provenance also needs source-level review. |

The other 43 notebooks were not assigned dummy bibliographies. Screening found
internal file links, implementation dates and numerical reference comparisons;
those are not automatically scholarly citations. The absence of a detected
bibliography declaration is **unclaimed**, not an all-clear citation audit.
The Kelvin notebook's laboratory webpage link remains an ordinary source link,
not a generated scholarly reference list.

## Follow-up source reconciliation (2026-09-14)

The remaining five rows were investigated, not silently filled with plausible
author/year matches. Their blockers are now more specific:

- **Analytical formulas:** the repository cross-reference and module docstrings
  supply report identifiers such as SA-02-28/RM-02-64 and SA-04-44/RM-04-68,
  but the nine source PDFs are described as lab-internal and are not included.
  The aggregate author list is not proof of each part's author list. Obtain the
  actual report front pages before creating nine parent records. The book/report
  and Ortner references also need their exact cited editions/publications.
- **Cohomology:** [Kotiuga's own publication list](https://people.bu.edu/prk/Publications.htm)
  identifies the 1987 cuts paper (J. Appl. Phys. 61(8), 3916--3918), and the
  [SIAM record](https://epubs.siam.org/doi/10.1137/130906556) identifies Pellikka
  et al. (2013). The bare Bossavit attribution still does not identify a work.
  `src/radia/esim_multiport.py` mentions 1998, and the
  [1998 Academic Press book](https://www.sciencedirect.com/book/monograph/9780121187101/computational-electromagnetism)
  is a candidate, not proof that it was the intended source. No candidate was
  substituted into the saved notebook as a verified original citation.
- **Mesh fusion:** [arXiv 2005.12020](https://arxiv.org/abs/2005.12020) lists
  Egger, Harutyunyan, Merkel and Schöps, unlike the saved code-cell attribution
  containing Loescher and Steinbach. [arXiv 2112.05572](https://arxiv.org/abs/2112.05572)
  has five authors, including Löscher and Merkel, and corresponds to a 2022
  journal publication (the preprint is 2021). The
  [composite-grid candidate](https://doi.org/10.1051/m2an:2003039) has **Anita
  Hansbo, Peter Hansbo and Mats G. Larson**, not just the abbreviated pair in
  the notebook. Its identity does not establish the notebook's specific
  corner-stabilization claim. That claim needs a matching passage/theorem;
  do not certify it merely by adding a bibliography entry.
- **Legacy mixed_galerkin:** Senior's *A note on impedance boundary conditions*
  (Canadian J. Phys. 40, 663--665, 1962) is a distinct candidate, not the
  1960 article currently stored under `Senior1962`. The saved notebook gives
  neither title nor DOI, so the intended work must be confirmed. Yuferev/Ida's
  book also needs the cited edition/date convention resolved: the parent uses
  the 2009 publication date while the book's copyright page says 2010. Keep
  legacy calculation cells unchanged until the bibliographic intent is clear.
- **PEEC:** `validation_test/peec_integration/ngsbem_peec_demo/compute_L_final.py`
  repeats the FastHenry 21.6 nH comparison as a printed constant. The sibling
  `verify_loop_peec_vs_ngsbem.py` gives a Grover-labelled formula but no edition,
  page or equation locator. FastHenry-like Radia examples are not independent
  upstream FastHenry run evidence. Locate the actual run/model and Grover
  passage before claiming those numerical references are sourced. Adding the
  general FastHenry or PRIMA papers would not validate the saved numbers or
  figure provenance.

These findings do not change solver code or invalidate/recompute saved output;
they limit what this bibliography migration can responsibly claim.

## Reproduction and acceptance

Use the [notebook finalization procedure](../../packages/radia-mcp/docs/operations/bibliography-finalization.md#notebook-references)
to generate from `metadata.radia.bibliography.keys` through the existing
`bibliography_make_bbl` tool and render the bbl with installed TeX4ht.
No separate MCP tool, parser, notebook-local `.bib`, or JSON sidecar was added.
The fast contract checks only notebooks that explicitly declare migration;
it must not be reported as resolving the five rows above.

No shared editable installation, live client or solver binary was changed.
