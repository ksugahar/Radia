# Notebook bibliography migration: initial audit

Base: `fbfecfa6f`, 2026-09-13. All 51 `docs/**/*.ipynb` files were parsed and
their cell sources screened for bibliography headings, citation syntax,
author/year references, arXiv/DOI links and source attributions. This is not a
claim to have checked every scientific statement for sufficient citation.
Initially there were no `.bib` or `.bbl` files below docs and no notebook
declarations using the canonical bibliography pipeline.

## Two migrated notebooks

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

## Six notebooks remain unresolved

| Notebook (below docs) | Citation work still required |
| --- | --- |
| `analytical_formulas/analytical_formulas.ipynb` | Identify and register the IEEJ Parts 1–9 series, Weissenburger/Christensen report, Montgomery book and Ortner formulation from the source cross-reference table. |
| `cohomology/tomega_wire.ipynb` | Identify Kotiuga 1987, the particular Bossavit source, and Pellikka 2013; an author surname alone does not identify a publication. |
| `mesh_fusion/mesh_fusion.ipynb` | Register the exact Becker/Hansbo/Stenberg, Hansbo/Larson, Buffa/Maday/Rapetti and Egger papers, including arXiv 2005.12020 and 2112.05572. |
| `mixed_galerkin/mixed_galerkin_results.ipynb` | Resolve legacy Senior 1962 versus parent key `Senior1962` (verified record year 1960), and Yuferev–Ida 2010 versus parent `YuferevIda2009`. Mitzner, Kameari and Quarteroni/Valli have candidate parent keys. Do not silently change publication identity to match a year. |
| `nodal_force/nodal_force.ipynb` | Verify Coulomb 1983, Henrotte/Hameyer 2004, and the eggshell attribution written as 2003 (parent `henrotte2004` is a 2004 COMPEL paper). |
| `peec_integration/peec_showcase.ipynb` | Resolve Grover and FastHenry reference identities; Dowell has parent `dowell1966eddy`. Named PRIMA/source figure provenance also needs source-level review. |

The other 43 notebooks were not assigned dummy bibliographies. Screening found
internal file links, implementation dates and numerical reference comparisons;
those are not automatically scholarly citations. The absence of a detected
bibliography declaration is **unclaimed**, not an all-clear citation audit.
The Kelvin notebook's laboratory webpage link remains an ordinary source link,
not a generated scholarly reference list.

## Reproduction and acceptance

Use the [notebook finalization procedure](../../packages/radia-mcp/docs/operations/bibliography-finalization.md#notebook-references)
to generate from `metadata.radia.bibliography.keys` through the existing
`bibliography_make_bbl` tool and render the bbl with installed TeX4ht.
No separate MCP tool, parser, notebook-local `.bib`, or JSON sidecar was added.
The fast contract checks only notebooks that explicitly declare migration;
it must not be reported as resolving the six rows above.

No shared editable installation, live client or solver binary was changed.
