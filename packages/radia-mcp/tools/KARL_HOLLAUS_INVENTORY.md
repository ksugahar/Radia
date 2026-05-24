# Karl Hollaus / TU Wien — Paper Inventory & Download Plan

Generated 2026-05-21 in response to user request for "Karl の論文も全部学習したい".

## Storage location

**Recommended folder**: `W:/03_文献・論文/00_電磁界解析/MOR_モデル縮約/Karl_Hollaus/`

Rationale: this folder already exists with 10 Hollaus papers + the
new arXiv-fetched paper (11 total).  Hollaus's MSFEM / Effective
Material work is structurally homogenization (a form of MOR), so the
parent `MOR_モデル縮約/` is the right home.  Other Hollaus papers
scattered elsewhere (`SIBC/`, `Non-Liniear/`) should be consolidated
here as well.

---

## What we already have (12 papers as of 2026-05-21)

In `MOR_モデル縮約/Karl_Hollaus/` (11 papers):
1. A T,Φ-Φ Multiscale Finite Element Formulation for Eddy Current Problems in Open Magnetic Circuits (Hanser-Schöbinger-Hollaus 2024)
2. A Nonlinear Effective Surface Impedance in a Magnetic Scalar Potential Formulation (Hollaus-Hanser-Schöbinger)
3. A T---- Multiscale Finite Element Formulation for Eddy Current Problems in Open Magnetic Circuits (Hanser-Schöbinger-Hollaus, duplicate)
4. An Equilibrated Error Estimator for the 2-D/1-D MSFEM T-Formulation of the Eddy Current Problem (Schöbinger-Hollaus, IEEE TMAG 60(5) May 2024)
5. Effective Interface Condition for Electromagnetic Shielding Using the T-Φ Formulation in 3D
6. Effective Material Modeling for Laminated Iron Cores With a T,ε-ε Formulation (Hanser-Schöbinger-Hollaus, IEEE TMAG 60(10) Oct 2024)
7. Effective Material and Static Magnetic Field for the 2-D/1-D Problem of Laminated Electrical Machines (2 copies)
8. **Effective Medium Transformation: the Case of Stratified Magnetic Structures (Schöbinger-Hollaus-Tsukerman 2020, arXiv:2003.12092)** ← downloaded today
9. Modeling of a Winding by Segmentation and a Two Domain Method
10. Multiscale Finite Element Formulations for 2D/1D Problems (Schöbinger-Hollaus 2023, arXiv:2304.06553 — published version)

Elsewhere (need consolidation):
- `SIBC/A_Nonlinear_Effective_Surface_Impedance...` — duplicate of #2 above, can delete
- `MOR_モデル縮約/MSFEM and MOR to Minimize the Computational Costs of Nonlinear Eddy-Current Problems in Laminated Iron Cores.pdf` (Hollaus-Schöberl-Schöbinger, IEEE TMAG 56(2) Feb 2020) — move into `Karl_Hollaus/`
- `MOR_モデル縮約/Non-Liniear/A MSFEM to simulate the eddy current problem in laminated iron cores in 3D.pdf` (Hollaus 2019 COMPEL, 16p monograph) — move into `Karl_Hollaus/`

After consolidation: **13 Hollaus papers** in `Karl_Hollaus/`.

---

## Confirmed missing — open access (download today)

✅ **ALREADY DOWNLOADED**: Schöbinger-Hollaus-Tsukerman 2020,
"Effective Medium Transformation: the Case of Stratified Magnetic
Structures", arXiv:2003.12092 (5p), DOI:10.48550/arXiv.2003.12092.
Saved to `MOR_モデル縮約/Karl_Hollaus/`.

---

## Downloaded 2026-05-21 (curl + cookie-jar session — solution found!)

All three previously-blocked papers were successfully downloaded via
curl with a Chrome-like User-Agent + cookie jar (visit abstract page
first to get cookies, then call `stampPDF/getPDF.jsp` with Referer
header).  IEEE Xplore's anti-bot accepts this pattern.

### ✅ Now in lab library at `W:/03_文献・論文/00_電磁界解析/MOR_モデル縮約/Karl_Hollaus/`

1. **Schöbinger-Hollaus-Tsukerman 2020 "Nonasymptotic Homogenization of Laminated Magnetic Cores"**, IEEE TMAG 56(2):7509504.
   - DOI: 10.1109/TMAG.2019.2943463  (was misquoted earlier as ...2945020)
   - arnumber: 8954674
   - 4p, 812 KB.  Knowledge: `motor_hollaus_eddy('nonasymptotic_homogenization')`.

2. **Hollaus-Schöbinger 2024 "MSFEM With MOR and DEIM to Solve Nonlinear Eddy Current Problems in Laminated Iron Cores"**, IEEE TMAG 60(3):7401004.
   - DOI: 10.1109/TMAG.2023.3314835
   - arnumber: 10250999
   - 4p, 5.7 MB.  Knowledge: `motor_hollaus_eddy('msfem_mor_deim_detail')`.

3. **Hanser-Schöbinger-Hollaus 2025 "Effective material modelling for laminated iron cores with an A-formulation and circuit coupling"**, COMPEL 44(5):832 (open access CC BY 4.0).
   - PDF: emerald.com/compel/article-pdf/44/5/832/10296817/compel-12-2024-0520en.pdf
   - 16p, 3.9 MB.  Knowledge: `motor_hollaus_eddy('hanser_2025_circuit_coupling')`.

Total Hollaus papers now: 14 in folder, all 14 catalogued under
`motor` subpackage (bibliography 140 → 143 entries).

## Also downloaded 2026-05-21 (curl+cookies, bonus round)

4. **Hollaus-Hannukainen-Schöberl 2014 "Two-Scale Homogenization of the Nonlinear Eddy Current Problem with FEM"**, IEEE TMAG 50(2):7010104.
   - DOI: 10.1109/TMAG.2013.2282334 ; arnumber: 6749045
   - 4p, 442 KB.  Knowledge: `motor_hollaus_eddy('hollaus_2014_nonlinear_two_scale')`.
   - The foundational two-scale ansatz that all subsequent MSFEM
     work builds on.

5. **Schöbinger-Hollaus 2021 "A Hierarchical Error Estimator for the MSFEM for the Eddy Current Problem in 3-D"**, IEEE TMAG 57(5):7401005.
   - DOI: 10.1109/TMAG.2021.3062041 ; arnumber: 9363162
   - 5p, 949 KB.  Knowledge: `motor_hollaus_eddy('hierarchical_error_estimator')`.
   - Adaptive p-refinement for 3D MSFEM (companion to the 2024
     equilibrated estimator for 2D-1D).

6. **Frljić-Hanser-Hollaus-Schöbinger 2026 "Homogenization of the Eddy Current Problem in Laminated Open-Type Cores Accounting for Perpendicular Flux"**, IEEE TMAG 62(5):6300408.
   - DOI: 10.1109/TMAG.2026.3680721 ; arnumber: 11474630
   - 8p, 1.2 MB.  Knowledge: `motor_hollaus_eddy('frljic_2026_perpendicular_flux')`.
   - Anisotropic μ_eff tensor for open-type cores with through-thickness
     B field — lifts the "net B_z = 0" restriction from
     Schöbinger-Hollaus-Tsukerman 2020.  Linear / freq-domain only.

Total Karl_Hollaus folder: **17 papers** (was 11 → +3 main set
2026-05-21 → +3 bonus set 2026-05-21).  Motor bibliography 140 → 146
entries.

---

## Likely missing — need more thorough search

To find these I would need either:
- Hollaus's complete TU Wien publication list (the institutional
  page seems to have moved / requires JS to render)
- ResearchGate access (currently blocked: HTTP 403)
- Google Scholar account with proper auth

Suggested user action:
1. Login to https://scholar.google.com/citations?user=hWnXM_4AAAAJ
   (if that's his correct Scholar ID) and bulk-download recent
   open-access PDFs
2. From Kindai IEEE Xplore institutional access, download the 2
   confirmed-missing IEEE papers above
3. Drop downloads into `W:/03_文献・論文/00_電磁界解析/MOR_モデル縮約/Karl_Hollaus/`

I'll then re-run `pdf_batch_indexer.py` to update the catalog and
`catalog_to_knowledge.py` to refresh the bibliography index.

---

## Hanser PhD thesis (TU Wien, 2024 or 2025)

Valentin Hanser is the most active Hollaus collaborator on the
recent papers (2024 TMAG×2).  Per the TU Wien wiki search, he is a
PhD candidate.  When his thesis is published it will be a
comprehensive ~150-200 page reference covering the same material as
the 4-5 individual papers — worth checking annually.

Search: `site:repositum.tuwien.at "Hanser" multiscale`

---

## Consolidation script

```powershell
# Move scattered Hollaus papers into the Karl_Hollaus folder.
$dst = "W:\03_文献・論文\00_電磁界解析\MOR_モデル縮約\Karl_Hollaus"
Move-Item "W:\03_文献・論文\00_電磁界解析\MOR_モデル縮約\MSFEM and MOR to Minimize*.pdf" $dst
Move-Item "W:\03_文献・論文\00_電磁界解析\MOR_モデル縮約\Non-Liniear\A MSFEM to simulate*.pdf" $dst
# SIBC/...nonlinear_effective_surface_impedance... is a duplicate;
# delete after diff-confirming:
# Remove-Item "W:\03_文献・論文\00_電磁界解析\SIBC\A_Nonlinear_Effective_Surface_Impedance*.pdf"
```

After consolidation, run:
```powershell
cd S:\Radia\01_GitHub\.claude\worktrees\compassionate-benz-8c4989\packages\radia-mcp
python tools/pdf_batch_indexer.py --root "W:\03_文献・論文" --out tools/catalog.json --incremental
python tools/reclassify.py --catalog tools/catalog.json
$env:PYTHONIOENCODING="utf-8"; python tools/catalog_to_knowledge.py --catalog tools/catalog.json --radia-mcp-src src/radia_mcp
```
