"""radia_mcp.motor — Rotating electric machine analysis.

Builds on the ONELAB/GetDP ElectricMachines reference template family
(im, im_3kW, pmsm, pmsm_cbmag, lomonova, wfsm_4p, srm, t30, shaded)
and on the Wakao group's auto-encoder + level-set topology optimization
work (Liu Xinyao's 2025 thesis, Kishi-Wakao-Murata-Makino-Takeuchi
2025 IEEJ).  Provides:

- ONELAB → NGSolve translation patterns for 2D rotating machines
- Moving-band air-gap coupling (sliding interface) recipes
- Motor-type vocabulary: PMSM, IM (cage), SRM, WFSM, shaded-pole
- Park/Clarke transformation + abc-dq0 voltage-driven simulation
- SynRM topology optimization (auto-encoder + LS method)
- Darwin-model time-domain solver (Kaimori-Mifune-Kameari-Wakao 2024)
- FEMM newbuild Lange-2009 incremental-permeability transient
  (David Meeker's implementation + sliding-band air-gap BC)
- Pointers from GetDP source-region groups → NGSolve material labels

Distilled from:
- public-safe curated corpus (Sabariego-Gyselinck-Geuzaine, ULiège-ULB)
- public-safe curated corpus (Liu Xinyao thesis, SynRM)
- Kishi-Wakao-Murata et al., IEEJ 2025 (autoencoder + LS, SynRM)
- Kaimori-Mifune-Kameari-Wakao, Darwin TD with Coulomb-type gauge
- https://www.femm.info/doku/doku.php?id=newbuild (David Meeker)
- D. Lange et al., IEEE Trans. Mag. 45(3):1258-1261, 2009
  (DOI: 10.1109/TMAG.2009.2012585)

Promoted from research/ to the radia-mcp public subpackage tree
2026-05-21.
"""
