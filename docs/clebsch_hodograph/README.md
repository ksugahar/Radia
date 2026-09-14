# Clebsch Hodograph Docs

Explore how coordinate-based magnet design changes bending and focusing,
and which field assumptions survive numerical checks. Start with the saved
results, equations and references in these Python notebooks:

- [Two-dimensional feasibility](hodograph_feasibility_2d.ipynb)
- [Bending and transverse response](hodograph_bending_sy.ipynb)
- [Excitation-invariant field construction](excitation_invariant_field.ipynb)
- [Edge focusing and particle tracking](edge_focusing_tracking.ipynb)

The Radia MCP electromagnet and NGSolve tool families own current operating
guidance; Python is driven through MCP by an LLM. These research demonstrations
do not certify every geometry or replace the formal Simulink UI.
The former source catalogs are [internal maintenance records](../../validation_test/documentation_maintenance/README.md),
not executable capability showcases.

Use the existing notes as the theory spine:

| Topic | Read |
| --- | --- |
| Coordinate-transform backbone | [HODOGRAPH_BACKBONE.md](HODOGRAPH_BACKBONE.md) |
| Magnet-design methodology | [DESIGN_METHODOLOGY.md](DESIGN_METHODOLOGY.md) |
| HDiv/VIM/Clebsch bridge | [HDIV_VIM_CLEBSCH_BRIDGE.md](HDIV_VIM_CLEBSCH_BRIDGE.md) |
| Symbolic differential-geometry WLS index | [DIFFERENTIAL_GEOMETRY_WLS.md](DIFFERENTIAL_GEOMETRY_WLS.md) |

The runnable companion set is
[`docs/clebsch_hodograph/demos/`](demos/). It is part
of the documented result chain because
`validation_test/feec/test_clebsch_hodograph_research.py` golden-locks those
research examples directly.
