# Clebsch Hodograph Docs

The public entry point is the executed, self-contained
[`public_demo.ipynb`](public_demo.ipynb). It reads its maintained companion
sources and figures directly; no notebook sidecar or migration ledger is part
of the public contract. Numerical evidence belongs under `validation_test/`
with checked JSON.

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
