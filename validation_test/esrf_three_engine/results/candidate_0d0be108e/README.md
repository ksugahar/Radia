# Radia 5.0.3 exact-wheel acceptance

Source: `0d0be108e6ddb84ddf90d14b7fce58a5e1bd4e03` (`v5.0.3`).
Tag build: [36426351986](https://github.com/ksugahar/Radia/actions/runs/36426351986).
Wheel SHA-256: `623054289d882352c6624806453f2e164983b93cb67e2af6a43a55f10f65287f`.

All four acceptance hosts passed the same wheel without rebuilding. Each
record verifies all 309 packaged Python sources and the primary native binary,
pip dependencies, six strong-field native HDiv nonlinear cases (TET/HEX/WEDGE,
orders 1 and 2, drive 300000), and the seven named publication regression cases.
The six-case driver is the existing
`../candidate_4d85e72cc/hdiv_newton_native_smoke_strong.py`.
The frozen WEDGE case uses the tracked JSON/NPZ in `../candidate_d0d0bc4b5/`.

The first isolated test input transfer omitted those frozen reference files.
After adding them, only the failed focused step and subsequent dependency
inventory were rerun; successful nonlinear results were retained. On mdx1,
the default 38-thread assembly also exhausted its localheap allocation.
The focused rerun explicitly uses four NGSolve threads, as recorded in its
command. Failed diagnostics remain in the private operational archive.
The focused suite selects the seven immutable publication-contract cases;
the newer batched indefinite-operator test is covered by numerical CI.

These are release correctness checks, not scaling measurements or a claim that
all ESRF models have completed. ESRF case 7, the case 3 direct-factorization
failure and FFAG validation remain separate work.

The same-source Simulink/MEX package was built separately and passed the full
application check on all four hosts through their existing official MATLAB
MCP sessions, including the nonlinear HDiv reactor and electromagnet topology
optimization. Package SHA-256:
`6ffc8ecb0aa20d0356e9bd9bec0c897468f4cef2fe6b54bb7d472c1d842e2a20`.
[Promotion run 36428106261](https://github.com/ksugahar/Radia/actions/runs/36428106261)
passed and published these exact bytes to PyPI. The matching wheel and
Simulink ZIP are attached to [v5.0.3](https://github.com/ksugahar/Radia/releases/tag/v5.0.3).
Artifact acceptance is complete. Deployed `release_quad done` remains a
separate gate; at publication LAB/mdx1/mdx2 were updated, while 100 was
waiting for existing MCP consumers to release the old binaries.
