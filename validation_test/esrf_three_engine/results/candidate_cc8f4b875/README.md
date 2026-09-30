# Radia 5.1.0 exact-wheel acceptance

These records gate publication of the Windows CPython 3.12 wheel built by
GitHub Actions run `36675230770` from source commit
`cc8f4b8752d34b79e8914a6bd47a0a4d98ce262b` and tag `v5.1.0`.

- Wheel SHA-256: `2f5633376a33582a8bffaa2a696e24f9da6abb3b801eb4c6bd388755e334d22b`
- Native extension SHA-256: `32486195d84799c7c41c1736533b645c6d0cb59be9ab1aa612273933c08d3732`
- Wheel Python sources verified: 313
- Targets: LAB, 100, mdx1, and mdx2

Each target used a fresh virtual environment and the same retained wheel.
`acceptance.json` records package identity and successful commands,
`full6.json` records all six tetrahedron/hexahedron/wedge order-1/order-2
strong-field nonlinear HDiv cases, and `focused.xml` records the eight
fail-loud PCG tests. The repository promotion verifier accepted every record.

An initial diagnostic run placed an unpacked source archive below the working
directory, which could shadow the installed package through Python's import
search path. That diagnostic output was discarded. The retained records are
from the clean rerun after moving the archive outside the importable tree; the
recorded import paths point to each fresh environment's `site-packages`.
