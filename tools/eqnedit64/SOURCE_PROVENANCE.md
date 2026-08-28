# Source provenance

Eqnedit64 was imported into the public Radia monorepo from the Sugahara Lab
development repository at commit `87bc7d7` on 2026-08-28.

This public component contains the independent 64-bit TeX editor source,
tests, documentation, and redistributable assets. It intentionally excludes
the Microsoft Eqnedit32 executable, MTEF conversion code, reverse-engineering
working material, private signing keys, and internal deployment targets.

The original repository, including its Git history and preserved legacy
comparison assets, remains in the laboratory archive. Those legacy assets are
not build inputs and are not required to use Eqnedit64.

The Web/JavaScript edition was imported from the laboratory homepage source on
2026-08-28 (homepage build `2026-08-27a`). From this Radia change onward,
`web/equation-editor.js` and `web/equation-editor.fragment.html` are canonical;
the homepage receives a verified publication copy rather than maintaining an
independent implementation.
