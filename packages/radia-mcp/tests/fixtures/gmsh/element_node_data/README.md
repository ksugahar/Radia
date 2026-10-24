These are synthetic, post-processing-only MSH 4.2 fixtures. They are never solver
inputs. `generate.py` writes them with the qualified native Gmsh API, opens them
again with Gmsh, and records `getModelData` results in `provenance.json` as the
independent oracle. Fixtures cover point and line elements, scalar and vector
data, and a line with an interpolation sample count different from its mesh
node count. Tags are supplied out of order to expose value-association errors.
Only CRLF-to-LF normalization is applied; both byte hashes are recorded.

MSH 4.2 places `elementType`, `numValuesPerElement` and the element count before
the block's element tags. There are `numComponents * numValuesPerElement`
numbers per element. The sample count is not a component-inclusive width and
can differ from connectivity when an interpolation scheme is declared.

Primary references inspected in the official Gmsh 5 development source snapshot
(`gmsh-git-source.tgz`, SHA256
`4aacb2c2f9deebdf8838a4238ac860be36a5a2c7ef9a43b210dd6677adca082a`):

- `src/post/PViewDataGModelIO_MSH.cpp:128-142` reads the count and multiplies it
  by the component count; lines 343-367 write the same structure.
- `doc/texinfo/gmsh.texi:5181-5186` gives the block grammar; lines 5334-5346
  describe sample multiplicity, including interpolation schemes.

See the [official source snapshot](https://gmsh.info/src/gmsh-git-source.tgz)
and [development format specification](https://gmsh.info/dev/doc/texinfo/gmsh.html#MSH-file-format).
The source snapshot hash identifies the inspected reference, not a dependency
pin or a claim that the wheel was rebuilt from that archive. The native build
and wheel provenance actually used are recorded in `provenance.json`.
