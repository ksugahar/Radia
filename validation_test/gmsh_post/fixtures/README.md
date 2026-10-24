# Post-processing mesh fixture

`unit_cube_post.msh` is an analytic Netgen unit cube exported by the existing
Radia MSH 4.2 post writer. It is a mesh-only display/format fixture, never a
solver input. Its provenance records the generating source SHA, command and
exact file digest. It contains eight vertices and six linear tetrahedra.

The standing durable-lane MSH audit consumes this file. It is independently
generated data; none of the deleted conductor solver fixtures were reused or
relabeled. Mesh-based analysis continues to consume checked `.vol` only.
