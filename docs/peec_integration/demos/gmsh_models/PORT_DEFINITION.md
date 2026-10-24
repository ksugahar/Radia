# PEEC ports on checked Netgen meshes

Solver port labels belong to the `.vol` label contract. Declare
`port_positive` and `port_negative` as boundary names in the mesher and include
them in the versioned contract's `required.boundaries` list. Gmsh physical
names are display metadata and do not define solver ports.

Run `../basic_peec/demo_port_extraction.py --mesh conductor.vol --contract
conductor.labels.json` to inspect the native Netgen boundary faces and derive
port edges and centroids. The other PEEC demos retain their documented
coordinate-based selection; the supplied mesh must represent that geometry.

SI units apply to `.vol` coordinates. Material parameters and PEEC
cross-sections remain explicit in each calculation.
