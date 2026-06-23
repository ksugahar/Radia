"""Bad Cubit script: intentionally triggers cubit lint rules.

Expected findings:
- missing-cubit-init (HIGH): imports cubit but no cubit.init()
- missing-mesh-command (CRITICAL): export without mesh command
- missing-block-registration (CRITICAL): export without block command
- geometry-block-2nd-order (HIGH): block adds geometry + element type tetra10
- wrong-connectivity-2nd-order (HIGH): get_connectivity with 2nd-order elements
- element-type-before-add (HIGH): element type set before elements added
- nodeset-sideset-usage (HIGH): nodeset/sideset with non-Exodus export
- deleted-api-usage (CRITICAL): usage of removed export_netgen API
- hardcoded-absolute-path (MODERATE): absolute Windows path
- missing-boundary-block (MODERATE): volume block without surface block + Netgen export
- wrong-file-extension (MODERATE): export with mismatched extension
- missing-block-names (LOW): blocks without names
- non-ascii-byte: NOT triggered here (would need actual non-ASCII bytes)
"""

import sys
import cubit

# hardcoded-absolute-path: absolute Windows path
sys.path.insert(0, r'D:\MyProject\custom_lib')

# Import cubit_mesh_export for export functions
import cubit_mesh_export as cme

# missing-cubit-init: no cubit.init() call anywhere

# geometry-block-2nd-order: block adds geometry then sets element type
cubit.cmd("block 1 add volume 1")
cubit.cmd("block 1 element type tetra10")

# element-type-before-add: element type set before elements added to block 2
cubit.cmd("block 2 element type tetra10")
cubit.cmd("block 2 add tet all")

# wrong-connectivity-2nd-order: get_connectivity with 2nd-order elements
nodes = cubit.get_connectivity("tet", 1)

# nodeset-sideset-usage: nodeset with non-Exodus export
cubit.cmd("nodeset 1 add node all")
cubit.cmd("sideset 1 add tri all")

# deleted-api-usage: usage of removed export_netgen API
cme.export_netgen(cubit, "output.vol")

# missing-mesh-command + missing-block-registration:
# export called but no mesh command and no proper block registration
cme.export_Gmsh_ver4(cubit, "output.vtk")
# wrong-file-extension: .vtk extension for Gmsh export

# missing-boundary-block + missing-block-names:
# volume block without surface block, no block names
cubit.cmd('export netgen "output.vol" order 2 overwrite')
