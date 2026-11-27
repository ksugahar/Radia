"""Bad Cubit script: intentionally triggers exporter-owned Cubit lint rules.

Expected findings:
- missing-cubit-init (HIGH): imports cubit but no cubit.init()
- missing-mesh-command (CRITICAL): export without mesh command
- geometry-block-2nd-order (HIGH): block adds geometry + element type tetra10
- wrong-connectivity-2nd-order (HIGH): get_connectivity with 2nd-order elements
- element-type-before-add (HIGH): element type set before elements added
- hardcoded-absolute-path (MODERATE): absolute Windows path
- missing-boundary-block (MODERATE): volume block without surface block + Netgen export
- missing-block-names (LOW): blocks without names
- non-ascii-byte: NOT triggered here (would need actual non-ASCII bytes)
"""

import sys
import cubit

# hardcoded-absolute-path: absolute Windows path
sys.path.insert(0, r'D:\MyProject\custom_lib')

# missing-cubit-init: no cubit.init() call anywhere

# geometry-block-2nd-order: block adds geometry then sets element type
cubit.cmd("block 1 add volume 1")
cubit.cmd("block 1 element type tetra10")

# element-type-before-add: element type set before elements added to block 2
cubit.cmd("block 2 element type tetra10")
cubit.cmd("block 2 add tet all")

# wrong-connectivity-2nd-order: get_connectivity with 2nd-order elements
nodes = cubit.get_connectivity("tet", 1)

# missing-mesh-command + missing-boundary-block + missing-block-names:
# export without a mesh command, volume block without surface block, no names
cubit.cmd('export netgen "output.vol" order 2 overwrite')
