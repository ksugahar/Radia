# cubit-mesh-export

High-order curved mesh export from [Coreform Cubit](https://coreform.com/products/coreform-cubit/) to [NGSolve](https://ngsolve.org)/[Netgen](https://github.com/NGSolve/netgen).

## Features

- **Arbitrary-order curving** (order 1-5) via ACIS geometry projection
- **Label preservation**: material (block), boundary (sideset), edge (BBND)
- **Companion JSON**: CAD reference values for Volume/Area/Length consistency checking
- **Standalone checker**: verify mesh quality without Cubit (`check-vol` CLI)

## Install

```bash
pip install cubit-mesh-export
```

## Usage

### Export (requires Cubit)

```python
import netgen       # must import before cubit
import ngsolve
import cubit

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd('open "model.cub5"')

from cubit_mesh_export import extract_curved_mesh
ng_mesh = extract_curved_mesh(cubit, order=3)
ng_mesh.Save("model.vol")
```

### Check (does NOT require Cubit)

```bash
check-vol model.vol
```

```python
from cubit_mesh_export.check import check_consistency
results = check_consistency("model.vol")
```

## Part of the Radia project

Source: [github.com/ksugahar/Radia](https://github.com/ksugahar/Radia)
