"""Graph Neural Networks for PDE / EM problems."""


# Authoritative topic enum for the dispatcher tool (wired into
# `<short>_topics()` via common.register_topics_tool).
TOPICS: dict[str, str] = {
    "overview": "(see knowledge file for details)",
    "physics_embedded": "(see knowledge file for details)",
    "equivariant_gnn": "(see knowledge file for details)",
    "all": "return \"\n\n\".join([OVERVIEW, PHYSICS_EMBEDDED, EQUIVARIANT_GNN])",
}

OVERVIEW = r"""
# Graph Neural Networks (GNN) Overview

A neural network architecture where each layer aggregates information
from neighbors in a graph. Generalizes CNN (regular grid) to
arbitrary connectivity.

## General message-passing formulation

For node `i` with feature `h_i` at layer `l`:

    h_i^{l+1} = phi(h_i^l, AGG_{j in N(i)} psi(h_i^l, h_j^l, e_ij))

where:
- N(i) = neighbors of i
- AGG = permutation-invariant aggregation (sum, mean, max)
- phi, psi = small MLPs
- e_ij = edge feature (length, angle, material contrast, ...)

After L layers, each node has aggregated info from its L-hop
neighborhood.

## Key variants

| Variant | Innovation | Use |
|---------|-----------|-----|
| **GCN** (Kipf-Welling 2017) | Normalized adjacency convolution | Citation networks, simple use |
| **GAT** (Velickovic 2018) | Attention per edge | Variable edge importance |
| **MPNN** (Gilmer 2017) | General message-passing | Chemistry, default for new problems |
| **GIN** (Xu 2019) | Sum aggregation = WL test | Maximum expressive power |
| **SchNet** (Schütt 2017) | Continuous filters in 3D space | Molecular dynamics, atoms |
| **EquivariantGNN** (E(n)-GNN) | SO(3) / E(3) equivariance | Physical systems |
| **MeshGraphNet** (Pfaff 2021) | Edges = mesh edges; rollout in time | PDE simulation surrogate |
| **GraphPDE** (e.g. Brandstetter 2022) | Encode PDE state as graph | Replace FEM solve with GNN |

## Why GNN for EM

Three natural connections:

1. **FEM mesh IS a graph**: nodes = vertices, edges = element edges,
   node features = field values. A GNN trained on FEM solutions
   becomes a fast surrogate.
2. **Topology design = graph**: choosing which voxels are iron vs
   air = node label prediction.
3. **Circuit netlists = graph**: PEEC, equivalent circuits, lumped
   models are all graph data.

## Position in radia-mcp universe

| Need | Use |
|------|-----|
| Mesh-based PDE surrogate | GNN (this MCP) |
| Continuous-coordinate PDE surrogate | PINN (`radia_mcp.pinn`) |
| Topology synthesis | GNN + topology_optimization |
| Circuit reduction | PRIMA / Cauer (`radia_mcp.mor`) |
| Graph-structured Bayesian inference | NumPyro + custom GNN |
"""


PHYSICS_EMBEDDED = r"""
# Physics-Embedded GNN PDE Solvers

`Physics-Embedded Neural Networks: Graph Neural PDE Solvers with
Mixed Boundary Conditions.pdf` (public-safe curated corpus)

## Idea

Hard-code conservation laws + boundary conditions into the GNN
architecture (not just into the loss like PINN does).

Two embedding strategies:

1. **Hard PDE residual via projection**: each GNN layer's output is
   projected onto the manifold of PDE-satisfying functions
2. **Equivariance**: build the GNN to respect physical symmetries
   (rotation, reflection, translation) by construction

## Mixed boundary conditions

Real EM problems mix:
- Dirichlet (fixed potential at electrodes)
- Neumann (zero normal flux at insulators)
- Robin (SIBC, impedance)
- Periodic (Kelvin, symmetry)
- Far-field decay (open boundary)

Standard PINN penalizes BC violations via loss terms; GNN can:
- Mask Dirichlet nodes (set them to known values, don't learn)
- Add Neumann/Robin as extra edge features
- Periodic via graph identification (twin nodes)
- Far-field via Kelvin transform pre-conditioning (cross-link:
  `radia_mcp.fem(topic="kelvin_transform")`)

## Where this beats PINN

- **Sharp material interfaces**: GNN handles by edge feature; PINN
  needs domain decomposition
- **Many similar geometries**: train once, run on new mesh in ms;
  PINN re-train per geometry
- **Adaptive mesh refinement**: GNN naturally adapts; PINN
  collocation needs re-sampling

## Where PINN beats GNN

- **No mesh available**: PINN samples; GNN needs mesh first
- **Inverse parameter ID**: PINN gradient w.r.t. params is clean
- **Differential operators inside loss**: PINN composes autograd
  naturally; GNN needs custom layers

## Skeleton (PyTorch + PyG)

```python
import torch
import torch.nn as nn
from torch_geometric.nn import MessagePassing
from torch_geometric.data import Data

class EMPDELayer(MessagePassing):
    def __init__(self, hidden_dim):
        super().__init__(aggr='mean')
        self.mlp_message = nn.Sequential(
            nn.Linear(2*hidden_dim + 1, hidden_dim),  # node_i + node_j + edge_feat
            nn.ReLU(), nn.Linear(hidden_dim, hidden_dim))
        self.mlp_update = nn.Sequential(
            nn.Linear(2*hidden_dim, hidden_dim),
            nn.ReLU(), nn.Linear(hidden_dim, hidden_dim))

    def forward(self, x, edge_index, edge_attr):
        out = self.propagate(edge_index, x=x, edge_attr=edge_attr)
        return self.mlp_update(torch.cat([x, out], dim=-1))

    def message(self, x_i, x_j, edge_attr):
        return self.mlp_message(torch.cat([x_i, x_j, edge_attr], dim=-1))


class EMPDESurrogate(nn.Module):
    def __init__(self, in_dim, hidden_dim, n_layers, out_dim):
        super().__init__()
        self.encoder = nn.Linear(in_dim, hidden_dim)
        self.layers = nn.ModuleList([
            EMPDELayer(hidden_dim) for _ in range(n_layers)])
        self.decoder = nn.Linear(hidden_dim, out_dim)

    def forward(self, data):
        x = self.encoder(data.x)
        for layer in self.layers:
            x = x + layer(x, data.edge_index, data.edge_attr)  # residual
        x_out = self.decoder(x)
        # HARD enforce Dirichlet:
        x_out[data.dirichlet_mask] = data.dirichlet_values
        return x_out
```

## Cross-references

- `radia_mcp.pinn` — soft physics constraints alternative
- `radia_mcp.fem` — for the mesh + Kelvin transform conventions
- `radia_mcp.topology_optimization` — GNN topology decoder
"""


EQUIVARIANT_GNN = r"""
# Equivariant GNN: Isometry-respecting architectures

`Isometric Transformation Invariant and Equivariant Graph
Convolutional Networks.pdf` (public-safe curated corpus)

## Definition

A function `f` is **equivariant** under group G if `f(g . x) = g . f(x)`
for all g in G. For 3D EM, the natural group is:
- **SO(3)** (rotations) — fields rotate with geometry
- **E(3)** = SO(3) + translations — full Euclidean group

A function is **invariant** if `f(g . x) = f(x)` (scalar output:
energy, loss).

## Why equivariance matters for EM

- A rotated coil produces a rotated B field (vector → vector
  equivariant)
- Energy / loss / inductance are rotation-invariant scalars
- A model that DOESN'T respect this needs to learn the symmetry
  from data (huge data inefficiency)

## Architecture families

| Architecture | Equivariance | Reference |
|--------------|--------------|-----------|
| **SchNet** | Translation + rotation invariance (scalar features only) | Schütt 2017 |
| **E(n)-GNN** | E(n) equivariance via radial features + cross product | Satorras-Hoogeboom 2021 |
| **NequIP** | E(3) equivariant up to tensor order L (for forces in MD) | Batzner 2022 |
| **MACE** | E(3) equivariant + many-body messages | Batatia 2022 |
| **eSCN** | Efficient SO(3) convolution | Passaro-Zitnick 2023 |

## Quick decision

| Need | Use |
|------|-----|
| Scalar prediction (energy, loss) | SchNet — simple |
| Vector prediction (force, B-field) | E(n)-GNN or NequIP |
| High accuracy on molecular forces | MACE / NequIP |
| Best speed | eSCN |

## E(n)-GNN skeleton

```python
import torch
import torch.nn as nn
from torch_geometric.nn import MessagePassing

class EquivariantLayer(MessagePassing):
    def __init__(self, hidden_dim):
        super().__init__(aggr='sum')
        self.mlp_h = nn.Sequential(
            nn.Linear(2*hidden_dim + 1, hidden_dim),
            nn.SiLU(), nn.Linear(hidden_dim, hidden_dim))
        self.mlp_x = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.SiLU(), nn.Linear(hidden_dim, 1))

    def forward(self, h, x, edge_index):
        # h: invariant features [N, hidden_dim]
        # x: equivariant coords [N, 3]
        row, col = edge_index
        r_ij = x[row] - x[col]
        d_ij_sq = (r_ij**2).sum(dim=-1, keepdim=True)
        # message based on invariant features + distance
        m_ij = self.mlp_h(torch.cat([h[row], h[col], d_ij_sq], dim=-1))
        # update h via message-passing (invariant)
        h_new = h + self.propagate(edge_index, m_ij=m_ij)
        # update x equivariantly: weighted sum of r_ij
        w_ij = self.mlp_x(m_ij)
        x_new = x + self.propagate(edge_index, r_ij=r_ij, w_ij=w_ij)
        return h_new, x_new

    def message(self, m_ij=None, r_ij=None, w_ij=None):
        if r_ij is None:
            return m_ij
        return w_ij * r_ij  # equivariant message
```

## EM application: rotation-equivariant motor torque predictor

```python
# Input: stator/rotor geometry as graph; rotor angle as edge feature
# Output: equivariant force/torque vector

class MotorTorqueGNN(nn.Module):
    def __init__(self, hidden_dim=128, n_layers=5):
        super().__init__()
        self.embed = nn.Embedding(N_materials, hidden_dim)
        self.layers = nn.ModuleList([
            EquivariantLayer(hidden_dim) for _ in range(n_layers)])
        self.torque_head = nn.Linear(hidden_dim, 1)  # invariant scalar

    def forward(self, data):
        h = self.embed(data.material_id)
        x = data.position
        for layer in self.layers:
            h, x = layer(h, x, data.edge_index)
        # Pool to graph-level scalar
        return self.torque_head(h.mean(dim=0))
```

Trained on a dataset of FEM-computed (geometry, torque) pairs, this
can replace the FEM solve in MCTS / Optuna loops at ~ms per call
(vs minutes for FEM).

## Cross-references

- `radia_mcp.differential_forms` — formal symmetry / Lie group theory
- `radia_mcp.motor` — application target
"""


def get_gnn_documentation(topic: str = "all") -> str:
    """Dispatch GNN topics.

    Topics:
        overview            - GNN families + when for EM
        physics_embedded    - Hard-constraint PDE GNN (vs PINN soft)
        equivariant_gnn     - E(n)/SO(3) equivariant for vector fields
        all                 - Everything
    """
    topic = topic.lower().strip()
    if topic in ("overview", "intro"):
        return OVERVIEW
    if topic in ("physics_embedded", "pde_gnn", "physics_gnn"):
        return PHYSICS_EMBEDDED
    if topic in ("equivariant_gnn", "equivariant", "e_n_gnn",
                 "isometry"):
        return EQUIVARIANT_GNN
    if topic == "all":
        return "\n\n".join([OVERVIEW, PHYSICS_EMBEDDED, EQUIVARIANT_GNN])
    return (f"Unknown topic '{topic}'. Available: overview, "
            "physics_embedded, equivariant_gnn, all.")
