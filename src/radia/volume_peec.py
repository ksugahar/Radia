"""Finite-volume pulse-current topology for explicitly supplied rectilinear metal.

Linear nonmagnetic MQS model: finite conductor volume, no surface impedance or
skin-depth substitution. Each oriented branch carries uniform J=I/area on its
rectangular dual support; orthogonal components may overlap, parallel supports
have disjoint interiors. Nodal KCL is integral conservation, not a claim of
pointwise H(div) conformity. Geometry/current convergence is required.

No CAD inference, staircase approximation of a supplied CAD, field solver,
material fitting, implicit backend or fallback lives in this module.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np


def _frozen(values, dtype):
    result = np.array(values, dtype=dtype, copy=True)
    result.setflags(write=False)
    return result


@dataclass(frozen=True)
class VolumePEECTopology:
    nodes_m: np.ndarray
    branches: np.ndarray
    support_boxes_m: np.ndarray
    axes: np.ndarray
    signs: np.ndarray
    areas_m2: np.ndarray
    resistance_ohm: np.ndarray
    support_cells: np.ndarray
    half_resistance_ohm: np.ndarray
    cell_indices: np.ndarray
    cell_boxes_m: np.ndarray
    sigma_S_per_m: np.ndarray
    port_plus: int
    port_minus: int

    def native_geometry(self):
        """Return Cartesian bar arrays matching each actual dual support.

        Width/height use the native PEEC frame: x/y bars have width along z,
        z bars along y; height is the other transverse Cartesian coordinate.
        Conductivity is the branch's SERIES-equivalent value, never fitted.
        """
        size = self.support_boxes_m[:, 1]-self.support_boxes_m[:, 0]
        n = len(self.branches)
        direction = np.zeros((n, 3))
        direction[np.arange(n), self.axes] = self.signs
        length = size[np.arange(n), self.axes]
        width = size[np.arange(n), np.where(self.axes == 2, 1, 2)]
        height = size[np.arange(n), np.where(self.axes == 0, 1, 0)]
        sigma = length/(self.resistance_ohm*self.areas_m2)
        return dict(centers=self.support_boxes_m.mean(axis=1),
                    directions=direction, lengths=length, widths=width,
                    heights=height, sigmas=sigma)

    def native_builder(self, *, quadrature_order):
        """Select the finite-volume native operator explicitly; no Zs added.

        Native import/build identity remains the caller's responsibility.
        Missing/new-API-mismatched artifacts raise rather than select filaments.
        """
        from . import peec_matrices
        builder = peec_matrices.PEECBuilder(volume_order=quadrature_order)
        geometry = self.native_geometry()
        for i in range(len(self.branches)):
            builder.add_segment(geometry['centers'][i], geometry['directions'][i],
                                geometry['lengths'][i], geometry['widths'][i],
                                geometry['heights'][i], geometry['sigmas'][i])
        return builder

    def incidence(self):
        """Return node-by-branch CSR: -1 at origin, +1 at destination."""
        from scipy.sparse import csr_matrix
        n = len(self.branches)
        return csr_matrix((np.r_[-np.ones(n), np.ones(n)],
                           (np.r_[self.branches[:, 0], self.branches[:, 1]],
                            np.tile(np.arange(n), 2))),
                          shape=(len(self.nodes_m), n))

    def mean_branch_power_W(self, current_peak_A):
        current = np.asarray(current_peak_A, dtype=complex)
        if current.shape != self.resistance_ohm.shape or not np.isfinite(current).all():
            raise ValueError('Finite complex peak current required for every branch')
        return .5*self.resistance_ohm*np.abs(current)**2

    def mean_cell_power_W(self, current_peak_A):
        """Integrate each branch half in its actual material cell, no heat fitting."""
        current = np.asarray(current_peak_A, dtype=complex)
        if current.shape != self.resistance_ohm.shape or not np.isfinite(current).all():
            raise ValueError('Finite complex peak current required for every branch')
        result = np.zeros(len(self.cell_indices))
        power = .5*self.half_resistance_ohm*np.abs(current[:, None])**2
        for half in (0, 1):
            active = self.support_cells[:, half] >= 0
            np.add.at(result, self.support_cells[active, half], power[active, half])
        return result


def build_volume_peec(grid_edges_m, metal_cells, conductivity_S_per_m, *,
                      source_faces, sink_faces):
    """Build one connected metal volume, with ONLY physical contact faces tied.

    grid_edges_m contains x/y/z edges; metal_cells is a boolean (nx,ny,nz)
    occupancy of those declared exact boxes. Material can be a positive scalar
    or a cell-shaped array (inactive cells are ignored). Contacts are explicit
    tuples (i,j,k,axis,side), side=-1/+1, on EXTERIOR metal faces. Each contact
    is equipotential; no interior plane or cell nodes are merged.

    The resistance of an internal branch is the SERIES integral over its two
    cell halves, preserving conductivity jumps. This topology does not supply
    the finite-volume Green inductance operator; callers must select it, and
    must never add half-space Zs to this finite-volume current model.
    """
    if len(grid_edges_m) != 3:
        raise ValueError('Explicit x/y/z SI grid edges required')
    edges = [np.asarray(e, dtype=float) for e in grid_edges_m]
    if any(e.ndim != 1 or len(e) < 2 or not np.isfinite(e).all()
           or np.any(np.diff(e) <= 0) for e in edges):
        raise ValueError('Finite strictly increasing SI edges required')
    mask = np.asarray(metal_cells)
    shape = tuple(len(e)-1 for e in edges)
    if mask.dtype != np.bool_ or mask.shape != shape or not mask.any():
        raise ValueError('Nonempty boolean metal occupancy must match the grid')
    sigma = np.asarray(conductivity_S_per_m, dtype=float)
    if sigma.ndim == 0:
        sigma = np.full(shape, float(sigma))
    if sigma.shape != shape or not np.isfinite(sigma[mask]).all() or np.any(sigma[mask] <= 0):
        raise ValueError('Every metal cell requires finite positive conductivity')
    cells = [tuple(map(int, c)) for c in np.argwhere(mask)]
    ids = {c: i for i, c in enumerate(cells)}
    widths = [np.diff(e) for e in edges]
    centers = [(e[:-1]+e[1:])/2 for e in edges]
    positions = [[centers[a][c[a]] for a in range(3)] for c in cells]
    boxes = [[[edges[a][c[a]] for a in range(3)],
              [edges[a][c[a]+1] for a in range(3)]] for c in cells]
    # No isolated islands are silently discarded or shorted to a port.
    visited = {cells[0]};stack = [cells[0]]
    while stack:
        c = stack.pop()
        for axis in range(3):
            for step in (-1, 1):
                q = list(c);q[axis] += step;q = tuple(q)
                if q in ids and q not in visited:
                    visited.add(q);stack.append(q)
    if len(visited) != len(cells):
        raise ValueError('One connected metal component required; floating islands unsupported')
    contacts = []
    used = set()
    for faces in (source_faces, sink_faces):
        parsed = []
        for raw in faces:
            if len(raw) != 5 or any(isinstance(x, (bool, np.bool_))
                                   or not isinstance(x, (int, np.integer)) for x in raw):
                raise ValueError('Contact face must be integer (i,j,k,axis,side)')
            i,j,k,axis,side = map(int, raw);c=(i,j,k)
            if c not in ids or axis not in (0,1,2) or side not in (-1,1):
                raise ValueError('Contact must reference a metal cell and oriented face')
            q=list(c);q[axis]+=side
            if tuple(q) in ids:
                raise ValueError('Contact face must be exterior, not an internal short')
            face=(*c,axis,side)
            if face in used:
                raise ValueError('A contact face may belong to exactly one port')
            used.add(face);parsed.append(face)
        if not parsed:
            raise ValueError('Both explicitly declared contacts must be nonempty')
        contacts.append(parsed)
    nodes = list(positions)
    for faces in contacts:
        points=[];weights=[]
        for i,j,k,a,side in faces:
            c=(i,j,k);p=list(positions[ids[c]]);p[a]=edges[a][c[a]+(side>0)]
            area=np.prod([widths[d][c[d]] for d in range(3) if d!=a])
            points.append(p);weights.append(area)
        nodes.append(np.average(points,axis=0,weights=weights))
    branches=[];supports=[];axes=[];signs=[];areas=[];r=[];members=[];halves=[]
    def append(start,end,box,axis,sign,parts):
        area=float(np.prod([box[1][a]-box[0][a] for a in range(3) if a!=axis]))
        branches.append((start,end));supports.append(box);axes.append(axis);signs.append(sign);areas.append(area)
        indices=[-1,-1];loss=[0.,0.]
        for slot,(cell,length) in enumerate(parts):
            indices[slot]=ids[cell];loss[slot]=float(length/(sigma[cell]*area))
        r.append(sum(loss));members.append(indices);halves.append(loss)
    for c in cells:
        for a in range(3):
            q=list(c);q[a]+=1;q=tuple(q)
            if q not in ids:continue
            box=[list(boxes[ids[c]][0]),list(boxes[ids[c]][1])]
            box[0][a]=positions[ids[c]][a];box[1][a]=positions[ids[q]][a]
            append(ids[c],ids[q],box,a,1,[(c,widths[a][c[a]]/2),(q,widths[a][q[a]]/2)])
    for port,faces in enumerate(contacts):
        for i,j,k,a,side in faces:
            c=(i,j,k);node=ids[c];port_node=len(cells)+port
            box=[list(boxes[node][0]),list(boxes[node][1])]
            if side==-1:box[1][a]=positions[node][a]
            else:box[0][a]=positions[node][a]
            start,end=(port_node,node) if port==0 else (node,port_node)
            append(start,end,box,a,-side if port==0 else side,[(c,widths[a][c[a]]/2)])
    return VolumePEECTopology(*[_frozen(v,t) for v,t in (
        (nodes,float),(branches,np.int64),(supports,float),(axes,np.int8),
        (signs,np.int8),(areas,float),(r,float),(members,np.int64),
        (halves,float),(cells,np.int64),(boxes,float),
        ([sigma[c] for c in cells],float))],len(cells),len(cells)+1)
