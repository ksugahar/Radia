"""Conservative longitudinal/surface-mesh PEEC topology for swept conductors.

This is the topology and reference-solve layer, not a surface-impedance or
partial-inductance kernel.  ``rings`` are CAD-derived, ordered points on each
cross-section perimeter.  Axial branches exist everywhere; circumferential
branches exist only at selected interior stations.  Every transition shares
the *individual* axial-lane nodes with the surface mesh.  No zero-impedance
fan or artificial equipotential plane is inserted at a transition.

The reference solve accepts independently assembled R and L.  It deliberately
does not reuse the series-bundle reduction, which would suppress all local
transverse currents.  A physically validated fin backend must supply surface
partial elements and SIBC consistently before this topology is used in IH.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class HybridSurfaceTopology:
    nodes: np.ndarray
    branches: np.ndarray
    branch_xyz: np.ndarray
    branch_kind: tuple[str, ...]
    port_plus: int
    port_minus: int
    n_stations: int
    n_lanes: int
    mesh_stations: tuple[int, ...]

    def incidence(self):
        """Node-by-branch incidence, with -1 at source and +1 at sink."""
        from scipy.sparse import csr_matrix

        m = len(self.branches)
        rows = np.r_[self.branches[:, 0], self.branches[:, 1]]
        cols = np.tile(np.arange(m), 2)
        data = np.r_[-np.ones(m), np.ones(m)]
        return csr_matrix((data, (rows, cols)),
                          shape=(len(self.nodes), m))

    def solve_reference(self, frequency, resistance, inductance,
                        port_current=1.0, emf=None):
        """Small dense MNA reference solve; return branch currents and KCL.

        ``resistance`` may be a positive diagonal vector or a symmetric
        matrix.  ``inductance`` is the full symmetric partial-inductance
        matrix.  This is for topology verification, not a production PEEC
        material model or a claim of physical fin accuracy.
        """
        m = len(self.branches)
        r = np.asarray(resistance, dtype=float)
        l = np.asarray(inductance, dtype=float)
        if r.shape == (m,):
            r = np.diag(r)
        if r.shape != (m, m) or l.shape != (m, m):
            raise ValueError("R and L must be branch-sized matrices")
        if not np.all(np.isfinite(r)) or not np.all(np.isfinite(l)):
            raise ValueError("R and L must be finite")
        if not np.allclose(r, r.T) or not np.allclose(l, l.T):
            raise ValueError("R and L must be symmetric")
        if np.min(np.linalg.eigvalsh(r)) <= 0:
            raise ValueError("R must be positive definite")
        if np.min(np.linalg.eigvalsh(l)) < -1e-12 * max(1.0, np.linalg.norm(l)):
            raise ValueError("L must be positive semidefinite")
        z = r + 2j * np.pi * float(frequency) * l
        a = self.incidence().toarray()
        keep = np.arange(len(self.nodes)) != self.port_minus
        ah = a[keep]
        rhs_kcl = np.zeros(len(self.nodes), dtype=complex)
        rhs_kcl[self.port_plus] = -complex(port_current)
        rhs_kcl[self.port_minus] = complex(port_current)
        drive = np.zeros(m, dtype=complex) if emf is None else np.asarray(
            emf, dtype=complex)
        if drive.shape != (m,):
            raise ValueError("emf must have one value per branch")
        mat = np.block([[z, ah.T], [ah, np.zeros((len(ah), len(ah)))]])
        sol = np.linalg.solve(mat, np.r_[drive, rhs_kcl[keep]])
        current = sol[:m]
        residual = a @ current - rhs_kcl
        return current, residual


def build_hybrid_surface_topology(rings, mesh_stations):
    """Connect CAD perimeter rings without shorting transition lanes.

    ``rings`` has shape (station, lane, xyz).  Each lane is continuous
    along the conductor.  Circumferential branches form closed rings only
    at ``mesh_stations``; those stations must be interior.  Both terminal
    cross-sections are equipotential *ports*, not interior transitions.
    """
    xyz = np.asarray(rings, dtype=float)
    if xyz.ndim != 3 or xyz.shape[2] != 3:
        raise ValueError("rings must have shape (stations, lanes, 3)")
    s_count, k_count, _ = xyz.shape
    if s_count < 3 or k_count < 3 or not np.all(np.isfinite(xyz)):
        raise ValueError("need >=3 finite stations and lanes")
    stations = tuple(sorted({int(i) for i in mesh_stations}))
    if any(i <= 0 or i >= s_count - 1 for i in stations):
        raise ValueError("mesh stations must be interior")
    if any(np.linalg.norm(xyz[s + 1, k] - xyz[s, k]) <= 0
           for s in range(s_count - 1) for k in range(k_count)):
        raise ValueError("longitudinal branches must have nonzero length")
    if any(np.linalg.norm(xyz[s, (k + 1) % k_count] - xyz[s, k]) <= 0
           for s in stations for k in range(k_count)):
        raise ValueError("mesh ring branches must have nonzero length")

    # Interior nodes are unique per lane; only real terminal planes are tied.
    center0 = np.mean(xyz[0], axis=0)
    center1 = np.mean(xyz[-1], axis=0)
    nodes = [center0, center1]
    node_id = np.empty((s_count, k_count), dtype=int)
    node_id[0, :] = 0
    node_id[-1, :] = 1
    for s in range(1, s_count - 1):
        for k in range(k_count):
            node_id[s, k] = len(nodes)
            nodes.append(xyz[s, k])
    branches = []
    branch_xyz = []
    kinds = []
    for s in range(s_count - 1):
        for k in range(k_count):
            branches.append((node_id[s, k], node_id[s + 1, k]))
            branch_xyz.append((xyz[s, k], xyz[s + 1, k]))
            kinds.append("longitudinal")
    for s in stations:
        for k in range(k_count):
            branches.append((node_id[s, k], node_id[s, (k + 1) % k_count]))
            branch_xyz.append((xyz[s, k], xyz[s, (k + 1) % k_count]))
            kinds.append("transverse")
    return HybridSurfaceTopology(
        np.asarray(nodes, dtype=float), np.asarray(branches, dtype=int),
        np.asarray(branch_xyz, dtype=float), tuple(kinds),
        0, 1, s_count, k_count, stations)


def rings_from_filament_paths(filament_paths):
    """Recover ordered CAD cross-section samples from equal-station paths.

    A different number of pieces per lane is not silently interpolated:
    mismatched stations would create nonphysical transverse connections.
    """
    if len(filament_paths) < 3:
        raise ValueError("need at least three perimeter paths")
    n_piece = len(filament_paths[0])
    if n_piece < 2 or any(len(path) != n_piece for path in filament_paths):
        raise ValueError("perimeter paths must share at least two pieces")
    xyz = np.empty((n_piece + 1, len(filament_paths), 3), dtype=float)
    for k, path in enumerate(filament_paths):
        xyz[0, k] = path[0][0]
        for s, (start, end) in enumerate(path):
            if not np.allclose(xyz[s, k], start, rtol=0.0, atol=1e-10):
                raise ValueError("perimeter path has a disconnected piece")
            xyz[s + 1, k] = end
    if not np.all(np.isfinite(xyz)):
        raise ValueError("perimeter paths must be finite")
    return xyz


def assemble_experimental_fin_peec(graph, rings, *, sigma, sheet_depth):
    """Assemble an exploratory thin-sheet PEEC with the existing C++ kernel.

    The local transverse width is the dual axial spacing and the axial
    width is the dual perimeter spacing.  Both have ``sheet_depth`` as the
    other cross-section dimension.  This supplies a real, conservative
    branched PEEC reference but is NOT an accepted SIBC discretization:
    the existing rectangular partial-element kernel has not been validated
    for a folded thin sheet, corners, or its internal-reactance split.
    Consequently this function is intentionally absent from IH production.
    """
    from radia.peec_matrices import PEECBuilder
    from radia.peec_topology import PEECCircuitSolver

    xyz = np.asarray(rings, dtype=float)
    if xyz.shape != (graph.n_stations, graph.n_lanes, 3):
        raise ValueError("rings do not match the graph")
    if sigma <= 0 or sheet_depth <= 0:
        raise ValueError("sigma and sheet_depth must be positive")
    s_count, k_count, _ = xyz.shape
    circumference_dual = np.empty((s_count, k_count))
    for k in range(k_count):
        circumference_dual[:, k] = 0.5 * (
            np.linalg.norm(xyz[:, k] - xyz[:, (k - 1) % k_count], axis=1)
            + np.linalg.norm(xyz[:, (k + 1) % k_count] - xyz[:, k], axis=1))
    axial_dual = np.empty((s_count, k_count))
    for s in range(1, s_count - 1):
        axial_dual[s] = 0.5 * (
            np.linalg.norm(xyz[s] - xyz[s - 1], axis=1)
            + np.linalg.norm(xyz[s + 1] - xyz[s], axis=1))
    axial_dual[0] = np.linalg.norm(xyz[1] - xyz[0], axis=1)
    axial_dual[-1] = np.linalg.norm(xyz[-1] - xyz[-2], axis=1)
    widths = []
    for s in range(s_count - 1):
        for k in range(k_count):
            widths.append(0.5 * (circumference_dual[s, k]
                                 + circumference_dual[s + 1, k]))
    for s in graph.mesh_stations:
        for k in range(k_count):
            widths.append(0.5 * (axial_dual[s, k]
                                 + axial_dual[s, (k + 1) % k_count]))
    if np.min(widths) <= 0:
        raise ValueError("zero dual width in CAD surface mesh")

    builder = PEECBuilder()
    for endpoints, width in zip(graph.branch_xyz, widths):
        p0, p1 = endpoints
        n0 = builder.add_node_at(*map(float, p0))
        n1 = builder.add_node_at(*map(float, p1))
        builder.add_connected_segment(n0, n1, float(width),
                                      float(sheet_depth), sigma=float(sigma))
    topo = builder.build_topology()
    if int(topo["n_loop"]) != len(graph.branches):
        raise RuntimeError("PEEC kernel dropped a surface branch")
    topo["segment_nodes"] = graph.branches.copy()
    topo["n_nodes"] = len(graph.nodes)
    topo["ports"] = [(graph.port_plus, graph.port_minus, 0)]
    return PEECCircuitSolver(topo), np.asarray(widths, dtype=float)


def build_hybrid_surface_topology_from_step(step_path, *, sigma,
                                            n_peri=16, mesh_stations=None,
                                            cad_units_per_meter=1.0):
    """Extract one STEP conductor's sampled perimeter and build its graph.

    The existing CAD extractor already validates solid coverage.  Only its
    direct surface/section paths are allowed: an equivalent-circle fallback
    would silently remove the fin.  ``unknown`` cross-sections are meshed at
    every interior station by default.  This is intentionally conservative;
    selective removal requires a validated AC error indicator later.

    Returns ``(graph, cad_metadata)``.  The extractor's PEEC solver is not
    reused, because its series-bundle topology has no transverse branches.
    """
    from radia.coil_from_cad import filaments_from_step

    topo = filaments_from_step(
        str(step_path), sigma=float(sigma), n_peri=int(n_peri),
        cad_units_per_meter=float(cad_units_per_meter), use_coil_builder=True)
    source = str(topo.get("source", ""))
    if source not in {"step_uv", "step_per_station", "step_section_planes"}:
        raise ValueError(
            "fin surface PEEC requires direct CAD perimeter extraction; "
            f"got {source!r}")
    paths = topo["filament_paths"]
    if len(paths) != n_peri:
        raise ValueError("first fin implementation requires one n_peri-lane conductor")
    rings = rings_from_filament_paths(paths)
    kind = str(topo.get("cross_section_kind", "unknown"))
    if mesh_stations is None:
        mesh_stations = (range(1, len(rings) - 1)
                         if kind not in {"rect", "circle"} else ())
    graph = build_hybrid_surface_topology(rings, mesh_stations)
    return graph, {"cad_source": source, "cross_section_kind": kind,
                   "mesh_stations": graph.mesh_stations,
                   "n_lanes": graph.n_lanes,
                   "n_stations": graph.n_stations}


def build_hybrid_surface_topology_from_straight_prism_step(
        step_path, *, n_peri=32, n_stations=9, cad_units_per_meter=1000.0):
    """Directly section a single z-extruded STEP solid, including a fin tip.

    This deliberately narrow fixture route does not infer a curved coil spine.
    Constant section area and a single face at each station are required.
    """
    from build123d import Plane, import_step, section

    from radia.coil_from_cad import _sample_face_perimeter_in_pt_frame

    if n_peri < 4 or n_stations < 3 or cad_units_per_meter <= 0:
        raise ValueError("invalid perimeter, station, or CAD-unit count")
    shape = import_step(str(step_path))
    solids = shape.solids()
    if len(solids) != 1:
        raise ValueError("straight-prism route requires exactly one solid")
    solid = solids[0]
    bbox = solid.bounding_box()
    z0, z1 = float(bbox.min.Z), float(bbox.max.Z)
    if z1 <= z0:
        raise ValueError("zero-length STEP solid")
    rings = []
    areas = []
    for z in np.linspace(z0, z1, n_stations):
        # Intersect strictly inside the end caps to avoid coincident-face
        # ambiguity in OpenCascade's section operation.
        probe_z = float(np.clip(z, z0 + 1e-7 * (z1 - z0),
                                z1 - 1e-7 * (z1 - z0)))
        cut = section(solid, section_by=Plane(origin=(0, 0, probe_z),
                                              z_dir=(0, 0, 1)))
        faces = cut.faces()
        if len(faces) != 1:
            raise ValueError("STEP section must have exactly one face")
        face = faces[0]
        areas.append(float(face.area))
        center = face.center()
        uv = _sample_face_perimeter_in_pt_frame(
            face, np.array([center.X, center.Y, center.Z]),
            np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]),
            int(n_peri))
        ring = np.column_stack((uv[:, 0] + center.X,
                                uv[:, 1] + center.Y,
                                np.full(n_peri, z))) / cad_units_per_meter
        rings.append(ring)
    if not np.allclose(areas, areas[0], rtol=1e-6, atol=1e-9):
        raise ValueError("STEP is not a constant-section straight prism")
    if not np.allclose(np.asarray(rings)[:, :, :2], rings[0][:, :2],
                       rtol=0, atol=1e-9):
        raise ValueError("STEP section changes along the extrusion axis")
    graph = build_hybrid_surface_topology(
        np.asarray(rings), range(1, n_stations - 1))
    return graph, {"cad_source": "step_straight_prism_sections",
                   "cross_section_kind": "unknown",
                   "section_area_m2": areas[0] / cad_units_per_meter**2,
                   "n_lanes": graph.n_lanes,
                   "n_stations": graph.n_stations}
