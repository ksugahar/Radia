"""
CoilGeometry — CAD first, analysis model second.

build123d CAD shape から解析モデルを自動生成する。
2つの用途に対応:

- DC（加速器電磁石）: to_wire_segments()    → Biot-Savart 線電流源
- AC（誘導加熱）    : to_surface_current()  → Biot-Savart 表面電流源

Usage:
    from build123d import *
    from radia.coil_geometry_unified import CoilGeometry

    coil_solid = sweep(Circle(5), helix_path)
    geo = CoilGeometry(coil_solid, current=1000, centerline=helix_path)

    # DC (accelerator): line current Biot-Savart
    segments, I = geo.to_wire_segments()

    # AC (IH): surface current Biot-Savart (path-length + skin effect)
    panels = geo.to_surface_current(freq=10e3)

    # Manufacturing
    geo.export_step("coil.step")
    geo.export_stl("coil.stl")
"""

import numpy as np


class CoilGeometry:
    """Unified coil representation derived from CAD geometry.

    The CAD solid (build123d Shape) is the single source of truth.

    DC (accelerator):
        Centerline → wire segments → Biot-Savart line current source.
        Coil is NOT meshed. Field is computed analytically.

    AC (IH):
        CAD surface → surface current panels → Biot-Savart source.
        Coil is NOT meshed. Current distribution accounts for:
        1. Path-length effect (循環電流, J ∝ 1/r at curved sections)
        2. Skin effect (表皮効果, exponential decay from surface)
        Proximity effect is NOT included (use FEM for precision).

    Attributes:
        solid: build123d Shape (the CAD truth)
        current: total current [A]
        freq: operating frequency [Hz] (0 for DC)
        sigma: conductor conductivity [S/m]
    """

    def __init__(self, solid, current, freq=0, sigma=5.8e7, centerline=None):
        """
        Args:
            solid: build123d Part/Solid/Compound
            current: coil current [A]
            freq: operating frequency [Hz], 0 for DC
            sigma: conductivity [S/m], default copper
            centerline: build123d Wire or Edge for coil center path.
                        Required for DC (to_wire_segments).
                        Not needed for AC (to_mesh_region).
        """
        self.solid = solid
        self.current = current
        self.freq = freq
        self.sigma = sigma
        self._centerline = centerline

    # ================================================================
    # Centerline
    # ================================================================

    @property
    def centerline_points(self):
        """Ordered centerline points as (N, 3) numpy array.

        Raises ValueError if no centerline was provided.
        """
        if self._centerline is None:
            raise ValueError(
                "No centerline. Provide centerline= in constructor.")

        n = max(100, int(self._centerline.length * 1000))
        positions = [self._centerline.position_at(i / n)
                     for i in range(n + 1)]
        return np.array([(p.X, p.Y, p.Z) for p in positions])

    # ================================================================
    # Closure check
    # ================================================================

    @property
    def is_closed(self):
        """True if the centerline forms a closed loop (gap < 1µm)."""
        return self.gap < 1e-6

    @property
    def gap(self):
        """Distance between start and end of centerline [m]."""
        pts = self.centerline_points
        return float(np.linalg.norm(pts[-1] - pts[0]))

    def close(self, tolerance=1e-6):
        """Verify the coil forms a closed loop.

        Raises ValueError if gap exceeds tolerance.
        Closure is achieved in build123d by constructing a closed Wire
        (Polyline(..., close=True), full-circle Edge, etc.).
        This method validates, not repairs.

        Returns:
            self
        """
        g = self.gap
        if g > tolerance:
            raise ValueError(
                f"Coil not closed: gap = {g:.2e} m. "
                f"Build a closed Wire in build123d.")
        return self

    # ================================================================
    # DC output: Biot-Savart wire segments (accelerator)
    # ================================================================

    def to_wire_segments(self, n=None):
        """Centerline wire segments for DC Biot-Savart.

        The coil is represented as a thin wire carrying total current I.
        Each segment is a straight line between consecutive centerline points.

        Args:
            n: number of segments (None = use full centerline resolution)

        Returns:
            (segments, current) where:
                segments: list of ((x1,y1,z1), (x2,y2,z2))
                current: coil current [A]
        """
        pts = self.centerline_points
        if n is not None and n < len(pts) - 1:
            indices = np.linspace(0, len(pts) - 1, n + 1, dtype=int)
            pts = pts[indices]

        segments = [(tuple(pts[i]), tuple(pts[i + 1]))
                    for i in range(len(pts) - 1)]
        return segments, self.current

    # ================================================================
    # AC output: surface current Biot-Savart (IH)
    # ================================================================

    def to_surface_current(self, freq=None, sigma=None):
        """Surface current Biot-Savart source for IH coils.

        Tessellates the conductor surface and assigns surface current
        density K [A/m] to each panel, accounting for:
        1. Path-length effect (循環電流): J ∝ 1/r at curved sections
        2. Skin effect (表皮効果): current concentrated on surface

        Proximity effect is NOT included. For precision validation,
        use FEM with the solid from to_mesh_region().

        Args:
            freq: frequency [Hz] (overrides self.freq)
            sigma: conductivity [S/m] (overrides self.sigma)

        Returns:
            dict with:
                vertices: (Nv, 3) array
                triangles: (Nt, 3) index array
                normals: (Nt, 3) outward normals
                centers: (Nt, 3) triangle centroids
                areas: (Nt,) triangle areas [m²]
                K: (Nt, 3) surface current density [A/m]
                skin_depth: δ [m]
        """
        freq = freq or self.freq
        sigma = sigma or self.sigma
        mu_0 = 4 * np.pi * 1e-7

        if freq <= 0:
            raise ValueError(
                "freq must be > 0 for surface current. "
                "Use to_wire_segments() for DC.")

        delta = np.sqrt(2 / (2 * np.pi * freq * mu_0 * sigma))

        # Tessellate solid surface
        verts_raw, tris_raw = self.solid.tessellate(0.001)
        verts = np.array([(v.X, v.Y, v.Z) for v in verts_raw])
        tris = np.array(tris_raw)

        centers, normals, areas = _triangle_properties(verts, tris)

        # Centerline tangents
        pts = self.centerline_points
        tangent_mids = (pts[:-1] + pts[1:]) / 2
        tangents = np.diff(pts, axis=0)
        t_norms = np.linalg.norm(tangents, axis=1, keepdims=True)
        tangents = tangents / np.maximum(t_norms, 1e-30)

        # Curvature at each centerline station
        curvatures = np.zeros(len(tangent_mids))
        curvature_dirs = np.zeros_like(tangent_mids)
        for i in range(1, len(tangents) - 1):
            dv = tangents[i] - tangents[i - 1]
            ds = np.linalg.norm(tangent_mids[i] - tangent_mids[i - 1])
            if ds > 1e-15:
                curvatures[i] = np.linalg.norm(dv) / ds
                if curvatures[i] > 1e-15:
                    curvature_dirs[i] = dv / np.linalg.norm(dv)

        # Compute K for each surface panel
        K = np.zeros_like(centers)
        for i, (c, n, a) in enumerate(zip(centers, normals, areas)):
            # Find nearest centerline point
            dists = np.linalg.norm(tangent_mids - c, axis=1)
            nearest = np.argmin(dists)
            t = tangents[nearest]

            # Project tangent onto surface plane
            t_surf = t - np.dot(t, n) * n
            t_norm = np.linalg.norm(t_surf)
            if t_norm < 1e-10:
                continue
            t_surf /= t_norm

            # 1. Path-length effect (1/r)
            path_factor = 1.0
            kappa = curvatures[nearest]
            if kappa > 1e-10:
                R_curv = 1.0 / kappa
                # Radial offset from centerline
                cl_point = tangent_mids[nearest]
                radial_vec = c - cl_point
                # Remove tangential component
                radial_vec -= np.dot(radial_vec, t) * t
                radial_offset = np.dot(radial_vec, curvature_dirs[nearest])
                r_local = R_curv + radial_offset
                if r_local > 1e-10:
                    path_factor = R_curv / r_local

            # K magnitude: I / (circumferential width) * path_factor
            # For a tube: circumference ≈ perimeter of cross-section
            # Approximate: I / (perimeter) where perimeter ~ sqrt(area_total/n_stations)
            # Better: total current distributed over surface with path-length weighting
            K[i] = t_surf * path_factor

        # Normalize so that total current integrates to I
        # ∫ K · dl = I  →  sum(|K| * panel_width_along_current) = I
        # Approximate: sum(|K_i| * sqrt(area_i)) ≈ I
        K_mags = np.linalg.norm(K, axis=1)
        integral = np.sum(K_mags * np.sqrt(areas))
        if integral > 1e-30:
            K *= self.current / integral

        return {
            'vertices': verts,
            'triangles': tris,
            'normals': normals,
            'centers': centers,
            'areas': areas,
            'K': K,
            'skin_depth': delta,
        }

    # ================================================================
    # FEM mesh region (for precision validation)
    # ================================================================

    def to_mesh_region(self, label="coil"):
        """Return CAD solid as labeled mesh region for FEM validation.

        For precision work when the surface current approximation is
        insufficient. Mesh with Netgen (high-order curved elements)
        and solve with full FEM or SIBC on the workpiece side.

        Args:
            label: region label for NGSolve

        Returns:
            build123d Shape with label set
        """
        self.solid.label = label
        return self.solid

    # ================================================================
    # Export
    # ================================================================

    def export_step(self, path):
        """Export to STEP (Cubit interchange, manufacturing)."""
        from build123d import export_step
        export_step(self.solid, str(path))

    def export_stl(self, path, tolerance=0.001):
        """Export to STL (3D printing)."""
        from build123d import export_stl
        export_stl(self.solid, str(path), tolerance=tolerance)

    def export_brep(self, path):
        """Export to BREP (Netgen, lossless OCCT transfer)."""
        from build123d import export_brep
        export_brep(self.solid, str(path))

    # ================================================================
    # Properties
    # ================================================================

    @property
    def skin_depth(self):
        """Skin depth δ [m] at operating frequency. inf for DC."""
        if self.freq <= 0:
            return float('inf')
        mu_0 = 4 * np.pi * 1e-7
        return np.sqrt(2 / (2 * np.pi * self.freq * mu_0 * self.sigma))

    @property
    def volume(self):
        """Coil volume [unit³]."""
        return self.solid.volume

    @property
    def area(self):
        """Coil surface area [unit²]."""
        return self.solid.area

    def info(self):
        """Summary dict of coil geometry and electrical properties."""
        edges = self.solid.edges()
        d = {
            'current': self.current,
            'freq': self.freq,
            'sigma': self.sigma,
            'skin_depth': self.skin_depth,
            'volume': self.volume,
            'area': self.area,
            'face_count': len(self.solid.faces()),
            'edge_count': len(edges),
            'min_edge_length': min(e.length for e in edges) if edges else None,
            'is_valid': self.solid.is_valid,
        }
        if self._centerline is not None:
            d['is_closed'] = self.is_closed
            d['gap'] = self.gap
        return d

    def __repr__(self):
        mode = "DC" if self.freq == 0 else f"AC {self.freq:.0f}Hz"
        return (f"CoilGeometry({mode}, I={self.current}A, "
                f"V={self.volume:.6g})")


# ================================================================
# Module-level helpers
# ================================================================

def _triangle_properties(verts, tris):
    """Compute centers, normals, areas for triangle mesh."""
    p0 = verts[tris[:, 0]]
    p1 = verts[tris[:, 1]]
    p2 = verts[tris[:, 2]]
    centers = (p0 + p1 + p2) / 3
    cross = np.cross(p1 - p0, p2 - p0)
    areas = np.linalg.norm(cross, axis=1) / 2
    normals = cross / (np.linalg.norm(cross, axis=1, keepdims=True) + 1e-30)
    return centers, normals, areas
