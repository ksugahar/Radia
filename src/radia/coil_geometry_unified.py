"""
CoilGeometry — CAD first, analysis model second.

build123d CAD shape から DC/AC 解析モデルを自動生成する。

Usage:
    from build123d import *
    from radia.coil_geometry_unified import CoilGeometry

    coil_solid = sweep(Circle(5), helix_path)
    geo = CoilGeometry(coil_solid, current=1000, centerline=helix_path)

    segments, I = geo.to_wire_segments()         # DC Biot-Savart
    filaments   = geo.to_ac_filaments(freq=10e3) # AC (path-length + skin)
    panels      = geo.to_surface_current()       # AC high-freq (surface current)
    solid       = geo.to_mesh_region("coil")     # FEM/SIBC
    geo.export_step("coil.step")                 # Manufacturing
    geo.export_stl("coil.stl")                   # 3D printing
"""

import numpy as np


class CoilGeometry:
    """Unified coil representation derived from CAD geometry.

    The CAD solid is the single source of truth. Analysis models are
    derived automatically based on the application:

    - DC (accelerator): to_wire_segments() → Biot-Savart
    - AC low-freq: to_ac_filaments() → filament model with skin effect
    - AC high-freq (IH): to_surface_current() → surface current sheet
    - FEM/SIBC: to_mesh_region() → solid for meshing
    - Manufacturing: export_step/stl/brep

    Attributes:
        solid: build123d Shape (the CAD truth)
        current: total current [A]
        freq: operating frequency [Hz] (0 for DC)
        sigma: conductor conductivity [S/m]
    """

    def __init__(self, solid, current, freq=0, sigma=5.8e7, centerline=None):
        self.solid = solid
        self.current = current
        self.freq = freq
        self.sigma = sigma
        self._centerline = centerline

    # ================================================================
    # Centerline
    # ================================================================

    # ================================================================
    # Centerline
    # ================================================================

    @property
    def centerline_points(self):
        """Ordered centerline points as (N, 3) numpy array."""
        if self._centerline is None:
            raise ValueError(
                "No centerline available. Provide centerline= in constructor.")

        n = max(100, int(self._centerline.length * 1000))
        positions = [self._centerline.position_at(i / n)
                     for i in range(n + 1)]
        return np.array([(p.X, p.Y, p.Z) for p in positions])

    # ================================================================
    # Geometry utilities
    # ================================================================

    @property
    def is_closed(self):
        """Check if the centerline forms a closed loop."""
        pts = self.centerline_points
        gap = np.linalg.norm(pts[-1] - pts[0])
        return gap < 1e-6

    @property
    def gap(self):
        """Distance between start and end of centerline [m]."""
        pts = self.centerline_points
        return float(np.linalg.norm(pts[-1] - pts[0]))

    def close(self, tolerance=1e-6):
        """Verify the coil forms a closed loop.

        Raises ValueError if the centerline gap exceeds tolerance.
        For build123d coils, closure is achieved by constructing a
        closed Wire (e.g., using close=True in Polyline, or a
        full-circle Edge). This method validates, not repairs.

        Args:
            tolerance: maximum allowed gap [m]

        Returns:
            self
        """
        g = self.gap
        if g > tolerance:
            raise ValueError(
                f"Coil not closed: gap = {g:.2e} m. "
                f"Ensure the centerline Wire is closed "
                f"(e.g., Polyline(..., close=True) or full-circle Edge).")
        return self

    # ================================================================
    # DC: Biot-Savart wire segments
    # ================================================================

    def to_wire_segments(self, n=None):
        """Centerline wire segments for DC Biot-Savart.

        Args:
            n: number of segments (None = use existing discretization)

        Returns:
            (segments, current) where segments is list of
            ((x1,y1,z1), (x2,y2,z2)) pairs and current is [A].
        """
        pts = self.centerline_points
        if n is not None and n < len(pts) - 1:
            indices = np.linspace(0, len(pts) - 1, n + 1, dtype=int)
            pts = pts[indices]

        segments = [(tuple(pts[i]), tuple(pts[i + 1]))
                    for i in range(len(pts) - 1)]
        return segments, self.current

    # ================================================================
    # AC: filament model (path-length + skin effect)
    # ================================================================

    def to_ac_filaments(self, freq=None, sigma=None, n_path=100,
                        n_cross=16):
        """AC Biot-Savart filament model with path-length and skin effects.

        Three effects modeled (in order of importance for IH):
        1. Path-length effect: J ∝ 1/r at curved sections (DC too)
        2. Skin effect: exponential decay from surface
        3. (Near-field/proximity: NOT included — use SIBC for this)

        Args:
            freq: frequency [Hz] (overrides self.freq)
            sigma: conductivity [S/m] (overrides self.sigma)
            n_path: stations along coil path
            n_cross: filaments per cross-section

        Returns:
            dict with keys:
                positions: (M, 3) filament midpoints
                directions: (M, 3) current direction unit vectors
                currents: (M,) complex currents [A]
                skin_depth: δ [m]
        """
        freq = freq or self.freq
        sigma = sigma or self.sigma
        mu_0 = 4 * np.pi * 1e-7

        if freq <= 0:
            raise ValueError("freq must be > 0 for AC. Use to_wire_segments() for DC.")

        delta = np.sqrt(2 / (2 * np.pi * freq * mu_0 * sigma))

        pts = self.centerline_points
        if len(pts) > n_path + 1:
            indices = np.linspace(0, len(pts) - 1, n_path + 1, dtype=int)
            pts = pts[indices]

        positions = []
        directions = []
        currents = []

        for i in range(len(pts) - 1):
            p1, p2 = pts[i], pts[i + 1]
            mid = (p1 + p2) / 2
            tangent = p2 - p1
            seg_len = np.linalg.norm(tangent)
            if seg_len < 1e-15:
                continue
            tangent_unit = tangent / seg_len

            # Local curvature estimate
            curvature = 0.0
            curvature_dir = np.zeros(3)
            if 0 < i < len(pts) - 2:
                v1 = pts[i] - pts[i - 1]
                v2 = pts[i + 1] - pts[i]
                l1 = np.linalg.norm(v1)
                l2 = np.linalg.norm(v2)
                if l1 > 1e-15 and l2 > 1e-15:
                    dv = v2 / l2 - v1 / l1
                    curvature = np.linalg.norm(dv) / (0.5 * (l1 + l2))
                    if curvature > 1e-15:
                        curvature_dir = dv / np.linalg.norm(dv)

            # Local frame perpendicular to tangent
            n1, n2 = _perpendicular_frame(tangent_unit)

            # Cross-section radius estimate
            cs_radius = self._estimate_cross_section_radius(mid, tangent_unit)

            # Distribute filaments on concentric ring
            for j in range(n_cross):
                angle = 2 * np.pi * j / n_cross
                r_frac = 0.7 * cs_radius  # 70% of radius
                offset = r_frac * (np.cos(angle) * n1 + np.sin(angle) * n2)
                fil_pos = mid + offset

                # 1. Path-length effect
                path_factor = 1.0
                if curvature > 1e-10:
                    R_curv = 1.0 / curvature
                    radial_offset = np.dot(offset, curvature_dir)
                    path_factor = R_curv / max(R_curv + radial_offset, 1e-10)

                # 2. Skin effect
                d_surface = cs_radius - np.linalg.norm(offset)
                skin_factor = np.exp(-(1 + 1j) * d_surface / delta)

                I_fil = (self.current / n_cross) * path_factor * skin_factor

                positions.append(fil_pos)
                directions.append(tangent_unit)
                currents.append(I_fil)

        return {
            'positions': np.array(positions),
            'directions': np.array(directions),
            'currents': np.array(currents),
            'skin_depth': delta,
        }

    # ================================================================
    # AC high-freq: surface current (IH)
    # ================================================================

    def to_surface_current(self):
        """Surface current model for high-frequency (δ << conductor size).

        Tessellates the conductor surface and assigns surface current
        density K [A/m] to each triangle. Current direction is the
        projection of the centerline tangent onto the surface.

        For fast optimization loops. Use SIBC for precision.

        Returns:
            dict with keys:
                vertices: (Nv, 3) array
                triangles: (Nt, 3) index array
                normals: (Nt, 3) outward normals
                centers: (Nt, 3) triangle centroids
                areas: (Nt,) triangle areas
                K: (Nt, 3) surface current density [A/m]
        """
        verts_raw, tris_raw = self.solid.tessellate(0.001)
        verts = np.array([(v.X, v.Y, v.Z) for v in verts_raw])
        tris = np.array(tris_raw)

        centers, normals, areas = _triangle_properties(verts, tris)

        # Centerline tangents for current direction
        pts = self.centerline_points
        tangent_mids = (pts[:-1] + pts[1:]) / 2
        tangents = np.diff(pts, axis=0)
        t_norms = np.linalg.norm(tangents, axis=1, keepdims=True)
        tangents = tangents / np.maximum(t_norms, 1e-30)

        # Circumference estimate for K magnitude
        # K = I / circumference (surface current density)
        try:
            circumference = 2 * np.sqrt(np.pi * self.solid.area /
                                         max(len(pts) - 1, 1))
        except Exception:
            circumference = 1.0

        K = np.zeros_like(centers)
        for i, (c, n) in enumerate(zip(centers, normals)):
            dists = np.linalg.norm(tangent_mids - c, axis=1)
            nearest = np.argmin(dists)
            t = tangents[nearest]

            # Project tangent onto surface (remove normal component)
            t_surf = t - np.dot(t, n) * n
            t_norm = np.linalg.norm(t_surf)
            if t_norm > 1e-10:
                t_surf /= t_norm
            K[i] = t_surf * self.current / max(circumference, 1e-10)

        return {
            'vertices': verts,
            'triangles': tris,
            'normals': normals,
            'centers': centers,
            'areas': areas,
            'K': K,
        }

    # ================================================================
    # FEM / SIBC mesh region
    # ================================================================

    def to_mesh_region(self, label="coil"):
        """Return CAD solid as labeled region for FEM/SIBC meshing.

        Use with Netgen (high-order curved elements) for SIBC accuracy.
        """
        self.solid.label = label
        return self.solid

    # ================================================================
    # Export
    # ================================================================

    def export_step(self, path):
        """Export to STEP (Cubit, interchange)."""
        from build123d import export_step
        export_step(self.solid, str(path))

    def export_stl(self, path, tolerance=0.001):
        """Export to STL (3D printing)."""
        from build123d import export_stl
        export_stl(self.solid, str(path), tolerance=tolerance)

    def export_brep(self, path):
        """Export to BREP (Netgen, lossless)."""
        from build123d import export_brep
        export_brep(self.solid, str(path))

    # ================================================================
    # Properties
    # ================================================================

    @property
    def skin_depth(self):
        """Skin depth δ [m]."""
        if self.freq <= 0:
            return float('inf')
        mu_0 = 4 * np.pi * 1e-7
        return np.sqrt(2 / (2 * np.pi * self.freq * mu_0 * self.sigma))

    @property
    def volume(self):
        return self.solid.volume

    @property
    def area(self):
        return self.solid.area

    def info(self):
        """Summary dict of coil properties."""
        edges = self.solid.edges()
        return {
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

    def __repr__(self):
        return (f"CoilGeometry(I={self.current}A, f={self.freq}Hz, "
                f"δ={self.skin_depth:.4g}m, V={self.volume:.6g})")

    # ================================================================
    # Internal helpers
    # ================================================================

    def _estimate_cross_section_radius(self, point, tangent):
        """Rough cross-section radius at a point along the centerline."""
        try:
            bb = self.solid.bounding_box()
            return min(bb.size) / 2
        except Exception:
            return 0.01


# ================================================================
# Module-level helpers
# ================================================================

def _perpendicular_frame(tangent):
    """Build two unit vectors perpendicular to tangent."""
    if abs(tangent[2]) < 0.9:
        up = np.array([0.0, 0.0, 1.0])
    else:
        up = np.array([1.0, 0.0, 0.0])
    n1 = np.cross(tangent, up)
    n1 /= np.linalg.norm(n1) + 1e-30
    n2 = np.cross(tangent, n1)
    n2 /= np.linalg.norm(n2) + 1e-30
    return n1, n2


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
