"""
Entity wrapper methods for GmshBuilder.

Wraps GMSH entity queries in a Gmsh-based functional pattern,
e.g., gb.volume_surfaces(tag) for entity queries.
"""

import math


class EntityWrapperMixin:
    """Mixin providing Gmsh-based entity wrapper queries.

    All methods operate on raw GMSH entity tags (not GmshBuilder volume IDs).
    This provides a flat, functional API mirroring entity query
    methods for volumes, surfaces, curves, and vertices.
    """

    # ------------------------------------------------------------------
    # Volume queries (14 methods)
    # ------------------------------------------------------------------

    def volume_surfaces(self, vol_tag):
        """Get surface tags bounding a volume.

        Parameters
        ----------
        vol_tag : int
            GMSH volume tag.

        Returns
        -------
        list of int
            Sorted list of surface tags bounding the volume.
        """
        import gmsh
        self._check_initialized()

        boundary = gmsh.model.getBoundary(
            [(3, vol_tag)], oriented=False)
        surfs = sorted(set(
            abs(dt[1]) for dt in boundary if dt[0] == 2))

        if self._verbose:
            print(f"[GmshBuilder] volume_surfaces vol={vol_tag} "
                  f"-> {len(surfs)} surfaces")
        return surfs

    def volume_curves(self, vol_tag):
        """Get curve tags bounding a volume.

        Retrieves all surfaces of the volume, then collects curves
        from each surface boundary.

        Parameters
        ----------
        vol_tag : int
            GMSH volume tag.

        Returns
        -------
        list of int
            Sorted list of unique curve tags.
        """
        import gmsh
        self._check_initialized()

        surf_boundary = gmsh.model.getBoundary(
            [(3, vol_tag)], combined=False, oriented=False, recursive=False)
        surf_tags = [abs(dt[1]) for dt in surf_boundary if dt[0] == 2]

        curves = set()
        for st in surf_tags:
            curve_boundary = gmsh.model.getBoundary(
                [(2, st)], oriented=False)
            for dt in curve_boundary:
                if dt[0] == 1:
                    curves.add(abs(dt[1]))

        result = sorted(curves)
        if self._verbose:
            print(f"[GmshBuilder] volume_curves vol={vol_tag} "
                  f"-> {len(result)} curves")
        return result

    def volume_vertices(self, vol_tag):
        """Get vertex (point) tags of a volume.

        Recursively traverses the boundary hierarchy down to
        dimension 0 (points).

        Parameters
        ----------
        vol_tag : int
            GMSH volume tag.

        Returns
        -------
        list of int
            Sorted list of unique point tags.
        """
        import gmsh
        self._check_initialized()

        boundary = gmsh.model.getBoundary(
            [(3, vol_tag)], oriented=False, recursive=True)
        points = sorted(set(
            abs(dt[1]) for dt in boundary if dt[0] == 0))

        if self._verbose:
            print(f"[GmshBuilder] volume_vertices vol={vol_tag} "
                  f"-> {len(points)} vertices")
        return points

    def volume_bounding_box(self, vol_tag):
        """Get axis-aligned bounding box of a volume.

        Parameters
        ----------
        vol_tag : int
            GMSH volume tag.

        Returns
        -------
        dict
            Keys: xmin, ymin, zmin, xmax, ymax, zmax.
        """
        import gmsh
        self._check_initialized()

        bb = gmsh.model.occ.getBoundingBox(3, vol_tag)

        if self._verbose:
            print(f"[GmshBuilder] volume_bounding_box vol={vol_tag}")
        return {
            'xmin': bb[0], 'ymin': bb[1], 'zmin': bb[2],
            'xmax': bb[3], 'ymax': bb[4], 'zmax': bb[5],
        }

    def volume_center(self, vol_tag):
        """Get center of mass of a volume.

        Parameters
        ----------
        vol_tag : int
            GMSH volume tag.

        Returns
        -------
        list of float
            [x, y, z] center of mass coordinates.
        """
        import gmsh
        self._check_initialized()

        x, y, z = gmsh.model.occ.getCenterOfMass(3, vol_tag)

        if self._verbose:
            print(f"[GmshBuilder] volume_center vol={vol_tag} "
                  f"-> [{x:.6g}, {y:.6g}, {z:.6g}]")
        return [x, y, z]

    def volume_measure(self, vol_tag):
        """Get volume (3D measure) of a volume entity.

        Parameters
        ----------
        vol_tag : int
            GMSH volume tag.

        Returns
        -------
        float
            Volume in cubic length units.
        """
        import gmsh
        self._check_initialized()

        vol = gmsh.model.occ.getMass(3, vol_tag)

        if self._verbose:
            print(f"[GmshBuilder] volume_measure vol={vol_tag} -> {vol:.6g}")
        return vol

    def volume_inertia(self, vol_tag):
        """Get matrix of inertia of a volume (uniform unit density).

        Parameters
        ----------
        vol_tag : int
            GMSH volume tag.

        Returns
        -------
        list of float
            9 elements [Ixx, Ixy, Ixz, Iyx, Iyy, Iyz, Izx, Izy, Izz]
            in row-major order.
        """
        import gmsh
        self._check_initialized()

        mat = gmsh.model.occ.getMatrixOfInertia(3, vol_tag)

        if self._verbose:
            print(f"[GmshBuilder] volume_inertia vol={vol_tag}")
        return list(mat)

    def volume_type(self, vol_tag):
        """Get geometric type of a volume entity.

        Parameters
        ----------
        vol_tag : int
            GMSH volume tag.

        Returns
        -------
        str
            Type name (e.g., 'Solid').
        """
        import gmsh
        self._check_initialized()

        vtype = gmsh.model.getType(3, vol_tag)

        if self._verbose:
            print(f"[GmshBuilder] volume_type vol={vol_tag} -> '{vtype}'")
        return vtype

    def volume_adjacent_volumes(self, vol_tag):
        """Get volume tags sharing a surface with the given volume.

        For each bounding surface of the volume, finds other volumes
        that share that surface via upward adjacency.

        Parameters
        ----------
        vol_tag : int
            GMSH volume tag.

        Returns
        -------
        list of int
            Sorted list of adjacent volume tags (excludes vol_tag itself).
        """
        import gmsh
        self._check_initialized()

        surfs = self.volume_surfaces(vol_tag)
        adjacent = set()
        for st in surfs:
            upward, _downward = gmsh.model.getAdjacencies(2, st)
            for adj_vol in upward:
                if adj_vol != vol_tag:
                    adjacent.add(adj_vol)

        result = sorted(adjacent)
        if self._verbose:
            print(f"[GmshBuilder] volume_adjacent_volumes vol={vol_tag} "
                  f"-> {result}")
        return result

    def volume_physical_groups(self, vol_tag):
        """Get physical group tags that contain the given volume.

        Parameters
        ----------
        vol_tag : int
            GMSH volume tag.

        Returns
        -------
        list of int
            Physical group tags for dimension 3 containing this volume.
        """
        import gmsh
        self._check_initialized()

        pg_tags = list(gmsh.model.getPhysicalGroupsForEntity(3, vol_tag))

        if self._verbose:
            print(f"[GmshBuilder] volume_physical_groups vol={vol_tag} "
                  f"-> {pg_tags}")
        return pg_tags

    def volume_is_inside(self, vol_tag, point):
        """Check if a point is inside a volume.

        Parameters
        ----------
        vol_tag : int
            GMSH volume tag.
        point : list of float
            [x, y, z] coordinates of the query point.

        Returns
        -------
        bool
            True if the point is inside the volume.
        """
        import gmsh
        self._check_initialized()

        result = gmsh.model.isInside(3, vol_tag, point)

        if self._verbose:
            print(f"[GmshBuilder] volume_is_inside vol={vol_tag} "
                  f"point={point} -> {result}")
        return bool(result)

    def volume_surface_area(self, vol_tag):
        """Get total surface area of a volume.

        Sums the area (2D measure) of each bounding surface.

        Parameters
        ----------
        vol_tag : int
            GMSH volume tag.

        Returns
        -------
        float
            Total surface area in square length units.
        """
        import gmsh
        self._check_initialized()

        surfs = self.volume_surfaces(vol_tag)
        total = 0.0
        for st in surfs:
            total += gmsh.model.occ.getMass(2, st)

        if self._verbose:
            print(f"[GmshBuilder] volume_surface_area vol={vol_tag} "
                  f"-> {total:.6g}")
        return total

    def volume_curve_count(self, vol_tag):
        """Count unique curves bounding a volume.

        Parameters
        ----------
        vol_tag : int
            GMSH volume tag.

        Returns
        -------
        int
            Number of unique curves.
        """
        self._check_initialized()

        curves = self.volume_curves(vol_tag)

        if self._verbose:
            print(f"[GmshBuilder] volume_curve_count vol={vol_tag} "
                  f"-> {len(curves)}")
        return len(curves)

    def volume_vertex_count(self, vol_tag):
        """Count unique vertices bounding a volume.

        Parameters
        ----------
        vol_tag : int
            GMSH volume tag.

        Returns
        -------
        int
            Number of unique vertices.
        """
        self._check_initialized()

        vertices = self.volume_vertices(vol_tag)

        if self._verbose:
            print(f"[GmshBuilder] volume_vertex_count vol={vol_tag} "
                  f"-> {len(vertices)}")
        return len(vertices)

    # ------------------------------------------------------------------
    # Surface queries (14 methods)
    # ------------------------------------------------------------------

    def surface_curves(self, surf_tag):
        """Get curve tags bounding a surface.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.

        Returns
        -------
        list of int
            Sorted list of curve tags.
        """
        import gmsh
        self._check_initialized()

        boundary = gmsh.model.getBoundary(
            [(2, surf_tag)], oriented=False)
        curves = sorted(set(
            abs(dt[1]) for dt in boundary if dt[0] == 1))

        if self._verbose:
            print(f"[GmshBuilder] surface_curves surf={surf_tag} "
                  f"-> {len(curves)} curves")
        return curves

    def surface_vertices(self, surf_tag):
        """Get vertex (point) tags of a surface.

        Recursively traverses boundary to dimension 0.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.

        Returns
        -------
        list of int
            Sorted list of unique point tags.
        """
        import gmsh
        self._check_initialized()

        boundary = gmsh.model.getBoundary(
            [(2, surf_tag)], oriented=False, recursive=True)
        points = sorted(set(
            abs(dt[1]) for dt in boundary if dt[0] == 0))

        if self._verbose:
            print(f"[GmshBuilder] surface_vertices surf={surf_tag} "
                  f"-> {len(points)} vertices")
        return points

    def surface_volumes(self, surf_tag):
        """Get volume tags adjacent to a surface (upward adjacency).

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.

        Returns
        -------
        list of int
            Sorted list of volume tags containing this surface.
        """
        import gmsh
        self._check_initialized()

        upward, _downward = gmsh.model.getAdjacencies(2, surf_tag)
        result = sorted(int(v) for v in upward)

        if self._verbose:
            print(f"[GmshBuilder] surface_volumes surf={surf_tag} "
                  f"-> {result}")
        return result

    def surface_bounding_box(self, surf_tag):
        """Get axis-aligned bounding box of a surface.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.

        Returns
        -------
        dict
            Keys: xmin, ymin, zmin, xmax, ymax, zmax.
        """
        import gmsh
        self._check_initialized()

        bb = gmsh.model.occ.getBoundingBox(2, surf_tag)

        if self._verbose:
            print(f"[GmshBuilder] surface_bounding_box surf={surf_tag}")
        return {
            'xmin': bb[0], 'ymin': bb[1], 'zmin': bb[2],
            'xmax': bb[3], 'ymax': bb[4], 'zmax': bb[5],
        }

    def surface_center(self, surf_tag):
        """Get center of mass of a surface.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.

        Returns
        -------
        list of float
            [x, y, z] center of mass coordinates.
        """
        import gmsh
        self._check_initialized()

        x, y, z = gmsh.model.occ.getCenterOfMass(2, surf_tag)

        if self._verbose:
            print(f"[GmshBuilder] surface_center surf={surf_tag} "
                  f"-> [{x:.6g}, {y:.6g}, {z:.6g}]")
        return [x, y, z]

    def surface_area(self, surf_tag):
        """Get area (2D measure) of a surface.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.

        Returns
        -------
        float
            Surface area in square length units.
        """
        import gmsh
        self._check_initialized()

        area = gmsh.model.occ.getMass(2, surf_tag)

        if self._verbose:
            print(f"[GmshBuilder] surface_area surf={surf_tag} -> {area:.6g}")
        return area

    def surface_type(self, surf_tag):
        """Get geometric type of a surface entity.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.

        Returns
        -------
        str
            Type name (e.g., 'Plane', 'Cylinder', 'BSpline Surface').
        """
        import gmsh
        self._check_initialized()

        stype = gmsh.model.getType(2, surf_tag)

        if self._verbose:
            print(f"[GmshBuilder] surface_type surf={surf_tag} -> '{stype}'")
        return stype

    def surface_normal_at(self, surf_tag, u=0.5, v=0.5):
        """Get outward unit normal at parametric coordinates on a surface.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.
        u : float, optional
            Parametric u coordinate. Default is 0.5.
        v : float, optional
            Parametric v coordinate. Default is 0.5.

        Returns
        -------
        list of float
            [nx, ny, nz] unit normal vector.
        """
        import gmsh
        self._check_initialized()

        # Normalize u,v from [0,1] to actual parametric bounds
        bounds = gmsh.model.getParametrizationBounds(2, surf_tag)
        u_actual = bounds[0][0] + u * (bounds[1][0] - bounds[0][0])
        v_actual = bounds[0][1] + v * (bounds[1][1] - bounds[0][1])

        normal = gmsh.model.getNormal(surf_tag, [u_actual, v_actual])

        if self._verbose:
            print(f"[GmshBuilder] surface_normal_at surf={surf_tag} "
                  f"u={u} v={v}")
        return list(normal)

    def surface_closest_point(self, surf_tag, point):
        """Find the closest point on a surface to a given point.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.
        point : list of float
            [x, y, z] query point coordinates.

        Returns
        -------
        dict
            'point': [x, y, z] closest point on surface,
            'parametric': [u, v] parametric coordinates.
        """
        import gmsh
        self._check_initialized()

        cp = gmsh.model.getClosestPoint(2, surf_tag, point)
        closest = list(cp[0])
        parametric = list(cp[1])

        # Normalize raw parametric coords to [0, 1]
        bounds = gmsh.model.getParametrizationBounds(2, surf_tag)
        u_range = bounds[1][0] - bounds[0][0]
        v_range = bounds[1][1] - bounds[0][1]
        u_norm = ((parametric[0] - bounds[0][0]) / u_range
                  if u_range != 0.0 else 0.0)
        v_norm = ((parametric[1] - bounds[0][1]) / v_range
                  if v_range != 0.0 else 0.0)

        if self._verbose:
            print(f"[GmshBuilder] surface_closest_point surf={surf_tag} "
                  f"query={point}")
        return {'point': closest, 'parametric': [u_norm, v_norm]}

    def surface_is_planar(self, surf_tag):
        """Check if a surface is planar.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.

        Returns
        -------
        bool
            True if the surface type contains 'Plane'.
        """
        self._check_initialized()

        stype = self.surface_type(surf_tag)
        return 'Plane' in stype

    def surface_is_cylindrical(self, surf_tag):
        """Check if a surface is cylindrical.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.

        Returns
        -------
        bool
            True if the surface type contains 'Cylinder'.
        """
        self._check_initialized()

        stype = self.surface_type(surf_tag)
        return 'Cylinder' in stype

    def surface_is_conical(self, surf_tag):
        """Check if a surface is conical.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.

        Returns
        -------
        bool
            True if the surface type contains 'Cone'.
        """
        self._check_initialized()

        stype = self.surface_type(surf_tag)
        return 'Cone' in stype

    def surface_is_spherical(self, surf_tag):
        """Check if a surface is spherical.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.

        Returns
        -------
        bool
            True if the surface type contains 'Sphere'.
        """
        self._check_initialized()

        stype = self.surface_type(surf_tag)
        return 'Sphere' in stype

    def surface_physical_groups(self, surf_tag):
        """Get physical group tags that contain the given surface.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.

        Returns
        -------
        list of int
            Physical group tags for dimension 2 containing this surface.
        """
        import gmsh
        self._check_initialized()

        pg_tags = list(gmsh.model.getPhysicalGroupsForEntity(2, surf_tag))

        if self._verbose:
            print(f"[GmshBuilder] surface_physical_groups surf={surf_tag} "
                  f"-> {pg_tags}")
        return pg_tags

    # ------------------------------------------------------------------
    # Curve queries (14 methods)
    # ------------------------------------------------------------------

    def curve_vertices(self, curve_tag):
        """Get vertex (point) tags at the endpoints of a curve.

        Parameters
        ----------
        curve_tag : int
            GMSH curve tag.

        Returns
        -------
        list of int
            Sorted list of point tags (typically 2 for open curves,
            0 or 1 for closed/periodic curves).
        """
        import gmsh
        self._check_initialized()

        boundary = gmsh.model.getBoundary(
            [(1, curve_tag)], oriented=False)
        points = sorted(set(
            abs(dt[1]) for dt in boundary if dt[0] == 0))

        if self._verbose:
            print(f"[GmshBuilder] curve_vertices curve={curve_tag} "
                  f"-> {points}")
        return points

    def curve_surfaces(self, curve_tag):
        """Get surface tags adjacent to a curve (upward adjacency).

        Parameters
        ----------
        curve_tag : int
            GMSH curve tag.

        Returns
        -------
        list of int
            Sorted list of surface tags containing this curve.
        """
        import gmsh
        self._check_initialized()

        upward, _downward = gmsh.model.getAdjacencies(1, curve_tag)
        result = sorted(int(s) for s in upward)

        if self._verbose:
            print(f"[GmshBuilder] curve_surfaces curve={curve_tag} "
                  f"-> {result}")
        return result

    def curve_bounding_box(self, curve_tag):
        """Get axis-aligned bounding box of a curve.

        Parameters
        ----------
        curve_tag : int
            GMSH curve tag.

        Returns
        -------
        dict
            Keys: xmin, ymin, zmin, xmax, ymax, zmax.
        """
        import gmsh
        self._check_initialized()

        bb = gmsh.model.occ.getBoundingBox(1, curve_tag)

        if self._verbose:
            print(f"[GmshBuilder] curve_bounding_box curve={curve_tag}")
        return {
            'xmin': bb[0], 'ymin': bb[1], 'zmin': bb[2],
            'xmax': bb[3], 'ymax': bb[4], 'zmax': bb[5],
        }

    def curve_center(self, curve_tag):
        """Get center of mass of a curve.

        Parameters
        ----------
        curve_tag : int
            GMSH curve tag.

        Returns
        -------
        list of float
            [x, y, z] center of mass coordinates.
        """
        import gmsh
        self._check_initialized()

        x, y, z = gmsh.model.occ.getCenterOfMass(1, curve_tag)

        if self._verbose:
            print(f"[GmshBuilder] curve_center curve={curve_tag} "
                  f"-> [{x:.6g}, {y:.6g}, {z:.6g}]")
        return [x, y, z]

    def curve_length(self, curve_tag):
        """Get length (1D measure) of a curve.

        Parameters
        ----------
        curve_tag : int
            GMSH curve tag.

        Returns
        -------
        float
            Curve length in length units.
        """
        import gmsh
        self._check_initialized()

        length = gmsh.model.occ.getMass(1, curve_tag)

        if self._verbose:
            print(f"[GmshBuilder] curve_length curve={curve_tag} "
                  f"-> {length:.6g}")
        return length

    def curve_type(self, curve_tag):
        """Get geometric type of a curve entity.

        Parameters
        ----------
        curve_tag : int
            GMSH curve tag.

        Returns
        -------
        str
            Type name (e.g., 'Line', 'Circle', 'BSpline Curve').
        """
        import gmsh
        self._check_initialized()

        ctype = gmsh.model.getType(1, curve_tag)

        if self._verbose:
            print(f"[GmshBuilder] curve_type curve={curve_tag} -> '{ctype}'")
        return ctype

    def curve_point_at(self, curve_tag, t):
        """Get 3D point on a curve at a parametric value.

        Parameters
        ----------
        curve_tag : int
            GMSH curve tag.
        t : float
            Parametric value. The interpretation depends on the
            curve parameterization (use getParametrizationBounds
            to determine the valid range).

        Returns
        -------
        list of float
            [x, y, z] point coordinates.
        """
        import gmsh
        self._check_initialized()

        # Normalize t from [0,1] to actual parametric bounds
        bounds = gmsh.model.getParametrizationBounds(1, curve_tag)
        t_actual = bounds[0][0] + t * (bounds[1][0] - bounds[0][0])

        point = gmsh.model.getValue(1, curve_tag, [t_actual])

        if self._verbose:
            print(f"[GmshBuilder] curve_point_at curve={curve_tag} t={t}")
        return list(point)

    def curve_tangent_at(self, curve_tag, t):
        """Get tangent vector on a curve at a parametric value.

        Parameters
        ----------
        curve_tag : int
            GMSH curve tag.
        t : float
            Parametric value.

        Returns
        -------
        list of float
            [tx, ty, tz] tangent vector (not necessarily unit length).
        """
        import gmsh
        self._check_initialized()

        # Normalize t from [0,1] to actual parametric bounds
        bounds = gmsh.model.getParametrizationBounds(1, curve_tag)
        t_actual = bounds[0][0] + t * (bounds[1][0] - bounds[0][0])

        deriv = gmsh.model.getDerivative(1, curve_tag, [t_actual])

        if self._verbose:
            print(f"[GmshBuilder] curve_tangent_at curve={curve_tag} t={t}")
        return list(deriv[:3])

    def curve_curvature_at(self, curve_tag, t):
        """Get curvature of a curve at a parametric value.

        Parameters
        ----------
        curve_tag : int
            GMSH curve tag.
        t : float
            Parametric value.

        Returns
        -------
        float
            Curvature value (1/radius for circular arcs).
        """
        import gmsh
        self._check_initialized()

        # Normalize t from [0,1] to actual parametric bounds
        bounds = gmsh.model.getParametrizationBounds(1, curve_tag)
        t_actual = bounds[0][0] + t * (bounds[1][0] - bounds[0][0])

        curv = gmsh.model.getCurvature(1, curve_tag, [t_actual])

        if self._verbose:
            print(f"[GmshBuilder] curve_curvature_at curve={curve_tag} "
                  f"t={t} -> {curv[0]}")
        return float(curv[0])

    def curve_closest_point(self, curve_tag, point):
        """Find the closest point on a curve to a given point.

        Parameters
        ----------
        curve_tag : int
            GMSH curve tag.
        point : list of float
            [x, y, z] query point coordinates.

        Returns
        -------
        dict
            'point': [x, y, z] closest point on curve,
            'parametric': [t] parametric coordinate.
        """
        import gmsh
        self._check_initialized()

        cp = gmsh.model.getClosestPoint(1, curve_tag, point)
        closest = list(cp[0])
        parametric = list(cp[1])

        # Normalize raw parametric coord to [0, 1]
        bounds = gmsh.model.getParametrizationBounds(1, curve_tag)
        t_range = bounds[1][0] - bounds[0][0]
        t_norm = ((parametric[0] - bounds[0][0]) / t_range
                  if t_range != 0.0 else 0.0)

        if self._verbose:
            print(f"[GmshBuilder] curve_closest_point curve={curve_tag} "
                  f"query={point}")
        return {'point': closest, 'parametric': [t_norm]}

    def curve_is_line(self, curve_tag):
        """Check if a curve is a straight line.

        Parameters
        ----------
        curve_tag : int
            GMSH curve tag.

        Returns
        -------
        bool
            True if the curve type contains 'Line'.
        """
        self._check_initialized()

        ctype = self.curve_type(curve_tag)
        return 'Line' in ctype

    def curve_is_circle(self, curve_tag):
        """Check if a curve is a circular arc.

        Parameters
        ----------
        curve_tag : int
            GMSH curve tag.

        Returns
        -------
        bool
            True if the curve type contains 'Circle'.
        """
        self._check_initialized()

        ctype = self.curve_type(curve_tag)
        return 'Circle' in ctype

    def curve_is_bspline(self, curve_tag):
        """Check if a curve is a B-spline.

        Parameters
        ----------
        curve_tag : int
            GMSH curve tag.

        Returns
        -------
        bool
            True if the curve type contains 'BSpline'.
        """
        self._check_initialized()

        ctype = self.curve_type(curve_tag)
        return 'BSpline' in ctype

    def curve_physical_groups(self, curve_tag):
        """Get physical group tags that contain the given curve.

        Parameters
        ----------
        curve_tag : int
            GMSH curve tag.

        Returns
        -------
        list of int
            Physical group tags for dimension 1 containing this curve.
        """
        import gmsh
        self._check_initialized()

        pg_tags = list(gmsh.model.getPhysicalGroupsForEntity(1, curve_tag))

        if self._verbose:
            print(f"[GmshBuilder] curve_physical_groups curve={curve_tag} "
                  f"-> {pg_tags}")
        return pg_tags

    # ------------------------------------------------------------------
    # Vertex queries (13 methods)
    # ------------------------------------------------------------------

    def vertex_coordinates(self, point_tag):
        """Get 3D coordinates of a vertex (point entity).

        Parameters
        ----------
        point_tag : int
            GMSH point tag.

        Returns
        -------
        list of float
            [x, y, z] coordinates.
        """
        import gmsh
        self._check_initialized()

        coords = gmsh.model.getValue(0, point_tag, [])

        if self._verbose:
            print(f"[GmshBuilder] vertex_coordinates point={point_tag} "
                  f"-> [{coords[0]:.6g}, {coords[1]:.6g}, {coords[2]:.6g}]")
        return list(coords)

    def vertex_curves(self, point_tag):
        """Get curve tags adjacent to a vertex (upward adjacency).

        Parameters
        ----------
        point_tag : int
            GMSH point tag.

        Returns
        -------
        list of int
            Sorted list of curve tags that have this vertex as an
            endpoint.
        """
        import gmsh
        self._check_initialized()

        upward, _downward = gmsh.model.getAdjacencies(0, point_tag)
        result = sorted(int(c) for c in upward)

        if self._verbose:
            print(f"[GmshBuilder] vertex_curves point={point_tag} "
                  f"-> {result}")
        return result

    def vertex_surfaces(self, point_tag):
        """Get surface tags reachable from a vertex via upward traversal.

        Traverses from vertex to curves, then from curves to surfaces.

        Parameters
        ----------
        point_tag : int
            GMSH point tag.

        Returns
        -------
        list of int
            Sorted list of unique surface tags.
        """
        import gmsh
        self._check_initialized()

        curves = self.vertex_curves(point_tag)
        surfaces = set()
        for ct in curves:
            upward, _downward = gmsh.model.getAdjacencies(1, ct)
            for st in upward:
                surfaces.add(int(st))

        result = sorted(surfaces)
        if self._verbose:
            print(f"[GmshBuilder] vertex_surfaces point={point_tag} "
                  f"-> {result}")
        return result

    def vertex_volumes(self, point_tag):
        """Get volume tags reachable from a vertex via upward traversal.

        Traverses from vertex to curves, curves to surfaces, then
        surfaces to volumes.

        Parameters
        ----------
        point_tag : int
            GMSH point tag.

        Returns
        -------
        list of int
            Sorted list of unique volume tags.
        """
        import gmsh
        self._check_initialized()

        surfaces = self.vertex_surfaces(point_tag)
        volumes = set()
        for st in surfaces:
            upward, _downward = gmsh.model.getAdjacencies(2, st)
            for vt in upward:
                volumes.add(int(vt))

        result = sorted(volumes)
        if self._verbose:
            print(f"[GmshBuilder] vertex_volumes point={point_tag} "
                  f"-> {result}")
        return result

    def vertex_physical_groups(self, point_tag):
        """Get physical group tags that contain the given vertex.

        Parameters
        ----------
        point_tag : int
            GMSH point tag.

        Returns
        -------
        list of int
            Physical group tags for dimension 0 containing this vertex.
        """
        import gmsh
        self._check_initialized()

        pg_tags = list(gmsh.model.getPhysicalGroupsForEntity(0, point_tag))

        if self._verbose:
            print(f"[GmshBuilder] vertex_physical_groups point={point_tag} "
                  f"-> {pg_tags}")
        return pg_tags

    def vertex_distance_to(self, point_tag, point):
        """Compute Euclidean distance from a vertex to a given point.

        Parameters
        ----------
        point_tag : int
            GMSH point tag.
        point : list of float
            [x, y, z] target point coordinates.

        Returns
        -------
        float
            Euclidean distance.
        """
        self._check_initialized()

        coords = self.vertex_coordinates(point_tag)
        dist = math.sqrt(sum(
            (a - b) ** 2 for a, b in zip(coords, point)))

        if self._verbose:
            print(f"[GmshBuilder] vertex_distance_to point={point_tag} "
                  f"target={point} -> {dist:.6g}")
        return dist

    def vertex_on_curve(self, point_tag, curve_tag):
        """Check if a vertex is an endpoint of a curve.

        Parameters
        ----------
        point_tag : int
            GMSH point tag.
        curve_tag : int
            GMSH curve tag.

        Returns
        -------
        bool
            True if the vertex is in the boundary of the curve.
        """
        import gmsh
        self._check_initialized()

        boundary = gmsh.model.getBoundary(
            [(1, curve_tag)], oriented=False)
        boundary_points = set(
            abs(dt[1]) for dt in boundary if dt[0] == 0)

        result = point_tag in boundary_points
        if self._verbose:
            print(f"[GmshBuilder] vertex_on_curve point={point_tag} "
                  f"curve={curve_tag} -> {result}")
        return result

    def vertex_on_surface(self, point_tag, surf_tag):
        """Check if a vertex is on the boundary of a surface.

        Recursively checks if the point tag appears in the
        surface boundary hierarchy.

        Parameters
        ----------
        point_tag : int
            GMSH point tag.
        surf_tag : int
            GMSH surface tag.

        Returns
        -------
        bool
            True if the vertex is in the recursive boundary of
            the surface.
        """
        import gmsh
        self._check_initialized()

        boundary = gmsh.model.getBoundary(
            [(2, surf_tag)], oriented=False, recursive=True)
        boundary_points = set(
            abs(dt[1]) for dt in boundary if dt[0] == 0)

        result = point_tag in boundary_points
        if self._verbose:
            print(f"[GmshBuilder] vertex_on_surface point={point_tag} "
                  f"surf={surf_tag} -> {result}")
        return result

    def vertex_on_volume(self, point_tag, vol_tag):
        """Check if a vertex is on the boundary of a volume.

        Recursively checks if the point tag appears in the
        volume boundary hierarchy.

        Parameters
        ----------
        point_tag : int
            GMSH point tag.
        vol_tag : int
            GMSH volume tag.

        Returns
        -------
        bool
            True if the vertex is in the recursive boundary of
            the volume.
        """
        import gmsh
        self._check_initialized()

        boundary = gmsh.model.getBoundary(
            [(3, vol_tag)], oriented=False, recursive=True)
        boundary_points = set(
            abs(dt[1]) for dt in boundary if dt[0] == 0)

        result = point_tag in boundary_points
        if self._verbose:
            print(f"[GmshBuilder] vertex_on_volume point={point_tag} "
                  f"vol={vol_tag} -> {result}")
        return result

    def vertex_bounding_box(self, point_tag):
        """Get bounding box of a vertex (degenerate, point = box).

        Parameters
        ----------
        point_tag : int
            GMSH point tag.

        Returns
        -------
        dict
            Keys: xmin, ymin, zmin, xmax, ymax, zmax (all equal
            to the vertex coordinates for a point entity).
        """
        import gmsh
        self._check_initialized()

        bb = gmsh.model.occ.getBoundingBox(0, point_tag)

        if self._verbose:
            print(f"[GmshBuilder] vertex_bounding_box point={point_tag}")
        return {
            'xmin': bb[0], 'ymin': bb[1], 'zmin': bb[2],
            'xmax': bb[3], 'ymax': bb[4], 'zmax': bb[5],
        }

    def get_vertex_count(self):
        """Get the total number of point entities in the model.

        Returns
        -------
        int
            Number of point (dim=0) entities.
        """
        import gmsh
        self._check_initialized()

        count = len(gmsh.model.getEntities(0))

        if self._verbose:
            print(f"[GmshBuilder] get_vertex_count -> {count}")
        return count

    def get_all_vertices(self):
        """Get all point entity tags in the model.

        Returns
        -------
        list of int
            Sorted list of all point tags.
        """
        import gmsh
        self._check_initialized()

        entities = gmsh.model.getEntities(0)
        result = sorted(e[1] for e in entities)

        if self._verbose:
            print(f"[GmshBuilder] get_all_vertices -> {len(result)} vertices")
        return result

    def vertices_in_bounding_box(self, xmin, ymin, zmin, xmax, ymax, zmax):
        """Get point tags within an axis-aligned bounding box.

        Parameters
        ----------
        xmin : float
            Minimum x coordinate.
        ymin : float
            Minimum y coordinate.
        zmin : float
            Minimum z coordinate.
        xmax : float
            Maximum x coordinate.
        ymax : float
            Maximum y coordinate.
        zmax : float
            Maximum z coordinate.

        Returns
        -------
        list of int
            Sorted list of point tags inside the bounding box.
        """
        import gmsh
        self._check_initialized()

        entities = gmsh.model.getEntitiesInBoundingBox(
            xmin, ymin, zmin, xmax, ymax, zmax, dim=0)
        result = sorted(e[1] for e in entities)

        if self._verbose:
            print(f"[GmshBuilder] vertices_in_bounding_box "
                  f"[{xmin},{ymin},{zmin}]-[{xmax},{ymax},{zmax}] "
                  f"-> {len(result)} vertices")
        return result
