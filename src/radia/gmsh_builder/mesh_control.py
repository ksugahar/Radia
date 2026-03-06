"""
Mesh control methods for GmshBuilder.

Gmsh equivalents: size, interval, scheme, bias, boundary layer, sizing fields.
"""


class MeshControlMixin:
    """Mixin providing mesh size/scheme/field control."""

    def set_mesh_size(self, size):
        """Set global mesh element size (Gmsh 'size auto').

        Parameters
        ----------
        size : float
            Target element size in meters.
        """
        self._mesh_size = size
        if self._verbose:
            print(f"[GmshBuilder] set_mesh_size = {size}")

    def set_mesh_size_at(self, vol_id, size):
        """Set mesh size at a specific volume's vertices.

        Parameters
        ----------
        vol_id : int
            Volume ID.
        size : float
            Target element size.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        for vt in self._volumes[vol_id]:
            bounds = gmsh.model.getBoundary([(3, vt)], recursive=True)
            points = [dt for dt in bounds if dt[0] == 0]
            if points:
                gmsh.model.mesh.setSize(points, size)

        if self._verbose:
            print(f"[GmshBuilder] set_mesh_size_at vol={vol_id} size={size}")

    def set_divisions(self, vol_id, nx, ny, nz):
        """Set transfinite divisions for a volume (Gmsh 'interval').

        Parameters
        ----------
        vol_id : int
            Volume ID.
        nx, ny, nz : int
            Divisions along x, y, z.
        """
        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        for tag in self._volumes[vol_id]:
            self._divisions[tag] = (nx, ny, nz)

        if self._verbose:
            print(f"[GmshBuilder] set_divisions vol={vol_id} "
                  f"nx={nx} ny={ny} nz={nz}")

    def set_divisions_all(self, nx, ny, nz):
        """Set transfinite divisions for all volumes.

        Parameters
        ----------
        nx, ny, nz : int
            Divisions along x, y, z.
        """
        self._global_divisions = (nx, ny, nz)
        if self._verbose:
            print(f"[GmshBuilder] set_divisions_all nx={nx} ny={ny} nz={nz}")

    def set_order(self, order):
        """Set mesh element order (Gmsh 'block element_type').

        Parameters
        ----------
        order : int
            1 (linear: Hex8/Tet4) or 2 (quadratic: Hex20/Tet10).
        """
        import gmsh
        self._check_initialized()

        if not self._meshed:
            raise RuntimeError("generate() must be called before set_order()")

        gmsh.model.mesh.setOrder(order)
        if self._verbose:
            print(f"[GmshBuilder] set_order({order})")

    def set_scheme(self, vol_id, scheme):
        """Set meshing scheme for a volume (Gmsh 'scheme').

        Parameters
        ----------
        vol_id : int
            Volume ID.
        scheme : str
            'map' (transfinite hex), 'sweep', 'tet', or 'auto'.
        """
        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        # Store scheme for use during generate
        if not hasattr(self, '_schemes'):
            self._schemes = {}
        self._schemes[vol_id] = scheme

        if self._verbose:
            print(f"[GmshBuilder] set_scheme vol={vol_id} scheme='{scheme}'")

    def set_algorithm(self, algo_2d=6, algo_3d=1):
        """Set mesh algorithm (algorithm selection).

        Parameters
        ----------
        algo_2d : int
            2D algorithm: 1=MeshAdapt, 5=Delaunay, 6=Frontal-Delaunay,
            8=Frontal-Delaunay-Quad, 9=Packing-Parallelograms.
        algo_3d : int
            3D algorithm: 1=Delaunay, 4=Frontal, 10=HXT.
        """
        import gmsh
        self._check_initialized()

        gmsh.option.setNumber('Mesh.Algorithm', algo_2d)
        gmsh.option.setNumber('Mesh.Algorithm3D', algo_3d)

        if self._verbose:
            print(f"[GmshBuilder] set_algorithm 2D={algo_2d} 3D={algo_3d}")

    def set_curve_divisions(self, curve_tag, n, mesh_type="Progression",
                             coef=1.0):
        """Set transfinite curve divisions (Gmsh 'curve interval').

        Parameters
        ----------
        curve_tag : int
            GMSH curve tag.
        n : int
            Number of divisions.
        mesh_type : str
            'Progression', 'Bump', 'Beta'.
        coef : float
            Grading coefficient.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.setTransfiniteCurve(curve_tag, n + 1,
                                              mesh_type, coef)
        if self._verbose:
            print(f"[GmshBuilder] set_curve_divisions curve={curve_tag} "
                  f"n={n} type={mesh_type} coef={coef}")

    def set_bias(self, curve_tag, n, ratio):
        """Set biased curve divisions (Gmsh 'curve bias').

        Parameters
        ----------
        curve_tag : int
            GMSH curve tag.
        n : int
            Number of divisions.
        ratio : float
            Size ratio (first/last element).
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.setTransfiniteCurve(
            curve_tag, n + 1, "Progression", ratio)

        if self._verbose:
            print(f"[GmshBuilder] set_bias curve={curve_tag} n={n} ratio={ratio}")

    def set_double_bias(self, curve_tag, n, ratio):
        """Set double-biased curve (Gmsh 'double bias'): small at both ends.

        Parameters
        ----------
        curve_tag : int
            GMSH curve tag.
        n : int
            Number of divisions.
        ratio : float
            Bump ratio.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.setTransfiniteCurve(
            curve_tag, n + 1, "Bump", ratio)

        if self._verbose:
            print(f"[GmshBuilder] set_double_bias curve={curve_tag} "
                  f"n={n} ratio={ratio}")

    def set_boundary_layer(self, surf_tags, thickness, n_layers, ratio=1.2):
        """Set boundary layer mesh (Gmsh 'boundary layer').

        Uses GMSH BoundaryLayer field.

        Parameters
        ----------
        surf_tags : list of int
            GMSH surface tags for BL.
        n_layers : int
            Number of layers.
        thickness : float
            Total BL thickness.
        ratio : float
            Growth ratio between layers.
        """
        import gmsh
        self._check_initialized()

        field = gmsh.model.mesh.field.add("BoundaryLayer")
        # BoundaryLayer operates on curves, not surfaces directly
        # Extract boundary curves from the given surfaces
        curve_tags = []
        for st in surf_tags:
            boundary = gmsh.model.getBoundary([(2, st)], oriented=False)
            curve_tags.extend([abs(dt[1]) for dt in boundary if dt[0] == 1])
        curve_tags = list(set(curve_tags))  # deduplicate
        gmsh.model.mesh.field.setNumbers(field, "CurvesList", curve_tags)
        gmsh.model.mesh.field.setNumber(field, "Thickness", thickness)
        gmsh.model.mesh.field.setNumber(field, "NbLayers", n_layers)
        gmsh.model.mesh.field.setNumber(field, "Ratio", ratio)
        gmsh.model.mesh.field.setAsBoundaryLayer(field)

        if self._verbose:
            print(f"[GmshBuilder] set_boundary_layer n={n_layers} "
                  f"t={thickness} ratio={ratio}")

    def add_field_distance(self, entity_tags, dim=0):
        """Add a Distance sizing field.

        Parameters
        ----------
        entity_tags : list of int
            Entity tags to measure distance from.
        dim : int
            Entity dimension (0=point, 1=curve, 2=surface).

        Returns
        -------
        int
            Field tag.
        """
        import gmsh
        self._check_initialized()

        field = gmsh.model.mesh.field.add("Distance")
        key = {0: "PointsList", 1: "CurvesList", 2: "SurfacesList"}[dim]
        gmsh.model.mesh.field.setNumbers(field, key, entity_tags)
        return field

    def add_field_threshold(self, distance_field, dist_min, dist_max,
                              size_min, size_max):
        """Add a Threshold sizing field.

        Parameters
        ----------
        distance_field : int
            Distance field tag.
        dist_min, dist_max : float
            Distance range.
        size_min, size_max : float
            Element size range.

        Returns
        -------
        int
            Field tag.
        """
        import gmsh
        self._check_initialized()

        field = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(field, "InField", distance_field)
        gmsh.model.mesh.field.setNumber(field, "DistMin", dist_min)
        gmsh.model.mesh.field.setNumber(field, "DistMax", dist_max)
        gmsh.model.mesh.field.setNumber(field, "SizeMin", size_min)
        gmsh.model.mesh.field.setNumber(field, "SizeMax", size_max)
        return field

    def add_field_box(self, xmin, xmax, ymin, ymax, zmin, zmax,
                       size_in, size_out):
        """Add a Box sizing field.

        Parameters
        ----------
        xmin, xmax, ymin, ymax, zmin, zmax : float
            Box extents.
        size_in : float
            Element size inside box.
        size_out : float
            Element size outside box.

        Returns
        -------
        int
            Field tag.
        """
        import gmsh
        self._check_initialized()

        field = gmsh.model.mesh.field.add("Box")
        gmsh.model.mesh.field.setNumber(field, "XMin", xmin)
        gmsh.model.mesh.field.setNumber(field, "XMax", xmax)
        gmsh.model.mesh.field.setNumber(field, "YMin", ymin)
        gmsh.model.mesh.field.setNumber(field, "YMax", ymax)
        gmsh.model.mesh.field.setNumber(field, "ZMin", zmin)
        gmsh.model.mesh.field.setNumber(field, "ZMax", zmax)
        gmsh.model.mesh.field.setNumber(field, "VIn", size_in)
        gmsh.model.mesh.field.setNumber(field, "VOut", size_out)
        return field

    def add_field_math(self, expression):
        """Add a MathEval sizing field.

        Parameters
        ----------
        expression : str
            Mathematical expression using x, y, z variables.
            Example: "0.001 + 0.01*sqrt(x*x+y*y)"

        Returns
        -------
        int
            Field tag.
        """
        import gmsh
        self._check_initialized()

        field = gmsh.model.mesh.field.add("MathEval")
        gmsh.model.mesh.field.setString(field, "F", expression)
        return field

    def add_field_min(self, field_tags):
        """Add a Min field (takes minimum of other fields).

        Parameters
        ----------
        field_tags : list of int
            Field tags to take minimum of.

        Returns
        -------
        int
            Field tag.
        """
        import gmsh
        self._check_initialized()

        field = gmsh.model.mesh.field.add("Min")
        gmsh.model.mesh.field.setNumbers(field, "FieldsList", field_tags)
        return field

    def set_background_field(self, field_tag):
        """Set a field as the background mesh size field.

        Parameters
        ----------
        field_tag : int
            Field tag to use as background.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.field.setAsBackgroundMesh(field_tag)

    def set_curvature_mesh(self, n_per_2pi):
        """Set curvature-based mesh sizing.

        Parameters
        ----------
        n_per_2pi : int
            Number of elements per 2*pi curvature.
        """
        import gmsh
        self._check_initialized()

        gmsh.option.setNumber('Mesh.MeshSizeFromCurvature', n_per_2pi)
        if self._verbose:
            print(f"[GmshBuilder] set_curvature_mesh n_per_2pi={n_per_2pi}")

    def set_recombine(self, surf_tag, angle=45.0):
        """Set recombine on a surface (tri -> quad).

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.
        angle : float, optional
            Maximum deviation angle (in degrees) from a right angle
            for element recombination.  Default is 45.0.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.setRecombine(2, surf_tag, angle)

    def set_recombine_all(self):
        """Set recombine on all surfaces."""
        import gmsh
        self._check_initialized()

        gmsh.option.setNumber('Mesh.RecombineAll', 1)
        if self._verbose:
            print("[GmshBuilder] set_recombine_all")

    def set_recombination_algorithm(self, algo):
        """Set the recombination algorithm for quad/hex meshing.

        Parameters
        ----------
        algo : int or str
            Recombination algorithm.  Accepts integer codes or string
            names:

            - 0 or 'simple'  : Simple (default)
            - 1 or 'blossom' : Blossom (best quality, slower)
            - 2 or 'full_quad' : Simple full-quad
            - 3 or 'blossom_full_quad' : Blossom full-quad
        """
        import gmsh
        self._check_initialized()

        algo_map = {
            'simple': 0, 'blossom': 1,
            'full_quad': 2, 'blossom_full_quad': 3,
        }
        if isinstance(algo, str):
            algo = algo_map.get(algo.lower())
            if algo is None:
                raise ValueError(
                    f"Unknown recombination algorithm. "
                    f"Use one of: {list(algo_map.keys())} or 0-3")

        gmsh.option.setNumber('Mesh.RecombinationAlgorithm', algo)
        if self._verbose:
            print(f"[GmshBuilder] set_recombination_algorithm = {algo}")

    def set_smoothing(self, n_steps):
        """Set number of smoothing steps.

        Parameters
        ----------
        n_steps : int
            Number of smoothing iterations.
        """
        import gmsh
        self._check_initialized()

        gmsh.option.setNumber('Mesh.Smoothing', n_steps)

    def set_optimization(self, method='Netgen'):
        """Set mesh optimization method.

        Parameters
        ----------
        method : str
            'Netgen', 'Laplace2D', 'HighOrder', 'HighOrderElastic',
            'HighOrderFastCurving', 'Relocate3D', 'UntangleMeshGeometry'.
        """
        import gmsh
        self._check_initialized()

        if not hasattr(self, '_optimization_method'):
            self._optimization_method = None
        self._optimization_method = method
        # Enable optimization
        gmsh.option.setNumber('Mesh.Optimize', 1)
        if method.lower() == 'netgen':
            gmsh.option.setNumber('Mesh.OptimizeNetgen', 1)
        else:
            gmsh.option.setNumber('Mesh.OptimizeNetgen', 0)
        if self._verbose:
            print(f"[GmshBuilder] set_optimization: method={method}")

    def set_min_max_size(self, min_size, max_size):
        """Set global min/max mesh element size.

        Parameters
        ----------
        min_size : float
            Minimum element size.
        max_size : float
            Maximum element size.
        """
        import gmsh
        self._check_initialized()

        gmsh.option.setNumber('Mesh.CharacteristicLengthMin', min_size)
        gmsh.option.setNumber('Mesh.CharacteristicLengthMax', max_size)

    def set_growth_rate(self, rate):
        """Set mesh growth rate via curvature-based sizing.

        Enables curvature-based mesh sizing and sets the
        ``Mesh.MeshSizeFromCurvature`` option, which controls the
        number of mesh elements per 2*pi radians of curvature.
        A smaller ``rate`` produces more elements (finer mesh) near
        curved surfaces.

        Parameters
        ----------
        rate : float
            Growth factor (e.g., 1.3 = 30% growth between neighbors).
            Internally mapped to ``MeshSizeFromCurvature`` as
            ``max(1, int(2*pi / rate))``.
        """
        import gmsh
        self._check_initialized()

        gmsh.option.setNumber('Mesh.CharacteristicLengthExtendFromBoundary', 1)
        gmsh.option.setNumber('Mesh.CharacteristicLengthFromCurvature', 1)
        elements_per_2pi = max(1, int(2 * 3.14159265 / rate))
        gmsh.option.setNumber('Mesh.MeshSizeFromCurvature', elements_per_2pi)
        self._growth_rate = rate
        if self._verbose:
            print(f"[GmshBuilder] set_growth_rate: rate={rate} "
                  f"(MeshSizeFromCurvature={elements_per_2pi})")

    def set_angle_threshold(self, angle_deg):
        """Set feature angle threshold.

        Parameters
        ----------
        angle_deg : float
            Angle in degrees below which edges are considered features.
        """
        import gmsh
        self._check_initialized()

        gmsh.option.setNumber('Mesh.AngleToleranceFacetOverlap', angle_deg)

    # ------------------------------------------------------------------
    # Additional sizing fields
    # ------------------------------------------------------------------

    def add_field_cylinder(self, xc, yc, zc, xa, ya, za, radius,
                           size_in, size_out):
        """Add a Cylinder sizing field.

        Parameters
        ----------
        xc, yc, zc : float
            Center coordinates of the cylinder.
        xa, ya, za : float
            Axis direction vector of the cylinder.
        radius : float
            Cylinder radius.
        size_in : float
            Element size inside the cylinder.
        size_out : float
            Element size outside the cylinder.

        Returns
        -------
        int
            Field tag.
        """
        import gmsh
        self._check_initialized()

        field = gmsh.model.mesh.field.add("Cylinder")
        gmsh.model.mesh.field.setNumber(field, "XCenter", xc)
        gmsh.model.mesh.field.setNumber(field, "YCenter", yc)
        gmsh.model.mesh.field.setNumber(field, "ZCenter", zc)
        gmsh.model.mesh.field.setNumber(field, "XAxis", xa)
        gmsh.model.mesh.field.setNumber(field, "YAxis", ya)
        gmsh.model.mesh.field.setNumber(field, "ZAxis", za)
        gmsh.model.mesh.field.setNumber(field, "Radius", radius)
        gmsh.model.mesh.field.setNumber(field, "VIn", size_in)
        gmsh.model.mesh.field.setNumber(field, "VOut", size_out)
        return field

    def add_field_ball(self, xc, yc, zc, radius, size_in, size_out):
        """Add a Ball sizing field.

        Parameters
        ----------
        xc, yc, zc : float
            Center coordinates of the ball.
        radius : float
            Ball radius.
        size_in : float
            Element size inside the ball.
        size_out : float
            Element size outside the ball.

        Returns
        -------
        int
            Field tag.
        """
        import gmsh
        self._check_initialized()

        field = gmsh.model.mesh.field.add("Ball")
        gmsh.model.mesh.field.setNumber(field, "XCenter", xc)
        gmsh.model.mesh.field.setNumber(field, "YCenter", yc)
        gmsh.model.mesh.field.setNumber(field, "ZCenter", zc)
        gmsh.model.mesh.field.setNumber(field, "Radius", radius)
        gmsh.model.mesh.field.setNumber(field, "VIn", size_in)
        gmsh.model.mesh.field.setNumber(field, "VOut", size_out)
        return field

    def add_field_max(self, field_tags):
        """Add a Max field (takes maximum of other fields).

        Parameters
        ----------
        field_tags : list of int
            Field tags to take maximum of.

        Returns
        -------
        int
            Field tag.
        """
        import gmsh
        self._check_initialized()

        field = gmsh.model.mesh.field.add("Max")
        gmsh.model.mesh.field.setNumbers(field, "FieldsList", field_tags)
        return field

    def add_field_mean(self, field_tag, delta=1.0):
        """Add a Mean field that smooths another field.

        The Mean field computes a moving average of the referenced field
        within a ball of radius Delta around each point.

        Parameters
        ----------
        field_tag : int
            Tag of the input field to smooth.
        delta : float
            Radius of the averaging ball (default 1.0).

        Returns
        -------
        int
            Field tag.
        """
        import gmsh
        self._check_initialized()

        field = gmsh.model.mesh.field.add("Mean")
        gmsh.model.mesh.field.setNumber(field, "InField", field_tag)
        gmsh.model.mesh.field.setNumber(field, "Delta", delta)
        return field

    def add_field_restrict(self, field_tag, vol_ids):
        """Restrict a sizing field to specific volumes.

        Parameters
        ----------
        field_tag : int
            Existing field tag to restrict.
        vol_ids : list of int
            Volume tags to restrict the field to.

        Returns
        -------
        int
            Field tag of the new Restrict field.
        """
        import gmsh
        self._check_initialized()

        field = gmsh.model.mesh.field.add("Restrict")
        gmsh.model.mesh.field.setNumber(field, "InField", field_tag)
        gmsh.model.mesh.field.setNumbers(field, "VolumesList", vol_ids)
        return field

    def add_field_attractor(self, entity_tags, dim=0, n_per_radius=10):
        """Add an Attractor sizing field with sampling.

        Parameters
        ----------
        entity_tags : list of int
            Entity tags to attract mesh refinement toward.
        dim : int
            Entity dimension (0=point, 1=curve, 2=surface).
        n_per_radius : int
            Number of sampling points per radius of curvature.

        Returns
        -------
        int
            Field tag.
        """
        import gmsh
        self._check_initialized()

        field = gmsh.model.mesh.field.add("Attractor")
        key = {0: "PointsList", 1: "CurvesList", 2: "SurfacesList"}[dim]
        gmsh.model.mesh.field.setNumbers(field, key, entity_tags)
        gmsh.model.mesh.field.setNumber(field, "NNodesByEdge", n_per_radius)
        return field

    def remove_field(self, field_tag):
        """Remove a sizing field.

        Parameters
        ----------
        field_tag : int
            Field tag to remove.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.field.remove(field_tag)

        if self._verbose:
            print(f"[GmshBuilder] remove_field tag={field_tag}")

    def remove_all_fields(self):
        """Remove all sizing fields."""
        import gmsh
        self._check_initialized()

        field_tags = gmsh.model.mesh.field.list()
        for tag in field_tags:
            gmsh.model.mesh.field.remove(tag)

        if self._verbose:
            print(f"[GmshBuilder] remove_all_fields ({len(field_tags)} removed)")

    def get_field_count(self):
        """Return the number of active sizing fields.

        Returns
        -------
        int
            Number of active fields.
        """
        import gmsh
        self._check_initialized()

        return len(gmsh.model.mesh.field.list())

    # ------------------------------------------------------------------
    # Transfinite and entity-level mesh size control
    # ------------------------------------------------------------------

    def set_transfinite_surface(self, surf_tag, arrangement="Left",
                                corners=None):
        """Set transfinite meshing on a surface directly.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.
        arrangement : str
            Corner arrangement: 'Left', 'Right', 'AlternateLeft',
            'AlternateRight'.
        corners : list of int or None
            Corner point tags. If None, GMSH selects automatically.
        """
        import gmsh
        self._check_initialized()

        if corners is None:
            corners = []
        gmsh.model.mesh.setTransfiniteSurface(surf_tag, arrangement, corners)

        if self._verbose:
            print(f"[GmshBuilder] set_transfinite_surface surf={surf_tag} "
                  f"arrangement={arrangement}")

    def set_transfinite_volume(self, vol_tag, corners=None):
        """Set transfinite meshing on a volume directly.

        Parameters
        ----------
        vol_tag : int
            GMSH volume tag.
        corners : list of int or None
            Corner point tags (5, 6, or 8 points). If None, GMSH selects
            automatically.
        """
        import gmsh
        self._check_initialized()

        if corners is None:
            corners = []
        gmsh.model.mesh.setTransfiniteVolume(vol_tag, corners)

        if self._verbose:
            print(f"[GmshBuilder] set_transfinite_volume vol={vol_tag}")

    def set_mesh_size_on_surface(self, surf_tag, size):
        """Set mesh size on all boundary points of a surface.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.
        size : float
            Target element size.
        """
        import gmsh
        self._check_initialized()

        bounds = gmsh.model.getBoundary([(2, surf_tag)], recursive=True)
        points = [dt for dt in bounds if dt[0] == 0]
        if points:
            gmsh.model.mesh.setSize(points, size)

        if self._verbose:
            print(f"[GmshBuilder] set_mesh_size_on_surface surf={surf_tag} "
                  f"size={size}")

    def set_mesh_size_on_curve(self, curve_tag, size):
        """Set mesh size on all boundary points of a curve.

        Parameters
        ----------
        curve_tag : int
            GMSH curve tag.
        size : float
            Target element size.
        """
        import gmsh
        self._check_initialized()

        bounds = gmsh.model.getBoundary([(1, curve_tag)], recursive=True)
        points = [dt for dt in bounds if dt[0] == 0]
        if points:
            gmsh.model.mesh.setSize(points, size)

        if self._verbose:
            print(f"[GmshBuilder] set_mesh_size_on_curve curve={curve_tag} "
                  f"size={size}")

    def set_mesh_size_on_point(self, point_tag, size):
        """Set mesh size at a specific point.

        Parameters
        ----------
        point_tag : int
            GMSH point tag.
        size : float
            Target element size.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.setSize([(0, point_tag)], size)

        if self._verbose:
            print(f"[GmshBuilder] set_mesh_size_on_point point={point_tag} "
                  f"size={size}")

    def set_mesh_size_factor(self, factor):
        """Set global mesh size factor (multiplier).

        All mesh sizes are multiplied by this factor.

        Parameters
        ----------
        factor : float
            Mesh size scaling factor (e.g., 0.5 halves all sizes).
        """
        import gmsh
        self._check_initialized()

        gmsh.option.setNumber('Mesh.CharacteristicLengthFactor', factor)

        if self._verbose:
            print(f"[GmshBuilder] set_mesh_size_factor = {factor}")

    def set_mesh_size_from_boundary(self, extend):
        """Enable or disable extending mesh size from boundary.

        Parameters
        ----------
        extend : bool
            If True, mesh sizes are extended from boundary into volume.
            If False, only sizing fields and global size apply.
        """
        import gmsh
        self._check_initialized()

        gmsh.option.setNumber(
            'Mesh.CharacteristicLengthExtendFromBoundary',
            1 if extend else 0)

        if self._verbose:
            print(f"[GmshBuilder] set_mesh_size_from_boundary = {extend}")

    # ------------------------------------------------------------------
    # Algorithm control
    # ------------------------------------------------------------------

    def set_algorithm_2d(self, algo):
        """Set only the 2D meshing algorithm.

        Parameters
        ----------
        algo : int
            2D algorithm: 1=MeshAdapt, 5=Delaunay, 6=Frontal-Delaunay,
            8=Frontal-Delaunay-Quad, 9=Packing-Parallelograms.
        """
        import gmsh
        self._check_initialized()

        gmsh.option.setNumber('Mesh.Algorithm', algo)

        if self._verbose:
            print(f"[GmshBuilder] set_algorithm_2d = {algo}")

    def set_algorithm_3d(self, algo):
        """Set only the 3D meshing algorithm.

        Parameters
        ----------
        algo : int
            3D algorithm: 1=Delaunay, 4=Frontal, 10=HXT.
        """
        import gmsh
        self._check_initialized()

        gmsh.option.setNumber('Mesh.Algorithm3D', algo)

        if self._verbose:
            print(f"[GmshBuilder] set_algorithm_3d = {algo}")

    def set_subdivision_algorithm(self, algo):
        """Set the subdivision algorithm.

        Parameters
        ----------
        algo : int
            0 = none (default), 1 = all quads, 2 = all hexas.
        """
        import gmsh
        self._check_initialized()

        gmsh.option.setNumber('Mesh.SubdivisionAlgorithm', algo)

        if self._verbose:
            print(f"[GmshBuilder] set_subdivision_algorithm = {algo}")

    def get_mesh_options(self):
        """Return a dict of current mesh-related GMSH options.

        Returns
        -------
        dict
            Dictionary with keys: Algorithm, Algorithm3D,
            CharacteristicLengthMin, CharacteristicLengthMax,
            CharacteristicLengthFactor, CharacteristicLengthExtendFromBoundary,
            MeshSizeFromCurvature, RecombineAll, Smoothing,
            SubdivisionAlgorithm, ElementOrder.
        """
        import gmsh
        self._check_initialized()

        option_keys = [
            'Mesh.Algorithm',
            'Mesh.Algorithm3D',
            'Mesh.CharacteristicLengthMin',
            'Mesh.CharacteristicLengthMax',
            'Mesh.CharacteristicLengthFactor',
            'Mesh.CharacteristicLengthExtendFromBoundary',
            'Mesh.MeshSizeFromCurvature',
            'Mesh.RecombineAll',
            'Mesh.Smoothing',
            'Mesh.SubdivisionAlgorithm',
            'Mesh.ElementOrder',
        ]

        result = {}
        for key in option_keys:
            short_key = key.replace('Mesh.', '')
            result[short_key] = gmsh.option.getNumber(key)

        return result

    # ------------------------------------------------------------------
    # Additional sizing fields (Frustum, Laplacian, Gradient, etc.)
    # ------------------------------------------------------------------

    def add_field_frustum(self, x1, y1, z1, x2, y2, z2,
                          r1_inner, r1_outer, r2_inner, r2_outer,
                          v1_inner, v1_outer, v2_inner, v2_outer):
        """Add a Frustum sizing field.

        A frustum is a truncated cone defined by two endpoints, inner/outer
        radii at each end, and element sizes at each end inside/outside.

        Parameters
        ----------
        x1, y1, z1 : float
            Coordinates of the first endpoint.
        x2, y2, z2 : float
            Coordinates of the second endpoint.
        r1_inner : float
            Inner radius at the first endpoint.
        r1_outer : float
            Outer radius at the first endpoint.
        r2_inner : float
            Inner radius at the second endpoint.
        r2_outer : float
            Outer radius at the second endpoint.
        v1_inner : float
            Element size at the inner boundary of the first endpoint.
        v1_outer : float
            Element size at the outer boundary of the first endpoint.
        v2_inner : float
            Element size at the inner boundary of the second endpoint.
        v2_outer : float
            Element size at the outer boundary of the second endpoint.

        Returns
        -------
        int
            Field tag.
        """
        import gmsh
        self._check_initialized()

        field = gmsh.model.mesh.field.add("Frustum")
        gmsh.model.mesh.field.setNumber(field, "X1", x1)
        gmsh.model.mesh.field.setNumber(field, "Y1", y1)
        gmsh.model.mesh.field.setNumber(field, "Z1", z1)
        gmsh.model.mesh.field.setNumber(field, "X2", x2)
        gmsh.model.mesh.field.setNumber(field, "Y2", y2)
        gmsh.model.mesh.field.setNumber(field, "Z2", z2)
        gmsh.model.mesh.field.setNumber(field, "InnerR1", r1_inner)
        gmsh.model.mesh.field.setNumber(field, "OuterR1", r1_outer)
        gmsh.model.mesh.field.setNumber(field, "InnerR2", r2_inner)
        gmsh.model.mesh.field.setNumber(field, "OuterR2", r2_outer)
        gmsh.model.mesh.field.setNumber(field, "V1_inner", v1_inner)
        gmsh.model.mesh.field.setNumber(field, "V1_outer", v1_outer)
        gmsh.model.mesh.field.setNumber(field, "V2_inner", v2_inner)
        gmsh.model.mesh.field.setNumber(field, "V2_outer", v2_outer)

        if self._verbose:
            print(f"[GmshBuilder] add_field_frustum tag={field} "
                  f"({x1},{y1},{z1})->({x2},{y2},{z2})")

        return field

    def add_field_laplacian(self, field_tag):
        """Add a Laplacian sizing field referencing another field.

        The Laplacian field smooths the referenced field by solving a
        Laplace equation on the mesh sizing.

        Parameters
        ----------
        field_tag : int
            Tag of the field to apply the Laplacian to.

        Returns
        -------
        int
            New field tag.
        """
        import gmsh
        self._check_initialized()

        new_field = gmsh.model.mesh.field.add("Laplacian")
        gmsh.model.mesh.field.setNumber(new_field, "InField", field_tag)

        if self._verbose:
            print(f"[GmshBuilder] add_field_laplacian tag={new_field} "
                  f"InField={field_tag}")

        return new_field

    def add_field_gradient(self, field_tag):
        """Add a Gradient sizing field referencing another field.

        The Gradient field computes the gradient magnitude of the
        referenced field and uses it as a sizing metric.

        Parameters
        ----------
        field_tag : int
            Tag of the field to compute the gradient of.

        Returns
        -------
        int
            New field tag.
        """
        import gmsh
        self._check_initialized()

        new_field = gmsh.model.mesh.field.add("Gradient")
        gmsh.model.mesh.field.setNumber(new_field, "InField", field_tag)

        if self._verbose:
            print(f"[GmshBuilder] add_field_gradient tag={new_field} "
                  f"InField={field_tag}")

        return new_field

    def add_field_curvature_field(self):
        """Add a Curvature sizing field.

        The Curvature field sets mesh element size based on the local
        surface curvature of the geometry.

        Returns
        -------
        int
            Field tag.
        """
        import gmsh
        self._check_initialized()

        field = gmsh.model.mesh.field.add("Curvature")

        if self._verbose:
            print(f"[GmshBuilder] add_field_curvature_field tag={field}")

        return field

    def add_field_octree(self, in_field):
        """Add an Octree sizing field if available.

        The Octree field provides octree-based adaptive mesh sizing
        driven by a referenced input field.
        This field type may not be available in all GMSH versions.

        Parameters
        ----------
        in_field : int
            Tag of the input field that drives the octree refinement.

        Returns
        -------
        int
            Field tag.

        Raises
        ------
        RuntimeError
            If the Octree field type is not available in the current
            GMSH version.
        """
        import gmsh
        self._check_initialized()

        try:
            field = gmsh.model.mesh.field.add("Octree")
        except Exception as exc:
            raise RuntimeError(
                "Octree field type is not available in this GMSH version"
            ) from exc

        gmsh.model.mesh.field.setNumber(field, "InField", in_field)

        if self._verbose:
            print(f"[GmshBuilder] add_field_octree tag={field} "
                  f"InField={in_field}")

        return field

    def add_field_external_process(self, command):
        """Add an ExternalProcess sizing field.

        The ExternalProcess field delegates mesh sizing to an external
        program specified by the command string. The external program
        reads coordinates from stdin and writes sizes to stdout.

        Parameters
        ----------
        command : str
            Command to execute for the external sizing process.

        Returns
        -------
        int
            Field tag.
        """
        import gmsh
        self._check_initialized()

        field = gmsh.model.mesh.field.add("ExternalProcess")
        gmsh.model.mesh.field.setString(field, "CommandLine", command)

        if self._verbose:
            print(f"[GmshBuilder] add_field_external_process tag={field} "
                  f"command='{command}'")

        return field

    def add_field_structured(self, filename):
        """Add a Structured sizing field from a file.

        The Structured field reads mesh size data from a structured
        grid file (e.g., a .pos file with a regular grid of values).

        Parameters
        ----------
        filename : str
            Path to the structured grid file.

        Returns
        -------
        int
            Field tag.
        """
        import gmsh
        self._check_initialized()

        field = gmsh.model.mesh.field.add("Structured")
        gmsh.model.mesh.field.setString(field, "FileName", filename)

        if self._verbose:
            print(f"[GmshBuilder] add_field_structured tag={field} "
                  f"file='{filename}'")

        return field

    def add_field_post_view(self, view_tag):
        """Add a PostView sizing field.

        The PostView field uses data from a post-processing view
        (e.g., an existing .pos file loaded in GMSH) as the mesh
        sizing function.

        Parameters
        ----------
        view_tag : int
            GMSH view tag to use as size data source.

        Returns
        -------
        int
            Field tag.
        """
        import gmsh
        self._check_initialized()

        field = gmsh.model.mesh.field.add("PostView")
        gmsh.model.mesh.field.setNumber(field, "ViewTag", view_tag)

        if self._verbose:
            print(f"[GmshBuilder] add_field_post_view tag={field} "
                  f"ViewTag={view_tag}")

        return field

    def add_field_param(self, field_tag, key, value):
        """Set a numeric or string parameter on an existing sizing field.

        Parameters
        ----------
        field_tag : int
            Tag of the existing field.
        key : str
            Parameter name (e.g., 'VIn', 'FileName').
        value : float, int, or str
            Parameter value. Numeric types use ``setNumber``, string
            types use ``setString``.
        """
        import gmsh
        self._check_initialized()

        if isinstance(value, str):
            gmsh.model.mesh.field.setString(field_tag, key, value)
        else:
            gmsh.model.mesh.field.setNumber(field_tag, key, value)

        if self._verbose:
            print(f"[GmshBuilder] add_field_param field={field_tag} "
                  f"{key}={value}")

    def get_field_type(self, field_tag):
        """Get the type name of an existing sizing field.

        Parameters
        ----------
        field_tag : int
            Field tag to query.

        Returns
        -------
        str
            Field type name (e.g., 'Box', 'Threshold', 'Distance').
        """
        import gmsh
        self._check_initialized()

        field_type = gmsh.model.mesh.field.getType(field_tag)

        if self._verbose:
            print(f"[GmshBuilder] get_field_type field={field_tag} "
                  f"-> '{field_type}'")

        return field_type

    # ------------------------------------------------------------------
    # Compound meshing
    # ------------------------------------------------------------------

    def set_compound(self, dim, tags):
        """Set compound meshing on entities of given dimension.

        Compound meshing treats multiple entities as a single entity
        for meshing purposes, producing a conformal mesh across them.

        Parameters
        ----------
        dim : int
            Entity dimension (1=curve, 2=surface, 3=volume).
        tags : list of int
            Entity tags to combine into a compound.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.setCompound(dim, tags)
        self._meshed = False

        if self._verbose:
            print(f"[GmshBuilder] set_compound dim={dim} tags={tags}")

    def set_compound_surfaces(self, surf_tags):
        """Set compound meshing on surfaces.

        Treats multiple surfaces as a single surface for meshing.

        Parameters
        ----------
        surf_tags : list of int
            Surface tags to combine into a compound surface.
        """
        self.set_compound(2, surf_tags)

    def set_compound_curves(self, curve_tags):
        """Set compound meshing on curves.

        Treats multiple curves as a single curve for meshing.

        Parameters
        ----------
        curve_tags : list of int
            Curve tags to combine into a compound curve.
        """
        self.set_compound(1, curve_tags)

    def set_transfinite_automatic(self, vol_tag=None, corner_angle=2.35,
                                  recombine=True):
        """Set automatic transfinite meshing constraints.

        Automatically applies transfinite constraints to surfaces and
        volumes that have a compatible topology (4-sided surfaces,
        6-faced volumes).

        Parameters
        ----------
        vol_tag : int or None
            Volume tag to apply to. If None, applies to all entities.
        corner_angle : float
            Maximum corner angle (in radians) for automatic corner
            detection. Default is 2.35 (~135 degrees).
        recombine : bool
            If True, recombine triangles into quads on transfinite
            surfaces.
        """
        import gmsh
        self._check_initialized()

        if vol_tag is not None:
            dimTags = [(3, vol_tag)]
        else:
            dimTags = []

        gmsh.model.mesh.setTransfiniteAutomatic(dimTags, corner_angle,
                                                recombine)
        self._meshed = False

        if self._verbose:
            target = f"vol={vol_tag}" if vol_tag is not None else "all"
            print(f"[GmshBuilder] set_transfinite_automatic {target} "
                  f"corner_angle={corner_angle} recombine={recombine}")

    def set_transfinite_automatic_all(self, corner_angle=2.35,
                                      recombine=True):
        """Apply automatic transfinite meshing to all volumes.

        Convenience method that calls ``set_transfinite_automatic``
        with ``vol_tag=None``.

        Parameters
        ----------
        corner_angle : float
            Maximum corner angle (in radians) for automatic corner
            detection. Default is 2.35 (~135 degrees).
        recombine : bool
            If True, recombine triangles into quads on transfinite
            surfaces.
        """
        self.set_transfinite_automatic(vol_tag=None,
                                       corner_angle=corner_angle,
                                       recombine=recombine)

    # ------------------------------------------------------------------
    # Per-entity algorithm control
    # ------------------------------------------------------------------

    def set_algorithm_on_surface(self, surf_tag, algo):
        """Set the meshing algorithm on a specific surface.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.
        algo : int
            Meshing algorithm: 1=MeshAdapt, 5=Delaunay,
            6=Frontal-Delaunay, 8=Frontal-Delaunay for Quads,
            9=Packing of Parallelograms.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.setAlgorithm(2, surf_tag, algo)
        self._meshed = False

        if self._verbose:
            print(f"[GmshBuilder] set_algorithm_on_surface "
                  f"surf={surf_tag} algo={algo}")

    def set_algorithm_on_volume(self, vol_tag, algo):
        """Set the 3D meshing algorithm.

        GMSH does not support per-volume algorithm selection via
        ``setAlgorithm`` for dim=3. This method sets the global
        ``Mesh.Algorithm3D`` option instead. The ``vol_tag`` parameter
        is accepted for interface consistency but is not used.

        Parameters
        ----------
        vol_tag : int
            GMSH volume tag (unused -- algorithm is set globally).
        algo : int
            3D meshing algorithm: 1=Delaunay, 4=Frontal, 10=HXT.
        """
        import gmsh
        self._check_initialized()

        gmsh.option.setNumber('Mesh.Algorithm3D', algo)
        self._meshed = False

        if self._verbose:
            print(f"[GmshBuilder] set_algorithm_on_volume "
                  f"algo={algo} (global Mesh.Algorithm3D)")

    def set_smoothing_on_surface(self, surf_tag, n_steps):
        """Set the number of smoothing steps on a specific surface.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.
        n_steps : int
            Number of smoothing iterations.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.setSmoothing(2, surf_tag, n_steps)

        if self._verbose:
            print(f"[GmshBuilder] set_smoothing_on_surface "
                  f"surf={surf_tag} n_steps={n_steps}")

    def set_smoothing_on_volume(self, vol_tag, n_steps):
        """Set the number of smoothing steps on a specific volume.

        Parameters
        ----------
        vol_tag : int
            GMSH volume tag.
        n_steps : int
            Number of smoothing iterations.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.setSmoothing(3, vol_tag, n_steps)

        if self._verbose:
            print(f"[GmshBuilder] set_smoothing_on_volume "
                  f"vol={vol_tag} n_steps={n_steps}")

    def set_reverse_mesh(self, dim, tag):
        """Reverse the mesh orientation on an entity.

        This reverses the element node ordering for the specified entity,
        effectively flipping the normal direction for surface elements.

        Parameters
        ----------
        dim : int
            Entity dimension (1=curve, 2=surface).
        tag : int
            Entity tag.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.setReverse(dim, tag)

        if self._verbose:
            print(f"[GmshBuilder] set_reverse_mesh dim={dim} tag={tag}")

    # ------------------------------------------------------------------
    # Size callback
    # ------------------------------------------------------------------

    def set_size_callback(self, callback):
        """Set a Python callback function for mesh sizing.

        The callback is called for each mesh vertex during meshing and
        should return the desired element size at that location.

        Parameters
        ----------
        callback : callable
            Function with signature
            ``callback(dim, tag, x, y, z, lc) -> float``
            where ``dim`` is the entity dimension, ``tag`` is the entity
            tag, ``(x, y, z)`` are the coordinates, and ``lc`` is the
            default characteristic length. Must return the desired
            element size as a float.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.setSizeCallback(callback)
        self._meshed = False

        if self._verbose:
            print(f"[GmshBuilder] set_size_callback -> {callback}")

    def remove_size_callback(self):
        """Remove a previously set size callback function.

        Restores default mesh sizing behavior.
        """
        import gmsh
        self._check_initialized()

        try:
            gmsh.model.mesh.setSizeCallback(None)
        except Exception:
            # Some GMSH versions may not accept None; pass a no-op
            # lambda that returns the default lc instead.
            gmsh.model.mesh.setSizeCallback(
                lambda dim, tag, x, y, z, lc: lc
            )

        if self._verbose:
            print("[GmshBuilder] remove_size_callback")

    def set_size_at_parametric_points(self, dim, tag, parametric_coords,
                                      sizes):
        """Set mesh sizes at specific parametric points on an entity.

        Parameters
        ----------
        dim : int
            Entity dimension (1=curve, 2=surface).
        tag : int
            Entity tag.
        parametric_coords : list of float
            Parametric coordinates on the entity. For curves, a flat
            list of u values. For surfaces, a flat list of (u, v) pairs.
        sizes : list of float
            Desired element sizes at each parametric point.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.setSizeAtParametricPoints(dim, tag,
                                                   parametric_coords,
                                                   sizes)

        if self._verbose:
            n_pts = len(sizes)
            print(f"[GmshBuilder] set_size_at_parametric_points "
                  f"dim={dim} tag={tag} n_points={n_pts}")

    # ------------------------------------------------------------------
    # Constraints and visibility
    # ------------------------------------------------------------------

    def set_mesh_visibility(self, dim, tag, visible):
        """Set mesh visibility for an entity.

        Controls whether the mesh of a specific entity is displayed
        in the GMSH GUI.

        Parameters
        ----------
        dim : int
            Entity dimension (0=point, 1=curve, 2=surface, 3=volume).
        tag : int
            Entity tag.
        visible : bool
            If True, the entity mesh is visible; if False, hidden.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.setVisibility([(dim, tag)], 1 if visible else 0)

        if self._verbose:
            state = "visible" if visible else "hidden"
            print(f"[GmshBuilder] set_mesh_visibility ({dim},{tag}) "
                  f"-> {state}")

    def remove_mesh_constraints(self, dim=-1, tag=-1):
        """Remove mesh constraints from entities.

        Removes transfinite, embedded, and other mesh constraints.
        If both dim and tag are -1, removes constraints from all
        entities.

        Parameters
        ----------
        dim : int
            Entity dimension. -1 for all dimensions.
        tag : int
            Entity tag. -1 for all entities of the given dimension.
        """
        import gmsh
        self._check_initialized()

        try:
            if dim == -1 and tag == -1:
                gmsh.model.mesh.removeConstraints()
            elif tag == -1:
                # Remove all constraints for given dimension
                entities = gmsh.model.getEntities(dim)
                if entities:
                    gmsh.model.mesh.removeConstraints(entities)
            elif dim == -1:
                # Remove constraints for specific tag in all dimensions
                for d in range(4):
                    try:
                        gmsh.model.mesh.removeConstraints([(d, tag)])
                    except Exception:
                        pass
            else:
                gmsh.model.mesh.removeConstraints([(dim, tag)])
        except AttributeError:
            # removeConstraints may not exist in older GMSH versions.
            # Attempt per-entity removal via clear.
            try:
                gmsh.model.mesh.clear([(dim, tag)] if dim >= 0 else [])
            except Exception:
                pass

        self._meshed = False

        if self._verbose:
            target = (f"({dim},{tag})" if dim >= 0
                      else "all entities")
            print(f"[GmshBuilder] remove_mesh_constraints {target}")

    def set_outward_orientation(self, tag):
        """Set outward orientation for surface mesh normals of a volume.

        Ensures that all surface mesh normals of the specified volume
        point outward.

        Parameters
        ----------
        tag : int
            Volume tag whose boundary surface normals are oriented
            outward.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.setOutwardOrientation(tag)

        if self._verbose:
            print(f"[GmshBuilder] set_outward_orientation vol={tag}")

    def remove_embedded(self, dim, tag):
        """Remove embedded entities from a higher-dimensional entity.

        Parameters
        ----------
        dim : int
            Dimension of the entity to remove embeddings from
            (2=surface, 3=volume).
        tag : int
            Tag of the entity to remove embeddings from.

        Raises
        ------
        RuntimeError
            If removeEmbedded is not available in the current GMSH
            version.
        """
        import gmsh
        self._check_initialized()

        try:
            gmsh.model.mesh.removeEmbedded([(dim, tag)])
        except (AttributeError, Exception) as exc:
            raise RuntimeError(
                "removeEmbedded is not available in this GMSH version"
            ) from exc

        self._meshed = False

        if self._verbose:
            print(f"[GmshBuilder] remove_embedded ({dim},{tag})")

    # ------------------------------------------------------------------
    # Get queries
    # ------------------------------------------------------------------

    def get_mesh_sizes(self):
        """Get mesh edge lengths for all elements in the current mesh.

        Iterates over all 3D elements, extracts unique edges, and
        computes their lengths. Requires a mesh to have been generated.

        Returns
        -------
        list of float
            Sorted list of all unique edge lengths in meters.

        Raises
        ------
        RuntimeError
            If no mesh has been generated yet.
        """
        import gmsh
        import numpy as np
        self._check_initialized()

        if not self._meshed:
            raise RuntimeError(
                "generate() must be called before get_mesh_sizes()")

        # Collect all node coordinates
        node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
        # Build tag -> coordinate mapping
        coord_map = {}
        for i, nt in enumerate(node_tags):
            coord_map[int(nt)] = np.array(node_coords[3*i:3*i+3])

        # Element-type-specific edge definitions (vertex index pairs)
        HEX_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0),
                      (4, 5), (5, 6), (6, 7), (7, 4),
                      (0, 4), (1, 5), (2, 6), (3, 7)]
        TET_EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
        TRI_EDGES = [(0, 1), (1, 2), (2, 0)]
        WEDGE_EDGES = [(0, 1), (1, 2), (2, 0),
                        (3, 4), (4, 5), (5, 3),
                        (0, 3), (1, 4), (2, 5)]
        PYRAMID_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0),
                          (0, 4), (1, 4), (2, 4), (3, 4)]
        QUAD_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0)]
        LINE_EDGES = [(0, 1)]

        # Collect unique edges from all element types
        edges = set()
        for dim in [1, 2, 3]:
            elem_types, elem_tags, elem_node_tags = (
                gmsh.model.mesh.getElements(dim))
            for et_idx, etype in enumerate(elem_types):
                props = gmsh.model.mesh.getElementProperties(etype)
                # props: (name, dim, order, numNodes, localNodeCoords,
                #         numPrimaryNodes)
                n_nodes = props[3]
                n_primary = props[5] if len(props) > 5 else n_nodes
                nodes_flat = elem_node_tags[et_idx]

                # Select proper edge topology based on primary node count
                if n_primary == 8:      # hex
                    edge_defs = HEX_EDGES
                elif n_primary == 4 and dim == 3:  # tet
                    edge_defs = TET_EDGES
                elif n_primary == 4 and dim == 2:  # quad
                    edge_defs = QUAD_EDGES
                elif n_primary == 6:    # wedge
                    edge_defs = WEDGE_EDGES
                elif n_primary == 5:    # pyramid
                    edge_defs = PYRAMID_EDGES
                elif n_primary == 3:    # triangle
                    edge_defs = TRI_EDGES
                elif n_primary == 2:    # line
                    edge_defs = LINE_EDGES
                else:
                    # Fallback: all pairs (for unknown element types)
                    edge_defs = [(i, j) for i in range(n_primary)
                                 for j in range(i + 1, n_primary)]

                for k in range(0, len(nodes_flat), n_nodes):
                    elem_nodes = [int(nodes_flat[k + m])
                                  for m in range(n_primary)]
                    for a_idx, b_idx in edge_defs:
                        a = elem_nodes[a_idx]
                        b = elem_nodes[b_idx]
                        edge = (min(a, b), max(a, b))
                        edges.add(edge)

        # Compute edge lengths
        lengths = []
        for a, b in edges:
            if a in coord_map and b in coord_map:
                lengths.append(float(np.linalg.norm(
                    coord_map[a] - coord_map[b])))

        lengths.sort()

        if self._verbose:
            if lengths:
                print(f"[GmshBuilder] get_mesh_sizes: {len(lengths)} edges, "
                      f"min={lengths[0]:.6g} max={lengths[-1]:.6g}")
            else:
                print("[GmshBuilder] get_mesh_sizes: no edges found")

        return lengths

    def get_embedded_entities(self, dim, tag):
        """Get embedded entities within a higher-dimensional entity.

        Parameters
        ----------
        dim : int
            Dimension of the host entity (2=surface, 3=volume).
        tag : int
            Tag of the host entity.

        Returns
        -------
        list of tuple
            List of (dimension, tag) tuples for each embedded entity.
            Returns an empty list if the API is not available.
        """
        import gmsh
        self._check_initialized()

        try:
            result = gmsh.model.mesh.getEmbedded(dim, tag)
            if self._verbose:
                print(f"[GmshBuilder] get_embedded_entities ({dim},{tag}) "
                      f"-> {len(result)} entities")
            return list(result)
        except (AttributeError, Exception):
            if self._verbose:
                print(f"[GmshBuilder] get_embedded_entities ({dim},{tag}) "
                      "-> not available, returning []")
            return []

    def get_field_list(self):
        """Get the list of all active sizing field tags.

        Returns
        -------
        list of int
            List of active field tags.
        """
        import gmsh
        self._check_initialized()

        fields = list(gmsh.model.mesh.field.list())

        if self._verbose:
            print(f"[GmshBuilder] get_field_list -> {fields}")

        return fields

    # ------------------------------------------------------------------
    # CFD convenience methods
    # ------------------------------------------------------------------

    def add_size_gradient(self, entity_tags, dim, lc_near, lc_far,
                          dist_near, dist_far, set_as_background=True):
        """Add distance-based graded mesh sizing near entities.

        Combines Distance + Threshold fields into a one-step refinement
        gradient.  This is the most common CFD sizing pattern: fine mesh
        near walls/boundaries, coarse mesh in the far field.

        Parameters
        ----------
        entity_tags : list of int
            GMSH entity tags to measure distance from.
        dim : int
            Entity dimension (0=point, 1=curve, 2=surface).
        lc_near : float
            Element size near the entities (within dist_near).
        lc_far : float
            Element size far from the entities (beyond dist_far).
        dist_near : float
            Distance at which lc_near applies.
        dist_far : float
            Distance at which lc_far applies.
        set_as_background : bool
            If True, set this field as the background mesh.

        Returns
        -------
        int
            Threshold field tag.
        """
        import gmsh
        self._check_initialized()

        dist_tag = gmsh.model.mesh.field.add("Distance")
        key = {0: "PointsList", 1: "CurvesList", 2: "SurfacesList"}[dim]
        gmsh.model.mesh.field.setNumbers(dist_tag, key, entity_tags)

        thresh_tag = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(thresh_tag, "InField", dist_tag)
        gmsh.model.mesh.field.setNumber(thresh_tag, "SizeMin", lc_near)
        gmsh.model.mesh.field.setNumber(thresh_tag, "SizeMax", lc_far)
        gmsh.model.mesh.field.setNumber(thresh_tag, "DistMin", dist_near)
        gmsh.model.mesh.field.setNumber(thresh_tag, "DistMax", dist_far)

        if set_as_background:
            gmsh.model.mesh.field.setAsBackgroundMesh(thresh_tag)

        if self._verbose:
            print(f"[GmshBuilder] add_size_gradient dim={dim} "
                  f"lc=[{lc_near},{lc_far}] dist=[{dist_near},{dist_far}] "
                  f"-> field={thresh_tag}")

        return thresh_tag

    def set_background_mesh_from_fields(self, field_tags, operator='min'):
        """Combine multiple sizing fields and set as background mesh.

        Creates a combining field (Min, Max, or Mean) from the given
        field tags and sets it as the active background mesh field.

        Parameters
        ----------
        field_tags : list of int
            Field tags to combine.
        operator : str
            Combining operator: 'min', 'max', or 'mean'.

        Returns
        -------
        int
            Combined field tag.

        Raises
        ------
        ValueError
            If operator is not 'min', 'max', or 'mean'.
        """
        import gmsh
        self._check_initialized()

        op_map = {'min': 'Min', 'max': 'Max', 'mean': 'Mean'}
        op_name = op_map.get(operator.lower())
        if op_name is None:
            raise ValueError(
                f"Unknown operator '{operator}'. Use 'min', 'max', or 'mean'.")

        combined = gmsh.model.mesh.field.add(op_name)
        gmsh.model.mesh.field.setNumbers(combined, "FieldsList", field_tags)
        gmsh.model.mesh.field.setAsBackgroundMesh(combined)

        if self._verbose:
            print(f"[GmshBuilder] set_background_mesh_from_fields "
                  f"op={op_name} fields={field_tags} -> tag={combined}")

        return combined

    def set_field_numbers(self, field_tag, key, values):
        """Set a list-valued parameter on an existing sizing field.

        This is needed for setting entity lists (e.g., CurvesList,
        SurfacesList, FieldsList) on fields. The existing
        ``add_field_param`` handles scalar (setNumber/setString) values
        only.

        Parameters
        ----------
        field_tag : int
            Tag of the existing field.
        key : str
            Parameter name (e.g., 'CurvesList', 'SurfacesList',
            'FieldsList').
        values : list of int or float
            List of values to set.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.field.setNumbers(field_tag, key, values)

        if self._verbose:
            print(f"[GmshBuilder] set_field_numbers field={field_tag} "
                  f"{key}={values}")

    def add_boundary_layer_field(self, entity_tags, dim,
                                 lc_min, lc_max, dist_min, dist_max,
                                 n_layers=None, ratio=None,
                                 fan_points=None,
                                 intersect_metrics=False):
        """Add a boundary layer sizing field (BoundaryLayer).

        Creates a BoundaryLayer field for viscous wall resolution in
        CFD meshes.  Unlike ``set_boundary_layer`` which immediately
        activates, this method returns the field tag for composition
        with other fields via ``set_background_mesh_from_fields``.

        Parameters
        ----------
        entity_tags : list of int
            GMSH entity tags (curves or surfaces) defining the wall.
        dim : int
            Entity dimension (1=curve, 2=surface).
        lc_min : float
            First layer element size (wall-adjacent).
        lc_max : float
            Outer layer element size (far from wall).
        dist_min : float
            Distance of first layer from wall.
        dist_max : float
            Total boundary layer thickness.
        n_layers : int or None
            Number of boundary layers.  If None, GMSH auto-computes.
        ratio : float or None
            Growth ratio between successive layers.
        fan_points : list of int or None
            Point tags for fan-shaped boundary layers at corners.
        intersect_metrics : bool
            If True, use IntersectMetrics for improved BL quality.

        Returns
        -------
        int
            BoundaryLayer field tag.
        """
        import gmsh
        self._check_initialized()

        field = gmsh.model.mesh.field.add("BoundaryLayer")

        key = {1: "CurvesList", 2: "SurfacesList"}[dim]
        gmsh.model.mesh.field.setNumbers(field, key, entity_tags)

        gmsh.model.mesh.field.setNumber(field, "hwall_n", lc_min)
        gmsh.model.mesh.field.setNumber(field, "hfar", lc_max)
        gmsh.model.mesh.field.setNumber(field, "thickness", dist_max)

        if ratio is not None:
            gmsh.model.mesh.field.setNumber(field, "ratio", ratio)

        if n_layers is not None:
            gmsh.model.mesh.field.setNumber(field, "NbLayers", n_layers)

        if fan_points is not None:
            gmsh.model.mesh.field.setNumbers(
                field, "FanPointsList", fan_points)

        if intersect_metrics:
            gmsh.model.mesh.field.setNumber(
                field, "IntersectMetrics", 1)

        if self._verbose:
            print(f"[GmshBuilder] add_boundary_layer_field dim={dim} "
                  f"lc=[{lc_min},{lc_max}] dist=[{dist_min},{dist_max}] "
                  f"-> field={field}")

        return field

    def add_field_auto_mesh_size(self):
        """Add an AutomaticMeshSizeField for geometry-based sizing.

        This field automatically computes element sizes based on
        local geometry (curvature, proximity, feature size).

        Returns
        -------
        int
            Field tag.

        Raises
        ------
        RuntimeError
            If the field type is not available (requires GMSH >= 4.11).
        """
        import gmsh
        self._check_initialized()

        try:
            field = gmsh.model.mesh.field.add("AutomaticMeshSizeField")
        except Exception as exc:
            raise RuntimeError(
                "AutomaticMeshSizeField not available. "
                "Requires GMSH >= 4.11.") from exc

        if self._verbose:
            print(f"[GmshBuilder] add_field_auto_mesh_size -> field={field}")

        return field

    def set_mesh_size_at_points_direct(self, point_tags, size):
        """Set mesh size at specific geometric points (batch).

        Unlike ``set_mesh_size_on_point`` which takes a single point,
        this method sets the same size on a list of point tags in a
        single GMSH call.

        Parameters
        ----------
        point_tags : list of int
            GMSH point tags.
        size : float
            Target element size.
        """
        import gmsh
        self._check_initialized()

        dim_tags = [(0, t) for t in point_tags]
        gmsh.model.mesh.setSize(dim_tags, size)

        if self._verbose:
            print(f"[GmshBuilder] set_mesh_size_at_points_direct "
                  f"{len(point_tags)} points, size={size}")

    def add_wake_refinement(self, curve_tags, direction, length,
                            lc_min, lc_max, width=None):
        """Add anisotropic refinement in a wake region downstream of curves.

        Creates a Box field aligned with the wake direction for
        downstream refinement commonly needed in airfoil and bluff-body
        CFD simulations.

        Parameters
        ----------
        curve_tags : list of int
            Curve tags at the trailing edge / downstream boundary.
        direction : list of float
            Wake direction vector [dx, dy, dz] (will be normalized).
        length : float
            Wake refinement length downstream.
        lc_min : float
            Element size in the wake core.
        lc_max : float
            Element size outside the wake region.
        width : float or None
            Wake width.  If None, estimated from curve bounding box.

        Returns
        -------
        int
            Box field tag for the wake refinement region.
        """
        import gmsh
        import math
        self._check_initialized()

        # Get bounding box of source curves
        xmin_all = float('inf')
        ymin_all = float('inf')
        zmin_all = float('inf')
        xmax_all = float('-inf')
        ymax_all = float('-inf')
        zmax_all = float('-inf')
        for ct in curve_tags:
            x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(1, ct)
            xmin_all = min(xmin_all, x0)
            ymin_all = min(ymin_all, y0)
            zmin_all = min(zmin_all, z0)
            xmax_all = max(xmax_all, x1)
            ymax_all = max(ymax_all, y1)
            zmax_all = max(zmax_all, z1)

        # Normalize direction
        dx, dy, dz = direction
        norm = math.sqrt(dx * dx + dy * dy + dz * dz)
        if norm < 1e-30:
            raise ValueError("Direction vector has zero length")
        dx, dy, dz = dx / norm, dy / norm, dz / norm

        # Compute wake width from bounding box if not given
        if width is None:
            bbox_size = max(xmax_all - xmin_all,
                            ymax_all - ymin_all,
                            zmax_all - zmin_all)
            width = bbox_size * 2.0

        # Compute box aligned with wake
        cx = (xmin_all + xmax_all) / 2 + dx * length / 2
        cy = (ymin_all + ymax_all) / 2 + dy * length / 2
        cz = (zmin_all + zmax_all) / 2 + dz * length / 2

        field = gmsh.model.mesh.field.add("Box")
        gmsh.model.mesh.field.setNumber(field, "VIn", lc_min)
        gmsh.model.mesh.field.setNumber(field, "VOut", lc_max)
        gmsh.model.mesh.field.setNumber(
            field, "XMin", cx - abs(dx) * length / 2 - width / 2)
        gmsh.model.mesh.field.setNumber(
            field, "XMax", cx + abs(dx) * length / 2 + width / 2)
        gmsh.model.mesh.field.setNumber(
            field, "YMin", cy - abs(dy) * length / 2 - width / 2)
        gmsh.model.mesh.field.setNumber(
            field, "YMax", cy + abs(dy) * length / 2 + width / 2)
        gmsh.model.mesh.field.setNumber(
            field, "ZMin", cz - abs(dz) * length / 2 - width / 2)
        gmsh.model.mesh.field.setNumber(
            field, "ZMax", cz + abs(dz) * length / 2 + width / 2)

        if self._verbose:
            print(f"[GmshBuilder] add_wake_refinement "
                  f"direction=[{dx:.2f},{dy:.2f},{dz:.2f}] "
                  f"length={length} lc=[{lc_min},{lc_max}] "
                  f"-> field={field}")

        return field

    def add_terrain_field(self, filename, format='structured'):
        """Add a terrain/elevation data file as a sizing field.

        Loads structured or post-view data and creates a sizing field.
        Commonly used in ocean, geophysical, and terrain meshing.

        Parameters
        ----------
        filename : str
            Path to the data file (.pos for structured grid,
            or GMSH view file for post_view format).
        format : str
            Data format: 'structured' (regular grid .pos file) or
            'post_view' (GMSH view data).

        Returns
        -------
        int
            Field tag.

        Raises
        ------
        ValueError
            If format is not 'structured' or 'post_view'.
        """
        import gmsh
        self._check_initialized()

        if format == 'structured':
            field = gmsh.model.mesh.field.add("Structured")
            gmsh.model.mesh.field.setString(field, "FileName", filename)
        elif format == 'post_view':
            gmsh.merge(filename)
            views = gmsh.view.getTags()
            if not views:
                raise RuntimeError(
                    f"No views found after merging '{filename}'")
            view_tag = views[-1]
            field = gmsh.model.mesh.field.add("PostView")
            gmsh.model.mesh.field.setNumber(field, "ViewTag", view_tag)
        else:
            raise ValueError(
                f"Unknown format '{format}'. Use 'structured' or 'post_view'.")

        if self._verbose:
            print(f"[GmshBuilder] add_terrain_field '{filename}' "
                  f"format={format} -> field={field}")

        return field
