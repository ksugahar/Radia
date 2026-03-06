"""
Geometry and topology query methods for GmshBuilder.

Gmsh equivalents: get_bounding_box, get_center, get_volume, list entities.
"""

import math


class QueryMixin:
    """Mixin providing geometry/topology queries."""

    def get_volumes(self):
        """Get all managed volume IDs.

        Returns
        -------
        list of int
        """
        return list(self._volumes.keys())

    def get_gmsh_tags(self, vol_id=None):
        """Get underlying GMSH volume tags.

        Parameters
        ----------
        vol_id : int, optional
            Specific volume. If None, return all.

        Returns
        -------
        list of int
        """
        if vol_id is not None:
            if vol_id not in self._volumes:
                raise KeyError(f"Volume {vol_id} not found")
            return list(self._volumes[vol_id])

        all_tags = []
        for tags in self._volumes.values():
            all_tags.extend(tags)
        return all_tags

    def get_mesh_info(self):
        """Get mesh statistics.

        Returns
        -------
        dict
        """
        import gmsh
        self._check_initialized()

        if not self._meshed:
            return {'meshed': False}

        n_hex = n_tet = n_wedge = 0

        node_tags, _, _ = gmsh.model.mesh.getNodes()
        n_nodes = len(node_tags)

        elem_types, _, _ = gmsh.model.mesh.getElements(3)
        for et in elem_types:
            name = gmsh.model.mesh.getElementProperties(et)[0]
            tags, _ = gmsh.model.mesh.getElementsByType(et)
            count = len(tags)
            if 'Hexahedron' in name:
                n_hex += count
            elif 'Tetrahedron' in name:
                n_tet += count
            elif 'Prism' in name or 'Wedge' in name:
                n_wedge += count

        return {
            'meshed': True,
            'n_nodes': n_nodes,
            'n_hex': n_hex,
            'n_tet': n_tet,
            'n_wedge': n_wedge,
            'n_total': n_hex + n_tet + n_wedge,
            'n_volumes': sum(len(t) for t in self._volumes.values()),
        }

    def get_bounding_box(self, vol_id):
        """Get bounding box of a volume (Gmsh 'get_bounding_box').

        Parameters
        ----------
        vol_id : int
            Volume ID.

        Returns
        -------
        dict
            xmin, ymin, zmin, xmax, ymax, zmax, dx, dy, dz.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        xmin, ymin, zmin = float('inf'), float('inf'), float('inf')
        xmax, ymax, zmax = float('-inf'), float('-inf'), float('-inf')
        for vt in self._volumes[vol_id]:
            bb = gmsh.model.occ.getBoundingBox(3, vt)
            xmin = min(xmin, bb[0])
            ymin = min(ymin, bb[1])
            zmin = min(zmin, bb[2])
            xmax = max(xmax, bb[3])
            ymax = max(ymax, bb[4])
            zmax = max(zmax, bb[5])

        return {
            'xmin': xmin, 'ymin': ymin, 'zmin': zmin,
            'xmax': xmax, 'ymax': ymax, 'zmax': zmax,
            'dx': xmax - xmin, 'dy': ymax - ymin, 'dz': zmax - zmin,
        }

    def get_center(self, vol_id):
        """Get center of mass of a volume (Gmsh 'get_center_point').

        For multi-tag volumes, computes the volume-weighted center of mass
        across all GMSH tags.

        Parameters
        ----------
        vol_id : int
            Volume ID.

        Returns
        -------
        list
            [x, y, z] center of mass.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        tags = self._volumes[vol_id]
        if len(tags) == 1:
            x, y, z = gmsh.model.occ.getCenterOfMass(3, tags[0])
        else:
            total_mass = 0
            cx, cy, cz = 0, 0, 0
            for t in tags:
                m = gmsh.model.occ.getMass(3, t)
                x, y, z = gmsh.model.occ.getCenterOfMass(3, t)
                cx += m * x; cy += m * y; cz += m * z
                total_mass += m
            if total_mass > 0:
                x, y, z = cx / total_mass, cy / total_mass, cz / total_mass
            else:
                x, y, z = gmsh.model.occ.getCenterOfMass(3, tags[0])
        return [x, y, z]

    def get_volume_measure(self, vol_id):
        """Get volume of a volume (Gmsh 'get_volume_volume').

        Parameters
        ----------
        vol_id : int
            Volume ID.

        Returns
        -------
        float
            Volume in cubic meters.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        total = 0
        for vt in self._volumes[vol_id]:
            total += gmsh.model.occ.getMass(3, vt)
        return total

    def get_surface_area(self, surf_tag):
        """Get area of a surface.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.

        Returns
        -------
        float
        """
        import gmsh
        self._check_initialized()

        return gmsh.model.occ.getMass(2, surf_tag)

    def get_curve_length(self, curve_tag):
        """Get length of a curve.

        Parameters
        ----------
        curve_tag : int
            GMSH curve tag.

        Returns
        -------
        float
        """
        import gmsh
        self._check_initialized()

        return gmsh.model.occ.getMass(1, curve_tag)

    def get_surfaces(self, vol_id):
        """Get surface tags of a volume (Gmsh 'vol.surfaces()').

        Parameters
        ----------
        vol_id : int
            Volume ID.

        Returns
        -------
        list of int
            Surface tags.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        surfs = set()
        for vt in self._volumes[vol_id]:
            boundary = gmsh.model.getBoundary([(3, vt)],
                                               oriented=False, recursive=False)
            for dt in boundary:
                if dt[0] == 2:
                    surfs.add(abs(dt[1]))

        return sorted(surfs)

    def get_curves(self, surf_tag):
        """Get curve tags of a surface (Gmsh 'surf.curves()').

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.

        Returns
        -------
        list of int
            Curve tags.
        """
        import gmsh
        self._check_initialized()

        boundary = gmsh.model.getBoundary([(2, surf_tag)],
                                           oriented=False, recursive=False)
        return sorted(abs(dt[1]) for dt in boundary if dt[0] == 1)

    def get_vertices_of(self, curve_tag):
        """Get vertex points of a curve.

        Parameters
        ----------
        curve_tag : int
            GMSH curve tag.

        Returns
        -------
        list of int
            Point tags.
        """
        import gmsh
        self._check_initialized()

        boundary = gmsh.model.getBoundary([(1, curve_tag)],
                                           oriented=False, recursive=False)
        return sorted(abs(dt[1]) for dt in boundary if dt[0] == 0)

    def get_surface_normal(self, surf_tag, uv=None):
        """Get surface normal at a point.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.
        uv : list, optional
            Parametric coordinates [u, v]. Default: center.

        Returns
        -------
        list
            [nx, ny, nz] unit normal.
        """
        import gmsh
        self._check_initialized()

        if uv is None:
            bounds = gmsh.model.getParametrizationBounds(2, surf_tag)
            u = (bounds[0][0] + bounds[1][0]) / 2
            v = (bounds[0][1] + bounds[1][1]) / 2
        else:
            u, v = uv

        normal = gmsh.model.getNormal(surf_tag, [u, v])
        return list(normal)

    def get_surface_type(self, surf_tag):
        """Get surface type (plane, cylinder, etc.).

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.

        Returns
        -------
        str
            Surface type name.
        """
        import gmsh
        self._check_initialized()

        return gmsh.model.getType(2, surf_tag)

    def is_planar(self, surf_tag):
        """Check if a surface is planar.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.

        Returns
        -------
        bool
        """
        stype = self.get_surface_type(surf_tag)
        return 'Plane' in stype

    def is_meshed(self, vol_id=None):
        """Check if geometry/volume is meshed.

        Parameters
        ----------
        vol_id : int, optional
            Volume ID. If None, checks global mesh state.

        Returns
        -------
        bool
        """
        return self._meshed

    def get_adjacent_volumes(self, vol_id):
        """Get volumes adjacent to a given volume (sharing a surface).

        Parameters
        ----------
        vol_id : int
            Volume ID.

        Returns
        -------
        list of int
            Adjacent volume IDs.
        """
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        my_surfs = set(self.get_surfaces(vol_id))
        adjacent = []

        for other_id in self._volumes:
            if other_id == vol_id:
                continue
            other_surfs = set(self.get_surfaces(other_id))
            if my_surfs & other_surfs:
                adjacent.append(other_id)

        return adjacent

    def get_shared_surfaces(self, vol_id1, vol_id2):
        """Get surfaces shared between two volumes.

        Parameters
        ----------
        vol_id1, vol_id2 : int
            Volume IDs.

        Returns
        -------
        list of int
            Shared surface tags.
        """
        s1 = set(self.get_surfaces(vol_id1))
        s2 = set(self.get_surfaces(vol_id2))
        return sorted(s1 & s2)

    def get_point_on_curve(self, curve_tag, t):
        """Get point on a curve at parameter t.

        Parameters
        ----------
        curve_tag : int
            GMSH curve tag.
        t : float
            Parameter value (0 to 1 normalized).

        Returns
        -------
        list
            [x, y, z].
        """
        import gmsh
        self._check_initialized()

        bounds = gmsh.model.getParametrizationBounds(1, curve_tag)
        t_min, t_max = bounds[0][0], bounds[1][0]
        t_actual = t_min + t * (t_max - t_min)

        point = gmsh.model.getValue(1, curve_tag, [t_actual])
        return list(point)

    def get_distance(self, vol_id1, vol_id2):
        """Get minimum distance between two volumes.

        Parameters
        ----------
        vol_id1, vol_id2 : int
            Volume IDs.

        Returns
        -------
        float
        """
        c1 = self.get_center(vol_id1)
        c2 = self.get_center(vol_id2)
        return math.sqrt(sum((a - b) ** 2 for a, b in zip(c1, c2)))

    def list_entities(self, dim=None):
        """List all entities (Gmsh 'list volume').

        Parameters
        ----------
        dim : int, optional
            Filter by dimension.

        Returns
        -------
        dict
            Entity counts and IDs by dimension.
        """
        import gmsh
        self._check_initialized()

        result = {}
        for d in ([dim] if dim is not None else [0, 1, 2, 3]):
            entities = gmsh.model.getEntities(d)
            result[d] = {
                'count': len(entities),
                'tags': [e[1] for e in entities],
            }

        return result

    def get_entity_count(self, dim=3):
        """Get count of entities of given dimension.

        Parameters
        ----------
        dim : int
            Entity dimension.

        Returns
        -------
        int
        """
        import gmsh
        self._check_initialized()

        return len(gmsh.model.getEntities(dim))

    def find_volumes_by_name(self, pattern):
        """Find volumes by name pattern.

        Parameters
        ----------
        pattern : str
            Name or substring to search for.

        Returns
        -------
        list of int
            Matching volume IDs.
        """
        result = []
        for vid, name in self._names.items():
            if pattern in name:
                result.append(vid)
        return result

    def set_name(self, vol_id, name):
        """Set a name for a volume (Gmsh 'set_entity_name').

        Parameters
        ----------
        vol_id : int
            Volume ID.
        name : str
            Name.
        """
        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        self._names[vol_id] = name
        if self._verbose:
            print(f"[GmshBuilder] set_name vol={vol_id} name='{name}'")

    def get_name(self, vol_id):
        """Get name of a volume.

        Parameters
        ----------
        vol_id : int
            Volume ID.

        Returns
        -------
        str or None
        """
        return self._names.get(vol_id, None)

    # ------------------------------------------------------------------
    # Extended query methods
    # ------------------------------------------------------------------

    def get_all_surfaces(self):
        """Get all surface tags in the model.

        Returns
        -------
        list of int
            Sorted list of all surface tags.
        """
        import gmsh
        self._check_initialized()

        entities = gmsh.model.getEntities(2)
        return sorted(e[1] for e in entities)

    def get_all_curves(self):
        """Get all curve tags in the model.

        Returns
        -------
        list of int
            Sorted list of all curve tags.
        """
        import gmsh
        self._check_initialized()

        entities = gmsh.model.getEntities(1)
        return sorted(e[1] for e in entities)

    def get_all_points(self):
        """Get all point tags in the model.

        Returns
        -------
        list of int
            Sorted list of all point tags.
        """
        import gmsh
        self._check_initialized()

        entities = gmsh.model.getEntities(0)
        return sorted(e[1] for e in entities)

    def get_bounding_box_all(self):
        """Get bounding box of all volumes combined.

        Returns
        -------
        dict
            xmin, ymin, zmin, xmax, ymax, zmax, dx, dy, dz.
        """
        import gmsh
        self._check_initialized()

        xmin, ymin, zmin = float('inf'), float('inf'), float('inf')
        xmax, ymax, zmax = float('-inf'), float('-inf'), float('-inf')

        for vol_tags in self._volumes.values():
            for vt in vol_tags:
                bb = gmsh.model.occ.getBoundingBox(3, vt)
                xmin = min(xmin, bb[0])
                ymin = min(ymin, bb[1])
                zmin = min(zmin, bb[2])
                xmax = max(xmax, bb[3])
                ymax = max(ymax, bb[4])
                zmax = max(zmax, bb[5])

        return {
            'xmin': xmin, 'ymin': ymin, 'zmin': zmin,
            'xmax': xmax, 'ymax': ymax, 'zmax': zmax,
            'dx': xmax - xmin, 'dy': ymax - ymin, 'dz': zmax - zmin,
        }

    def get_center_of_mass(self, vol_id):
        """Get center of mass of a volume (uniform density).

        Uses gmsh.model.occ.getCenterOfMass to compute the geometric
        centroid of the volume. For multi-tag volumes, computes the
        volume-weighted center of mass across all GMSH tags.

        Parameters
        ----------
        vol_id : int
            Volume ID.

        Returns
        -------
        list
            [x, y, z] center of mass coordinates.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        tags = self._volumes[vol_id]
        if len(tags) == 1:
            x, y, z = gmsh.model.occ.getCenterOfMass(3, tags[0])
        else:
            total_mass = 0
            cx, cy, cz = 0, 0, 0
            for t in tags:
                m = gmsh.model.occ.getMass(3, t)
                x, y, z = gmsh.model.occ.getCenterOfMass(3, t)
                cx += m * x; cy += m * y; cz += m * z
                total_mass += m
            if total_mass > 0:
                x, y, z = cx / total_mass, cy / total_mass, cz / total_mass
            else:
                x, y, z = gmsh.model.occ.getCenterOfMass(3, tags[0])
        return [x, y, z]

    def get_mass(self, vol_id, density=1.0):
        """Get mass of a volume (volume times density).

        Parameters
        ----------
        vol_id : int
            Volume ID.
        density : float, optional
            Material density in kg/m^3. Default is 1.0.

        Returns
        -------
        float
            Mass in kg (for density in kg/m^3 and volume in m^3).
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        total_volume = 0.0
        for vt in self._volumes[vol_id]:
            total_volume += gmsh.model.occ.getMass(3, vt)
        return total_volume * density

    def get_inertia(self, vol_id):
        """Get moment of inertia matrix of a volume.

        Uses gmsh.model.occ.getMass on the volume entity. Note that
        gmsh.model.occ.getMass(3, tag) returns the volume (3D measure).
        The inertia matrix is computed from the volume geometry assuming
        uniform unit density.

        Parameters
        ----------
        vol_id : int
            Volume ID.

        Returns
        -------
        list
            [Ixx, Ixy, Ixz, Iyx, Iyy, Iyz, Izx, Izy, Izz] inertia
            matrix elements (row-major, 3x3).
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        tag = self._volumes[vol_id][0]
        mat = gmsh.model.occ.getMatrixOfInertia(3, tag)
        return list(mat)

    def get_surface_count(self, vol_id):
        """Count surfaces of a volume.

        Parameters
        ----------
        vol_id : int
            Volume ID.

        Returns
        -------
        int
            Number of surfaces bounding the volume.
        """
        return len(self.get_surfaces(vol_id))

    def get_curve_count(self, surf_tag):
        """Count curves of a surface.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.

        Returns
        -------
        int
            Number of curves bounding the surface.
        """
        return len(self.get_curves(surf_tag))

    def get_parametric_bounds(self, dim, tag):
        """Get parametric bounds of an entity.

        Parameters
        ----------
        dim : int
            Entity dimension (1 for curve, 2 for surface).
        tag : int
            Entity tag.

        Returns
        -------
        dict
            For curves: {'t_min': float, 't_max': float}.
            For surfaces: {'u_min': float, 'u_max': float,
                           'v_min': float, 'v_max': float}.
        """
        import gmsh
        self._check_initialized()

        bounds = gmsh.model.getParametrizationBounds(dim, tag)
        if dim == 1:
            return {'t_min': bounds[0][0], 't_max': bounds[1][0]}
        elif dim == 2:
            return {
                'u_min': bounds[0][0], 'v_min': bounds[0][1],
                'u_max': bounds[1][0], 'v_max': bounds[1][1],
            }
        else:
            raise ValueError(f"Parametric bounds not available for dim={dim}")

    def point_on_surface(self, surf_tag, u, v):
        """Get 3D point from parametric coordinates on a surface.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.
        u : float
            Parametric u coordinate.
        v : float
            Parametric v coordinate.

        Returns
        -------
        list
            [x, y, z] coordinates.
        """
        import gmsh
        self._check_initialized()

        bounds = gmsh.model.getParametrizationBounds(2, surf_tag)
        u_actual = bounds[0][0] + u * (bounds[1][0] - bounds[0][0])
        v_actual = bounds[0][1] + v * (bounds[1][1] - bounds[0][1])
        point = gmsh.model.getValue(2, surf_tag, [u_actual, v_actual])
        return list(point)

    def is_curve_closed(self, curve_tag):
        """Check if a curve forms a closed loop.

        A curve is closed when its start and end vertices coincide
        or when it has no boundary vertices (periodic curve).

        Parameters
        ----------
        curve_tag : int
            GMSH curve tag.

        Returns
        -------
        bool
            True if the curve is closed.
        """
        import gmsh
        self._check_initialized()

        boundary = gmsh.model.getBoundary([(1, curve_tag)],
                                           oriented=False, recursive=False)
        points = [abs(dt[1]) for dt in boundary if dt[0] == 0]

        if len(points) == 0:
            return True
        if len(points) == 2 and points[0] == points[1]:
            return True
        return False

    def get_curve_tangent(self, curve_tag, t):
        """Get tangent vector at parameter t on a curve.

        The tangent is computed by evaluating the curve derivative
        via gmsh.model.getDerivative.

        Parameters
        ----------
        curve_tag : int
            GMSH curve tag.
        t : float
            Normalized parameter (0 to 1).

        Returns
        -------
        list
            [tx, ty, tz] tangent vector (not necessarily unit length).
        """
        import gmsh
        self._check_initialized()

        bounds = gmsh.model.getParametrizationBounds(1, curve_tag)
        t_min, t_max = bounds[0][0], bounds[1][0]
        t_actual = t_min + t * (t_max - t_min)

        deriv = gmsh.model.getDerivative(1, curve_tag, [t_actual])
        return list(deriv[:3])

    def get_surface_curvature(self, surf_tag, u, v):
        """Get principal curvatures at a parametric point on a surface.

        Uses gmsh.model.getCurvature to obtain the mean curvature,
        and gmsh.model.getPrincipalCurvatures for the two principal
        curvature values and directions.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.
        u : float
            Parametric u coordinate.
        v : float
            Parametric v coordinate.

        Returns
        -------
        dict
            'curvatures': [k1, k2] principal curvatures,
            'directions': [[d1x, d1y, d1z], [d2x, d2y, d2z]]
            principal directions.
        """
        import gmsh
        self._check_initialized()

        bounds = gmsh.model.getParametrizationBounds(2, surf_tag)
        u_actual = bounds[0][0] + u * (bounds[1][0] - bounds[0][0])
        v_actual = bounds[0][1] + v * (bounds[1][1] - bounds[0][1])
        curvatureMax, curvatureMin, directionMax, directionMin = \
            gmsh.model.getPrincipalCurvatures(surf_tag, [u_actual, v_actual])
        return {
            'curvatures': [float(curvatureMax[0]), float(curvatureMin[0])],
            'directions': [list(directionMax[:3]), list(directionMin[:3])],
        }

    def get_volume_topology(self, vol_id):
        """Get full topology tree of a volume.

        Returns a nested structure: the volume contains surfaces,
        each surface contains curves, each curve has point endpoints.

        Parameters
        ----------
        vol_id : int
            Volume ID.

        Returns
        -------
        dict
            {'vol_id': int, 'surfaces': [{'tag': int, 'curves':
            [{'tag': int, 'points': [int, ...]}]}]}.
        """
        self._check_initialized()

        surfs = self.get_surfaces(vol_id)
        surface_list = []
        for s in surfs:
            curves = self.get_curves(s)
            curve_list = []
            for c in curves:
                pts = self.get_vertices_of(c)
                curve_list.append({'tag': c, 'points': pts})
            surface_list.append({'tag': s, 'curves': curve_list})

        return {'vol_id': vol_id, 'surfaces': surface_list}

    def get_model_summary(self):
        """Return summary dict with entity counts per dimension.

        Returns
        -------
        dict
            {'points': int, 'curves': int, 'surfaces': int,
             'volumes': int, 'managed_volumes': int}.
        """
        import gmsh
        self._check_initialized()

        summary = {}
        dim_names = {0: 'points', 1: 'curves', 2: 'surfaces', 3: 'volumes'}
        for d, name in dim_names.items():
            summary[name] = len(gmsh.model.getEntities(d))
        summary['managed_volumes'] = len(self._volumes)
        return summary

    def get_entities(self, dim):
        """Get all entity tags of a given dimension.

        Parameters
        ----------
        dim : int
            Entity dimension (0=points, 1=curves, 2=surfaces, 3=volumes).

        Returns
        -------
        list of int
            Sorted list of entity tags.
        """
        import gmsh
        self._check_initialized()

        entities = gmsh.model.getEntities(dim)
        return sorted(e[1] for e in entities)

    def find_closest_entity(self, point, dim=3):
        """Find the closest entity of given dimension to a point.

        For dim=3, searches managed volumes by center-of-mass distance.
        For dim=2/1, uses gmsh.model.getClosestPoint for exact distance.
        For dim=0, computes distance to each point entity.

        Parameters
        ----------
        point : list
            [x, y, z] coordinates of the query point.
        dim : int, optional
            Entity dimension to search (default 3).

        Returns
        -------
        dict
            {'tag': int, 'distance': float, 'closest_point': [x, y, z]}.
        """
        import gmsh
        self._check_initialized()

        px, py, pz = point
        best_tag = None
        best_dist = float('inf')
        best_pt = None

        if dim == 3:
            for vol_id in self._volumes:
                center = self.get_center(vol_id)
                d = math.sqrt(sum((a - b) ** 2
                                  for a, b in zip(point, center)))
                if d < best_dist:
                    best_dist = d
                    best_tag = vol_id
                    best_pt = center

        elif dim == 2 or dim == 1:
            entities = gmsh.model.getEntities(dim)
            for _, tag in entities:
                cp = gmsh.model.getClosestPoint(dim, tag, point)
                closest = list(cp[0])
                d = math.sqrt(sum((a - b) ** 2
                                  for a, b in zip(point, closest)))
                if d < best_dist:
                    best_dist = d
                    best_tag = tag
                    best_pt = closest

        elif dim == 0:
            entities = gmsh.model.getEntities(0)
            for _, tag in entities:
                coords = gmsh.model.getValue(0, tag, [])
                ep = list(coords)
                d = math.sqrt(sum((a - b) ** 2
                                  for a, b in zip(point, ep)))
                if d < best_dist:
                    best_dist = d
                    best_tag = tag
                    best_pt = ep

        else:
            raise ValueError(f"Invalid dimension: {dim}")

        return {'tag': best_tag, 'distance': best_dist,
                'closest_point': best_pt}

    def get_named_entities(self):
        """Return dict of all named entities.

        Includes both volume names managed by GmshBuilder and any
        physical groups defined in the GMSH model.

        Returns
        -------
        dict
            {name: (dim, tag)} for all named entities.
        """
        import gmsh
        self._check_initialized()

        result = {}

        # Managed volume names
        for vol_id, name in self._names.items():
            if name:
                result[name] = (3, vol_id)

        # GMSH physical groups
        groups = gmsh.model.getPhysicalGroups()
        for dim, tag in groups:
            name = gmsh.model.getPhysicalName(dim, tag)
            if name:
                result[name] = (dim, tag)

        return result

    def measure_distance_point_to_surface(self, point, surf_tag):
        """Compute minimum distance from a point to a surface.

        Uses gmsh.model.getClosestPoint to find the closest point
        on the surface and returns the Euclidean distance.

        Parameters
        ----------
        point : list
            [x, y, z] coordinates of the query point.
        surf_tag : int
            GMSH surface tag.

        Returns
        -------
        dict
            {'distance': float, 'closest_point': [x, y, z],
             'parametric': [u, v]}.
        """
        import gmsh
        self._check_initialized()

        cp = gmsh.model.getClosestPoint(2, surf_tag, point)
        closest = list(cp[0])
        parametric = list(cp[1])
        d = math.sqrt(sum((a - b) ** 2 for a, b in zip(point, closest)))

        return {
            'distance': d,
            'closest_point': closest,
            'parametric': parametric,
        }

    # ------------------------------------------------------------------
    # Entity property queries
    # ------------------------------------------------------------------

    def get_surface_derivative(self, surf_tag, u=0.5, v=0.5):
        """Get first-order partial derivatives on a surface.

        Computes the partial derivatives of the surface parametrization
        at the given (u, v) coordinates.

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
            [du_x, du_y, du_z, dv_x, dv_y, dv_z] partial derivatives.
        """
        import gmsh
        self._check_initialized()

        bounds = gmsh.model.getParametrizationBounds(2, surf_tag)
        u_actual = bounds[0][0] + u * (bounds[1][0] - bounds[0][0])
        v_actual = bounds[0][1] + v * (bounds[1][1] - bounds[0][1])
        deriv = gmsh.model.getDerivative(2, surf_tag, [u_actual, v_actual])
        result = list(deriv[:6])

        if self._verbose:
            print(f"[GmshBuilder] get_surface_derivative surf={surf_tag} "
                  f"u={u} v={v} -> {result}")
        return result

    def get_second_derivative(self, dim, tag, params):
        """Get second-order derivatives of a parametric entity.

        Uses gmsh.model.getSecondDerivative if available in the
        current GMSH version.

        Parameters
        ----------
        dim : int
            Entity dimension (1 for curve, 2 for surface).
        tag : int
            Entity tag.
        params : list of float
            Parametric coordinates (e.g. [t] for curve, [u, v] for surface).

        Returns
        -------
        list of float
            Second derivative values. For a curve (dim=1), returns
            [d2x, d2y, d2z]. For a surface (dim=2), returns second
            derivatives [d2u_x, d2u_y, d2u_z, d2v_x, d2v_y, d2v_z,
            d2uv_x, d2uv_y, d2uv_z].

        Raises
        ------
        RuntimeError
            If gmsh.model.getSecondDerivative is not available.
        """
        import gmsh
        self._check_initialized()

        try:
            result = list(gmsh.model.getSecondDerivative(dim, tag, params))
        except (AttributeError, Exception) as e:
            raise RuntimeError(
                f"getSecondDerivative not available in this GMSH version: {e}"
            )

        if self._verbose:
            print(f"[GmshBuilder] get_second_derivative dim={dim} tag={tag} "
                  f"params={params} -> {result}")
        return result

    def get_parametrization(self, dim, tag, point):
        """Get parametric coordinates of a point on an entity.

        Given a 3D point on (or near) an entity, returns the corresponding
        parametric coordinates.

        Parameters
        ----------
        dim : int
            Entity dimension (1 for curve, 2 for surface).
        tag : int
            Entity tag.
        point : list of float
            [x, y, z] coordinates of the point.

        Returns
        -------
        list of float
            Parametric coordinates. For curves: [t]. For surfaces: [u, v].
        """
        import gmsh
        self._check_initialized()

        result = list(gmsh.model.getParametrization(dim, tag, point))

        if self._verbose:
            print(f"[GmshBuilder] get_parametrization dim={dim} tag={tag} "
                  f"point={point} -> {result}")
        return result

    def reparametrize_on_surface(self, dim, tag, surf_tag, params):
        """Reparametrize an entity on a surface.

        Given parametric coordinates on an entity (curve or point) that
        lies on a surface, returns the corresponding parametric
        coordinates on the surface.

        Parameters
        ----------
        dim : int
            Entity dimension (0 or 1).
        tag : int
            Entity tag.
        surf_tag : int
            Target surface tag.
        params : list of float
            Parametric coordinates on the entity.

        Returns
        -------
        list of float
            Parametric coordinates [u, v] on the surface.

        Raises
        ------
        RuntimeError
            If gmsh.model.reparametrizeOnSurface is not available.
        """
        import gmsh
        self._check_initialized()

        try:
            result = list(
                gmsh.model.reparametrizeOnSurface(dim, tag, params, surf_tag)
            )
        except (AttributeError, Exception) as e:
            raise RuntimeError(
                f"reparametrizeOnSurface not available in this GMSH "
                f"version: {e}"
            )

        if self._verbose:
            print(f"[GmshBuilder] reparametrize_on_surface dim={dim} "
                  f"tag={tag} surf={surf_tag} params={params} -> {result}")
        return result

    def get_surface_loops_query(self, vol_tag):
        """Get surface loops of a volume.

        Returns the surface loops that form the boundary of the given
        volume using the OCC kernel.

        Parameters
        ----------
        vol_tag : int
            GMSH volume tag.

        Returns
        -------
        list of list of int
            Each inner list contains the surface tags forming one loop.

        Raises
        ------
        RuntimeError
            If gmsh.model.occ.getSurfaceLoops is not available.
        """
        import gmsh
        self._check_initialized()

        try:
            result = gmsh.model.occ.getSurfaceLoops(vol_tag)
            if isinstance(result, tuple) and len(result) == 2:
                loop_tags, surf_arrays = result
                loops = [list(arr) for arr in surf_arrays]
            else:
                loops = [list(r) for r in result]
        except (AttributeError, Exception) as e:
            raise RuntimeError(
                f"getSurfaceLoops not available in this GMSH version: {e}"
            )

        if self._verbose:
            print(f"[GmshBuilder] get_surface_loops_query vol={vol_tag} "
                  f"-> {len(loops)} loops")
        return loops

    def get_curve_loops_query(self, surf_tag):
        """Get curve loops of a surface.

        Returns the curve loops that form the boundary of the given
        surface using the OCC kernel.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.

        Returns
        -------
        list of list of int
            Each inner list contains the curve tags forming one loop.

        Raises
        ------
        RuntimeError
            If gmsh.model.occ.getCurveLoops is not available.
        """
        import gmsh
        self._check_initialized()

        try:
            result = gmsh.model.occ.getCurveLoops(surf_tag)
            if isinstance(result, tuple) and len(result) == 2:
                loop_tags, curve_arrays = result
                loops = [list(arr) for arr in curve_arrays]
            else:
                loops = [list(r) for r in result]
        except (AttributeError, Exception) as e:
            raise RuntimeError(
                f"getCurveLoops not available in this GMSH version: {e}"
            )

        if self._verbose:
            print(f"[GmshBuilder] get_curve_loops_query surf={surf_tag} "
                  f"-> {len(loops)} loops")
        return loops

    def get_point_coordinates(self, point_tag):
        """Get 3D coordinates of a point entity.

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
        result = list(coords[:3])

        if self._verbose:
            print(f"[GmshBuilder] get_point_coordinates point={point_tag} "
                  f"-> {result}")
        return result

    def get_entity_partitions(self, dim, tag):
        """Get partition tags of an entity.

        Returns the partition(s) to which the entity belongs when
        the model has been partitioned.

        Parameters
        ----------
        dim : int
            Entity dimension.
        tag : int
            Entity tag.

        Returns
        -------
        list of int
            Partition tags.

        Raises
        ------
        RuntimeError
            If gmsh.model.getPartitions is not available.
        """
        import gmsh
        self._check_initialized()

        try:
            result = list(gmsh.model.getPartitions(dim, tag))
        except (AttributeError, Exception) as e:
            raise RuntimeError(
                f"getPartitions not available in this GMSH version: {e}"
            )

        if self._verbose:
            print(f"[GmshBuilder] get_entity_partitions dim={dim} tag={tag} "
                  f"-> {result}")
        return result

    # ------------------------------------------------------------------
    # Distance / proximity queries
    # ------------------------------------------------------------------

    def distance_point_to_curve(self, point, curve_tag):
        """Compute minimum distance from a point to a curve.

        Uses gmsh.model.getClosestPoint to find the closest point
        on the curve and returns distance information.

        Parameters
        ----------
        point : list of float
            [x, y, z] coordinates of the query point.
        curve_tag : int
            GMSH curve tag.

        Returns
        -------
        dict
            'distance': float, Euclidean distance to closest point.
            'closest_point': [x, y, z], coordinates of closest point.
            'parametric': float, parametric coordinate on the curve.
        """
        import gmsh
        self._check_initialized()

        cp = gmsh.model.getClosestPoint(1, curve_tag, point)
        closest = list(cp[0][:3])
        parametric = cp[1][0] if hasattr(cp[1], '__len__') else float(cp[1])
        d = math.sqrt(sum((a - b) ** 2 for a, b in zip(point, closest)))

        result = {
            'distance': d,
            'closest_point': closest,
            'parametric': parametric,
        }

        if self._verbose:
            print(f"[GmshBuilder] distance_point_to_curve curve={curve_tag} "
                  f"-> dist={d:.6e}")
        return result

    def distance_point_to_entity(self, point, dim, tag):
        """Compute Euclidean distance from a point to an entity.

        Uses gmsh.model.getClosestPoint for curves (dim=1) and
        surfaces (dim=2). For points (dim=0), computes distance
        directly. For volumes (dim=3), uses center of mass as
        approximation.

        Parameters
        ----------
        point : list of float
            [x, y, z] coordinates of the query point.
        dim : int
            Entity dimension.
        tag : int
            Entity tag.

        Returns
        -------
        float
            Euclidean distance from the point to the entity.
        """
        import gmsh
        self._check_initialized()

        if dim == 0:
            coords = gmsh.model.getValue(0, tag, [])
            closest = list(coords[:3])
        elif dim == 1 or dim == 2:
            cp = gmsh.model.getClosestPoint(dim, tag, point)
            closest = list(cp[0][:3])
        elif dim == 3:
            x, y, z = gmsh.model.occ.getCenterOfMass(3, tag)
            closest = [x, y, z]
        else:
            raise ValueError(f"Invalid dimension: {dim}")

        d = math.sqrt(sum((a - b) ** 2 for a, b in zip(point, closest)))

        if self._verbose:
            print(f"[GmshBuilder] distance_point_to_entity dim={dim} "
                  f"tag={tag} -> dist={d:.6e}")
        return d

    def distance_between_volumes(self, vol_id1, vol_id2):
        """Compute minimum distance between two volumes.

        Attempts to use gmsh.model.occ.getDistance for exact minimum
        distance. Falls back to center-of-mass distance if the API
        is not available.

        Parameters
        ----------
        vol_id1 : int
            First volume ID.
        vol_id2 : int
            Second volume ID.

        Returns
        -------
        float
            Minimum distance between the two volumes.
        """
        import gmsh
        self._check_initialized()

        if vol_id1 not in self._volumes:
            raise KeyError(f"Volume {vol_id1} not found")
        if vol_id2 not in self._volumes:
            raise KeyError(f"Volume {vol_id2} not found")

        tag1 = self._volumes[vol_id1][0]
        tag2 = self._volumes[vol_id2][0]

        try:
            d = gmsh.model.occ.getDistance(3, tag1, 3, tag2)
            if isinstance(d, (list, tuple)):
                d = d[0]
            d = float(d)
        except (AttributeError, Exception):
            c1 = self.get_center(vol_id1)
            c2 = self.get_center(vol_id2)
            d = math.sqrt(sum((a - b) ** 2 for a, b in zip(c1, c2)))

        if self._verbose:
            print(f"[GmshBuilder] distance_between_volumes "
                  f"vol1={vol_id1} vol2={vol_id2} -> dist={d:.6e}")
        return d

    def closest_point_on_curve(self, curve_tag, point):
        """Find the closest point on a curve to a given point.

        Parameters
        ----------
        curve_tag : int
            GMSH curve tag.
        point : list of float
            [x, y, z] coordinates of the query point.

        Returns
        -------
        dict
            'point': [x, y, z], coordinates of closest point on curve.
            'parametric': float, parametric coordinate on the curve.
            'distance': float, Euclidean distance.
        """
        import gmsh
        self._check_initialized()

        cp = gmsh.model.getClosestPoint(1, curve_tag, point)
        closest = list(cp[0][:3])
        parametric = cp[1][0] if hasattr(cp[1], '__len__') else float(cp[1])
        d = math.sqrt(sum((a - b) ** 2 for a, b in zip(point, closest)))

        result = {
            'point': closest,
            'parametric': parametric,
            'distance': d,
        }

        if self._verbose:
            print(f"[GmshBuilder] closest_point_on_curve curve={curve_tag} "
                  f"-> dist={d:.6e}")
        return result

    def closest_point_on_entity(self, dim, tag, point):
        """Find the closest point on an entity to a given point.

        Supports curves (dim=1) and surfaces (dim=2).

        Parameters
        ----------
        dim : int
            Entity dimension (1 or 2).
        tag : int
            Entity tag.
        point : list of float
            [x, y, z] coordinates of the query point.

        Returns
        -------
        dict
            'point': [x, y, z], coordinates of closest point.
            'parametric': list of float, parametric coordinates.
            'distance': float, Euclidean distance.
        """
        import gmsh
        self._check_initialized()

        if dim not in (1, 2):
            raise ValueError(
                f"getClosestPoint supports dim=1 or dim=2, got {dim}"
            )

        cp = gmsh.model.getClosestPoint(dim, tag, point)
        closest = list(cp[0][:3])
        parametric = list(cp[1])
        d = math.sqrt(sum((a - b) ** 2 for a, b in zip(point, closest)))

        result = {
            'point': closest,
            'parametric': parametric,
            'distance': d,
        }

        if self._verbose:
            print(f"[GmshBuilder] closest_point_on_entity dim={dim} "
                  f"tag={tag} -> dist={d:.6e}")
        return result

    def is_point_inside_volume(self, vol_tag, point):
        """Check if a point is inside a volume.

        Uses gmsh.model.isInside to test containment.

        Parameters
        ----------
        vol_tag : int
            GMSH volume tag.
        point : list of float
            [x, y, z] coordinates.

        Returns
        -------
        bool
            True if the point is inside the volume.
        """
        import gmsh
        self._check_initialized()

        result = gmsh.model.isInside(3, vol_tag, point)

        if self._verbose:
            print(f"[GmshBuilder] is_point_inside_volume vol={vol_tag} "
                  f"point={point} -> {result}")
        return bool(result)

    # ------------------------------------------------------------------
    # Volume aggregates
    # ------------------------------------------------------------------

    def get_total_volume(self):
        """Get total volume of all volume entities in the model.

        Sums gmsh.model.occ.getMass(3, tag) for all 3D entities.

        Returns
        -------
        float
            Total volume in cubic length units.
        """
        import gmsh
        self._check_initialized()

        entities = gmsh.model.getEntities(3)
        total = 0.0
        for _, tag in entities:
            total += gmsh.model.occ.getMass(3, tag)

        if self._verbose:
            print(f"[GmshBuilder] get_total_volume -> {total:.6e} "
                  f"({len(entities)} volumes)")
        return total

    def get_total_surface_area(self):
        """Get total surface area of all surface entities in the model.

        Sums gmsh.model.occ.getMass(2, tag) for all 2D entities.

        Returns
        -------
        float
            Total surface area in square length units.
        """
        import gmsh
        self._check_initialized()

        entities = gmsh.model.getEntities(2)
        total = 0.0
        for _, tag in entities:
            total += gmsh.model.occ.getMass(2, tag)

        if self._verbose:
            print(f"[GmshBuilder] get_total_surface_area -> {total:.6e} "
                  f"({len(entities)} surfaces)")
        return total

    def get_total_curve_length(self):
        """Get total curve length of all curve entities in the model.

        Sums gmsh.model.occ.getMass(1, tag) for all 1D entities.

        Returns
        -------
        float
            Total curve length in length units.
        """
        import gmsh
        self._check_initialized()

        entities = gmsh.model.getEntities(1)
        total = 0.0
        for _, tag in entities:
            total += gmsh.model.occ.getMass(1, tag)

        if self._verbose:
            print(f"[GmshBuilder] get_total_curve_length -> {total:.6e} "
                  f"({len(entities)} curves)")
        return total

    def get_volumes_in_bounding_box(self, xmin, ymin, zmin,
                                    xmax, ymax, zmax):
        """Get volume tags within a bounding box.

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
            Volume tags within the bounding box.
        """
        import gmsh
        self._check_initialized()

        entities = gmsh.model.getEntitiesInBoundingBox(
            xmin, ymin, zmin, xmax, ymax, zmax, 3
        )
        result = sorted(e[1] for e in entities)

        if self._verbose:
            print(f"[GmshBuilder] get_volumes_in_bounding_box "
                  f"[{xmin},{ymin},{zmin}]-[{xmax},{ymax},{zmax}] "
                  f"-> {len(result)} volumes")
        return result

    def get_surfaces_in_bounding_box(self, xmin, ymin, zmin,
                                     xmax, ymax, zmax):
        """Get surface tags within a bounding box.

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
            Surface tags within the bounding box.
        """
        import gmsh
        self._check_initialized()

        entities = gmsh.model.getEntitiesInBoundingBox(
            xmin, ymin, zmin, xmax, ymax, zmax, 2
        )
        result = sorted(e[1] for e in entities)

        if self._verbose:
            print(f"[GmshBuilder] get_surfaces_in_bounding_box "
                  f"[{xmin},{ymin},{zmin}]-[{xmax},{ymax},{zmax}] "
                  f"-> {len(result)} surfaces")
        return result

    def get_curves_in_bounding_box(self, xmin, ymin, zmin,
                                   xmax, ymax, zmax):
        """Get curve tags within a bounding box.

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
            Curve tags within the bounding box.
        """
        import gmsh
        self._check_initialized()

        entities = gmsh.model.getEntitiesInBoundingBox(
            xmin, ymin, zmin, xmax, ymax, zmax, 1
        )
        result = sorted(e[1] for e in entities)

        if self._verbose:
            print(f"[GmshBuilder] get_curves_in_bounding_box "
                  f"[{xmin},{ymin},{zmin}]-[{xmax},{ymax},{zmax}] "
                  f"-> {len(result)} curves")
        return result

    # ------------------------------------------------------------------
    # Model attribute system
    # ------------------------------------------------------------------

    def set_attribute(self, dim, tag, name, values):
        """Set an attribute on an entity.

        Uses gmsh.model.setAttribute to attach named string data
        to a model entity.

        Parameters
        ----------
        dim : int
            Entity dimension.
        tag : int
            Entity tag.
        name : str
            Attribute name.
        values : list
            Attribute values (will be converted to strings).

        Raises
        ------
        RuntimeError
            If gmsh.model.setAttribute is not available.
        """
        import gmsh
        self._check_initialized()

        try:
            attr_name = f"{dim}_{tag}_{name}"
            gmsh.model.setAttribute(attr_name, [str(v) for v in values])
        except (AttributeError, Exception) as e:
            raise RuntimeError(
                f"setAttribute not available in this GMSH version: {e}"
            )

        if self._verbose:
            print(f"[GmshBuilder] set_attribute dim={dim} tag={tag} "
                  f"name='{name}' values={values}")

    def get_attribute(self, dim, tag, name):
        """Get an attribute from an entity.

        Uses gmsh.model.getAttribute to retrieve named string data
        from a model entity.

        Parameters
        ----------
        dim : int
            Entity dimension.
        tag : int
            Entity tag.
        name : str
            Attribute name.

        Returns
        -------
        list of str
            Attribute values.

        Raises
        ------
        RuntimeError
            If gmsh.model.getAttribute is not available.
        """
        import gmsh
        self._check_initialized()

        try:
            attr_name = f"{dim}_{tag}_{name}"
            result = list(gmsh.model.getAttribute(attr_name))
        except (AttributeError, Exception) as e:
            raise RuntimeError(
                f"getAttribute not available in this GMSH version: {e}"
            )

        if self._verbose:
            print(f"[GmshBuilder] get_attribute dim={dim} tag={tag} "
                  f"name='{name}' -> {result}")
        return result

    def get_attribute_names(self, dim, tag):
        """Get all attribute names on an entity.

        Parameters
        ----------
        dim : int
            Entity dimension.
        tag : int
            Entity tag.

        Returns
        -------
        list of str
            Attribute names.

        Raises
        ------
        RuntimeError
            If gmsh.model.getAttributeNames is not available.
        """
        import gmsh
        self._check_initialized()

        try:
            all_names = gmsh.model.getAttributeNames()
            prefix = f"{dim}_{tag}_"
            result = [n[len(prefix):] for n in all_names
                      if n.startswith(prefix)]
        except (AttributeError, Exception) as e:
            raise RuntimeError(
                f"getAttributeNames not available in this GMSH "
                f"version: {e}"
            )

        if self._verbose:
            print(f"[GmshBuilder] get_attribute_names dim={dim} tag={tag} "
                  f"-> {result}")
        return result

    def remove_attribute(self, dim, tag, name):
        """Remove an attribute from an entity.

        Parameters
        ----------
        dim : int
            Entity dimension.
        tag : int
            Entity tag.
        name : str
            Attribute name to remove.

        Raises
        ------
        RuntimeError
            If gmsh.model.removeAttribute is not available.
        """
        import gmsh
        self._check_initialized()

        try:
            attr_name = f"{dim}_{tag}_{name}"
            gmsh.model.removeAttribute(attr_name)
        except (AttributeError, Exception) as e:
            raise RuntimeError(
                f"removeAttribute not available in this GMSH version: {e}"
            )

        if self._verbose:
            print(f"[GmshBuilder] remove_attribute dim={dim} tag={tag} "
                  f"name='{name}'")

    def set_model_filename(self, filename):
        """Set the filename associated with the current model.

        Parameters
        ----------
        filename : str
            File path to associate with the model.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.setFileName(filename)

        if self._verbose:
            print(f"[GmshBuilder] set_model_filename -> '{filename}'")

    def get_model_filename(self):
        """Get the filename associated with the current model.

        Returns
        -------
        str
            File path associated with the model.
        """
        import gmsh
        self._check_initialized()

        result = gmsh.model.getFileName()

        if self._verbose:
            print(f"[GmshBuilder] get_model_filename -> '{result}'")
        return result

    # ------------------------------------------------------------------
    # Topology traversal
    # ------------------------------------------------------------------

    def get_upward_adjacencies(self, dim, tag):
        """Get upward-adjacent entities (dim+1) of an entity.

        For example, for a surface (dim=2), returns the volumes
        (dim=3) that contain it.

        Parameters
        ----------
        dim : int
            Entity dimension (0, 1, or 2).
        tag : int
            Entity tag.

        Returns
        -------
        list of int
            Tags of adjacent entities of dimension dim+1.
        """
        import gmsh
        self._check_initialized()

        upward, _ = gmsh.model.getAdjacencies(dim, tag)
        result = list(upward)

        if self._verbose:
            print(f"[GmshBuilder] get_upward_adjacencies dim={dim} "
                  f"tag={tag} -> {result}")
        return result

    def get_downward_adjacencies(self, dim, tag):
        """Get downward-adjacent entities (dim-1) of an entity.

        For example, for a volume (dim=3), returns the surfaces
        (dim=2) that bound it.

        Parameters
        ----------
        dim : int
            Entity dimension (1, 2, or 3).
        tag : int
            Entity tag.

        Returns
        -------
        list of int
            Tags of adjacent entities of dimension dim-1.
        """
        import gmsh
        self._check_initialized()

        _, downward = gmsh.model.getAdjacencies(dim, tag)
        result = list(downward)

        if self._verbose:
            print(f"[GmshBuilder] get_downward_adjacencies dim={dim} "
                  f"tag={tag} -> {result}")
        return result

    def get_parent(self, dim, tag):
        """Get the parent entity of an entity.

        Returns the parent entity in the model hierarchy (e.g., the
        volume from which a surface was derived via boundary
        operations).

        Parameters
        ----------
        dim : int
            Entity dimension.
        tag : int
            Entity tag.

        Returns
        -------
        tuple or None
            (parent_dim, parent_tag) if a parent exists, None otherwise.

        Raises
        ------
        RuntimeError
            If gmsh.model.getParent is not available.
        """
        import gmsh
        self._check_initialized()

        try:
            parent_dim, parent_tag = gmsh.model.getParent(dim, tag)
            if parent_tag < 0:
                result = None
            else:
                result = (parent_dim, parent_tag)
        except (AttributeError, Exception) as e:
            raise RuntimeError(
                f"getParent not available in this GMSH version: {e}"
            )

        if self._verbose:
            print(f"[GmshBuilder] get_parent dim={dim} tag={tag} "
                  f"-> {result}")
        return result

    def get_number_of_partitions(self):
        """Get the number of partitions in the model.

        Returns
        -------
        int
            Number of partitions (0 if model is not partitioned).
        """
        import gmsh
        self._check_initialized()

        result = gmsh.model.getNumberOfPartitions()

        if self._verbose:
            print(f"[GmshBuilder] get_number_of_partitions -> {result}")
        return int(result)

    def get_changed_entities(self, dim=-1):
        """Get entities changed by the last geometry operation.

        Returns entity (dim, tag) tuples modified by the most recent
        OCC boolean, heal, or synchronize operation.

        Parameters
        ----------
        dim : int
            Entity dimension to filter by.  -1 returns all dimensions.

        Returns
        -------
        list of tuple
            List of (dim, tag) tuples for changed entities.
            Returns empty list if the API is not available.
        """
        import gmsh
        self._check_initialized()

        try:
            result = gmsh.model.getChangedEntities(dim)
            result = [tuple(dt) for dt in result]
        except (AttributeError, Exception):
            result = []

        if self._verbose:
            print(f"[GmshBuilder] get_changed_entities dim={dim} "
                  f"-> {len(result)} entities")
        return result
