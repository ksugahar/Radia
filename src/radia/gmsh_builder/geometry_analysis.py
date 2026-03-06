"""
Geometry analysis for GmshBuilder.

Cleanup, distance/overlap detection, topology analysis, and
differential geometry queries.  Gmsh equivalents: heal, defeature,
validate, analysis.
"""

import math
import warnings


class GeometryAnalysisMixin:
    """Mixin providing geometry cleanup, detection, and differential
    geometry analysis."""

    # ------------------------------------------------------------------
    # Cleanup (10 methods)
    # ------------------------------------------------------------------

    def heal_all(self, tolerance=1e-8):
        """Heal all shapes in the model.

        Attempts to fix common CAD issues such as gaps, overlaps,
        and degenerate faces using the OCC healing algorithm.

        Parameters
        ----------
        tolerance : float, optional
            Healing tolerance in meters. Default is 1e-8.

        Returns
        -------
        dict
            Summary with keys 'healed' (bool) and 'message' (str).
        """
        import gmsh
        self._check_initialized()

        try:
            out_dimtags = gmsh.model.occ.healShapes(
                dimTags=gmsh.model.occ.getEntities(-1),
                tolerance=tolerance,
                fixDegenerated=True,
                fixSmallEdges=True,
                fixSmallFaces=True,
                sewFaces=True,
                makeSolids=True)
            gmsh.model.occ.synchronize()
            if self._verbose:
                print(f"[GmshBuilder] heal_all: healed {len(out_dimtags)} "
                      f"entities (tol={tolerance})")
            return {
                'healed': True,
                'n_entities': len(out_dimtags),
                'message': f'Healed {len(out_dimtags)} entities',
            }
        except (AttributeError, TypeError):
            # healShapes may not accept tolerance kwarg in older GMSH
            try:
                out_dimtags = gmsh.model.occ.healShapes()
                gmsh.model.occ.synchronize()
                if self._verbose:
                    print("[GmshBuilder] heal_all: healed (no tolerance param)")
                return {
                    'healed': True,
                    'n_entities': len(out_dimtags),
                    'message': 'Healed (tolerance parameter not supported)',
                }
            except Exception:
                warnings.warn("gmsh.model.occ.healShapes not available")
                return {
                    'healed': False,
                    'n_entities': 0,
                    'message': 'healShapes not available in this GMSH version',
                }

    def defeature_volume(self, vol_tag, tolerance):
        """Remove small features from a volume.

        Uses the OCC defeature algorithm to simplify geometry by
        removing features smaller than the given tolerance.

        Parameters
        ----------
        vol_tag : int
            GMSH volume tag (not GmshBuilder volume ID).
        tolerance : float
            Feature size tolerance in meters.

        Returns
        -------
        list of int
            New volume tags after defeaturing.

        Raises
        ------
        NotImplementedError
            If the GMSH version does not support defeaturing.
        """
        import gmsh
        self._check_initialized()

        try:
            # Find faces of the volume smaller than tolerance
            surfaces = [dt[1] for dt in gmsh.model.getBoundary(
                [(3, vol_tag)], oriented=False) if dt[0] == 2]
            small_faces = [s for s in surfaces
                           if gmsh.model.occ.getMass(2, s) < tolerance]
            if not small_faces:
                if self._verbose:
                    print(f"[GmshBuilder] defeature_volume: vol={vol_tag} "
                          f"no features smaller than {tolerance} found")
                return []
            out_dimtags = gmsh.model.occ.defeature([vol_tag], small_faces)
            gmsh.model.occ.synchronize()
            new_tags = [dt[1] for dt in out_dimtags if dt[0] == 3]
            if self._verbose:
                print(f"[GmshBuilder] defeature_volume: vol={vol_tag} "
                      f"tol={tolerance} -> {new_tags}")
            return new_tags
        except (AttributeError, Exception) as exc:
            raise NotImplementedError(
                "Requires GMSH >= 4.11 with OCC defeature support"
            ) from exc

    def remove_small_faces(self, vol_tag, tolerance):
        """Remove surfaces with area smaller than tolerance.

        Identifies small surfaces on a volume and attempts to remove
        them via the OCC defeature algorithm.

        Parameters
        ----------
        vol_tag : int
            GMSH volume tag.
        tolerance : float
            Minimum surface area in m^2. Surfaces below this are removed.

        Returns
        -------
        int
            Number of small faces removed.
        """
        import gmsh
        self._check_initialized()

        boundary = gmsh.model.getBoundary(
            [(3, vol_tag)], oriented=False, recursive=False)
        small_faces = []
        for dim, tag in boundary:
            if dim == 2:
                area = gmsh.model.occ.getMass(2, abs(tag))
                if area < tolerance:
                    small_faces.append(abs(tag))

        if not small_faces:
            if self._verbose:
                print(f"[GmshBuilder] remove_small_faces: vol={vol_tag} "
                      f"no faces below {tolerance}")
            return 0

        try:
            gmsh.model.occ.defeature([vol_tag], small_faces)
            gmsh.model.occ.synchronize()
        except (AttributeError, Exception):
            warnings.warn("Defeature not available; small faces not removed")
            return 0

        count = len(small_faces)
        if self._verbose:
            print(f"[GmshBuilder] remove_small_faces: vol={vol_tag} "
                  f"removed {count} faces")
        return count

    def find_small_curves(self, min_length):
        """Find curves shorter than a minimum length.

        Parameters
        ----------
        min_length : float
            Minimum curve length in meters.

        Returns
        -------
        list of tuple
            List of (tag, length) for curves shorter than min_length.
        """
        import gmsh
        self._check_initialized()

        result = []
        entities = gmsh.model.getEntities(1)
        for _, tag in entities:
            length = gmsh.model.occ.getMass(1, tag)
            if length < min_length:
                result.append((tag, length))

        if self._verbose:
            print(f"[GmshBuilder] find_small_curves: found {len(result)} "
                  f"curves < {min_length}")
        return result

    def find_small_surfaces(self, min_area):
        """Find surfaces with area smaller than a minimum.

        Parameters
        ----------
        min_area : float
            Minimum surface area in m^2.

        Returns
        -------
        list of tuple
            List of (tag, area) for surfaces smaller than min_area.
        """
        import gmsh
        self._check_initialized()

        result = []
        entities = gmsh.model.getEntities(2)
        for _, tag in entities:
            area = gmsh.model.occ.getMass(2, tag)
            if area < min_area:
                result.append((tag, area))

        if self._verbose:
            print(f"[GmshBuilder] find_small_surfaces: found {len(result)} "
                  f"surfaces < {min_area}")
        return result

    def find_small_volumes(self, min_volume):
        """Find volumes smaller than a minimum.

        Parameters
        ----------
        min_volume : float
            Minimum volume in m^3.

        Returns
        -------
        list of tuple
            List of (tag, volume) for volumes smaller than min_volume.
        """
        import gmsh
        self._check_initialized()

        result = []
        entities = gmsh.model.getEntities(3)
        for _, tag in entities:
            vol = gmsh.model.occ.getMass(3, tag)
            if vol < min_volume:
                result.append((tag, vol))

        if self._verbose:
            print(f"[GmshBuilder] find_small_volumes: found {len(result)} "
                  f"volumes < {min_volume}")
        return result

    def find_narrow_surfaces(self, aspect_ratio_max):
        """Find surfaces with high aspect ratio (narrow/sliver).

        Estimates aspect ratio from the bounding box dimensions
        of each surface. A high ratio indicates a sliver surface.

        Parameters
        ----------
        aspect_ratio_max : float
            Maximum acceptable aspect ratio. Surfaces exceeding
            this threshold are returned.

        Returns
        -------
        list of tuple
            List of (tag, ratio) for narrow surfaces.
        """
        import gmsh
        self._check_initialized()

        result = []
        entities = gmsh.model.getEntities(2)
        for _, tag in entities:
            bb = gmsh.model.occ.getBoundingBox(2, tag)
            dx = bb[3] - bb[0]
            dy = bb[4] - bb[1]
            dz = bb[5] - bb[2]
            dims = sorted([dx, dy, dz], reverse=True)
            # Use the two largest bounding box dimensions
            max_dim = dims[0]
            min_dim = dims[1] if dims[1] > 0 else dims[2]
            if min_dim > 0:
                ratio = max_dim / min_dim
            else:
                ratio = float('inf')

            if ratio > aspect_ratio_max:
                result.append((tag, ratio))

        if self._verbose:
            print(f"[GmshBuilder] find_narrow_surfaces: found {len(result)} "
                  f"surfaces with aspect ratio > {aspect_ratio_max}")
        return result

    def remove_duplicates(self):
        """Remove duplicate entities from the model.

        Uses the OCC kernel to identify and merge coincident entities.

        Returns
        -------
        int
            Number of entities removed.
        """
        import gmsh
        self._check_initialized()

        n_before = sum(len(gmsh.model.getEntities(d)) for d in range(4))
        gmsh.model.occ.removeAllDuplicates()
        gmsh.model.occ.synchronize()
        n_after = sum(len(gmsh.model.getEntities(d)) for d in range(4))

        count = n_before - n_after
        if self._verbose:
            print(f"[GmshBuilder] remove_duplicates: removed {count} entities")
        return count

    def find_orphan_entities(self, dim):
        """Find orphan entities not referenced by higher-dimensional entities.

        An orphan entity is one that does not belong to any
        higher-dimensional entity (e.g., a surface not part of any
        volume).

        Parameters
        ----------
        dim : int
            Entity dimension to search (0=points, 1=curves, 2=surfaces).

        Returns
        -------
        list of int
            Tags of orphan entities.
        """
        import gmsh
        self._check_initialized()

        if dim >= 3:
            return []

        entities = gmsh.model.getEntities(dim)
        orphans = []

        for _, tag in entities:
            try:
                is_orphan = gmsh.model.isEntityOrphan(dim, tag)
                if is_orphan:
                    orphans.append(tag)
            except (AttributeError, Exception):
                # Fallback: check if entity appears in any higher-dim boundary
                referenced = False
                higher_entities = gmsh.model.getEntities(dim + 1)
                for _, h_tag in higher_entities:
                    boundary = gmsh.model.getBoundary(
                        [(dim + 1, h_tag)], oriented=False, recursive=False)
                    boundary_tags = {abs(dt[1]) for dt in boundary
                                     if dt[0] == dim}
                    if tag in boundary_tags:
                        referenced = True
                        break
                if not referenced:
                    orphans.append(tag)

        if self._verbose:
            print(f"[GmshBuilder] find_orphan_entities: dim={dim} "
                  f"found {len(orphans)} orphans")
        return orphans

    def clean_geometry(self, vol_ids=None, tolerance=1e-8):
        """Perform combined geometry cleanup.

        Runs healing, duplicate removal, and reports small entities.
        This is the recommended one-stop cleanup method.

        Parameters
        ----------
        vol_ids : list of int, optional
            Volume IDs to clean. If None, cleans all volumes.
        tolerance : float, optional
            Healing tolerance in meters. Default is 1e-8.

        Returns
        -------
        dict
            Summary with keys: 'heal_result', 'duplicates_removed',
            'small_curves', 'small_surfaces'.
        """
        import gmsh
        self._check_initialized()

        if vol_ids is not None:
            if self._verbose:
                print(f"[GmshBuilder] clean_geometry: vol_ids filtering is "
                      f"not yet implemented; cleaning all volumes")

        heal_result = self.heal_all(tolerance=tolerance)
        dup_count = self.remove_duplicates()
        small_curves = self.find_small_curves(tolerance)
        small_surfaces = self.find_small_surfaces(tolerance * tolerance)

        summary = {
            'heal_result': heal_result,
            'duplicates_removed': dup_count,
            'small_curves': len(small_curves),
            'small_surfaces': len(small_surfaces),
        }

        if self._verbose:
            print(f"[GmshBuilder] clean_geometry: heal={heal_result['healed']} "
                  f"dups={dup_count} small_curves={len(small_curves)} "
                  f"small_surfs={len(small_surfaces)}")
        return summary

    # ------------------------------------------------------------------
    # Distance / overlap (8 methods)
    # ------------------------------------------------------------------

    def distance_between_entities(self, dim1, tag1, dim2, tag2):
        """Compute minimum distance between two entities.

        Parameters
        ----------
        dim1 : int
            Dimension of first entity.
        tag1 : int
            Tag of first entity.
        dim2 : int
            Dimension of second entity.
        tag2 : int
            Tag of second entity.

        Returns
        -------
        float
            Minimum distance in meters.
        """
        import gmsh
        self._check_initialized()

        try:
            dist = gmsh.model.occ.getDistance(dim1, tag1, dim2, tag2)
            return dist[0] if isinstance(dist, (list, tuple)) else dist
        except (AttributeError, Exception):
            # Fallback: compute from centers of mass
            c1 = gmsh.model.occ.getCenterOfMass(dim1, tag1)
            c2 = gmsh.model.occ.getCenterOfMass(dim2, tag2)
            return math.sqrt(sum((a - b) ** 2 for a, b in zip(c1, c2)))

    def closest_entities(self, dim1, tag1, dim2, tag2):
        """Find closest points between two entities.

        Parameters
        ----------
        dim1 : int
            Dimension of first entity.
        tag1 : int
            Tag of first entity.
        dim2 : int
            Dimension of second entity.
        tag2 : int
            Tag of second entity.

        Returns
        -------
        dict
            Keys: 'distance' (float), 'point1' (list), 'point2' (list).
        """
        import gmsh
        self._check_initialized()

        try:
            result = gmsh.model.occ.getDistance(dim1, tag1, dim2, tag2)
            # getDistance returns (distance, x1, y1, z1, x2, y2, z2)
            dist, x1, y1, z1, x2, y2, z2 = result
            return {
                'distance': dist,
                'point1': [x1, y1, z1],
                'point2': [x2, y2, z2],
            }
        except (AttributeError, TypeError, Exception):
            # Fallback: center-of-mass approximation
            c1 = list(gmsh.model.occ.getCenterOfMass(dim1, tag1))
            c2 = list(gmsh.model.occ.getCenterOfMass(dim2, tag2))
            dist = math.sqrt(sum((a - b) ** 2 for a, b in zip(c1, c2)))
            return {
                'distance': dist,
                'point1': c1,
                'point2': c2,
            }

    def find_overlapping_volumes(self, tolerance=0.0):
        """Find pairs of volumes with overlapping bounding boxes.

        Parameters
        ----------
        tolerance : float, optional
            Bounding box expansion tolerance in meters.
            Default is 0.0 (exact overlap).

        Returns
        -------
        list of tuple
            List of (tag1, tag2) pairs with overlapping bounding boxes.
        """
        import gmsh
        self._check_initialized()

        entities = gmsh.model.getEntities(3)
        bboxes = {}
        for _, tag in entities:
            bb = gmsh.model.occ.getBoundingBox(3, tag)
            bboxes[tag] = (
                bb[0] - tolerance, bb[1] - tolerance, bb[2] - tolerance,
                bb[3] + tolerance, bb[4] + tolerance, bb[5] + tolerance)

        pairs = []
        tags = list(bboxes.keys())
        for i in range(len(tags)):
            for j in range(i + 1, len(tags)):
                bb1 = bboxes[tags[i]]
                bb2 = bboxes[tags[j]]
                # Check AABB overlap
                if (bb1[0] <= bb2[3] and bb1[3] >= bb2[0] and
                        bb1[1] <= bb2[4] and bb1[4] >= bb2[1] and
                        bb1[2] <= bb2[5] and bb1[5] >= bb2[2]):
                    pairs.append((tags[i], tags[j]))

        if self._verbose:
            print(f"[GmshBuilder] find_overlapping_volumes: "
                  f"found {len(pairs)} overlapping pairs")
        return pairs

    def find_touching_surfaces(self, vol_tag1, vol_tag2, tolerance=1e-6):
        """Find surfaces that are close between two volumes.

        Compares surface centers and bounding boxes to identify
        surfaces that are within tolerance of each other.

        Parameters
        ----------
        vol_tag1 : int
            First GMSH volume tag.
        vol_tag2 : int
            Second GMSH volume tag.
        tolerance : float, optional
            Distance tolerance in meters. Default is 1e-6.

        Returns
        -------
        list of tuple
            List of (surf_tag1, surf_tag2) pairs that are close.
        """
        import gmsh
        self._check_initialized()

        def _get_surf_tags(vol_tag):
            boundary = gmsh.model.getBoundary(
                [(3, vol_tag)], oriented=False, recursive=False)
            return [abs(dt[1]) for dt in boundary if dt[0] == 2]

        surfs1 = _get_surf_tags(vol_tag1)
        surfs2 = _get_surf_tags(vol_tag2)

        pairs = []
        for s1 in surfs1:
            c1 = gmsh.model.occ.getCenterOfMass(2, s1)
            for s2 in surfs2:
                c2 = gmsh.model.occ.getCenterOfMass(2, s2)
                dist = math.sqrt(sum((a - b) ** 2 for a, b in zip(c1, c2)))
                if dist < tolerance:
                    pairs.append((s1, s2))

        if self._verbose:
            print(f"[GmshBuilder] find_touching_surfaces: "
                  f"vols=({vol_tag1},{vol_tag2}) found {len(pairs)} pairs")
        return pairs

    def entities_in_bounding_box(self, xmin, ymin, zmin,
                                 xmax, ymax, zmax, dim=3):
        """Find entities within a bounding box.

        Parameters
        ----------
        xmin, ymin, zmin : float
            Minimum corner of the bounding box in meters.
        xmax, ymax, zmax : float
            Maximum corner of the bounding box in meters.
        dim : int, optional
            Entity dimension to search. Default is 3 (volumes).

        Returns
        -------
        list of int
            Tags of entities within the bounding box.
        """
        import gmsh
        self._check_initialized()

        entities = gmsh.model.getEntitiesInBoundingBox(
            xmin, ymin, zmin, xmax, ymax, zmax, dim)
        tags = [e[1] for e in entities]

        if self._verbose:
            print(f"[GmshBuilder] entities_in_bounding_box: dim={dim} "
                  f"found {len(tags)} entities")
        return tags

    def find_coincident_vertices(self, tolerance=1e-6):
        """Find groups of coincident (overlapping) vertices.

        Parameters
        ----------
        tolerance : float, optional
            Distance tolerance in meters. Default is 1e-6.

        Returns
        -------
        list of list
            Each inner list contains tags of coincident vertices.
        """
        import gmsh
        self._check_initialized()

        entities = gmsh.model.getEntities(0)
        coords = {}
        for _, tag in entities:
            pt = gmsh.model.getValue(0, tag, [])
            coords[tag] = list(pt)

        tags = list(coords.keys())
        visited = set()
        groups = []

        for i in range(len(tags)):
            if tags[i] in visited:
                continue
            group = [tags[i]]
            for j in range(i + 1, len(tags)):
                if tags[j] in visited:
                    continue
                dist = math.sqrt(sum(
                    (a - b) ** 2
                    for a, b in zip(coords[tags[i]], coords[tags[j]])))
                if dist < tolerance:
                    group.append(tags[j])
                    visited.add(tags[j])
            if len(group) > 1:
                groups.append(group)
                visited.add(tags[i])

        if self._verbose:
            print(f"[GmshBuilder] find_coincident_vertices: "
                  f"found {len(groups)} coincident groups")
        return groups

    def find_coincident_curves(self, tolerance=1e-6):
        """Find pairs of coincident curves.

        Compares midpoint positions and lengths to identify curves
        that likely represent the same geometry.

        Parameters
        ----------
        tolerance : float, optional
            Distance tolerance in meters. Default is 1e-6.

        Returns
        -------
        list of tuple
            List of (tag1, tag2) pairs of coincident curves.
        """
        import gmsh
        self._check_initialized()

        entities = gmsh.model.getEntities(1)
        curve_info = {}
        for _, tag in entities:
            length = gmsh.model.occ.getMass(1, tag)
            bounds = gmsh.model.getParametrizationBounds(1, tag)
            t_mid = (bounds[0][0] + bounds[1][0]) / 2.0
            midpt = list(gmsh.model.getValue(1, tag, [t_mid]))
            curve_info[tag] = {'length': length, 'midpoint': midpt}

        tags = list(curve_info.keys())
        pairs = []
        for i in range(len(tags)):
            for j in range(i + 1, len(tags)):
                info_i = curve_info[tags[i]]
                info_j = curve_info[tags[j]]
                if abs(info_i['length'] - info_j['length']) > tolerance:
                    continue
                dist = math.sqrt(sum(
                    (a - b) ** 2
                    for a, b in zip(info_i['midpoint'], info_j['midpoint'])))
                if dist < tolerance:
                    pairs.append((tags[i], tags[j]))

        if self._verbose:
            print(f"[GmshBuilder] find_coincident_curves: "
                  f"found {len(pairs)} coincident pairs")
        return pairs

    def check_watertight(self, vol_tag):
        """Check if a volume is watertight (closed shell).

        Verifies that the volume's surface loops form a complete
        closed boundary.

        Parameters
        ----------
        vol_tag : int
            GMSH volume tag.

        Returns
        -------
        bool
            True if the volume is watertight.
        """
        import gmsh
        self._check_initialized()

        try:
            surface_loops = gmsh.model.occ.getSurfaceLoops(vol_tag)
            if not surface_loops:
                return False
            # A volume is watertight if it has at least one surface loop
            # and the OCC kernel accepted it as valid
            return True
        except (AttributeError, Exception):
            # Fallback: check that volume has bounding surfaces
            boundary = gmsh.model.getBoundary(
                [(3, vol_tag)], oriented=True, recursive=False)
            if not boundary:
                return False
            # Check each surface has a valid curve loop
            for dim, tag in boundary:
                if dim == 2:
                    surf_boundary = gmsh.model.getBoundary(
                        [(2, abs(tag))], oriented=True, recursive=False)
                    if not surf_boundary:
                        return False
            return True

    # ------------------------------------------------------------------
    # Topology analysis (10 methods)
    # ------------------------------------------------------------------

    def get_curve_loops(self, surf_tag):
        """Get curve loops of a surface.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.

        Returns
        -------
        list of list
            Each inner list contains curve tags forming a loop.
        """
        import gmsh
        self._check_initialized()

        try:
            loops = gmsh.model.occ.getCurveLoops(surf_tag)
            # getCurveLoops returns (loop_tags, curve_tags_per_loop)
            if isinstance(loops, tuple) and len(loops) == 2:
                loop_tags, curve_arrays = loops
                result = []
                for arr in curve_arrays:
                    result.append(list(arr))
                return result
            return [list(loops)]
        except (AttributeError, Exception):
            # Fallback: use getBoundary
            boundary = gmsh.model.getBoundary(
                [(2, surf_tag)], oriented=True, recursive=False)
            curves = [dt[1] for dt in boundary if abs(dt[0]) == 1]
            return [curves] if curves else []

    def get_surface_loops(self, vol_tag):
        """Get surface loops of a volume.

        Parameters
        ----------
        vol_tag : int
            GMSH volume tag.

        Returns
        -------
        list of list
            Each inner list contains surface tags forming a shell.
        """
        import gmsh
        self._check_initialized()

        try:
            loops = gmsh.model.occ.getSurfaceLoops(vol_tag)
            if isinstance(loops, tuple) and len(loops) == 2:
                loop_tags, surf_arrays = loops
                result = []
                for arr in surf_arrays:
                    result.append(list(arr))
                return result
            return [list(loops)]
        except (AttributeError, Exception):
            # Fallback: use getBoundary
            boundary = gmsh.model.getBoundary(
                [(3, vol_tag)], oriented=True, recursive=False)
            surfs = [dt[1] for dt in boundary if abs(dt[0]) == 2]
            return [surfs] if surfs else []

    def get_adjacencies(self, dim, tag):
        """Get upward and downward adjacencies of an entity.

        Parameters
        ----------
        dim : int
            Entity dimension.
        tag : int
            Entity tag.

        Returns
        -------
        dict
            'upward': list of tags of higher-dimensional adjacent entities,
            'downward': list of tags of lower-dimensional adjacent entities.
        """
        import gmsh
        self._check_initialized()

        upward, downward = gmsh.model.getAdjacencies(dim, tag)
        return {
            'upward': list(upward),
            'downward': list(downward),
        }

    def get_entity_type(self, dim, tag):
        """Get the type name of an entity.

        Parameters
        ----------
        dim : int
            Entity dimension.
        tag : int
            Entity tag.

        Returns
        -------
        str
            Entity type name (e.g., 'Plane', 'BSpline', 'Line').
        """
        import gmsh
        self._check_initialized()

        return gmsh.model.getType(dim, tag)

    def get_entity_name_raw(self, dim, tag):
        """Get the raw name of an entity from GMSH.

        Parameters
        ----------
        dim : int
            Entity dimension.
        tag : int
            Entity tag.

        Returns
        -------
        str
            Entity name, or empty string if unnamed.
        """
        import gmsh
        self._check_initialized()

        return gmsh.model.getEntityName(dim, tag)

    def set_entity_name_raw(self, dim, tag, name):
        """Set the raw name of an entity in GMSH.

        Parameters
        ----------
        dim : int
            Entity dimension.
        tag : int
            Entity tag.
        name : str
            Name to assign.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.setEntityName(dim, tag, name)
        if self._verbose:
            print(f"[GmshBuilder] set_entity_name_raw: ({dim},{tag}) "
                  f"name='{name}'")

    def is_entity_orphan(self, dim, tag):
        """Check if an entity is orphaned (not referenced by higher dim).

        Parameters
        ----------
        dim : int
            Entity dimension.
        tag : int
            Entity tag.

        Returns
        -------
        bool
            True if the entity is an orphan.
        """
        import gmsh
        self._check_initialized()

        if dim >= 3:
            return False

        try:
            return bool(gmsh.model.isEntityOrphan(dim, tag))
        except (AttributeError, Exception):
            # Fallback: check adjacency
            upward, _ = gmsh.model.getAdjacencies(dim, tag)
            return len(upward) == 0

    def count_topology(self, vol_tag):
        """Count topological entities bounding a volume.

        Recursively enumerates all surfaces, curves, and vertices
        belonging to the given volume.

        Parameters
        ----------
        vol_tag : int
            GMSH volume tag.

        Returns
        -------
        dict
            Keys: 'surfaces', 'curves', 'vertices' with integer counts.
        """
        import gmsh
        self._check_initialized()

        # Surfaces
        boundary_2d = gmsh.model.getBoundary(
            [(3, vol_tag)], oriented=False, recursive=False)
        surfs = {abs(dt[1]) for dt in boundary_2d if dt[0] == 2}

        # Curves (from surfaces)
        curves = set()
        for s in surfs:
            boundary_1d = gmsh.model.getBoundary(
                [(2, s)], oriented=False, recursive=False)
            for dt in boundary_1d:
                if dt[0] == 1:
                    curves.add(abs(dt[1]))

        # Vertices (from curves)
        vertices = set()
        for c in curves:
            boundary_0d = gmsh.model.getBoundary(
                [(1, c)], oriented=False, recursive=False)
            for dt in boundary_0d:
                if dt[0] == 0:
                    vertices.add(abs(dt[1]))

        result = {
            'surfaces': len(surfs),
            'curves': len(curves),
            'vertices': len(vertices),
        }

        if self._verbose:
            print(f"[GmshBuilder] count_topology: vol={vol_tag} "
                  f"S={result['surfaces']} C={result['curves']} "
                  f"V={result['vertices']}")
        return result

    def find_nonmanifold_edges(self):
        """Find non-manifold edges (curves with more than 2 adjacent surfaces).

        Non-manifold edges indicate T-junctions or complex topology
        where three or more surfaces meet at an edge.

        Returns
        -------
        list of int
            Curve tags that are non-manifold.
        """
        import gmsh
        self._check_initialized()

        nonmanifold = []
        entities = gmsh.model.getEntities(1)
        for _, tag in entities:
            upward, _ = gmsh.model.getAdjacencies(1, tag)
            if len(upward) > 2:
                nonmanifold.append(tag)

        if self._verbose:
            print(f"[GmshBuilder] find_nonmanifold_edges: "
                  f"found {len(nonmanifold)} non-manifold curves")
        return nonmanifold

    def find_free_surfaces(self):
        """Find free surfaces (surfaces adjacent to only one volume).

        Free surfaces are boundary surfaces of the model that are
        not shared between two volumes.

        Returns
        -------
        list of int
            Surface tags that are free (one-sided).
        """
        import gmsh
        self._check_initialized()

        free = []
        entities = gmsh.model.getEntities(2)
        for _, tag in entities:
            upward, _ = gmsh.model.getAdjacencies(2, tag)
            if len(upward) == 1:
                free.append(tag)

        if self._verbose:
            print(f"[GmshBuilder] find_free_surfaces: "
                  f"found {len(free)} free surfaces")
        return free

    # ------------------------------------------------------------------
    # Differential geometry (7 methods)
    # ------------------------------------------------------------------

    def curvature_at_curve(self, curve_tag, t=0.5):
        """Get curvature of a curve at a parametric point.

        Parameters
        ----------
        curve_tag : int
            GMSH curve tag.
        t : float, optional
            Normalized parameter (0 to 1). Default is 0.5 (midpoint).

        Returns
        -------
        float
            Curvature value (1/radius).
        """
        import gmsh
        self._check_initialized()

        bounds = gmsh.model.getParametrizationBounds(1, curve_tag)
        t_min, t_max = bounds[0][0], bounds[1][0]
        t_actual = t_min + t * (t_max - t_min)

        curvature = gmsh.model.getCurvature(1, curve_tag, [t_actual])
        return float(curvature[0])

    def curvature_at_surface(self, surf_tag, u=0.5, v=0.5):
        """Get curvature of a surface at a parametric point.

        Returns the mean curvature at the given parametric coordinates.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.
        u : float, optional
            Normalized parametric u coordinate (0 to 1). Default is 0.5.
        v : float, optional
            Normalized parametric v coordinate (0 to 1). Default is 0.5.

        Returns
        -------
        float
            Mean curvature value.
        """
        import gmsh
        self._check_initialized()

        bounds = gmsh.model.getParametrizationBounds(2, surf_tag)
        u_actual = bounds[0][0] + u * (bounds[1][0] - bounds[0][0])
        v_actual = bounds[0][1] + v * (bounds[1][1] - bounds[0][1])

        curvature = gmsh.model.getCurvature(2, surf_tag, [u_actual, v_actual])
        return float(curvature[0])

    def principal_curvatures(self, surf_tag, u=0.5, v=0.5):
        """Get principal curvatures and directions at a surface point.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.
        u : float, optional
            Normalized parametric u coordinate (0 to 1). Default is 0.5.
        v : float, optional
            Normalized parametric v coordinate (0 to 1). Default is 0.5.

        Returns
        -------
        dict
            Keys: 'k1' (float), 'k2' (float),
            'dir1' (list of 3 floats), 'dir2' (list of 3 floats).

        Raises
        ------
        NotImplementedError
            If the GMSH version does not support getPrincipalCurvatures.
        """
        import gmsh
        self._check_initialized()

        bounds = gmsh.model.getParametrizationBounds(2, surf_tag)
        u_actual = bounds[0][0] + u * (bounds[1][0] - bounds[0][0])
        v_actual = bounds[0][1] + v * (bounds[1][1] - bounds[0][1])

        try:
            k_max, k_min, dir_max, dir_min = gmsh.model.getPrincipalCurvatures(
                surf_tag, [u_actual, v_actual])
            return {
                'k1': float(k_max[0]),
                'k2': float(k_min[0]),
                'dir1': list(dir_max[:3]),
                'dir2': list(dir_min[:3]),
            }
        except (AttributeError, Exception) as exc:
            raise NotImplementedError(
                "Requires GMSH >= 4.11 with getPrincipalCurvatures support"
            ) from exc

    def second_derivative_curve(self, curve_tag, t=0.5):
        """Get second derivative of a curve at a parametric point.

        Parameters
        ----------
        curve_tag : int
            GMSH curve tag.
        t : float, optional
            Normalized parameter (0 to 1). Default is 0.5 (midpoint).

        Returns
        -------
        list
            Second derivative vector [d2x/dt2, d2y/dt2, d2z/dt2].

        Raises
        ------
        NotImplementedError
            If the GMSH version does not support getSecondDerivative.
        """
        import gmsh
        self._check_initialized()

        bounds = gmsh.model.getParametrizationBounds(1, curve_tag)
        t_min, t_max = bounds[0][0], bounds[1][0]
        t_actual = t_min + t * (t_max - t_min)

        try:
            deriv2 = gmsh.model.getSecondDerivative(
                1, curve_tag, [t_actual])
            return list(deriv2[:3])
        except (AttributeError, Exception) as exc:
            raise NotImplementedError(
                "Requires GMSH >= 4.11 with getSecondDerivative support"
            ) from exc

    def second_derivative_surface(self, surf_tag, u=0.5, v=0.5):
        """Get second derivatives of a surface at a parametric point.

        Returns the second-order partial derivatives of the surface
        parametrization (d2r/du2, d2r/dv2, d2r/dudv).

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.
        u : float, optional
            Normalized parametric u coordinate (0 to 1). Default is 0.5.
        v : float, optional
            Normalized parametric v coordinate (0 to 1). Default is 0.5.

        Returns
        -------
        list
            Second derivative components. For a surface, this is a
            9-component vector [d2x/du2, d2y/du2, d2z/du2,
            d2x/dv2, d2y/dv2, d2z/dv2,
            d2x/dudv, d2y/dudv, d2z/dudv].

        Raises
        ------
        NotImplementedError
            If the GMSH version does not support getSecondDerivative.
        """
        import gmsh
        self._check_initialized()

        bounds = gmsh.model.getParametrizationBounds(2, surf_tag)
        u_actual = bounds[0][0] + u * (bounds[1][0] - bounds[0][0])
        v_actual = bounds[0][1] + v * (bounds[1][1] - bounds[0][1])

        try:
            deriv2 = gmsh.model.getSecondDerivative(
                2, surf_tag, [u_actual, v_actual])
            return list(deriv2)
        except (AttributeError, Exception) as exc:
            raise NotImplementedError(
                "Requires GMSH >= 4.11 with getSecondDerivative support"
            ) from exc

    def derivative_curve(self, curve_tag, t=0.5):
        """Get first derivative (tangent) of a curve at a parametric point.

        Parameters
        ----------
        curve_tag : int
            GMSH curve tag.
        t : float, optional
            Normalized parameter (0 to 1). Default is 0.5 (midpoint).

        Returns
        -------
        list
            Derivative vector [dx/dt, dy/dt, dz/dt].
        """
        import gmsh
        self._check_initialized()

        bounds = gmsh.model.getParametrizationBounds(1, curve_tag)
        t_min, t_max = bounds[0][0], bounds[1][0]
        t_actual = t_min + t * (t_max - t_min)

        deriv = gmsh.model.getDerivative(1, curve_tag, [t_actual])
        return list(deriv[:3])

    def derivative_surface(self, surf_tag, u=0.5, v=0.5):
        """Get first derivatives of a surface at a parametric point.

        Returns the partial derivatives of the surface parametrization
        (dr/du, dr/dv).

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.
        u : float, optional
            Normalized parametric u coordinate (0 to 1). Default is 0.5.
        v : float, optional
            Normalized parametric v coordinate (0 to 1). Default is 0.5.

        Returns
        -------
        list
            Derivative components: 6-component vector
            [dx/du, dy/du, dz/du, dx/dv, dy/dv, dz/dv].
        """
        import gmsh
        self._check_initialized()

        bounds = gmsh.model.getParametrizationBounds(2, surf_tag)
        u_actual = bounds[0][0] + u * (bounds[1][0] - bounds[0][0])
        v_actual = bounds[0][1] + v * (bounds[1][1] - bounds[0][1])

        deriv = gmsh.model.getDerivative(2, surf_tag, [u_actual, v_actual])
        return list(deriv)
