"""
Mesh data access for GmshBuilder.

Gmsh equivalents: get_node_count, get_connectivity, get_nodal_coordinates.
"""

import numpy as np


class MeshAccessMixin:
    """Mixin providing mesh node/element data access."""

    def get_nodes(self):
        """Get all node tags.

        Returns
        -------
        numpy array
            1D array of node tags.
        """
        import gmsh
        self._check_initialized()

        tags, _, _ = gmsh.model.mesh.getNodes()
        return np.array(tags, dtype=int)

    def get_node_count(self):
        """Get total node count.

        Returns
        -------
        int
        """
        import gmsh
        self._check_initialized()

        tags, _, _ = gmsh.model.mesh.getNodes()
        return len(tags)

    def get_node_coordinates(self, node_tag):
        """Get coordinates of a specific node.

        Parameters
        ----------
        node_tag : int
            GMSH node tag.

        Returns
        -------
        list
            [x, y, z].
        """
        import gmsh
        self._check_initialized()

        coords, _, _, _ = gmsh.model.mesh.getNode(node_tag)
        return list(coords)

    def get_all_node_coordinates(self):
        """Get all node coordinates as numpy array.

        Returns
        -------
        numpy array
            (n_nodes, 3) array.
        """
        import gmsh
        self._check_initialized()

        tags, coords, _ = gmsh.model.mesh.getNodes()
        return np.array(coords).reshape(-1, 3)

    def get_elements(self, dim=3):
        """Get all elements of a given dimension.

        Parameters
        ----------
        dim : int
            Element dimension: 2 (surface), 3 (volume).

        Returns
        -------
        list of dict
            Each dict has 'tag', 'type', 'type_name', 'nodes'.
        """
        import gmsh
        self._check_initialized()

        elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(
            dim)
        elements = []
        for i, et in enumerate(elem_types):
            name, _, _, n_nodes, _, _ = gmsh.model.mesh.getElementProperties(
                et)
            tags = elem_tags[i]
            nodes = elem_node_tags[i]
            for j, tag in enumerate(tags):
                start = j * n_nodes
                elem_nodes = list(nodes[start:start + n_nodes])
                elements.append({
                    'tag': int(tag),
                    'type': int(et),
                    'type_name': name,
                    'nodes': [int(n) for n in elem_nodes],
                })

        return elements

    def get_element_nodes(self, elem_tag):
        """Get node tags of a specific element.

        Parameters
        ----------
        elem_tag : int
            GMSH element tag.

        Returns
        -------
        list of int
            Node tags.
        """
        import gmsh
        self._check_initialized()

        et, nodes, _, _ = gmsh.model.mesh.getElement(elem_tag)
        return [int(n) for n in nodes]

    def get_element_type(self, elem_tag):
        """Get element type name.

        Parameters
        ----------
        elem_tag : int
            GMSH element tag.

        Returns
        -------
        str
            Element type name (e.g., 'Hexahedron 8').
        """
        import gmsh
        self._check_initialized()

        et, _, _, _ = gmsh.model.mesh.getElement(elem_tag)
        name = gmsh.model.mesh.getElementProperties(et)[0]
        return name

    def get_hex_count(self):
        """Count hexahedral elements."""
        return self._count_elements_by_name('Hexahedron')

    def get_tet_count(self):
        """Count tetrahedral elements."""
        return self._count_elements_by_name('Tetrahedron')

    def get_tri_count(self):
        """Count triangular elements (surface)."""
        return self._count_elements_by_name('Triangle', dim=2)

    def get_quad_count(self):
        """Count quadrilateral elements (surface)."""
        return self._count_elements_by_name('Quadrilateral', dim=2)

    def get_wedge_count(self):
        """Count wedge/prism elements."""
        count = self._count_elements_by_name('Prism')
        count += self._count_elements_by_name('Wedge')
        return count

    def get_element_count(self, dim=3):
        """Count total elements of given dimension."""
        import gmsh
        self._check_initialized()

        elem_types, elem_tags, _ = gmsh.model.mesh.getElements(dim)
        return sum(len(tags) for tags in elem_tags)

    def _count_elements_by_name(self, name_pattern, dim=3):
        """Count elements matching a type name pattern."""
        import gmsh
        self._check_initialized()

        elem_types, elem_tags, _ = gmsh.model.mesh.getElements(dim)
        count = 0
        for i, et in enumerate(elem_types):
            name = gmsh.model.mesh.getElementProperties(et)[0]
            if name_pattern in name:
                count += len(elem_tags[i])
        return count

    def get_volume_elements(self, vol_id):
        """Get elements belonging to a specific volume.

        Parameters
        ----------
        vol_id : int
            Volume ID.

        Returns
        -------
        list of dict
            Element info dicts.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        elements = []
        for vt in self._volumes[vol_id]:
            elem_types, elem_tags, elem_nodes = gmsh.model.mesh.getElements(
                3, vt)
            for i, et in enumerate(elem_types):
                name, _, _, n_nodes, _, _ = \
                    gmsh.model.mesh.getElementProperties(et)
                tags = elem_tags[i]
                nodes = elem_nodes[i]
                for j, tag in enumerate(tags):
                    start = j * n_nodes
                    elements.append({
                        'tag': int(tag),
                        'type': int(et),
                        'type_name': name,
                        'nodes': [int(n) for n in nodes[start:start + n_nodes]],
                    })

        return elements

    def get_surface_elements(self, surf_tag):
        """Get elements on a specific surface.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag.

        Returns
        -------
        list of dict
            Element info dicts.
        """
        import gmsh
        self._check_initialized()

        elements = []
        elem_types, elem_tags, elem_nodes = gmsh.model.mesh.getElements(
            2, surf_tag)
        for i, et in enumerate(elem_types):
            name, _, _, n_nodes, _, _ = \
                gmsh.model.mesh.getElementProperties(et)
            tags = elem_tags[i]
            nodes = elem_nodes[i]
            for j, tag in enumerate(tags):
                start = j * n_nodes
                elements.append({
                    'tag': int(tag),
                    'type': int(et),
                    'type_name': name,
                    'nodes': [int(n) for n in nodes[start:start + n_nodes]],
                })

        return elements

    def get_closest_node(self, x, y, z):
        """Find the closest node to a point.

        Parameters
        ----------
        x, y, z : float
            Query point.

        Returns
        -------
        dict
            tag, coordinates, distance.
        """
        import gmsh
        self._check_initialized()

        tags, coords, _ = gmsh.model.mesh.getNodes()
        all_coords = np.array(coords).reshape(-1, 3)
        point = np.array([x, y, z])
        distances = np.linalg.norm(all_coords - point, axis=1)
        idx = np.argmin(distances)

        return {
            'tag': int(tags[idx]),
            'coordinates': all_coords[idx].tolist(),
            'distance': float(distances[idx]),
        }

    def get_nodes_in_box(self, bbox):
        """Get nodes within a bounding box.

        Parameters
        ----------
        bbox : list
            [xmin, ymin, zmin, xmax, ymax, zmax].

        Returns
        -------
        numpy array
            Node tags within the box.
        """
        import gmsh
        self._check_initialized()

        tags, coords, _ = gmsh.model.mesh.getNodes()
        all_coords = np.array(coords).reshape(-1, 3)

        mask = ((all_coords[:, 0] >= bbox[0]) &
                (all_coords[:, 0] <= bbox[3]) &
                (all_coords[:, 1] >= bbox[1]) &
                (all_coords[:, 1] <= bbox[4]) &
                (all_coords[:, 2] >= bbox[2]) &
                (all_coords[:, 2] <= bbox[5]))

        return np.array(tags, dtype=int)[mask]

    def get_boundary_nodes(self):
        """Get nodes on the mesh boundary.

        Returns
        -------
        numpy array
            Boundary node tags.
        """
        import gmsh
        self._check_initialized()

        # Get boundary surfaces, then their nodes
        boundary_nodes = set()
        elem_types, elem_tags, elem_nodes = gmsh.model.mesh.getElements(2)
        for i, et in enumerate(elem_types):
            _, _, _, n_nodes, _, _ = gmsh.model.mesh.getElementProperties(et)
            nodes = elem_nodes[i]
            for n in nodes:
                boundary_nodes.add(int(n))

        return np.array(sorted(boundary_nodes), dtype=int)

    def get_element_centroid(self, elem_tag):
        """Get centroid of an element.

        Parameters
        ----------
        elem_tag : int
            GMSH element tag.

        Returns
        -------
        list
            [x, y, z] centroid.
        """
        import gmsh
        self._check_initialized()

        _, nodes, _, _ = gmsh.model.mesh.getElement(elem_tag)
        coords = []
        for n in nodes:
            node_data = gmsh.model.mesh.getNode(int(n))
            coords.append(node_data[0])

        coords = np.array(coords)
        return coords.mean(axis=0).tolist()

    def get_elements_as_numpy(self, dim=3):
        """Get element connectivity as numpy arrays.

        Returns
        -------
        dict
            'node_coords': (n_nodes, 3) coords,
            'hex': (n_hex, 8) connectivity,
            'tet': (n_tet, 4) connectivity,
            'wedge': (n_wedge, 6) connectivity.
        """
        import gmsh
        self._check_initialized()

        # Nodes
        tags, coords, _ = gmsh.model.mesh.getNodes()
        node_coords = np.array(coords).reshape(-1, 3)
        tag_to_idx = {int(t): i for i, t in enumerate(tags)}

        result = {'node_coords': node_coords}

        elem_types, elem_tags, elem_nodes = gmsh.model.mesh.getElements(dim)
        for i, et in enumerate(elem_types):
            name, _, _, n_nodes, _, _ = \
                gmsh.model.mesh.getElementProperties(et)
            nodes = elem_nodes[i]
            n_elems = len(elem_tags[i])

            conn = np.array(nodes, dtype=int).reshape(n_elems, n_nodes)
            # Map to 0-based indices
            for r in range(n_elems):
                for c in range(n_nodes):
                    conn[r, c] = tag_to_idx.get(conn[r, c], conn[r, c])

            if 'Hexahedron' in name:
                result['hex'] = conn
            elif 'Tetrahedron' in name:
                result['tet'] = conn
            elif 'Prism' in name or 'Wedge' in name:
                result['wedge'] = conn
            elif 'Triangle' in name:
                result['tri'] = conn
            elif 'Quadrilateral' in name:
                result['quad'] = conn

        return result

    def get_node_coordinates_batch(self, node_tags):
        """Get coordinates for multiple nodes at once.

        Parameters
        ----------
        node_tags : list or array of int
            Node tags to query.

        Returns
        -------
        numpy array
            (n, 3) array of coordinates.
        """
        import gmsh
        self._check_initialized()

        coords = np.empty((len(node_tags), 3), dtype=float)
        for i, tag in enumerate(node_tags):
            node_data = gmsh.model.mesh.getNode(int(tag))
            coords[i] = node_data[0]
        return coords

    def get_elements_by_type(self, elem_type_name, dim=3):
        """Get elements of a specific type name.

        Parameters
        ----------
        elem_type_name : str
            Element type name (e.g. 'Hexahedron 8', 'Tetrahedron 4').
        dim : int
            Element dimension (default 3).

        Returns
        -------
        list of dict
            Each dict has 'tag', 'type', 'type_name', 'nodes'.
        """
        import gmsh
        self._check_initialized()

        elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(
            dim)
        elements = []
        for i, et in enumerate(elem_types):
            name, _, _, n_nodes, _, _ = gmsh.model.mesh.getElementProperties(
                et)
            if name != elem_type_name:
                continue
            tags = elem_tags[i]
            nodes = elem_node_tags[i]
            for j, tag in enumerate(tags):
                start = j * n_nodes
                elements.append({
                    'tag': int(tag),
                    'type': int(et),
                    'type_name': name,
                    'nodes': [int(n) for n in nodes[start:start + n_nodes]],
                })

        return elements

    def get_pyramid_count(self):
        """Count pyramid elements."""
        return self._count_elements_by_name('Pyramid')

    def get_prism_count(self):
        """Count prism elements (alias for wedge/prism).

        Returns
        -------
        int
            Number of prism/wedge elements.
        """
        return self.get_wedge_count()

    def get_element_volume(self, elem_tag):
        """Compute volume of a 3D element from its node coordinates.

        Supports tetrahedra (4 nodes), wedges/prisms (6 nodes),
        hexahedra (8 nodes), and pyramids (5 nodes).

        Parameters
        ----------
        elem_tag : int
            GMSH element tag.

        Returns
        -------
        float
            Element volume (always positive).
        """
        import gmsh
        self._check_initialized()

        et, nodes, _, _ = gmsh.model.mesh.getElement(elem_tag)
        name = gmsh.model.mesh.getElementProperties(et)[0]
        n_nodes = len(nodes)

        # Gather vertex coordinates
        verts = np.empty((n_nodes, 3), dtype=float)
        for i, n in enumerate(nodes):
            c, _, _, _ = gmsh.model.mesh.getNode(int(n))
            verts[i] = c

        if n_nodes == 4:
            # Tetrahedron: V = |det([v1-v0, v2-v0, v3-v0])| / 6
            mat = np.array([verts[1] - verts[0],
                            verts[2] - verts[0],
                            verts[3] - verts[0]])
            return abs(np.linalg.det(mat)) / 6.0

        elif n_nodes == 5:
            # Pyramid: decompose into 2 tetrahedra
            # Base: v0,v1,v2,v3; apex: v4
            # Tet1: v0,v1,v2,v4; Tet2: v0,v2,v3,v4
            def _tet_vol(a, b, c, d):
                mat = np.array([b - a, c - a, d - a])
                return abs(np.linalg.det(mat)) / 6.0
            vol = (_tet_vol(verts[0], verts[1], verts[2], verts[4]) +
                   _tet_vol(verts[0], verts[2], verts[3], verts[4]))
            return vol

        elif n_nodes == 6:
            # Wedge/prism: decompose into 3 tetrahedra
            # Tri faces: v0,v1,v2 and v3,v4,v5
            def _tet_vol(a, b, c, d):
                mat = np.array([b - a, c - a, d - a])
                return abs(np.linalg.det(mat)) / 6.0
            vol = (_tet_vol(verts[0], verts[1], verts[2], verts[3]) +
                   _tet_vol(verts[1], verts[2], verts[3], verts[4]) +
                   _tet_vol(verts[2], verts[3], verts[4], verts[5]))
            return vol

        elif n_nodes == 8:
            # Hexahedron: decompose into 5 tetrahedra
            def _tet_vol(a, b, c, d):
                mat = np.array([b - a, c - a, d - a])
                return abs(np.linalg.det(mat)) / 6.0
            vol = (_tet_vol(verts[0], verts[1], verts[3], verts[4]) +
                   _tet_vol(verts[1], verts[2], verts[3], verts[6]) +
                   _tet_vol(verts[1], verts[4], verts[5], verts[6]) +
                   _tet_vol(verts[3], verts[4], verts[6], verts[7]) +
                   _tet_vol(verts[1], verts[3], verts[4], verts[6]))
            return vol

        else:
            raise ValueError(
                f"Unsupported element type '{name}' with {n_nodes} nodes "
                f"for volume computation")

    def get_element_face_count(self, elem_tag):
        """Return number of faces for an element.

        Parameters
        ----------
        elem_tag : int
            GMSH element tag.

        Returns
        -------
        int
            Number of faces: 4 (tet), 5 (wedge/pyramid), 6 (hex).
        """
        import gmsh
        self._check_initialized()

        et, nodes, _, _ = gmsh.model.mesh.getElement(elem_tag)
        n_nodes = len(nodes)

        face_count_map = {4: 4, 5: 5, 6: 5, 8: 6}
        if n_nodes not in face_count_map:
            name = gmsh.model.mesh.getElementProperties(et)[0]
            raise ValueError(
                f"Unsupported element type '{name}' with {n_nodes} nodes")
        return face_count_map[n_nodes]

    def get_nodes_on_surface(self, surf_tag):
        """Get all node tags on a specific surface.

        Parameters
        ----------
        surf_tag : int
            GMSH surface entity tag.

        Returns
        -------
        numpy array
            1D array of node tags on the surface.
        """
        import gmsh
        self._check_initialized()

        tags, _, _ = gmsh.model.mesh.getNodes(dim=2, tag=surf_tag)
        return np.array(tags, dtype=int)

    def get_nodes_on_curve(self, curve_tag):
        """Get all node tags on a specific curve.

        Parameters
        ----------
        curve_tag : int
            GMSH curve entity tag.

        Returns
        -------
        numpy array
            1D array of node tags on the curve.
        """
        import gmsh
        self._check_initialized()

        tags, _, _ = gmsh.model.mesh.getNodes(dim=1, tag=curve_tag)
        return np.array(tags, dtype=int)

    def get_internal_nodes(self):
        """Get nodes NOT on the boundary (interior nodes only).

        Returns
        -------
        numpy array
            1D array of interior node tags.
        """
        import gmsh
        self._check_initialized()

        all_tags, _, _ = gmsh.model.mesh.getNodes()
        all_set = set(int(t) for t in all_tags)
        boundary_set = set(int(t) for t in self.get_boundary_nodes())
        internal = sorted(all_set - boundary_set)
        return np.array(internal, dtype=int)

    def get_node_to_element_map(self):
        """Build inverse connectivity: node tag -> list of element tags.

        Returns
        -------
        dict
            {node_tag: [elem_tag, ...]} mapping.
        """
        import gmsh
        self._check_initialized()

        node_to_elems = {}

        # Process both 2D and 3D elements
        for dim in (2, 3):
            elem_types, elem_tags, elem_node_tags = \
                gmsh.model.mesh.getElements(dim)
            for i, et in enumerate(elem_types):
                _, _, _, n_nodes, _, _ = \
                    gmsh.model.mesh.getElementProperties(et)
                tags = elem_tags[i]
                nodes = elem_node_tags[i]
                for j, tag in enumerate(tags):
                    start = j * n_nodes
                    elem_nodes = nodes[start:start + n_nodes]
                    tag_int = int(tag)
                    for n in elem_nodes:
                        n_int = int(n)
                        if n_int not in node_to_elems:
                            node_to_elems[n_int] = []
                        node_to_elems[n_int].append(tag_int)

        return node_to_elems

    def get_element_neighbors(self, elem_tag):
        """Find elements sharing at least one face with the given element.

        Parameters
        ----------
        elem_tag : int
            GMSH element tag.

        Returns
        -------
        list of int
            Tags of neighboring elements.
        """
        import gmsh
        self._check_initialized()

        # Get nodes of the target element
        _, target_nodes, _, _ = gmsh.model.mesh.getElement(elem_tag)
        target_nodes = set(int(n) for n in target_nodes)

        # Get face definitions for the target element
        target_faces = self._get_element_faces(elem_tag)

        # Build node-to-element map for lookup
        node_to_elems = self.get_node_to_element_map()

        # Candidate elements: share at least one node with target
        candidates = set()
        for n in target_nodes:
            for e in node_to_elems.get(n, []):
                if e != elem_tag:
                    candidates.add(e)

        # Check face sharing: two elements are neighbors if they share
        # a complete face (3 nodes for tri face, 4 nodes for quad face)
        neighbors = []
        for cand in candidates:
            cand_faces = self._get_element_faces(cand)
            for tf in target_faces:
                tf_set = set(tf)
                for cf in cand_faces:
                    if tf_set == set(cf):
                        neighbors.append(cand)
                        break
                else:
                    continue
                break

        return neighbors

    def _get_element_faces(self, elem_tag):
        """Get face node lists for an element.

        Parameters
        ----------
        elem_tag : int
            GMSH element tag.

        Returns
        -------
        list of list of int
            Each sub-list is node tags forming one face.
        """
        import gmsh

        et, nodes, _, _ = gmsh.model.mesh.getElement(elem_tag)
        nodes = [int(n) for n in nodes]
        n = len(nodes)

        # Use element type to get dimension and disambiguate Quad4 vs Tet4
        _, elem_dim, _, _, _, _ = gmsh.model.mesh.getElementProperties(et)

        if n == 4 and elem_dim == 2:
            # Quadrilateral (2D): 4 edges
            return [
                [nodes[0], nodes[1]],
                [nodes[1], nodes[2]],
                [nodes[2], nodes[3]],
                [nodes[3], nodes[0]],
            ]
        elif n == 4:
            # Tetrahedron (3D): 4 triangular faces
            return [
                [nodes[0], nodes[1], nodes[2]],
                [nodes[0], nodes[1], nodes[3]],
                [nodes[0], nodes[2], nodes[3]],
                [nodes[1], nodes[2], nodes[3]],
            ]
        elif n == 5:
            # Pyramid: 1 quad base + 4 tri faces
            return [
                [nodes[0], nodes[1], nodes[2], nodes[3]],
                [nodes[0], nodes[1], nodes[4]],
                [nodes[1], nodes[2], nodes[4]],
                [nodes[2], nodes[3], nodes[4]],
                [nodes[0], nodes[3], nodes[4]],
            ]
        elif n == 6:
            # Wedge/prism: 2 tri + 3 quad faces
            return [
                [nodes[0], nodes[1], nodes[2]],
                [nodes[3], nodes[4], nodes[5]],
                [nodes[0], nodes[1], nodes[4], nodes[3]],
                [nodes[1], nodes[2], nodes[5], nodes[4]],
                [nodes[0], nodes[2], nodes[5], nodes[3]],
            ]
        elif n == 8:
            # Hexahedron: 6 quad faces
            return [
                [nodes[0], nodes[1], nodes[2], nodes[3]],
                [nodes[4], nodes[5], nodes[6], nodes[7]],
                [nodes[0], nodes[1], nodes[5], nodes[4]],
                [nodes[1], nodes[2], nodes[6], nodes[5]],
                [nodes[2], nodes[3], nodes[7], nodes[6]],
                [nodes[0], nodes[3], nodes[7], nodes[4]],
            ]
        else:
            return []

    def get_edge_nodes(self):
        """Get unique edge node pairs from all 3D elements.

        Returns
        -------
        numpy array
            (n_edges, 2) array of node tag pairs (sorted per edge).
        """
        import gmsh
        self._check_initialized()

        # Edge definitions by number of element nodes
        edge_defs = {
            4: [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],
            5: [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 4),
                (2, 4), (3, 4)],
            6: [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3),
                (0, 3), (1, 4), (2, 5)],
            8: [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6),
                (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)],
        }

        edges = set()
        elem_types, elem_tags, elem_node_tags = \
            gmsh.model.mesh.getElements(3)

        for i, et in enumerate(elem_types):
            _, _, _, n_nodes, _, _ = \
                gmsh.model.mesh.getElementProperties(et)
            nodes = elem_node_tags[i]
            n_elems = len(elem_tags[i])
            local_edges = edge_defs.get(n_nodes, [])

            for j in range(n_elems):
                start = j * n_nodes
                elem_nodes = [int(nodes[start + k]) for k in range(n_nodes)]
                for e0, e1 in local_edges:
                    edge = tuple(sorted((elem_nodes[e0], elem_nodes[e1])))
                    edges.add(edge)

        if not edges:
            return np.empty((0, 2), dtype=int)
        return np.array(sorted(edges), dtype=int)

    def get_face_nodes(self, dim=3):
        """Get unique face node lists from elements of given dimension.

        For 3D elements, returns faces (triangles and quads).
        For 2D elements, returns edges.

        Parameters
        ----------
        dim : int
            Element dimension (default 3).

        Returns
        -------
        list of list of int
            Each sub-list is sorted node tags forming one face.
        """
        import gmsh
        self._check_initialized()

        faces = set()
        elem_types, elem_tags, elem_node_tags = \
            gmsh.model.mesh.getElements(dim)

        for i, et in enumerate(elem_types):
            _, e_dim, _, n_nodes, _, _ = \
                gmsh.model.mesh.getElementProperties(et)
            tags = elem_tags[i]
            nodes = elem_node_tags[i]

            for j in range(len(tags)):
                start = j * n_nodes
                elem_nodes = [int(nodes[start + k]) for k in range(n_nodes)]
                elem_faces = self._local_faces(elem_nodes, n_nodes, e_dim)
                for f in elem_faces:
                    faces.add(tuple(sorted(f)))

        return [list(f) for f in sorted(faces)]

    def _local_faces(self, nodes, n_nodes, elem_dim=3):
        """Return face node lists for a single element given its nodes.

        Parameters
        ----------
        nodes : list of int
            Element node tags.
        n_nodes : int
            Number of nodes in element.
        elem_dim : int
            Element dimension (2 or 3). Used to disambiguate Quad4
            (2D, 4 nodes -> edges) from Tet4 (3D, 4 nodes -> tri faces).

        Returns
        -------
        list of list of int
        """
        if n_nodes == 4 and elem_dim == 2:
            # Quadrilateral: 4 edges
            return [
                [nodes[0], nodes[1]],
                [nodes[1], nodes[2]],
                [nodes[2], nodes[3]],
                [nodes[3], nodes[0]],
            ]
        elif n_nodes == 4:
            # Tetrahedron (3D): 4 triangular faces
            return [
                [nodes[0], nodes[1], nodes[2]],
                [nodes[0], nodes[1], nodes[3]],
                [nodes[0], nodes[2], nodes[3]],
                [nodes[1], nodes[2], nodes[3]],
            ]
        elif n_nodes == 5:
            return [
                [nodes[0], nodes[1], nodes[2], nodes[3]],
                [nodes[0], nodes[1], nodes[4]],
                [nodes[1], nodes[2], nodes[4]],
                [nodes[2], nodes[3], nodes[4]],
                [nodes[0], nodes[3], nodes[4]],
            ]
        elif n_nodes == 6:
            return [
                [nodes[0], nodes[1], nodes[2]],
                [nodes[3], nodes[4], nodes[5]],
                [nodes[0], nodes[1], nodes[4], nodes[3]],
                [nodes[1], nodes[2], nodes[5], nodes[4]],
                [nodes[0], nodes[2], nodes[5], nodes[3]],
            ]
        elif n_nodes == 8:
            return [
                [nodes[0], nodes[1], nodes[2], nodes[3]],
                [nodes[4], nodes[5], nodes[6], nodes[7]],
                [nodes[0], nodes[1], nodes[5], nodes[4]],
                [nodes[1], nodes[2], nodes[6], nodes[5]],
                [nodes[2], nodes[3], nodes[7], nodes[6]],
                [nodes[0], nodes[3], nodes[7], nodes[4]],
            ]
        elif n_nodes == 3:
            # Triangle (2D): edges are faces
            return [
                [nodes[0], nodes[1]],
                [nodes[1], nodes[2]],
                [nodes[0], nodes[2]],
            ]
        else:
            return []

    def get_connectivity_matrix(self, dim=3):
        """Build sparse element connectivity matrix (CSR format).

        Entry (i, j) is 1 if elements i and j share a face.
        Uses 0-based indexing corresponding to the element order
        returned by get_elements(dim).

        Parameters
        ----------
        dim : int
            Element dimension (default 3).

        Returns
        -------
        scipy.sparse.csr_matrix
            (n_elem, n_elem) symmetric connectivity matrix.
        """
        import gmsh
        from scipy import sparse
        self._check_initialized()

        elements = self.get_elements(dim)
        n_elem = len(elements)
        if n_elem == 0:
            return sparse.csr_matrix((0, 0), dtype=int)

        # Map element tag to index
        tag_to_idx = {e['tag']: i for i, e in enumerate(elements)}

        # Build face-to-element map
        face_to_elems = {}
        for i, e in enumerate(elements):
            nodes = e['nodes']
            n_nodes = len(nodes)
            _, e_dim, _, _, _, _ = gmsh.model.mesh.getElementProperties(
                e['type'])
            faces = self._local_faces(nodes, n_nodes, e_dim)
            for f in faces:
                key = tuple(sorted(f))
                if key not in face_to_elems:
                    face_to_elems[key] = []
                face_to_elems[key].append(i)

        # Build adjacency
        rows = []
        cols = []
        for face, elems in face_to_elems.items():
            for a in range(len(elems)):
                for b in range(a + 1, len(elems)):
                    rows.append(elems[a])
                    cols.append(elems[b])
                    rows.append(elems[b])
                    cols.append(elems[a])

        data = np.ones(len(rows), dtype=int)
        return sparse.csr_matrix(
            (data, (np.array(rows), np.array(cols))),
            shape=(n_elem, n_elem))

    def get_elements_in_sphere(self, center, radius, dim=3):
        """Get elements whose centroid lies within a sphere.

        Parameters
        ----------
        center : list or array
            [x, y, z] sphere center.
        radius : float
            Sphere radius.
        dim : int
            Element dimension (default 3).

        Returns
        -------
        list of dict
            Element info dicts for elements inside the sphere.
        """
        import gmsh
        self._check_initialized()

        center = np.array(center, dtype=float)
        elements = self.get_elements(dim)
        result = []

        for e in elements:
            # Compute centroid from node coordinates
            coords = np.empty((len(e['nodes']), 3), dtype=float)
            for i, n in enumerate(e['nodes']):
                node_data = gmsh.model.mesh.getNode(n)
                coords[i] = node_data[0]
            centroid = coords.mean(axis=0)

            if np.linalg.norm(centroid - center) <= radius:
                result.append(e)

        return result

    # ------------------------------------------------------------------
    # Element data access (12 methods)
    # ------------------------------------------------------------------

    def get_element_edge_nodes(self, elem_type, tag=-1, primary=True):
        """Get edge nodes for elements of a given type.

        Parameters
        ----------
        elem_type : int
            GMSH element type (e.g. 4 for tet4, 5 for hex8).
        tag : int, optional
            Entity tag filter. -1 for all entities.
        primary : bool, optional
            If True, return only primary (corner) nodes.

        Returns
        -------
        numpy.ndarray
            1D array of node tags forming element edges, flattened.
        """
        import gmsh
        self._check_initialized()

        self.log(f"Getting element edge nodes for type {elem_type}, "
                 f"tag={tag}, primary={primary}")
        edge_nodes = gmsh.model.mesh.getElementEdgeNodes(
            elem_type, tag, primary)
        return np.array(edge_nodes, dtype=int)

    def get_element_face_nodes(self, elem_type, face_type,
                               tag=-1, primary=True):
        """Get face nodes for elements of a given type.

        Parameters
        ----------
        elem_type : int
            GMSH element type (e.g. 4 for tet4, 5 for hex8).
        face_type : int
            Number of corners per face (3 for triangular, 4 for quad).
        tag : int, optional
            Entity tag filter. -1 for all entities.
        primary : bool, optional
            If True, return only primary (corner) nodes.

        Returns
        -------
        numpy.ndarray
            1D array of node tags forming element faces, flattened.
        """
        import gmsh
        self._check_initialized()

        self.log(f"Getting element face nodes for type {elem_type}, "
                 f"face_type={face_type}, tag={tag}, primary={primary}")
        face_nodes = gmsh.model.mesh.getElementFaceNodes(
            elem_type, face_type, tag, primary)
        return np.array(face_nodes, dtype=int)

    def get_element_by_coordinates(self, x, y, z, dim=3):
        """Find the element containing the given point.

        Parameters
        ----------
        x : float
            X coordinate.
        y : float
            Y coordinate.
        z : float
            Z coordinate.
        dim : int, optional
            Element dimension to search (default 3).

        Returns
        -------
        dict
            Dictionary with keys 'element_tag', 'element_type',
            'node_tags', 'u', 'v', 'w' (parametric coordinates).
        """
        import gmsh
        self._check_initialized()

        self.log(f"Finding element at coordinates ({x}, {y}, {z}), dim={dim}")
        elem_tag, elem_type, node_tags, u, v, w = \
            gmsh.model.mesh.getElementByCoordinates(x, y, z, dim)
        return {
            'element_tag': int(elem_tag),
            'element_type': int(elem_type),
            'node_tags': [int(n) for n in node_tags],
            'u': float(u),
            'v': float(v),
            'w': float(w),
        }

    def get_elements_by_coordinates(self, x, y, z, dim=3):
        """Find elements at given coordinates.

        Note: Uses single-point lookup. For batch lookups, call
        get_element_by_coordinates() in a loop.

        Parameters
        ----------
        x : float
            X coordinate.
        y : float
            Y coordinate.
        z : float
            Z coordinate.
        dim : int, optional
            Element dimension to search (default 3).

        Returns
        -------
        dict or None
            Dict with 'element_tag', 'element_type', 'node_tags',
            'u', 'v', 'w', or None if no element found.
        """
        import gmsh
        self._check_initialized()

        self.log(f"Finding element at coordinates ({x}, {y}, {z}), "
                 f"dim={dim}")
        try:
            result = gmsh.model.mesh.getElementByCoordinates(x, y, z, dim)
            return {
                'element_tag': int(result[0]),
                'element_type': int(result[1]),
                'node_tags': [int(n) for n in result[2]],
                'u': float(result[3]),
                'v': float(result[4]),
                'w': float(result[5]),
            }
        except Exception as e:
            if self._verbose:
                print(f"[GmshBuilder] get_elements_by_coordinates: {e}")
            return None

    def get_element_properties(self, elem_type):
        """Get properties of a GMSH element type.

        Parameters
        ----------
        elem_type : int
            GMSH element type integer (e.g. 4 for tet4, 5 for hex8).

        Returns
        -------
        dict
            Dictionary with keys 'name', 'dim', 'order', 'n_nodes',
            'local_coords', 'n_primary_nodes'.
        """
        import gmsh
        self._check_initialized()

        self.log(f"Getting properties for element type {elem_type}")
        name, dim, order, n_nodes, local_coords, n_primary_nodes = \
            gmsh.model.mesh.getElementProperties(elem_type)
        return {
            'name': name,
            'dim': int(dim),
            'order': int(order),
            'n_nodes': int(n_nodes),
            'local_coords': np.array(local_coords, dtype=float),
            'n_primary_nodes': int(n_primary_nodes),
        }

    def get_max_element_tag(self):
        """Get the maximum element tag in the mesh.

        Returns
        -------
        int
            Maximum element tag.
        """
        import gmsh
        self._check_initialized()

        max_tag = gmsh.model.mesh.getMaxElementTag()
        self.log(f"Max element tag: {max_tag}")
        return int(max_tag)

    def get_max_node_tag(self):
        """Get the maximum node tag in the mesh.

        Returns
        -------
        int
            Maximum node tag.
        """
        import gmsh
        self._check_initialized()

        max_tag = gmsh.model.mesh.getMaxNodeTag()
        self.log(f"Max node tag: {max_tag}")
        return int(max_tag)

    def get_element_types(self, dim=-1, tag=-1):
        """Get all element types present in the mesh.

        Parameters
        ----------
        dim : int, optional
            Dimension filter. -1 for all dimensions.
        tag : int, optional
            Entity tag filter. -1 for all entities.

        Returns
        -------
        list of int
            GMSH element type integers.
        """
        import gmsh
        self._check_initialized()

        self.log(f"Getting element types for dim={dim}, tag={tag}")
        elem_types = gmsh.model.mesh.getElementTypes(dim, tag)
        result = [int(t) for t in elem_types]
        self.log(f"Found {len(result)} element types: {result}")
        return result

    def get_elements_by_type_raw(self, elem_type, tag=-1):
        """Get elements of a given type as raw arrays.

        Parameters
        ----------
        elem_type : int
            GMSH element type.
        tag : int, optional
            Entity tag filter. -1 for all entities.

        Returns
        -------
        tuple
            (element_tags, node_tags) as numpy arrays.
        """
        import gmsh
        self._check_initialized()

        self.log(f"Getting raw elements by type {elem_type}, tag={tag}")
        elem_tags, node_tags = gmsh.model.mesh.getElementsByType(
            elem_type, tag)
        self.log(f"Found {len(elem_tags)} elements of type {elem_type}")
        return (np.array(elem_tags, dtype=int),
                np.array(node_tags, dtype=int))

    def get_local_coords_in_element(self, elem_tag, x, y, z):
        """Get local parametric coordinates within an element.

        Uses ``getElementByCoordinates`` to find parametric coordinates
        at the given physical point.  If the point lies in a different
        element than *elem_tag*, the parametric coordinates in that
        element are still returned (with a verbose warning).

        Parameters
        ----------
        elem_tag : int
            GMSH element tag.
        x : float
            X coordinate of point.
        y : float
            Y coordinate of point.
        z : float
            Z coordinate of point.

        Returns
        -------
        list
            Local coordinates [u, v, w].
        """
        import gmsh
        self._check_initialized()

        self.log(f"Getting local coordinates in element {elem_tag} "
                 f"at ({x}, {y}, {z})")
        try:
            # Use getElementByCoordinates to find parametric coords
            tag, etype, nodes, u, v, w = gmsh.model.mesh.getElementByCoordinates(x, y, z)
            if tag != elem_tag:
                # Point found in different element
                if self._verbose:
                    print(f"[GmshBuilder] get_local_coords_in_element: "
                          f"point in element {tag}, not {elem_tag}")
            return [float(u), float(v), float(w)]
        except Exception as e:
            raise RuntimeError(f"Cannot find local coordinates: {e}")

    def get_nodes_for_physical_group(self, dim, tag):
        """Get nodes belonging to a physical group.

        Parameters
        ----------
        dim : int
            Dimension of the physical group.
        tag : int
            Physical group tag.

        Returns
        -------
        tuple
            (node_tags, coords) where node_tags is a 1D numpy array
            and coords is a (n, 3) numpy array.
        """
        import gmsh
        self._check_initialized()

        self.log(f"Getting nodes for physical group dim={dim}, tag={tag}")
        node_tags, coords = gmsh.model.mesh.getNodesForPhysicalGroup(
            dim, tag)
        node_tags_arr = np.array(node_tags, dtype=int)
        coords_arr = np.array(coords, dtype=float).reshape(-1, 3)
        self.log(f"Found {len(node_tags_arr)} nodes in physical "
                 f"group ({dim}, {tag})")
        return (node_tags_arr, coords_arr)

    def get_nodes_by_element_type(self, elem_type, tag=-1):
        """Get nodes associated with elements of a given type.

        Parameters
        ----------
        elem_type : int
            GMSH element type.
        tag : int, optional
            Entity tag filter. -1 for all entities.

        Returns
        -------
        tuple
            (node_tags, coords, parametric_coords) where node_tags is
            a 1D int array, coords is a (n, 3) float array, and
            parametric_coords is a 1D float array.
        """
        import gmsh
        self._check_initialized()

        self.log(f"Getting nodes by element type {elem_type}, tag={tag}")
        node_tags, coords, parametric_coords = \
            gmsh.model.mesh.getNodesByElementType(elem_type, tag)
        self.log(f"Found {len(node_tags)} nodes for element type {elem_type}")
        return (np.array(node_tags, dtype=int),
                np.array(coords, dtype=float).reshape(-1, 3),
                np.array(parametric_coords, dtype=float))

    # ------------------------------------------------------------------
    # Jacobian / basis functions (8 methods)
    # ------------------------------------------------------------------

    def get_jacobian(self, elem_type, local_coords, tag=-1):
        """Get Jacobian matrices at specified local coordinates.

        Parameters
        ----------
        elem_type : int
            GMSH element type.
        local_coords : list or array
            Local (parametric) coordinates, flattened [u1, v1, w1, ...].
        tag : int, optional
            Entity tag filter. -1 for all entities.

        Returns
        -------
        dict
            'jacobians': flat numpy array of 3x3 Jacobians,
            'determinants': numpy array of Jacobian determinants,
            'coords': flat numpy array of physical coordinates.
        """
        import gmsh
        self._check_initialized()

        self.log(f"Computing Jacobians for element type {elem_type}, "
                 f"tag={tag}")
        jacobians, determinants, coords = gmsh.model.mesh.getJacobians(
            elem_type, local_coords, tag)
        return {
            'jacobians': np.array(jacobians, dtype=float),
            'determinants': np.array(determinants, dtype=float),
            'coords': np.array(coords, dtype=float),
        }

    def get_jacobians(self, elem_type, local_coords, tag=-1):
        """Get Jacobian matrices (alias for get_jacobian).

        Parameters
        ----------
        elem_type : int
            GMSH element type.
        local_coords : list or array
            Local (parametric) coordinates, flattened [u1, v1, w1, ...].
        tag : int, optional
            Entity tag filter. -1 for all entities.

        Returns
        -------
        dict
            'jacobians': flat numpy array of 3x3 Jacobians,
            'determinants': numpy array of Jacobian determinants,
            'coords': flat numpy array of physical coordinates.
        """
        return self.get_jacobian(elem_type, local_coords, tag)

    def get_basis_functions(self, elem_type, local_coords,
                            function_type='Lagrange'):
        """Get basis function values at specified local coordinates.

        Parameters
        ----------
        elem_type : int
            GMSH element type.
        local_coords : list or array
            Local coordinates, flattened [u1, v1, w1, ...].
        function_type : str, optional
            Basis function type (default 'Lagrange'). Other options
            include 'GradLagrange', etc.

        Returns
        -------
        dict
            'n_components': int, number of components per function,
            'functions': numpy array of basis function values,
            'n_orient': int, number of orientations.
        """
        import gmsh
        self._check_initialized()

        self.log(f"Getting {function_type} basis functions for element "
                 f"type {elem_type}")
        n_components, functions, n_orient = \
            gmsh.model.mesh.getBasisFunctions(
                elem_type, local_coords, function_type)
        return {
            'n_components': int(n_components),
            'functions': np.array(functions, dtype=float),
            'n_orient': int(n_orient),
        }

    def get_integration_points(self, elem_type, order):
        """Get Gauss integration points and weights.

        Parameters
        ----------
        elem_type : int
            GMSH element type.
        order : int
            Quadrature order.

        Returns
        -------
        dict
            'local_coords': numpy array of local coordinates (flattened),
            'weights': numpy array of integration weights.
        """
        import gmsh
        self._check_initialized()

        self.log(f"Getting Gauss{order} integration points for element "
                 f"type {elem_type}")
        local_coords, weights = gmsh.model.mesh.getIntegrationPoints(
            elem_type, 'Gauss' + str(order))
        return {
            'local_coords': np.array(local_coords, dtype=float),
            'weights': np.array(weights, dtype=float),
        }

    def get_element_barycenters(self, elem_type, tag=-1, fast=True,
                                primary=True):
        """Get barycenters (centroids) of all elements of a given type.

        Parameters
        ----------
        elem_type : int
            GMSH element type.
        tag : int, optional
            Entity tag filter. -1 for all entities.
        fast : bool, optional
            Use fast computation (default True).
        primary : bool, optional
            Use only primary nodes (default True).

        Returns
        -------
        numpy.ndarray
            (n, 3) array of barycenter coordinates.
        """
        import gmsh
        self._check_initialized()

        self.log(f"Getting barycenters for element type {elem_type}, "
                 f"tag={tag}, fast={fast}")
        barycenters = gmsh.model.mesh.getBarycenters(
            elem_type, tag, fast, primary)
        result = np.array(barycenters, dtype=float).reshape(-1, 3)
        self.log(f"Got {len(result)} barycenters")
        return result

    def get_element_sizes(self, dim=-1, tag=-1):
        """Compute element sizes (characteristic lengths).

        For each element, the size is computed as the maximum edge length
        among all edges of the element.

        Parameters
        ----------
        dim : int, optional
            Element dimension. -1 for all dimensions.
        tag : int, optional
            Entity tag filter. -1 for all entities.

        Returns
        -------
        numpy.ndarray
            1D array of element sizes (one per element).
        """
        import gmsh
        self._check_initialized()

        self.log(f"Computing element sizes for dim={dim}, tag={tag}")

        # Edge definitions keyed by element type name to avoid
        # Quad4/Tet4 ambiguity (both have 4 nodes).
        edge_defs_by_name = {
            'Line': [(0, 1)],
            'Triangle': [(0, 1), (1, 2), (2, 0)],
            'Quadrilateral': [(0, 1), (1, 2), (2, 3), (3, 0)],
            'Quadrangle': [(0, 1), (1, 2), (2, 3), (3, 0)],
            'Tetrahedron': [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],
            'Pyramid': [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 4),
                        (2, 4), (3, 4)],
            'Prism': [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3),
                      (0, 3), (1, 4), (2, 5)],
            'Hexahedron': [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6),
                           (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)],
        }

        elem_types, elem_tags, elem_node_tags = \
            gmsh.model.mesh.getElements(dim, tag)

        sizes = []
        for i, et in enumerate(elem_types):
            name, _, _, n_nodes, _, _ = \
                gmsh.model.mesh.getElementProperties(et)
            nodes = elem_node_tags[i]
            n_elems = len(elem_tags[i])
            # Match element name to edge definitions
            local_edges = []
            for key, val in edge_defs_by_name.items():
                if key in name:
                    local_edges = val
                    break

            for j in range(n_elems):
                start = j * n_nodes
                elem_nodes = nodes[start:start + n_nodes]
                # Get coordinates for all nodes of this element
                coords = np.empty((n_nodes, 3), dtype=float)
                for k in range(n_nodes):
                    c, _, _, _ = gmsh.model.mesh.getNode(int(elem_nodes[k]))
                    coords[k] = c

                # Compute max edge length
                max_len = 0.0
                for e0, e1 in local_edges:
                    edge_len = np.linalg.norm(coords[e1] - coords[e0])
                    if edge_len > max_len:
                        max_len = edge_len
                sizes.append(max_len)

        result = np.array(sizes, dtype=float)
        self.log(f"Computed {len(result)} element sizes")
        return result

    def get_element_type_for_name(self, name):
        """Get GMSH element type code for a given name.

        Parses element names like ``'Hexahedron 8'`` into the family
        string and finite-element order expected by
        ``gmsh.model.mesh.getElementType(family, order)``.

        Parameters
        ----------
        name : str
            Element name, e.g. 'Hexahedron 8', 'Tetrahedron 4',
            'Triangle 3', 'Line 2'.

        Returns
        -------
        int
            GMSH element type integer.
        """
        import gmsh
        self._check_initialized()

        self.log(f"Getting element type for name '{name}'")

        # Map common element names to (family, order) or
        # (family, order, serendip) tuples
        name_map = {
            'Hexahedron 8': ('Hexahedron', 1),
            'Hexahedron 20': ('Hexahedron', 2, True),
            'Hexahedron 27': ('Hexahedron', 2),
            'Tetrahedron 4': ('Tetrahedron', 1),
            'Tetrahedron 10': ('Tetrahedron', 2),
            'Triangle 3': ('Triangle', 1),
            'Triangle 6': ('Triangle', 2),
            'Quadrangle 4': ('Quadrangle', 1),
            'Quadrangle 8': ('Quadrangle', 2, True),
            'Quadrangle 9': ('Quadrangle', 2),
            'Quadrilateral 4': ('Quadrangle', 1),
            'Quadrilateral 8': ('Quadrangle', 2, True),
            'Quadrilateral 9': ('Quadrangle', 2),
            'Line 2': ('Line', 1),
            'Line 3': ('Line', 2),
            'Prism 6': ('Prism', 1),
            'Prism 15': ('Prism', 2, True),
            'Prism 18': ('Prism', 2),
            'Pyramid 5': ('Pyramid', 1),
            'Pyramid 13': ('Pyramid', 2, True),
            'Pyramid 14': ('Pyramid', 2),
            'Point': ('Point', 1),
        }

        serendip = False
        if name in name_map:
            entry = name_map[name]
            family = entry[0]
            order = entry[1]
            if len(entry) > 2:
                serendip = entry[2]
        else:
            # Try splitting on last space to extract family name
            parts = name.rsplit(' ', 1)
            family = parts[0]
            # Node count -> order lookup for common element families
            node_order_table = {
                'Line': {2: 1, 3: 2, 4: 3},
                'Triangle': {3: 1, 6: 2, 10: 3},
                'Quadrangle': {4: 1, 9: 2, 8: (2, True)},
                'Tetrahedron': {4: 1, 10: 2, 20: 3},
                'Hexahedron': {8: 1, 27: 2, 20: (2, True)},
                'Prism': {6: 1, 18: 2, 15: (2, True)},
                'Pyramid': {5: 1, 14: 2, 13: (2, True)},
            }
            order = 1
            if len(parts) > 1 and parts[1].isdigit():
                n_count = int(parts[1])
                family_orders = node_order_table.get(family, {})
                if n_count in family_orders:
                    val = family_orders[n_count]
                    if isinstance(val, tuple):
                        order = val[0]
                        serendip = val[1]
                    else:
                        order = val
                else:
                    # Try common orders until one works
                    for try_order in (1, 2, 3):
                        try:
                            test_type = gmsh.model.mesh.getElementType(
                                family, try_order, serendip=serendip)
                            if test_type > 0:
                                order = try_order
                                break
                        except Exception:
                            continue

        elem_type = gmsh.model.mesh.getElementType(
            family, order, serendip=serendip)
        self.log(f"Element type for '{name}': {elem_type}")
        return int(elem_type)

    def get_nodes_per_element(self, elem_type):
        """Get the number of nodes for an element type.

        Parameters
        ----------
        elem_type : int
            GMSH element type.

        Returns
        -------
        int
            Number of nodes per element.
        """
        import gmsh
        self._check_initialized()

        self.log(f"Getting node count for element type {elem_type}")
        _, _, _, n_nodes, _, _ = gmsh.model.mesh.getElementProperties(
            elem_type)
        self.log(f"Element type {elem_type} has {n_nodes} nodes")
        return int(n_nodes)

    # ------------------------------------------------------------------
    # Node manipulation (5 methods)
    # ------------------------------------------------------------------

    def set_node_coordinates(self, node_tag, x, y, z):
        """Set the coordinates of a specific node.

        Parameters
        ----------
        node_tag : int
            GMSH node tag.
        x : float
            New X coordinate.
        y : float
            New Y coordinate.
        z : float
            New Z coordinate.
        """
        import gmsh
        self._check_initialized()

        self.log(f"Setting node {node_tag} coordinates to "
                 f"({x}, {y}, {z})")
        gmsh.model.mesh.setNode(node_tag, [x, y, z], [])

    def add_nodes_manual(self, dim, tag, coords, tags=None):
        """Add nodes manually to the mesh.

        Parameters
        ----------
        dim : int
            Dimension of the entity to add nodes to.
        tag : int
            Entity tag.
        coords : list or array
            Flat list of coordinates [x1, y1, z1, x2, y2, z2, ...].
        tags : list of int, optional
            Node tags. If empty or None, tags are assigned automatically.

        Returns
        -------
        list of int
            Assigned node tags.
        """
        import gmsh
        self._check_initialized()

        if tags is None:
            tags = []

        n_nodes = len(coords) // 3
        self.log(f"Adding {n_nodes} nodes manually to entity "
                 f"({dim}, {tag})")
        gmsh.model.mesh.addNodes(dim, tag, tags, coords)
        self.log(f"Added {n_nodes} nodes")

        # Return the tags that were assigned
        if tags:
            return list(tags)
        # If no tags specified, retrieve the latest nodes
        all_tags, _, _ = gmsh.model.mesh.getNodes(dim, tag)
        return [int(t) for t in all_tags[-n_nodes:]]

    def add_elements_manual(self, dim, tag, elem_type, node_tags,
                            elem_tags=None):
        """Add elements manually to the mesh.

        Parameters
        ----------
        dim : int
            Dimension of the entity (unused by GMSH but kept for API
            consistency; elements are added by entity tag).
        tag : int
            Entity tag to add elements to.
        elem_type : int
            GMSH element type.
        node_tags : list of int
            Flat list of node tags for all elements.
        elem_tags : list of int, optional
            Element tags. If empty or None, tags are assigned
            automatically.

        Returns
        -------
        list of int
            Assigned element tags.
        """
        import gmsh
        self._check_initialized()

        if elem_tags is None:
            elem_tags = []

        _, _, _, n_nodes, _, _ = gmsh.model.mesh.getElementProperties(
            elem_type)
        n_elems = len(node_tags) // n_nodes
        self.log(f"Adding {n_elems} elements of type {elem_type} "
                 f"to entity tag {tag}")
        gmsh.model.mesh.addElementsByType(tag, elem_type, elem_tags,
                                          node_tags)
        self.log(f"Added {n_elems} elements")

        if elem_tags:
            return list(elem_tags)
        # Retrieve elements to get assigned tags
        all_tags, _ = gmsh.model.mesh.getElementsByType(elem_type, tag)
        return [int(t) for t in all_tags[-n_elems:]]

    def remove_duplicate_nodes(self):
        """Remove duplicate nodes from the mesh."""
        import gmsh
        self._check_initialized()

        self.log("Removing duplicate nodes")
        gmsh.model.mesh.removeDuplicateNodes()
        self.log("Duplicate nodes removed")

    def remove_duplicate_elements(self):
        """Remove duplicate elements from the mesh."""
        import gmsh
        self._check_initialized()

        self.log("Removing duplicate elements")
        gmsh.model.mesh.removeDuplicateElements()
        self.log("Duplicate elements removed")

    # ------------------------------------------------------------------
    # Periodic mesh (5 methods)
    # ------------------------------------------------------------------

    def set_periodic_mesh(self, dim, tags, source_tags, affine_transform):
        """Set periodic mesh constraints.

        Parameters
        ----------
        dim : int
            Dimension of the entities.
        tags : list of int
            Target entity tags (slave surfaces).
        source_tags : list of int
            Source entity tags (master surfaces).
        affine_transform : list of float
            16-element affine transformation matrix (4x4, row-major).
            Maps source coordinates to target coordinates.
        """
        import gmsh
        self._check_initialized()

        self.log(f"Setting periodic mesh: dim={dim}, "
                 f"tags={tags}, source_tags={source_tags}")
        gmsh.model.mesh.setPeriodic(dim, tags, source_tags,
                                    affine_transform)
        self.log("Periodic mesh constraints set")

    def get_periodic_nodes(self, dim, tag):
        """Get periodic node correspondences for an entity.

        Parameters
        ----------
        dim : int
            Dimension of the entity.
        tag : int
            Entity tag (slave surface).

        Returns
        -------
        dict
            'tag_master': int, master entity tag,
            'node_tags': numpy array of slave node tags,
            'node_tags_master': numpy array of corresponding master
            node tags,
            'affine': list of 16 floats (4x4 affine transform).
        """
        import gmsh
        self._check_initialized()

        self.log(f"Getting periodic nodes for entity ({dim}, {tag})")
        tag_master, node_tags, node_tags_master, affine = \
            gmsh.model.mesh.getPeriodicNodes(dim, tag)
        self.log(f"Found {len(node_tags)} periodic node pairs, "
                 f"master entity: {tag_master}")
        return {
            'tag_master': int(tag_master),
            'node_tags': np.array(node_tags, dtype=int),
            'node_tags_master': np.array(node_tags_master, dtype=int),
            'affine': list(affine),
        }

    def get_periodic_keys(self, elem_type, function_type='Lagrange', tag=-1):
        """Get periodic keys for DOF mapping.

        Parameters
        ----------
        elem_type : int
            GMSH element type.
        function_type : str, optional
            Function space type (default 'Lagrange').
        tag : int, optional
            Entity tag (-1 for all entities).

        Returns
        -------
        dict
            Result from gmsh.model.mesh.getPeriodicKeys containing
            'tag_master', 'type_keys', 'type_keys_master',
            'entity_keys', 'entity_keys_master', 'coord', 'coord_master'.
        """
        import gmsh
        self._check_initialized()

        self.log(f"Getting periodic keys for element type {elem_type}, "
                 f"function_type='{function_type}', tag={tag}")
        (tagMaster, typeKeys, typeKeysMaster, entityKeys,
         entityKeysMaster, coord, coordMaster) = \
            gmsh.model.mesh.getPeriodicKeys(elem_type, function_type, tag)
        return {
            'tag_master': int(tagMaster),
            'type_keys': list(typeKeys),
            'type_keys_master': list(typeKeysMaster),
            'entity_keys': list(entityKeys),
            'entity_keys_master': list(entityKeysMaster),
            'coord': list(coord),
            'coord_master': list(coordMaster),
        }

    def get_duplicate_nodes(self):
        """Get duplicate node tags in the mesh.

        Returns
        -------
        list of int
            List of duplicate node tags. Returns empty list if the
            API is not available.
        """
        import gmsh
        self._check_initialized()

        self.log("Getting duplicate nodes")
        try:
            result = gmsh.model.mesh.getDuplicateNodes()
            tags = list(result)
            self.log(f"Found {len(tags)} duplicate node tags")
            return tags
        except AttributeError:
            self.log("getDuplicateNodes not available in this GMSH version")
            return []
        except Exception as e:
            self.log(f"getDuplicateNodes raised: {e}")
            return []

    def relocate_nodes(self, dim=-1, tag=-1):
        """Relocate nodes to their associated geometry.

        Parameters
        ----------
        dim : int, optional
            Dimension filter. -1 for all dimensions.
        tag : int, optional
            Entity tag filter. -1 for all entities.
        """
        import gmsh
        self._check_initialized()

        self.log(f"Relocating nodes for dim={dim}, tag={tag}")
        gmsh.model.mesh.relocateNodes(dim, tag)
        self.log("Nodes relocated")

    # ------------------------------------------------------------------
    # Volume-specific mesh access (5 methods)
    # ------------------------------------------------------------------

    def get_volume_node_count(self, vol_id):
        """Get the number of nodes in a specific volume.

        Parameters
        ----------
        vol_id : int
            Volume ID (user-assigned).

        Returns
        -------
        int
            Number of nodes in the volume.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        self.log(f"Getting node count for volume {vol_id}")
        node_set = set()
        for vt in self._volumes[vol_id]:
            tags, _, _ = gmsh.model.mesh.getNodes(3, vt)
            for t in tags:
                node_set.add(int(t))
        count = len(node_set)
        self.log(f"Volume {vol_id} has {count} nodes")
        return count

    def get_volume_element_count(self, vol_id):
        """Get the total number of elements in a specific volume.

        Parameters
        ----------
        vol_id : int
            Volume ID (user-assigned).

        Returns
        -------
        int
            Total number of elements in the volume.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        self.log(f"Getting element count for volume {vol_id}")
        count = 0
        for vt in self._volumes[vol_id]:
            _, elem_tags, _ = gmsh.model.mesh.getElements(3, vt)
            for tags in elem_tags:
                count += len(tags)
        self.log(f"Volume {vol_id} has {count} elements")
        return count

    def get_volume_hex_count(self, vol_id):
        """Get the number of hexahedral elements in a specific volume.

        Parameters
        ----------
        vol_id : int
            Volume ID (user-assigned).

        Returns
        -------
        int
            Number of hexahedral (type 5) elements.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        self.log(f"Getting hex count for volume {vol_id}")
        count = 0
        for vt in self._volumes[vol_id]:
            elem_types, elem_tags, _ = gmsh.model.mesh.getElements(3, vt)
            for i, et in enumerate(elem_types):
                if int(et) == 5:  # Hex8
                    count += len(elem_tags[i])
        self.log(f"Volume {vol_id} has {count} hex elements")
        return count

    def get_volume_tet_count(self, vol_id):
        """Get the number of tetrahedral elements in a specific volume.

        Parameters
        ----------
        vol_id : int
            Volume ID (user-assigned).

        Returns
        -------
        int
            Number of tetrahedral (type 4) elements.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        self.log(f"Getting tet count for volume {vol_id}")
        count = 0
        for vt in self._volumes[vol_id]:
            elem_types, elem_tags, _ = gmsh.model.mesh.getElements(3, vt)
            for i, et in enumerate(elem_types):
                if int(et) == 4:  # Tet4
                    count += len(elem_tags[i])
        self.log(f"Volume {vol_id} has {count} tet elements")
        return count

    def get_surface_node_count(self, surf_tag):
        """Get the number of nodes on a specific surface.

        Parameters
        ----------
        surf_tag : int
            GMSH surface entity tag.

        Returns
        -------
        int
            Number of nodes on the surface.
        """
        import gmsh
        self._check_initialized()

        self.log(f"Getting node count for surface {surf_tag}")
        tags, _, _ = gmsh.model.mesh.getNodes(2, surf_tag)
        count = len(tags)
        self.log(f"Surface {surf_tag} has {count} nodes")
        return count

    # ------------------------------------------------------------------
    # Edge / face global access (5 methods)
    # ------------------------------------------------------------------

    def create_edges(self):
        """Create unique edges for the mesh.

        This must be called before get_all_edges() or
        get_edges_for_element().
        """
        import gmsh
        self._check_initialized()

        self.log("Creating edges for the mesh")
        gmsh.model.mesh.createEdges()
        self.log("Edges created")

    def create_faces(self, dim=-1, tag=-1):
        """Create unique faces for the mesh.

        Parameters
        ----------
        dim : int, optional
            Entity dimension (-1 for all entities).
        tag : int, optional
            Entity tag (-1 for all entities).
        """
        import gmsh
        self._check_initialized()

        if dim == -1 and tag == -1:
            self.log("Creating faces for the mesh (all entities)")
            gmsh.model.mesh.createFaces()
        else:
            self.log(f"Creating faces for entity ({dim}, {tag})")
            gmsh.model.mesh.createFaces([(dim, tag)])
        self.log("Faces created")

    def get_all_edges(self):
        """Get all unique edges in the mesh.

        Must call create_edges() first.

        Returns
        -------
        tuple
            (edge_tags, edge_nodes) where edge_tags is a 1D numpy
            array of edge tags and edge_nodes is a 1D numpy array
            of node tags (pairs, flattened).
        """
        import gmsh
        self._check_initialized()

        self.log("Getting all edges")
        edge_tags, edge_nodes = gmsh.model.mesh.getAllEdges()
        self.log(f"Found {len(edge_tags)} edges")
        return (np.array(edge_tags, dtype=int),
                np.array(edge_nodes, dtype=int))

    def get_all_faces(self, face_type=3):
        """Get all unique faces in the mesh.

        Must call create_faces() first.

        Parameters
        ----------
        face_type : int, optional
            Face type: 3 for triangles, 4 for quads (default 3).

        Returns
        -------
        dict
            Dictionary with 'tags' (list of face tags) and
            'nodes' (list of node tags, flattened).
        """
        import gmsh
        self._check_initialized()

        self.log(f"Getting all faces (face_type={face_type})")
        face_tags, face_nodes = gmsh.model.mesh.getAllFaces(face_type)
        self.log(f"Found {len(face_tags)} faces")
        return {'tags': list(face_tags), 'nodes': list(face_nodes)}

    def get_edges_for_element(self, elem_tag):
        """Get edge tags for a specific element.

        Must call create_edges() first.

        Parameters
        ----------
        elem_tag : int
            GMSH element tag.

        Returns
        -------
        list of int
            Edge tags belonging to the element.
        """
        import gmsh
        self._check_initialized()

        self.log(f"Getting edges for element {elem_tag}")
        try:
            # Get element type and nodes
            elem_type, node_tags_arr, _, _ = gmsh.model.mesh.getElement(
                elem_tag)
            node_tags_list = list(node_tags_arr)

            # Standard edge connectivity per element type
            # (local node index pairs, 0-based)
            edge_connectivity = {
                2: [(0, 1), (1, 2), (2, 0)],  # Tri3
                3: [(0, 1), (1, 2), (2, 3), (3, 0)],  # Quad4
                4: [(0, 1), (0, 2), (0, 3),
                    (1, 2), (1, 3), (2, 3)],  # Tet4
                5: [(0, 1), (1, 2), (2, 3), (3, 0),
                    (4, 5), (5, 6), (6, 7), (7, 4),
                    (0, 4), (1, 5), (2, 6), (3, 7)],  # Hex8
            }

            local_edges = edge_connectivity.get(int(elem_type))
            if local_edges is None:
                self.log(f"Unknown element type {elem_type} for edge lookup")
                return []

            # Build node pairs for defined edges only
            edge_nodes = []
            for i, j in local_edges:
                edge_nodes.extend([node_tags_list[i], node_tags_list[j]])

            try:
                edge_tags, _ = gmsh.model.mesh.getEdges(edge_nodes)
                result = [int(t) for t in edge_tags]
                self.log(f"Element {elem_tag} has {len(result)} edges")
                return result
            except Exception:
                return []
        except Exception as e:
            self.log(f"getEdges raised: {e}")
            raise
