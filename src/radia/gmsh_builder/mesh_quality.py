"""
Mesh quality evaluation for GmshBuilder.

Gmsh equivalents: quality vol, quality histogram, scaled jacobian.
"""

import numpy as np


class MeshQualityMixin:
    """Mixin providing mesh quality evaluation."""

    def _get_quality_data(self, metric='minSICN'):
        """Get quality data for all 3D elements.

        Parameters
        ----------
        metric : str
            GMSH quality metric: 'minSICN' (scaled Jacobian),
            'minSIGE' (inscribed/circumscribed ratio),
            'gamma' (aspect ratio), 'minAniso' (anisotropy).

        Returns
        -------
        numpy array
            Quality values for each element.
        """
        import gmsh
        self._check_initialized()

        if not self._meshed:
            raise RuntimeError("Mesh not generated yet")

        elem_types, _, _ = gmsh.model.mesh.getElements(3)
        all_qualities = []
        for et in elem_types:
            tags, _ = gmsh.model.mesh.getElementsByType(et)
            if len(tags) > 0:
                qualities = gmsh.model.mesh.getElementQualities(
                    tags, metric)
                all_qualities.extend(qualities)

        return np.array(all_qualities)

    def get_quality(self, metric='minSICN'):
        """Get quality array for all elements (Gmsh 'quality vol').

        Parameters
        ----------
        metric : str
            'minSICN' (scaled Jacobian), 'gamma' (aspect ratio),
            'minSIGE', 'minAniso'.

        Returns
        -------
        numpy array
            Quality values.
        """
        return self._get_quality_data(metric)

    def get_quality_stats(self, metric='minSICN'):
        """Get quality statistics (Gmsh 'quality stats').

        Returns
        -------
        dict
            min, max, mean, std, n_elements, n_negative.
        """
        q = self._get_quality_data(metric)
        if len(q) == 0:
            return {'min': 0, 'max': 0, 'mean': 0, 'std': 0,
                    'n_elements': 0, 'n_negative': 0}

        return {
            'min': float(np.min(q)),
            'max': float(np.max(q)),
            'mean': float(np.mean(q)),
            'std': float(np.std(q)),
            'n_elements': len(q),
            'n_negative': int(np.sum(q < 0)),
        }

    def get_quality_histogram(self, metric='minSICN', n_bins=10):
        """Get quality histogram (Gmsh 'quality histogram').

        Returns
        -------
        dict
            bin_edges, counts, percentages.
        """
        q = self._get_quality_data(metric)
        if len(q) == 0:
            return {'bin_edges': [], 'counts': [], 'percentages': []}

        counts, bin_edges = np.histogram(q, bins=n_bins)
        percentages = counts / len(q) * 100

        return {
            'bin_edges': bin_edges.tolist(),
            'counts': counts.tolist(),
            'percentages': percentages.tolist(),
        }

    def get_scaled_jacobian(self):
        """Get scaled Jacobian for all elements.

        Returns
        -------
        numpy array
        """
        return self._get_quality_data('minSICN')

    def get_aspect_ratio(self):
        """Get aspect ratio for all elements.

        Returns
        -------
        numpy array
        """
        return self._get_quality_data('gamma')

    def get_min_quality(self, metric='minSICN'):
        """Get minimum quality value.

        Returns
        -------
        float
        """
        q = self._get_quality_data(metric)
        return float(np.min(q)) if len(q) > 0 else 0.0

    def get_max_quality(self, metric='minSICN'):
        """Get maximum quality value.

        Returns
        -------
        float
        """
        q = self._get_quality_data(metric)
        return float(np.max(q)) if len(q) > 0 else 0.0

    def get_worst_elements(self, metric='minSICN', n=10):
        """Get worst elements by quality.

        Parameters
        ----------
        metric : str
            Quality metric.
        n : int
            Number of worst elements to return.

        Returns
        -------
        list of dict
            Element info with quality values.
        """
        import gmsh
        self._check_initialized()

        elem_types, _, _ = gmsh.model.mesh.getElements(3)
        all_info = []
        for et in elem_types:
            tags, _ = gmsh.model.mesh.getElementsByType(et)
            if len(tags) > 0:
                name = gmsh.model.mesh.getElementProperties(et)[0]
                qualities = gmsh.model.mesh.getElementQualities(tags, metric)
                for i, tag in enumerate(tags):
                    all_info.append({
                        'tag': int(tag),
                        'type': name,
                        'quality': float(qualities[i]),
                    })

        all_info.sort(key=lambda x: x['quality'])
        return all_info[:n]

    def check_inverted_elements(self):
        """Check for inverted (negative Jacobian) elements.

        Returns
        -------
        dict
            n_inverted, n_total, has_inverted.
        """
        q = self._get_quality_data('minSICN')
        inverted = np.where(q < 0)[0]

        return {
            'n_inverted': len(inverted),
            'n_total': len(q),
            'has_inverted': len(inverted) > 0,
        }

    def improve_quality(self, method=""):
        """Improve mesh quality (Gmsh 'smooth').

        Parameters
        ----------
        method : str
            Optimization method.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.optimize(method)
        if self._verbose:
            stats = self.get_quality_stats()
            print(f"[GmshBuilder] improve_quality: min={stats['min']:.4f} "
                  f"mean={stats['mean']:.4f}")

    def fix_negative_jacobians(self):
        """Attempt to fix negative Jacobian elements (untangle).

        Returns
        -------
        int
            Number of remaining inverted elements.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.optimize("UntangleMeshGeometry")
        check = self.check_inverted_elements()

        if self._verbose:
            print(f"[GmshBuilder] fix_negative_jacobians: "
                  f"{check['n_inverted']} remaining")
        return check['n_inverted']

    def get_element_quality(self, elem_tag, metric='minSICN'):
        """Get quality value for a single element.

        Parameters
        ----------
        elem_tag : int
            GMSH element tag.
        metric : str
            Quality metric.

        Returns
        -------
        float
            Quality value.
        """
        import gmsh
        self._check_initialized()

        qualities = gmsh.model.mesh.getElementQualities([elem_tag], metric)
        return float(qualities[0])

    def get_quality_by_type(self, metric='minSICN', elem_type='hex'):
        """Get quality for a specific element type only.

        Parameters
        ----------
        metric : str
            Quality metric: 'minSICN', 'minSIGE', 'gamma', 'minAniso'.
        elem_type : str
            Element type: 'hex', 'tet', 'wedge', 'prism', 'pyramid'.

        Returns
        -------
        numpy array
            Quality values for elements of the specified type.
        """
        import gmsh
        self._check_initialized()

        if not self._meshed:
            raise RuntimeError("Mesh not generated yet")

        type_map = {
            'hex': 'Hexahedron',
            'tet': 'Tetrahedron',
            'wedge': 'Prism',
            'prism': 'Prism',
            'pyramid': 'Pyramid',
        }
        target_name = type_map.get(elem_type, elem_type)

        elem_types, _, _ = gmsh.model.mesh.getElements(3)
        all_qualities = []
        for et in elem_types:
            name = gmsh.model.mesh.getElementProperties(et)[0]
            if target_name in name:
                tags, _ = gmsh.model.mesh.getElementsByType(et)
                if len(tags) > 0:
                    qualities = gmsh.model.mesh.getElementQualities(
                        tags, metric)
                    all_qualities.extend(qualities)

        return np.array(all_qualities)

    def get_quality_range(self, metric='minSICN'):
        """Return (min_val, max_val) quality range.

        Parameters
        ----------
        metric : str
            Quality metric.

        Returns
        -------
        tuple
            (min_val, max_val).
        """
        q = self._get_quality_data(metric)
        if len(q) == 0:
            return (0.0, 0.0)
        return (float(np.min(q)), float(np.max(q)))

    def get_quality_percentile(self, metric='minSICN', percentile=5):
        """Return quality value at the given percentile.

        Parameters
        ----------
        metric : str
            Quality metric.
        percentile : float
            Percentile (0-100). E.g. 5 means 5th percentile (worst 5%).

        Returns
        -------
        float
            Quality value at the given percentile.
        """
        q = self._get_quality_data(metric)
        if len(q) == 0:
            return 0.0
        return float(np.percentile(q, percentile))

    def count_elements_below_threshold(self, metric='minSICN', threshold=0.1):
        """Count elements with quality below a threshold.

        Parameters
        ----------
        metric : str
            Quality metric.
        threshold : float
            Quality threshold.

        Returns
        -------
        int
            Number of elements below threshold.
        """
        q = self._get_quality_data(metric)
        if len(q) == 0:
            return 0
        return int(np.sum(q < threshold))

    def get_element_quality_map(self, metric='minSICN'):
        """Return dict mapping element tag to quality value for all 3D elements.

        Parameters
        ----------
        metric : str
            Quality metric.

        Returns
        -------
        dict
            {element_tag: quality_value} for all 3D elements.
        """
        import gmsh
        self._check_initialized()

        if not self._meshed:
            raise RuntimeError("Mesh not generated yet")

        quality_map = {}
        elem_types, _, _ = gmsh.model.mesh.getElements(3)
        for et in elem_types:
            tags, _ = gmsh.model.mesh.getElementsByType(et)
            if len(tags) > 0:
                qualities = gmsh.model.mesh.getElementQualities(tags, metric)
                for i, tag in enumerate(tags):
                    quality_map[int(tag)] = float(qualities[i])

        return quality_map

    def get_volume_quality(self, vol_id, metric='minSICN'):
        """Get quality statistics for a specific volume.

        Parameters
        ----------
        vol_id : int
            Volume ID (user-facing ID from create_box, etc.).
        metric : str
            Quality metric.

        Returns
        -------
        dict
            min, max, mean, std, n_elements, n_negative.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        all_qualities = []
        for vt in self._volumes[vol_id]:
            elem_types, elem_tags, _ = gmsh.model.mesh.getElements(3, vt)
            for i, et in enumerate(elem_types):
                tags = elem_tags[i]
                if len(tags) > 0:
                    qualities = gmsh.model.mesh.getElementQualities(
                        list(tags), metric)
                    all_qualities.extend(qualities)

        q = np.array(all_qualities)
        if len(q) == 0:
            return {'min': 0, 'max': 0, 'mean': 0, 'std': 0,
                    'n_elements': 0, 'n_negative': 0}

        return {
            'min': float(np.min(q)),
            'max': float(np.max(q)),
            'mean': float(np.mean(q)),
            'std': float(np.std(q)),
            'n_elements': len(q),
            'n_negative': int(np.sum(q < 0)),
        }

    def get_quality_summary(self):
        """Return comprehensive dict with multiple quality metrics.

        Computes statistics for minSICN (scaled Jacobian), gamma
        (aspect ratio), minSIGE (shape regularity), and minAniso
        (anisotropy).

        Returns
        -------
        dict
            Keys are metric names, values are stat dicts
            (min, max, mean, std, n_elements, n_negative).
        """
        self._check_initialized()

        metrics = ['minSICN', 'gamma', 'minSIGE', 'minAniso']
        summary = {}
        for metric in metrics:
            try:
                summary[metric] = self.get_quality_stats(metric)
            except Exception:
                summary[metric] = None

        return summary

    def print_quality_report(self, metric='minSICN'):
        """Print a formatted quality report to the console.

        Parameters
        ----------
        metric : str
            Quality metric.
        """
        self._check_initialized()

        stats = self.get_quality_stats(metric)
        hist = self.get_quality_histogram(metric, n_bins=10)

        print(f"=== Mesh Quality Report (metric: {metric}) ===")
        print(f"  Elements:   {stats['n_elements']}")
        print(f"  Min:        {stats['min']:.6f}")
        print(f"  Max:        {stats['max']:.6f}")
        print(f"  Mean:       {stats['mean']:.6f}")
        print(f"  Std Dev:    {stats['std']:.6f}")
        print(f"  Negative:   {stats['n_negative']}")
        print()
        print("  Histogram:")
        edges = hist['bin_edges']
        counts = hist['counts']
        pcts = hist['percentages']
        for i in range(len(counts)):
            bar = '#' * max(1, int(pcts[i] / 2))
            print(f"    [{edges[i]:+.3f}, {edges[i+1]:+.3f}]: "
                  f"{counts[i]:5d} ({pcts[i]:5.1f}%) {bar}")
        print("=" * 48)

    def get_distortion(self):
        """Get mesh distortion (element shape deviation from ideal).

        Uses the 'minSIGE' metric which measures the ratio of
        inscribed to circumscribed sphere radii, normalized to [0,1].
        A value of 1 means ideal element shape.

        Returns
        -------
        dict
            min, max, mean distortion values and per-element array.
        """
        self._check_initialized()

        q = self._get_quality_data('minSIGE')
        if len(q) == 0:
            return {'min': 0.0, 'max': 0.0, 'mean': 0.0,
                    'values': np.array([])}

        return {
            'min': float(np.min(q)),
            'max': float(np.max(q)),
            'mean': float(np.mean(q)),
            'values': q,
        }

    def get_edge_length_range(self):
        """Return (min_edge_length, max_edge_length) across all 3D elements.

        Iterates through all 3D elements, computes edge lengths from
        node coordinates, and returns the global minimum and maximum.

        Returns
        -------
        tuple
            (min_edge_length, max_edge_length).
        """
        import gmsh
        self._check_initialized()

        if not self._meshed:
            raise RuntimeError("Mesh not generated yet")

        # Edge definitions per element type (pairs of local node indices)
        # GMSH uses 0-based local indexing within element connectivity
        edge_defs = {
            'Tetrahedron': [
                (0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],
            'Hexahedron': [
                (0, 1), (1, 2), (2, 3), (3, 0),
                (4, 5), (5, 6), (6, 7), (7, 4),
                (0, 4), (1, 5), (2, 6), (3, 7)],
            'Prism': [
                (0, 1), (1, 2), (2, 0),
                (3, 4), (4, 5), (5, 3),
                (0, 3), (1, 4), (2, 5)],
            'Pyramid': [
                (0, 1), (1, 2), (2, 3), (3, 0),
                (0, 4), (1, 4), (2, 4), (3, 4)],
        }

        # Build node coordinate lookup
        node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
        coords_array = np.array(node_coords).reshape(-1, 3)
        tag_to_idx = {int(t): i for i, t in enumerate(node_tags)}

        min_len = float('inf')
        max_len = 0.0

        elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(3)
        for i, et in enumerate(elem_types):
            name, _, _, n_nodes, _, _ = gmsh.model.mesh.getElementProperties(
                et)
            nodes = elem_node_tags[i]
            n_elems = len(elem_tags[i])

            # Find matching edge definition
            edges = None
            for key, val in edge_defs.items():
                if key in name:
                    edges = val
                    break
            if edges is None:
                continue

            for j in range(n_elems):
                start = j * n_nodes
                elem_nodes = nodes[start:start + n_nodes]
                for n0, n1 in edges:
                    idx0 = tag_to_idx[int(elem_nodes[n0])]
                    idx1 = tag_to_idx[int(elem_nodes[n1])]
                    diff = coords_array[idx0] - coords_array[idx1]
                    length = np.sqrt(np.dot(diff, diff))
                    if length < min_len:
                        min_len = length
                    if length > max_len:
                        max_len = length

        if min_len == float('inf'):
            return (0.0, 0.0)
        return (float(min_len), float(max_len))

    def get_volume_ratio_range(self):
        """Return (min_vol, max_vol, ratio) for 3D element volumes.

        Computes the volume of each 3D element from node coordinates
        and returns the minimum, maximum, and their ratio (max/min).

        Returns
        -------
        tuple
            (min_vol, max_vol, ratio). Ratio is max_vol/min_vol.
        """
        import gmsh
        self._check_initialized()

        if not self._meshed:
            raise RuntimeError("Mesh not generated yet")

        # Build node coordinate lookup
        node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
        coords_array = np.array(node_coords).reshape(-1, 3)
        tag_to_idx = {int(t): i for i, t in enumerate(node_tags)}

        volumes = []
        elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(3)
        for i, et in enumerate(elem_types):
            name, _, _, n_nodes, _, _ = gmsh.model.mesh.getElementProperties(
                et)
            nodes = elem_node_tags[i]
            n_elems = len(elem_tags[i])

            for j in range(n_elems):
                start = j * n_nodes
                elem_nodes = nodes[start:start + n_nodes]
                # Get vertex coordinates
                verts = np.array([
                    coords_array[tag_to_idx[int(elem_nodes[k])]]
                    for k in range(n_nodes)])

                if 'Tetrahedron' in name and n_nodes >= 4:
                    # V = |det([v1-v0, v2-v0, v3-v0])| / 6
                    mat = np.column_stack([
                        verts[1] - verts[0],
                        verts[2] - verts[0],
                        verts[3] - verts[0]])
                    vol = abs(np.linalg.det(mat)) / 6.0
                    volumes.append(vol)
                elif 'Hexahedron' in name and n_nodes >= 8:
                    # Decompose hex into 5 tets, sum volumes
                    tet_indices = [
                        (0, 1, 3, 4), (1, 2, 3, 6),
                        (4, 5, 1, 6), (4, 7, 3, 6),
                        (1, 3, 4, 6)]
                    vol = 0.0
                    for ti in tet_indices:
                        mat = np.column_stack([
                            verts[ti[1]] - verts[ti[0]],
                            verts[ti[2]] - verts[ti[0]],
                            verts[ti[3]] - verts[ti[0]]])
                        vol += abs(np.linalg.det(mat)) / 6.0
                    volumes.append(vol)
                elif ('Prism' in name or 'Wedge' in name) and n_nodes >= 6:
                    # Decompose wedge into 3 tets
                    tet_indices = [
                        (0, 1, 2, 3), (1, 2, 3, 4),
                        (2, 3, 4, 5)]
                    vol = 0.0
                    for ti in tet_indices:
                        mat = np.column_stack([
                            verts[ti[1]] - verts[ti[0]],
                            verts[ti[2]] - verts[ti[0]],
                            verts[ti[3]] - verts[ti[0]]])
                        vol += abs(np.linalg.det(mat)) / 6.0
                    volumes.append(vol)
                elif 'Pyramid' in name and n_nodes >= 5:
                    # Decompose pyramid into 2 tets
                    tet_indices = [(0, 1, 3, 4), (1, 2, 3, 4)]
                    vol = 0.0
                    for ti in tet_indices:
                        mat = np.column_stack([
                            verts[ti[1]] - verts[ti[0]],
                            verts[ti[2]] - verts[ti[0]],
                            verts[ti[3]] - verts[ti[0]]])
                        vol += abs(np.linalg.det(mat)) / 6.0
                    volumes.append(vol)

        if len(volumes) == 0:
            return (0.0, 0.0, 0.0)

        min_vol = float(min(volumes))
        max_vol = float(max(volumes))
        ratio = max_vol / min_vol if min_vol > 0 else float('inf')
        return (min_vol, max_vol, ratio)

    def smooth_mesh(self, n_iterations=3, method=''):
        """Smooth the mesh by running optimization multiple times.

        Alias for repeated optimize_mesh calls. Equivalent to
        Gmsh 'smooth vol all'.

        Parameters
        ----------
        n_iterations : int
            Number of smoothing iterations.
        method : str
            GMSH optimization method ('' for default Laplacian).
        """
        import gmsh
        self._check_initialized()

        for i in range(n_iterations):
            gmsh.model.mesh.optimize(method)

        if self._verbose:
            stats = self.get_quality_stats()
            print(f"[GmshBuilder] smooth_mesh ({n_iterations} iterations): "
                  f"min={stats['min']:.4f} mean={stats['mean']:.4f}")

    def untangle_mesh(self):
        """Fix tangled (inverted) elements.

        Uses GMSH 'UntangleMeshGeometry' optimizer to resolve
        inverted elements. Equivalent to Gmsh 'smooth vol untangle'.

        Returns
        -------
        int
            Number of remaining inverted elements after untangling.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.optimize("UntangleMeshGeometry")
        check = self.check_inverted_elements()

        if self._verbose:
            print(f"[GmshBuilder] untangle_mesh: "
                  f"{check['n_inverted']} inverted remaining")
        return check['n_inverted']

    def optimize_high_order(self):
        """Optimize high-order element quality.

        Uses GMSH 'HighOrderOptimize' to improve the placement
        of high-order nodes. Should be called after generating
        a high-order mesh with set_order().
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.optimize("HighOrder")

        if self._verbose:
            stats = self.get_quality_stats()
            print(f"[GmshBuilder] optimize_high_order: "
                  f"min={stats['min']:.4f} mean={stats['mean']:.4f}")

    def get_face_quality(self, metric='minSICN'):
        """Get quality of 2D (surface) elements.

        Parameters
        ----------
        metric : str
            Quality metric: 'minSICN', 'minSIGE', 'gamma', 'minAniso'.

        Returns
        -------
        numpy array
            Quality values for all 2D elements.
        """
        import gmsh
        self._check_initialized()

        if not self._meshed:
            raise RuntimeError("Mesh not generated yet")

        elem_types, _, _ = gmsh.model.mesh.getElements(2)
        all_qualities = []
        for et in elem_types:
            tags, _ = gmsh.model.mesh.getElementsByType(et)
            if len(tags) > 0:
                qualities = gmsh.model.mesh.getElementQualities(
                    tags, metric)
                all_qualities.extend(qualities)

        return np.array(all_qualities)
