"""
Mesh generation methods for GmshBuilder.

Gmsh equivalents: mesh vol, mesh surface, delete mesh, refine, smooth.
"""

import math

from .utils import apply_transfinite, compute_divisions_from_size


class MeshGenerateMixin:
    """Mixin providing mesh generation and optimization."""

    def generate(self, element_type='hex'):
        """Generate 3D mesh (Gmsh 'mesh vol all').

        Parameters
        ----------
        element_type : str
            'hex' for hexahedral (transfinite), 'tet' for tetrahedral.
        """
        import gmsh
        self._check_initialized()

        if element_type == 'hex':
            all_tags = []
            for tags in self._volumes.values():
                all_tags.extend(tags)

            for vol_tag in all_tags:
                if vol_tag in self._divisions:
                    nx, ny, nz = self._divisions[vol_tag]
                elif self._global_divisions is not None:
                    nx, ny, nz = self._global_divisions
                elif self._mesh_size is not None:
                    nx, ny, nz = compute_divisions_from_size(
                        vol_tag, self._mesh_size)
                else:
                    nx, ny, nz = compute_divisions_from_size(vol_tag, 0.01)

                apply_transfinite(vol_tag, nx, ny, nz)

            gmsh.option.setNumber('Mesh.Algorithm3D', 1)

        gmsh.model.mesh.generate(3)
        self._meshed = True

        if self._verbose:
            elem_types, _, _ = gmsh.model.mesh.getElements(3)
            n_hex = n_tet = 0
            for et in elem_types:
                name = gmsh.model.mesh.getElementProperties(et)[0]
                tags, _ = gmsh.model.mesh.getElementsByType(et)
                if 'Hexahedron' in name:
                    n_hex += len(tags)
                elif 'Tetrahedron' in name:
                    n_tet += len(tags)
            print(f"[GmshBuilder] generate({element_type}): "
                  f"{n_hex} hex, {n_tet} tet")

    def generate_surface(self, mesh_size=None):
        """Generate 2D surface mesh (Gmsh 'mesh surface', PEEC BEM).

        Parameters
        ----------
        mesh_size : float, optional
            Target element size.
        """
        import gmsh
        self._check_initialized()

        size = mesh_size or self._mesh_size or 0.01
        gmsh.option.setNumber('Mesh.CharacteristicLengthMax', size)
        gmsh.option.setNumber('Mesh.CharacteristicLengthMin', size * 0.5)
        gmsh.model.mesh.generate(2)
        self._meshed = True

        if self._verbose:
            n_tri = n_quad = 0
            elem_types, _, _ = gmsh.model.mesh.getElements(2)
            for et in elem_types:
                name = gmsh.model.mesh.getElementProperties(et)[0]
                tags, _ = gmsh.model.mesh.getElementsByType(et)
                if 'Triangle' in name:
                    n_tri += len(tags)
                elif 'Quadrilateral' in name:
                    n_quad += len(tags)
            print(f"[GmshBuilder] generate_surface: {n_tri} tri, {n_quad} quad")

    def generate_curve(self):
        """Generate 1D curve mesh (Gmsh 'mesh curve').

        Meshes all curves in the model.  GMSH's ``generate(1)`` always
        meshes all curves; selective curve meshing is not supported.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.generate(1)
        if self._verbose:
            print("[GmshBuilder] generate_curve")

    def generate_volume(self, vol_ids, element_type='hex'):
        """Generate mesh for specific volumes only (Gmsh 'mesh vol').

        Parameters
        ----------
        vol_ids : list of int
            Volume IDs to mesh.
        element_type : str
            'hex' or 'tet'.
        """
        import gmsh
        self._check_initialized()

        if element_type == 'hex':
            for vid in vol_ids:
                if vid not in self._volumes:
                    raise KeyError(f"Volume {vid} not found")
                for vol_tag in self._volumes[vid]:
                    if vol_tag in self._divisions:
                        nx, ny, nz = self._divisions[vol_tag]
                    elif self._global_divisions is not None:
                        nx, ny, nz = self._global_divisions
                    elif self._mesh_size is not None:
                        nx, ny, nz = compute_divisions_from_size(
                            vol_tag, self._mesh_size)
                    else:
                        nx, ny, nz = compute_divisions_from_size(
                            vol_tag, 0.01)
                    apply_transfinite(vol_tag, nx, ny, nz)

            gmsh.option.setNumber('Mesh.Algorithm3D', 1)

        gmsh.model.mesh.generate(3)
        self._meshed = True

        if self._verbose:
            print(f"[GmshBuilder] generate_volume {vol_ids}")

    def delete_mesh(self, vol_id=None):
        """Delete mesh (Gmsh 'delete mesh vol').

        Parameters
        ----------
        vol_id : int, optional
            Volume to clear mesh. If None, clears all.
        """
        import gmsh
        self._check_initialized()

        if vol_id is not None:
            if vol_id not in self._volumes:
                raise KeyError(f"Volume {vol_id} not found")
            for vt in self._volumes[vol_id]:
                gmsh.model.mesh.clear([(3, vt)])
        else:
            gmsh.model.mesh.clear()
            self._meshed = False

        if self._verbose:
            print(f"[GmshBuilder] delete_mesh vol={vol_id}")

    def delete_all_mesh(self):
        """Delete all mesh (Gmsh 'delete mesh')."""
        self.delete_mesh(None)

    def refine(self, n_refine=1):
        """Refine mesh uniformly (Gmsh 'refine vol').

        Parameters
        ----------
        n_refine : int
            Number of refinement levels.
        """
        import gmsh
        self._check_initialized()

        for _ in range(n_refine):
            gmsh.model.mesh.refine()

        if self._verbose:
            print(f"[GmshBuilder] refine n={n_refine}")

    def optimize_mesh(self, method=""):
        """Optimize mesh (Gmsh 'smooth vol').

        Parameters
        ----------
        method : str
            Optimization method: '' (default), 'Netgen', 'HighOrder',
            'HighOrderElastic', 'HighOrderFastCurving',
            'Laplace2D', 'Relocate3D', 'UntangleMeshGeometry'.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.optimize(method)
        if self._verbose:
            print(f"[GmshBuilder] optimize_mesh method='{method}'")

    def smooth_laplacian(self, n_iterations=3):
        """Laplacian smoothing (Gmsh 'smooth').

        Parameters
        ----------
        n_iterations : int
            Number of smoothing passes.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.optimize("Laplace2D", niter=n_iterations)

        if self._verbose:
            print(f"[GmshBuilder] smooth_laplacian n={n_iterations}")

    def convert_to_second_order(self):
        """Convert mesh to 2nd order elements (Gmsh 'modify mesh order 2')."""
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.setOrder(2)
        if self._verbose:
            print("[GmshBuilder] convert_to_second_order")

    def recombine_mesh(self):
        """Recombine triangles to quadrilaterals (Gmsh 'recombine')."""
        import gmsh
        self._check_initialized()

        gmsh.option.setNumber('Mesh.RecombineAll', 1)
        gmsh.model.mesh.recombine()
        if self._verbose:
            print("[GmshBuilder] recombine_mesh")

    def partition_mesh(self, n_parts):
        """Partition mesh for parallel processing.

        Parameters
        ----------
        n_parts : int
            Number of partitions.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.partition(n_parts)
        if self._verbose:
            print(f"[GmshBuilder] partition_mesh n={n_parts}")

    def generate_tet(self, mesh_size=None):
        """Generate tetrahedral mesh.

        Convenience method that sets the 3D algorithm to Delaunay and
        generates a pure tetrahedral mesh.

        Parameters
        ----------
        mesh_size : float, optional
            Target element size. If None, uses the instance mesh_size
            or a default of 0.01.
        """
        import gmsh
        self._check_initialized()

        size = mesh_size or self._mesh_size or 0.01
        gmsh.option.setNumber('Mesh.CharacteristicLengthMax', size)
        gmsh.option.setNumber('Mesh.CharacteristicLengthMin', size * 0.5)
        gmsh.option.setNumber('Mesh.Algorithm3D', 1)  # Delaunay
        gmsh.option.setNumber('Mesh.RecombineAll', 0)
        gmsh.model.mesh.generate(3)
        self._meshed = True

        if self._verbose:
            n_tet = 0
            elem_types, _, _ = gmsh.model.mesh.getElements(3)
            for et in elem_types:
                name = gmsh.model.mesh.getElementProperties(et)[0]
                tags, _ = gmsh.model.mesh.getElementsByType(et)
                if 'Tetrahedron' in name:
                    n_tet += len(tags)
            print(f"[GmshBuilder] generate_tet: {n_tet} tet")

    def generate_hex(self, mesh_size=None):
        """Generate hexahedral mesh.

        Convenience method that applies transfinite constraints to all
        volumes and generates a structured hexahedral mesh.

        Parameters
        ----------
        mesh_size : float, optional
            Target element size used to compute divisions. If None, uses
            the instance mesh_size or a default of 0.01.
        """
        import gmsh
        self._check_initialized()

        all_tags = []
        for tags in self._volumes.values():
            all_tags.extend(tags)

        for vol_tag in all_tags:
            if vol_tag in self._divisions:
                nx, ny, nz = self._divisions[vol_tag]
            elif self._global_divisions is not None:
                nx, ny, nz = self._global_divisions
            else:
                size = mesh_size or self._mesh_size or 0.01
                nx, ny, nz = compute_divisions_from_size(vol_tag, size)

            apply_transfinite(vol_tag, nx, ny, nz)

        gmsh.option.setNumber('Mesh.Algorithm3D', 1)
        gmsh.model.mesh.generate(3)
        self._meshed = True

        if self._verbose:
            n_hex = 0
            elem_types, _, _ = gmsh.model.mesh.getElements(3)
            for et in elem_types:
                name = gmsh.model.mesh.getElementProperties(et)[0]
                tags, _ = gmsh.model.mesh.getElementsByType(et)
                if 'Hexahedron' in name:
                    n_hex += len(tags)
            print(f"[GmshBuilder] generate_hex: {n_hex} hex")

    def generate_mixed(self):
        """Generate mesh allowing mixed element types.

        Produces a mesh that may contain hexahedra, tetrahedra, and wedge
        (prism) elements. Disables recombination so the mesher is free to
        choose the best element type per region.
        """
        import gmsh
        self._check_initialized()

        gmsh.option.setNumber('Mesh.RecombineAll', 0)
        gmsh.option.setNumber('Mesh.Algorithm3D', 1)
        gmsh.model.mesh.generate(3)
        self._meshed = True

        if self._verbose:
            n_hex = n_tet = n_wedge = 0
            elem_types, _, _ = gmsh.model.mesh.getElements(3)
            for et in elem_types:
                name = gmsh.model.mesh.getElementProperties(et)[0]
                tags, _ = gmsh.model.mesh.getElementsByType(et)
                if 'Hexahedron' in name:
                    n_hex += len(tags)
                elif 'Tetrahedron' in name:
                    n_tet += len(tags)
                elif 'Prism' in name or 'Wedge' in name:
                    n_wedge += len(tags)
            print(f"[GmshBuilder] generate_mixed: "
                  f"{n_hex} hex, {n_tet} tet, {n_wedge} wedge")

    def generate_surface_quad(self, mesh_size=None):
        """Generate quadrilateral surface mesh.

        Enables recombination on all surfaces and generates a 2D mesh
        composed primarily of quadrilateral elements.

        Parameters
        ----------
        mesh_size : float, optional
            Target element size. If None, uses the instance mesh_size
            or a default of 0.01.
        """
        import gmsh
        self._check_initialized()

        size = mesh_size or self._mesh_size or 0.01
        gmsh.option.setNumber('Mesh.CharacteristicLengthMax', size)
        gmsh.option.setNumber('Mesh.CharacteristicLengthMin', size * 0.5)
        gmsh.option.setNumber('Mesh.RecombineAll', 1)
        gmsh.model.mesh.generate(2)
        self._meshed = True

        if self._verbose:
            n_tri = n_quad = 0
            elem_types, _, _ = gmsh.model.mesh.getElements(2)
            for et in elem_types:
                name = gmsh.model.mesh.getElementProperties(et)[0]
                tags, _ = gmsh.model.mesh.getElementsByType(et)
                if 'Triangle' in name:
                    n_tri += len(tags)
                elif 'Quadrilateral' in name:
                    n_quad += len(tags)
            print(f"[GmshBuilder] generate_surface_quad: "
                  f"{n_quad} quad, {n_tri} tri")

    def generate_and_optimize(self, element_type='hex', method=''):
        """Generate mesh then optimize in one call.

        Parameters
        ----------
        element_type : str
            'hex' for hexahedral, 'tet' for tetrahedral.
        method : str
            Optimization method passed to gmsh.model.mesh.optimize().
            '' (default), 'Netgen', 'HighOrder', 'Laplace2D',
            'Relocate3D', 'UntangleMeshGeometry'.
        """
        import gmsh
        self._check_initialized()

        self.generate(element_type=element_type)
        gmsh.model.mesh.optimize(method)

        if self._verbose:
            print(f"[GmshBuilder] generate_and_optimize "
                  f"type='{element_type}' method='{method}'")

    def convert_to_first_order(self):
        """Convert mesh back to 1st order elements.

        Reverses a previous convert_to_second_order() call, reducing
        all elements to their linear (1st order) equivalents.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.setOrder(1)
        if self._verbose:
            print("[GmshBuilder] convert_to_first_order")

    def renumber_nodes(self):
        """Renumber mesh nodes to be contiguous starting from 1.

        Useful after mesh modification operations that may leave gaps
        in node numbering.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.renumberNodes()
        if self._verbose:
            print("[GmshBuilder] renumber_nodes")

    def renumber_elements(self):
        """Renumber mesh elements to be contiguous starting from 1.

        Useful after mesh modification operations that may leave gaps
        in element numbering.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.renumberElements()
        if self._verbose:
            print("[GmshBuilder] renumber_elements")

    def classify_surfaces(self):
        """Classify mesh surfaces.

        Reclassifies surface mesh elements based on the underlying
        geometry. This is useful after importing an unstructured mesh
        to associate mesh faces with geometric surfaces.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.classifySurfaces(
            math.pi,  # angle threshold (radians)
            True,      # boundary edges
            False,     # curved surfaces
            math.pi    # curveAngle threshold
        )
        if self._verbose:
            print("[GmshBuilder] classify_surfaces")

    def create_geometry_from_mesh(self):
        """Create geometry from the current mesh.

        Constructs geometric entities (points, curves, surfaces, volumes)
        from the existing mesh. Typically called after classify_surfaces().
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.createGeometry()
        if self._verbose:
            print("[GmshBuilder] create_geometry_from_mesh")

    def embed_point_in_volume(self, point_tag, vol_tag):
        """Embed a point in a volume for local mesh refinement.

        The embedded point forces the mesh generator to place a node
        at the point location, enabling local size control.

        Parameters
        ----------
        point_tag : int
            Tag of the point entity to embed.
        vol_tag : int
            Tag of the volume entity to embed the point in.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.embed(0, [point_tag], 3, vol_tag)
        if self._verbose:
            print(f"[GmshBuilder] embed_point_in_volume "
                  f"point={point_tag} vol={vol_tag}")

    def embed_curve_in_surface(self, curve_tag, surf_tag):
        """Embed a curve in a surface.

        The embedded curve forces the mesh generator to conform to
        the curve, creating edges along it in the surface mesh.

        Parameters
        ----------
        curve_tag : int
            Tag of the curve entity to embed.
        surf_tag : int
            Tag of the surface entity to embed the curve in.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.embed(1, [curve_tag], 2, surf_tag)
        if self._verbose:
            print(f"[GmshBuilder] embed_curve_in_surface "
                  f"curve={curve_tag} surf={surf_tag}")

    def embed_surface_in_volume(self, surf_tag, vol_tag):
        """Embed a surface in a volume.

        The embedded surface forces the mesh generator to conform to
        the surface, creating faces along it in the volume mesh.

        Parameters
        ----------
        surf_tag : int
            Tag of the surface entity to embed.
        vol_tag : int
            Tag of the volume entity to embed the surface in.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.embed(2, [surf_tag], 3, vol_tag)
        if self._verbose:
            print(f"[GmshBuilder] embed_surface_in_volume "
                  f"surf={surf_tag} vol={vol_tag}")

    def set_mesh_coherence(self, value=True):
        """Enable/disable coherent meshing (merge nodes at shared boundaries).

        When enabled (default), GMSH merges duplicate nodes at shared
        boundaries. This is the standard behavior for multi-body assemblies.

        Parameters
        ----------
        value : bool
            True to enable coherence, False to disable.
        """
        import gmsh
        self._check_initialized()

        # GMSH handles coherence via fragment() boolean operation
        # Store the preference for use during mesh generation
        if not hasattr(self, '_mesh_coherence'):
            self._mesh_coherence = True
        self._mesh_coherence = bool(value)
        if self._verbose:
            print(f"[GmshBuilder] set_mesh_coherence: {value}")

    def get_mesh_statistics(self):
        """Return dictionary with mesh statistics.

        Returns
        -------
        dict
            Dictionary containing:
            - 'n_nodes': total number of mesh nodes
            - 'n_elements': total number of mesh elements
            - 'element_counts': dict mapping element type name to count
            - 'quality_min': minimum element quality (0-1)
            - 'quality_max': maximum element quality (0-1)
        """
        import gmsh
        self._check_initialized()

        # Node count
        node_tags, _, _ = gmsh.model.mesh.getNodes()
        n_nodes = len(node_tags)

        # Element counts by type (all dimensions)
        element_counts = {}
        n_elements_total = 0
        for dim in range(4):
            elem_types, elem_tags, _ = gmsh.model.mesh.getElements(dim)
            for i, et in enumerate(elem_types):
                name = gmsh.model.mesh.getElementProperties(et)[0]
                count = len(elem_tags[i])
                element_counts[name] = element_counts.get(name, 0) + count
                n_elements_total += count

        # Quality range (for highest-dimension elements)
        quality_min = float('inf')
        quality_max = float('-inf')
        for dim in [3, 2, 1]:
            elem_types_dim, _, _ = gmsh.model.mesh.getElements(dim)
            if elem_types_dim:
                for et in elem_types_dim:
                    qualities = gmsh.model.mesh.getElementQualities(
                        gmsh.model.mesh.getElementsByType(et)[0].tolist(),
                        qualityName='minSJ'
                    )
                    if len(qualities) > 0:
                        quality_min = min(quality_min, min(qualities))
                        quality_max = max(quality_max, max(qualities))
                break  # Only check highest dimension present

        # Fallback when no elements were found for quality computation
        if quality_min == float('inf'):
            quality_min = 0.0
        if quality_max == float('-inf'):
            quality_max = 0.0

        stats = {
            'n_nodes': n_nodes,
            'n_elements': n_elements_total,
            'element_counts': element_counts,
            'quality_min': quality_min,
            'quality_max': quality_max,
        }

        if self._verbose:
            print(f"[GmshBuilder] get_mesh_statistics: "
                  f"{n_nodes} nodes, {n_elements_total} elements")
            for name, count in element_counts.items():
                print(f"  {name}: {count}")

        return stats

    # ------------------------------------------------------------------
    # Advanced generation methods
    # ------------------------------------------------------------------

    def generate_adaptive(self, field_tag=None, n_iter=3,
                          quality_threshold=0.3):
        """Iterative mesh refinement with quality check.

        Generates the mesh and evaluates element quality.  If any element
        quality falls below *quality_threshold*, a **uniform** refinement
        step (``gmsh.model.mesh.refine()``) is applied and the quality is
        re-evaluated.  This cycle repeats up to *n_iter* iterations or
        until all elements satisfy the quality threshold.

        .. note::

           This method uses uniform refinement, not truly adaptive
           (local-only) refinement.  Every ``refine()`` call subdivides
           **all** elements, which increases the total element count
           rapidly.  Use a small *n_iter* (1--3) to avoid excessive mesh
           sizes.

        Parameters
        ----------
        field_tag : int or None, optional
            GMSH background field tag that drives element sizing.  If
            ``None``, the current mesh size settings are used.
        n_iter : int, optional
            Maximum number of generate-refine iterations (default 3).
        quality_threshold : float, optional
            Minimum acceptable element quality in the range ``[0, 1]``.
            Elements below this value trigger refinement (default 0.3).
        """
        import gmsh
        self._check_initialized()

        if field_tag is not None:
            gmsh.model.mesh.field.setAsBackgroundMesh(field_tag)

        for iteration in range(n_iter):
            gmsh.model.mesh.generate(3)

            # Evaluate quality of 3-D elements
            elem_types, _, _ = gmsh.model.mesh.getElements(3)
            needs_refinement = False
            for et in elem_types:
                tags_arr, _ = gmsh.model.mesh.getElementsByType(et)
                if len(tags_arr) == 0:
                    continue
                qualities = gmsh.model.mesh.getElementQualities(
                    tags_arr.tolist(), qualityName='minSJ')
                if any(q < quality_threshold for q in qualities):
                    needs_refinement = True
                    break

            if self._verbose:
                print(f"[GmshBuilder] generate_adaptive iter={iteration + 1}/"
                      f"{n_iter} refine={needs_refinement}")

            if not needs_refinement:
                break

            gmsh.model.mesh.refine()

        self._meshed = True

    def generate_boundary_layer_mesh(self, surf_tags=None, thickness=0.001,
                                     n_layers=5, ratio=1.2):
        """Generate mesh with boundary layer refinement near surfaces.

        Sets up a GMSH ``BoundaryLayer`` mesh size field on the specified
        surfaces, then generates the full 3-D mesh.

        Parameters
        ----------
        surf_tags : list of int or None, optional
            Surface tags on which to create boundary layers.  If ``None``,
            all model surfaces are used.
        thickness : float, optional
            Total boundary layer thickness in model units (default 0.001).
        n_layers : int, optional
            Number of layers in the boundary layer (default 5).
        ratio : float, optional
            Growth ratio between successive layers (default 1.2).
        """
        import gmsh
        self._check_initialized()

        if surf_tags is None:
            surfs = gmsh.model.getEntities(2)
            surf_tags = [s[1] for s in surfs]

        # Compute first layer height from geometric series
        if ratio != 1.0:
            first_height = thickness * (1.0 - ratio) / (1.0 - ratio ** n_layers)
        else:
            first_height = thickness / n_layers

        # Extract boundary curves from the surfaces -- the BoundaryLayer
        # field requires CurvesList, not SurfacesList.
        curve_tags = []
        for st in surf_tags:
            boundary = gmsh.model.getBoundary([(2, st)], oriented=False)
            curve_tags.extend([abs(dt[1]) for dt in boundary if dt[0] == 1])
        curve_tags = list(set(curve_tags))

        field_id = gmsh.model.mesh.field.add("BoundaryLayer")
        gmsh.model.mesh.field.setNumbers(field_id, "CurvesList", curve_tags)
        gmsh.model.mesh.field.setNumber(field_id, "Thickness", thickness)
        gmsh.model.mesh.field.setNumber(field_id, "NbLayers", n_layers)
        gmsh.model.mesh.field.setNumber(field_id, "Ratio", ratio)
        gmsh.model.mesh.field.setNumber(field_id, "Size", first_height)
        gmsh.model.mesh.field.setAsBoundaryLayer(field_id)

        gmsh.model.mesh.generate(3)
        self._meshed = True

        if self._verbose:
            print(f"[GmshBuilder] generate_boundary_layer_mesh "
                  f"surfs={surf_tags} thickness={thickness} "
                  f"n_layers={n_layers} ratio={ratio}")

    def generate_prism_layers(self, surf_tags, n_layers, total_thickness,
                              ratio=1.2):
        """Generate prismatic (wedge) boundary layers on specified surfaces.

        Combines a ``BoundaryLayer`` field with hexahedral sweep settings to
        produce structured prism layers adjacent to the surfaces.

        Parameters
        ----------
        surf_tags : list of int
            Surface tags on which to grow prism layers.
        n_layers : int
            Number of prism layers.
        total_thickness : float
            Total thickness of the prism layer region.
        ratio : float, optional
            Growth ratio between successive layers (default 1.2).
        """
        import gmsh
        self._check_initialized()

        if ratio != 1.0:
            first_height = (total_thickness * (1.0 - ratio)
                            / (1.0 - ratio ** n_layers))
        else:
            first_height = total_thickness / n_layers

        # Extract boundary curves from the surfaces -- the BoundaryLayer
        # field requires CurvesList, not SurfacesList.
        curve_tags = []
        for st in surf_tags:
            boundary = gmsh.model.getBoundary([(2, st)], oriented=False)
            curve_tags.extend([abs(dt[1]) for dt in boundary if dt[0] == 1])
        curve_tags = list(set(curve_tags))

        field_id = gmsh.model.mesh.field.add("BoundaryLayer")
        gmsh.model.mesh.field.setNumbers(field_id, "CurvesList", curve_tags)
        gmsh.model.mesh.field.setNumber(field_id, "Thickness", total_thickness)
        gmsh.model.mesh.field.setNumber(field_id, "NbLayers", n_layers)
        gmsh.model.mesh.field.setNumber(field_id, "Ratio", ratio)
        gmsh.model.mesh.field.setNumber(field_id, "Size", first_height)
        gmsh.model.mesh.field.setNumber(field_id, "Quads", 1)
        gmsh.model.mesh.field.setAsBoundaryLayer(field_id)

        gmsh.option.setNumber('Mesh.Algorithm3D', 1)
        gmsh.model.mesh.generate(3)
        self._meshed = True

        if self._verbose:
            print(f"[GmshBuilder] generate_prism_layers "
                  f"surfs={surf_tags} n_layers={n_layers} "
                  f"thickness={total_thickness} ratio={ratio}")

    def generate_surface_on_entity(self, surf_tag):
        """Generate 2D mesh on a specific surface entity.

        Parameters
        ----------
        surf_tag : int
            Tag of the surface entity to mesh.
        """
        import gmsh
        self._check_initialized()

        gmsh.option.setNumber('Mesh.MeshOnlyVisible', 1)
        # Hide everything, then show only the target surface
        entities = gmsh.model.getEntities()
        for dim, tag in entities:
            gmsh.model.setVisibility([(dim, tag)], False)
        gmsh.model.setVisibility([(2, surf_tag)], True, recursive=True)

        gmsh.model.mesh.generate(2)

        # Restore visibility
        for dim, tag in entities:
            gmsh.model.setVisibility([(dim, tag)], True)
        gmsh.option.setNumber('Mesh.MeshOnlyVisible', 0)

        self._meshed = True

        if self._verbose:
            print(f"[GmshBuilder] generate_surface_on_entity "
                  f"surf_tag={surf_tag}")

    def generate_curve_on_entity(self, curve_tag):
        """Generate 1D mesh on a specific curve entity.

        Parameters
        ----------
        curve_tag : int
            Tag of the curve entity to mesh.
        """
        import gmsh
        self._check_initialized()

        gmsh.option.setNumber('Mesh.MeshOnlyVisible', 1)
        entities = gmsh.model.getEntities()
        for dim, tag in entities:
            gmsh.model.setVisibility([(dim, tag)], False)
        gmsh.model.setVisibility([(1, curve_tag)], True, recursive=True)

        gmsh.model.mesh.generate(1)

        # Restore visibility
        for dim, tag in entities:
            gmsh.model.setVisibility([(dim, tag)], True)
        gmsh.option.setNumber('Mesh.MeshOnlyVisible', 0)

        if self._verbose:
            print(f"[GmshBuilder] generate_curve_on_entity "
                  f"curve_tag={curve_tag}")

    def split_quadrangles(self):
        """Split all quadrilateral elements into triangles.

        Uses ``gmsh.model.mesh.splitQuadrangles()`` to convert each
        quadrilateral into two triangles.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.splitQuadrangles()

        if self._verbose:
            print("[GmshBuilder] split_quadrangles")

    def reverse_elements(self, dim=-1, tag=-1):
        """Reverse element orientation.

        Attempts to use ``gmsh.model.mesh.reverseElements()`` if available,
        otherwise falls back to ``gmsh.model.mesh.reverse()``.

        Parameters
        ----------
        dim : int, optional
            Dimension of the entities to reverse (-1 for all, default -1).
        tag : int, optional
            Tag of the entity to reverse (-1 for all, default -1).
        """
        import gmsh
        self._check_initialized()

        try:
            # Collect element tags for the specified dim/tag
            all_tags = []
            if dim == -1 and tag == -1:
                for d in range(4):
                    for ent in gmsh.model.getEntities(d):
                        elem_types, elem_tags_list, _ = (
                            gmsh.model.mesh.getElements(ent[0], ent[1]))
                        for tags in elem_tags_list:
                            all_tags.extend(tags)
            else:
                elem_types, elem_tags_list, _ = (
                    gmsh.model.mesh.getElements(dim, tag))
                for tags in elem_tags_list:
                    all_tags.extend(tags)
            if all_tags:
                gmsh.model.mesh.reverseElements(all_tags)
        except AttributeError:
            if dim == -1 and tag == -1:
                # Reverse all entities of all dimensions
                for d in range(4):
                    for ent in gmsh.model.getEntities(d):
                        gmsh.model.mesh.reverse([(ent[0], ent[1])])
            else:
                gmsh.model.mesh.reverse([(dim, tag)])

        if self._verbose:
            print(f"[GmshBuilder] reverse_elements dim={dim} tag={tag}")

    # ------------------------------------------------------------------
    # Reordering / topology methods
    # ------------------------------------------------------------------

    def reorder_elements(self, method='RCMK', dim=-1, tag=-1):
        """Reorder mesh elements for improved memory locality.

        .. note::

            GMSH's ``reorderElements(elementType, tag, ordering)`` requires
            an explicit ordering vector and does not accept a method name.
            Use :meth:`compute_renumbering` to obtain a permutation first,
            then apply it manually.  This convenience wrapper is therefore
            not implemented.

        Parameters
        ----------
        method : str, optional
            Desired reordering algorithm name (default ``'RCMK'``).
        dim : int, optional
            Dimension of entities to reorder (-1 for all, default -1).
        tag : int, optional
            Tag of the entity to reorder (-1 for all, default -1).

        Raises
        ------
        NotImplementedError
            Always raised. GMSH does not provide a simple
            reorder-by-method-name API.
        """
        raise NotImplementedError(
            "GMSH's reorderElements(elementType, tag, ordering) requires an "
            "explicit ordering vector and does not accept a method name. "
            "Use compute_renumbering() to obtain a permutation, then apply "
            "it via gmsh.model.mesh.reorderElements() directly.")

    def compute_renumbering(self, method='RCMK'):
        """Compute element renumbering permutation.

        Uses ``gmsh.model.mesh.computeRenumbering(method, elementTags)``
        if available.

        Parameters
        ----------
        method : str, optional
            Renumbering algorithm (default ``'RCMK'``).

        Returns
        -------
        tuple
            Permutation data (oldTags, newTags) returned by GMSH.
        """
        import gmsh
        self._check_initialized()

        try:
            result = gmsh.model.mesh.computeRenumbering(
                method=method, elementTags=[])
        except AttributeError:
            raise NotImplementedError(
                "gmsh.model.mesh.computeRenumbering() is not available. "
                "Upgrade GMSH to a version that supports this API.")

        if self._verbose:
            print(f"[GmshBuilder] compute_renumbering method='{method}'")

        return result

    def create_topology(self, make_simply_connected=True):
        """Create topology from mesh.

        Builds topological entities (volumes, surfaces, curves, points)
        from the current mesh.

        Parameters
        ----------
        make_simply_connected : bool, optional
            If ``True``, ensure the resulting topology is simply connected
            (default ``True``).
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.createTopology(make_simply_connected)

        if self._verbose:
            print(f"[GmshBuilder] create_topology "
                  f"make_simply_connected={make_simply_connected}")

    def rebuild_node_cache(self):
        """Rebuild the internal node cache.

        Forces GMSH to rebuild its internal node lookup data structures.
        Useful after manual mesh modifications that bypass the cache.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.rebuildNodeCache(True)

        if self._verbose:
            print("[GmshBuilder] rebuild_node_cache")

    def rebuild_element_cache(self):
        """Rebuild the internal element cache.

        Forces GMSH to rebuild its internal element lookup data structures.
        Useful after manual mesh modifications that bypass the cache.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.rebuildElementCache(True)

        if self._verbose:
            print("[GmshBuilder] rebuild_element_cache")

    # ------------------------------------------------------------------
    # Homology methods
    # ------------------------------------------------------------------

    def add_homology_request(self, type_str='Homology', domain_tags=None,
                             subdomain_tags=None, dims=None):
        """Add a homology computation request.

        Registers a request for homology or cohomology computation on the
        specified domain.

        Parameters
        ----------
        type_str : str, optional
            Type of computation: ``'Homology'`` or ``'Cohomology'``
            (default ``'Homology'``).
        domain_tags : list of int, optional
            Volume tags defining the domain.  Empty list means the whole
            model (default ``[]``).
        subdomain_tags : list of int, optional
            Volume tags defining the subdomain (default ``[]``).
        dims : list of int, optional
            Dimensions to compute (default ``[]`` meaning all).
        """
        import gmsh
        self._check_initialized()

        if domain_tags is None:
            domain_tags = []
        if subdomain_tags is None:
            subdomain_tags = []
        if dims is None:
            dims = []

        gmsh.model.mesh.addHomologyRequest(
            type_str, domain_tags, subdomain_tags, dims)

        if self._verbose:
            print(f"[GmshBuilder] add_homology_request type='{type_str}' "
                  f"domain={domain_tags} subdomain={subdomain_tags} "
                  f"dims={dims}")

    def compute_homology(self):
        """Compute homology or cohomology.

        Processes all previously registered homology requests (see
        ``add_homology_request``).
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.computeHomology()

        if self._verbose:
            print("[GmshBuilder] compute_homology")

    def compute_cross_field(self):
        """Compute a cross field on the current mesh.

        Uses ``gmsh.model.mesh.computeCrossField()`` if available.

        Returns
        -------
        list of int
            View tags of the computed cross field (typically 3 tags).
        """
        import gmsh
        self._check_initialized()

        try:
            view_tags = gmsh.model.mesh.computeCrossField()
        except AttributeError:
            raise NotImplementedError(
                "gmsh.model.mesh.computeCrossField() is not available. "
                "Upgrade GMSH to a version that supports this API.")

        if self._verbose:
            print(f"[GmshBuilder] compute_cross_field -> "
                  f"view_tags={view_tags}")

        return view_tags

    def clear_homology_requests(self):
        """Clear all previously registered homology requests.

        Uses ``gmsh.model.mesh.clearHomologyRequests()`` if available.
        """
        import gmsh
        self._check_initialized()

        try:
            gmsh.model.mesh.clearHomologyRequests()
        except AttributeError:
            raise NotImplementedError(
                "gmsh.model.mesh.clearHomologyRequests() is not available. "
                "Upgrade GMSH to a version that supports this API.")

        if self._verbose:
            print("[GmshBuilder] clear_homology_requests")

    # ------------------------------------------------------------------
    # Partition methods
    # ------------------------------------------------------------------

    def unpartition_mesh(self):
        """Remove all mesh partitions.

        Merges all partitions back into a single mesh.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.unpartition()

        if self._verbose:
            print("[GmshBuilder] unpartition_mesh")

    def get_ghost_elements(self, dim, tag):
        """Get ghost elements for a partition entity.

        Ghost elements are copies of elements owned by another partition
        that share nodes with the queried partition.

        Parameters
        ----------
        dim : int
            Dimension of the partition entity.
        tag : int
            Tag of the partition entity.

        Returns
        -------
        list
            Ghost element data returned by GMSH.
        """
        import gmsh
        self._check_initialized()

        try:
            result = gmsh.model.mesh.getGhostElements(dim, tag)
        except AttributeError:
            raise NotImplementedError(
                "gmsh.model.mesh.getGhostElements() is not available. "
                "Upgrade GMSH to a version that supports this API.")

        if self._verbose:
            print(f"[GmshBuilder] get_ghost_elements dim={dim} tag={tag}")

        return list(result)

    def get_partition_count(self):
        """Get the number of mesh partitions.

        Returns
        -------
        int
            Number of partitions in the current model.
        """
        import gmsh
        self._check_initialized()

        n_parts = gmsh.model.getNumberOfPartitions()

        if self._verbose:
            print(f"[GmshBuilder] get_partition_count -> {n_parts}")

        return n_parts

    # ------------------------------------------------------------------
    # STL import methods
    # ------------------------------------------------------------------

    def import_stl_mesh(self, filename, angle=40.0):
        """Import an STL mesh file and create geometry from it.

        Merges the STL file into the current model, classifies the
        surface mesh into geometric surfaces based on dihedral angle,
        creates geometric entities, and builds topology.

        Parameters
        ----------
        filename : str
            Path to the STL file.
        angle : float, optional
            Dihedral angle threshold (in degrees) for surface
            classification (default 40.0).
        """
        import gmsh
        self._check_initialized()

        gmsh.merge(filename)

        angle_rad = angle * math.pi / 180.0
        curve_angle_rad = 180.0 * math.pi / 180.0

        gmsh.model.mesh.classifySurfaces(
            angle_rad,
            True,   # boundary edges
            False,  # curved surfaces
            curve_angle_rad
        )
        gmsh.model.mesh.createGeometry()
        gmsh.model.mesh.createTopology()

        self._meshed = True

        if self._verbose:
            print(f"[GmshBuilder] import_stl_mesh file='{filename}' "
                  f"angle={angle}")

    def reclassify_nodes(self):
        """Reclassify mesh nodes on their associated geometric entities.

        Uses ``gmsh.model.mesh.reclassifyNodes()`` to update node
        classification after geometry or mesh modifications.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.reclassifyNodes()

        if self._verbose:
            print("[GmshBuilder] reclassify_nodes")

    # ------------------------------------------------------------------
    # CFD convenience methods
    # ------------------------------------------------------------------

    def classify_surfaces_parametric(self, angle=40.0, boundary=True,
                                     for_reparametrize=False,
                                     curve_angle=180.0):
        """Classify mesh surfaces with configurable parameters.

        Reclassifies surface mesh elements based on dihedral angle.
        Essential for STL and terrain meshing workflows where the
        angle threshold controls feature detection.

        Parameters
        ----------
        angle : float
            Dihedral angle threshold in degrees for surface
            classification (default 40.0).
        boundary : bool
            If True, create boundary edges between classified
            surfaces (default True).
        for_reparametrize : bool
            If True, prepare surfaces for reparametrization
            (default False).
        curve_angle : float
            Angle threshold in degrees for curve classification
            (default 180.0).
        """
        import gmsh
        import math
        self._check_initialized()

        angle_rad = angle * math.pi / 180.0
        curve_angle_rad = curve_angle * math.pi / 180.0

        gmsh.model.mesh.classifySurfaces(
            angle_rad, boundary, for_reparametrize, curve_angle_rad)

        if self._verbose:
            print(f"[GmshBuilder] classify_surfaces_parametric "
                  f"angle={angle} boundary={boundary} "
                  f"curve_angle={curve_angle}")

    def embed_entities(self, dim, tags, in_dim, in_tag):
        """Embed entities of given dimension into a higher-dimensional entity.

        General-purpose embedding method.  The embedded entities force
        the mesh generator to place nodes/edges/faces at their
        locations.

        Parameters
        ----------
        dim : int
            Dimension of entities to embed (0=point, 1=curve,
            2=surface).
        tags : list of int
            Tags of entities to embed.
        in_dim : int
            Dimension of the host entity (must be > dim).
        in_tag : int
            Tag of the host entity.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.embed(dim, tags, in_dim, in_tag)

        if self._verbose:
            print(f"[GmshBuilder] embed_entities dim={dim} "
                  f"tags={tags} in ({in_dim},{in_tag})")

    def adapt_mesh_iterative(self, field_tag, n_iterations=3,
                             size_factor=0.5):
        """Iterative mesh adaptation driven by a sizing field.

        Generates the mesh, evaluates the sizing field, updates mesh
        sizes, and regenerates.  Unlike ``generate_adaptive`` which
        refines uniformly, this method regenerates with locally
        adapted sizes driven by the specified field.

        Parameters
        ----------
        field_tag : int
            GMSH background field tag that drives element sizing.
        n_iterations : int
            Maximum number of generate/adapt iterations (default 3).
        size_factor : float
            Scaling factor applied to the Mesh.MeshSizeFactor option
            each iteration (default 0.5).  Values < 1 increase
            refinement.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.mesh.field.setAsBackgroundMesh(field_tag)

        for i in range(n_iterations):
            gmsh.model.mesh.generate(3)

            if i < n_iterations - 1:
                current_factor = gmsh.option.getNumber(
                    'Mesh.MeshSizeFactor')
                gmsh.option.setNumber(
                    'Mesh.MeshSizeFactor',
                    current_factor * size_factor)
                gmsh.model.mesh.clear()

        self._meshed = True

        if self._verbose:
            print(f"[GmshBuilder] adapt_mesh_iterative "
                  f"field={field_tag} iterations={n_iterations} "
                  f"factor={size_factor}")
