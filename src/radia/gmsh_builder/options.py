"""
GMSH global options wrapper for GmshBuilder.

Gmsh equivalents: set/get options for geometry, mesh behavior, colors, visibility.
"""


class OptionsMixin:
    """Mixin providing GMSH global option access."""

    def set_option_number(self, name, value):
        """Set a GMSH numeric option.

        Parameters
        ----------
        name : str
            GMSH option string (e.g. 'Mesh.Algorithm', 'General.Terminal').
        value : float
            Numeric value to set.
        """
        import gmsh
        self._check_initialized()

        gmsh.option.setNumber(name, value)

        if self._verbose:
            print(f"[GmshBuilder] set_option_number '{name}' = {value}")

    def get_option_number(self, name):
        """Get a GMSH numeric option.

        Parameters
        ----------
        name : str
            GMSH option string (e.g. 'Mesh.Algorithm').

        Returns
        -------
        float
            Current value of the option.
        """
        import gmsh
        self._check_initialized()

        value = gmsh.option.getNumber(name)

        if self._verbose:
            print(f"[GmshBuilder] get_option_number '{name}' = {value}")

        return value

    def set_option_string(self, name, value):
        """Set a GMSH string option.

        Parameters
        ----------
        name : str
            GMSH option string (e.g. 'Mesh.OptimizeNetgen').
        value : str
            String value to set.
        """
        import gmsh
        self._check_initialized()

        gmsh.option.setString(name, value)

        if self._verbose:
            print(f"[GmshBuilder] set_option_string '{name}' = '{value}'")

    def get_option_string(self, name):
        """Get a GMSH string option.

        Parameters
        ----------
        name : str
            GMSH option string.

        Returns
        -------
        str
            Current value of the option.
        """
        import gmsh
        self._check_initialized()

        value = gmsh.option.getString(name)

        if self._verbose:
            print(f"[GmshBuilder] get_option_string '{name}' = '{value}'")

        return value

    def restore_default_options(self):
        """Restore all GMSH options to their default values."""
        import gmsh
        self._check_initialized()

        gmsh.option.restoreDefaults()

        if self._verbose:
            print("[GmshBuilder] restore_default_options")

    def get_all_mesh_options(self):
        """Read common Mesh.* options and return as a dictionary.

        Returns
        -------
        dict
            Dictionary mapping short option names (without 'Mesh.' prefix)
            to their current numeric values. Includes Algorithm,
            Algorithm3D, CharacteristicLengthMin, CharacteristicLengthMax,
            CharacteristicLengthFactor, ElementOrder, MinimumCurveNodes,
            MeshSizeMin, MeshSizeMax, MeshSizeFromCurvature,
            MeshSizeExtendFromBoundary, Optimize, OptimizeNetgen,
            Smoothing, RecombinationAlgorithm, RecombineAll,
            SubdivisionAlgorithm, and others.
        """
        import gmsh
        self._check_initialized()

        option_keys = [
            'Mesh.Algorithm',
            'Mesh.Algorithm3D',
            'Mesh.CharacteristicLengthMin',
            'Mesh.CharacteristicLengthMax',
            'Mesh.CharacteristicLengthFactor',
            'Mesh.ElementOrder',
            'Mesh.MinimumCurveNodes',
            'Mesh.MeshSizeMin',
            'Mesh.MeshSizeMax',
            'Mesh.MeshSizeFromCurvature',
            'Mesh.MeshSizeExtendFromBoundary',
            'Mesh.Optimize',
            'Mesh.OptimizeNetgen',
            'Mesh.Smoothing',
            'Mesh.RecombinationAlgorithm',
            'Mesh.RecombineAll',
            'Mesh.SubdivisionAlgorithm',
            'Mesh.SaveAll',
            'Mesh.SaveGroupsOfElements',
            'Mesh.SaveGroupsOfNodes',
            'Mesh.Binary',
            'Mesh.RandomFactor',
            'Mesh.MeshOnlyVisible',
        ]

        result = {}
        for key in option_keys:
            short_key = key.replace('Mesh.', '')
            result[short_key] = gmsh.option.getNumber(key)

        if self._verbose:
            print(f"[GmshBuilder] get_all_mesh_options ({len(result)} options)")

        return result

    # ------------------------------------------------------------------
    # Geometry options
    # ------------------------------------------------------------------

    def set_geometry_tolerance(self, tolerance):
        """Set the geometry tolerance for OCC operations.

        Parameters
        ----------
        tolerance : float
            Geometry tolerance value (e.g. 1e-8).
        """
        import gmsh
        self._check_initialized()

        gmsh.option.setNumber('Geometry.Tolerance', tolerance)

        if self._verbose:
            print(f"[GmshBuilder] set_geometry_tolerance = {tolerance}")

    def get_geometry_tolerance(self):
        """Get the current geometry tolerance.

        Returns
        -------
        float
            Current geometry tolerance value.
        """
        import gmsh
        self._check_initialized()

        value = gmsh.option.getNumber('Geometry.Tolerance')

        if self._verbose:
            print(f"[GmshBuilder] get_geometry_tolerance = {value}")

        return value

    def set_occ_scaling(self, factor):
        """Set OCC geometry scaling factor.

        Parameters
        ----------
        factor : float
            Scaling factor applied to OCC geometry import (e.g. 0.001
            to convert mm to meters).
        """
        import gmsh
        self._check_initialized()

        gmsh.option.setNumber('Geometry.OCCScaling', factor)

        if self._verbose:
            print(f"[GmshBuilder] set_occ_scaling = {factor}")

    def set_occ_parallel(self, enabled):
        """Enable or disable parallel OCC geometry operations.

        Parameters
        ----------
        enabled : bool
            If True, OCC operations run in parallel where possible.
        """
        import gmsh
        self._check_initialized()

        gmsh.option.setNumber('Geometry.OCCParallel', 1 if enabled else 0)

        if self._verbose:
            print(f"[GmshBuilder] set_occ_parallel = {enabled}")

    # ------------------------------------------------------------------
    # Mesh behavior options
    # ------------------------------------------------------------------

    def set_mesh_only_visible(self, enabled):
        """Mesh only visible entities.

        Parameters
        ----------
        enabled : bool
            If True, only entities that are visible will be meshed.
        """
        import gmsh
        self._check_initialized()

        gmsh.option.setNumber('Mesh.MeshOnlyVisible', 1 if enabled else 0)

        if self._verbose:
            print(f"[GmshBuilder] set_mesh_only_visible = {enabled}")

    def set_mesh_save_all(self, enabled):
        """Save all elements regardless of physical group membership.

        Parameters
        ----------
        enabled : bool
            If True, all mesh elements are saved on export, not just
            those belonging to physical groups.
        """
        import gmsh
        self._check_initialized()

        gmsh.option.setNumber('Mesh.SaveAll', 1 if enabled else 0)

        if self._verbose:
            print(f"[GmshBuilder] set_mesh_save_all = {enabled}")

    def set_mesh_save_groups_of_elements(self, value):
        """Control which element groups are saved on export.

        Parameters
        ----------
        value : int
            0 = none, 1 = physical groups only, 2 = all elementary
            entity groups, 3 = all (physical + elementary).
        """
        import gmsh
        self._check_initialized()

        gmsh.option.setNumber('Mesh.SaveGroupsOfElements', value)

        if self._verbose:
            print(f"[GmshBuilder] set_mesh_save_groups_of_elements = {value}")

    def set_mesh_save_groups_of_nodes(self, value):
        """Control which node groups are saved on export.

        Parameters
        ----------
        value : int
            0 = none, 1 = physical groups only, 2 = all elementary
            entity groups, 3 = all (physical + elementary).
        """
        import gmsh
        self._check_initialized()

        gmsh.option.setNumber('Mesh.SaveGroupsOfNodes', value)

        if self._verbose:
            print(f"[GmshBuilder] set_mesh_save_groups_of_nodes = {value}")

    def set_mesh_binary(self, enabled):
        """Enable or disable binary mesh file output.

        Parameters
        ----------
        enabled : bool
            If True, mesh files are written in binary format.
        """
        import gmsh
        self._check_initialized()

        gmsh.option.setNumber('Mesh.Binary', 1 if enabled else 0)

        if self._verbose:
            print(f"[GmshBuilder] set_mesh_binary = {enabled}")

    def set_mesh_random_factor(self, factor):
        """Set the random factor for initial Delaunay point placement.

        Controls the amount of randomization applied to the initial
        point positions in the Delaunay triangulation. Useful for
        avoiding degenerate configurations.

        Parameters
        ----------
        factor : float
            Random perturbation factor (e.g. 1e-9). Larger values
            increase randomization of initial point positions.
        """
        import gmsh
        self._check_initialized()

        gmsh.option.setNumber('Mesh.RandomFactor', factor)

        if self._verbose:
            print(f"[GmshBuilder] set_mesh_random_factor = {factor}")

    # ------------------------------------------------------------------
    # Color and visibility options
    # ------------------------------------------------------------------

    def set_entity_color(self, dim, tag, r, g, b, a=255):
        """Set the color of a geometric entity.

        Parameters
        ----------
        dim : int
            Entity dimension (0=point, 1=curve, 2=surface, 3=volume).
        tag : int
            Entity tag.
        r : int
            Red component (0-255).
        g : int
            Green component (0-255).
        b : int
            Blue component (0-255).
        a : int, optional
            Alpha component (0-255). Default is 255 (fully opaque).
        """
        import gmsh
        self._check_initialized()

        gmsh.model.setColor([(dim, tag)], r, g, b, a)

        if self._verbose:
            print(f"[GmshBuilder] set_entity_color dim={dim} tag={tag} "
                  f"rgba=({r},{g},{b},{a})")

    def get_entity_color(self, dim, tag):
        """Get the color of a geometric entity.

        Parameters
        ----------
        dim : int
            Entity dimension (0=point, 1=curve, 2=surface, 3=volume).
        tag : int
            Entity tag.

        Returns
        -------
        tuple of (int, int, int, int)
            Color as (r, g, b, a) with each component in 0-255.
        """
        import gmsh
        self._check_initialized()

        r, g, b, a = gmsh.model.getColor(dim, tag)

        if self._verbose:
            print(f"[GmshBuilder] get_entity_color dim={dim} tag={tag} "
                  f"rgba=({r},{g},{b},{a})")

        return (r, g, b, a)

    def set_entity_visibility(self, dim, tag, visible):
        """Set the visibility of a geometric entity.

        Parameters
        ----------
        dim : int
            Entity dimension (0=point, 1=curve, 2=surface, 3=volume).
        tag : int
            Entity tag.
        visible : bool
            If True, the entity is visible; if False, it is hidden.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.setVisibility([(dim, tag)], int(visible))

        if self._verbose:
            print(f"[GmshBuilder] set_entity_visibility dim={dim} tag={tag} "
                  f"visible={visible}")

    def get_entity_visibility(self, dim, tag):
        """Get the visibility of a geometric entity.

        Parameters
        ----------
        dim : int
            Entity dimension (0=point, 1=curve, 2=surface, 3=volume).
        tag : int
            Entity tag.

        Returns
        -------
        bool
            True if the entity is visible, False otherwise.
        """
        import gmsh
        self._check_initialized()

        value = gmsh.model.getVisibility(dim, tag)

        if self._verbose:
            print(f"[GmshBuilder] get_entity_visibility dim={dim} tag={tag} "
                  f"= {bool(value)}")

        return bool(value)
