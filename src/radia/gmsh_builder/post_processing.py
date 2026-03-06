"""
Post-processing visualization for GmshBuilder.

Wraps GMSH's gmsh.view API for creating, populating, querying,
and configuring post-processing views (scalar/vector/tensor data
on nodes or elements).
"""

import numpy as np


class PostProcessingMixin:
    """Mixin providing post-processing view management."""

    # ------------------------------------------------------------------
    # View lifecycle
    # ------------------------------------------------------------------

    def add_view(self, name=''):
        """Create a new post-processing view.

        Parameters
        ----------
        name : str
            Optional view name.

        Returns
        -------
        int
            View tag.
        """
        import gmsh
        self._check_initialized()

        tag = gmsh.view.add(name)

        if self._verbose:
            print(f"[GmshBuilder] add_view name='{name}' -> tag={tag}")
        return tag

    def remove_view(self, view_tag):
        """Remove a post-processing view.

        Parameters
        ----------
        view_tag : int
            View tag returned by ``add_view``.
        """
        import gmsh
        self._check_initialized()

        gmsh.view.remove(view_tag)

        if self._verbose:
            print(f"[GmshBuilder] remove_view tag={view_tag}")

    def get_all_views(self):
        """Return tags of all existing views.

        Returns
        -------
        list of int
            View tags.
        """
        import gmsh
        self._check_initialized()

        return list(gmsh.view.getTags())

    def get_view_count(self):
        """Return the number of existing views.

        Returns
        -------
        int
            Number of views.
        """
        import gmsh
        self._check_initialized()

        return len(gmsh.view.getTags())

    def clear_all_views(self):
        """Remove all post-processing views."""
        import gmsh
        self._check_initialized()

        for tag in list(gmsh.view.getTags()):
            gmsh.view.remove(tag)

        if self._verbose:
            print("[GmshBuilder] clear_all_views")

    def write_view(self, view_tag, filename):
        """Write a view to file (.pos, .msh, .vtk, etc.).

        Parameters
        ----------
        view_tag : int
            View tag.
        filename : str
            Output file path.  Extension determines format.
        """
        import gmsh
        self._check_initialized()

        gmsh.view.write(view_tag, filename)

        if self._verbose:
            print(f"[GmshBuilder] write_view tag={view_tag} -> {filename}")

    def combine_views(self, what='steps', how='all', remove_original=True):
        """Combine existing views.

        Parameters
        ----------
        what : str
            Combination type: ``'elements'`` or ``'steps'``.
        how : str
            How to combine: ``'all'``, ``'visible'``, or ``'name'``.
        remove_original : bool
            If True, remove the original views after combining.
        """
        import gmsh
        self._check_initialized()

        gmsh.view.combine(what, how, 1 if remove_original else 0)

        if self._verbose:
            print(f"[GmshBuilder] combine_views what='{what}' "
                  f"remove_original={remove_original}")

    def add_alias_view(self, view_tag, copy_options=False):
        """Create an alias (shallow copy) of an existing view.

        Parameters
        ----------
        view_tag : int
            Source view tag.
        copy_options : bool
            If True, also copy display options.

        Returns
        -------
        int
            New view tag.
        """
        import gmsh
        self._check_initialized()

        new_tag = gmsh.view.addAlias(view_tag, copy_options)

        if self._verbose:
            print(f"[GmshBuilder] add_alias_view src={view_tag} -> {new_tag}")
        return new_tag

    # ------------------------------------------------------------------
    # Add data
    # ------------------------------------------------------------------

    def add_scalar_node_data(self, view_tag, step, time, node_tags, values):
        """Add scalar data at mesh nodes.

        Parameters
        ----------
        view_tag : int
            View tag.
        step : int
            Time step index.
        time : float
            Physical time value.
        node_tags : list of int
            GMSH node tags.
        values : list of float
            One scalar per node.
        """
        import gmsh
        self._check_initialized()

        data = [[v] for v in values]
        gmsh.view.addModelData(
            view_tag, step, self._model_name, 'NodeData',
            node_tags, data, time, numComponents=1)

        if self._verbose:
            print(f"[GmshBuilder] add_scalar_node_data tag={view_tag} "
                  f"step={step} n_nodes={len(node_tags)}")

    def add_vector_node_data(self, view_tag, step, time, node_tags, values):
        """Add 3-component vector data at mesh nodes.

        Parameters
        ----------
        view_tag : int
            View tag.
        step : int
            Time step index.
        time : float
            Physical time value.
        node_tags : list of int
            GMSH node tags.
        values : list of list
            Each inner list is ``[vx, vy, vz]``.
        """
        import gmsh
        self._check_initialized()

        data = [list(v) for v in values]
        gmsh.view.addModelData(
            view_tag, step, self._model_name, 'NodeData',
            node_tags, data, time, numComponents=3)

        if self._verbose:
            print(f"[GmshBuilder] add_vector_node_data tag={view_tag} "
                  f"step={step} n_nodes={len(node_tags)}")

    def add_tensor_node_data(self, view_tag, step, time, node_tags, values):
        """Add 9-component tensor data at mesh nodes.

        Parameters
        ----------
        view_tag : int
            View tag.
        step : int
            Time step index.
        time : float
            Physical time value.
        node_tags : list of int
            GMSH node tags.
        values : list of list
            Each inner list has 9 components (row-major 3x3 tensor).
        """
        import gmsh
        self._check_initialized()

        data = [list(v) for v in values]
        gmsh.view.addModelData(
            view_tag, step, self._model_name, 'NodeData',
            node_tags, data, time, numComponents=9)

        if self._verbose:
            print(f"[GmshBuilder] add_tensor_node_data tag={view_tag} "
                  f"step={step} n_nodes={len(node_tags)}")

    def add_scalar_element_data(self, view_tag, step, time,
                                elem_tags, values):
        """Add scalar data on mesh elements.

        Parameters
        ----------
        view_tag : int
            View tag.
        step : int
            Time step index.
        time : float
            Physical time value.
        elem_tags : list of int
            GMSH element tags.
        values : list of float
            One scalar per element.
        """
        import gmsh
        self._check_initialized()

        data = [[v] for v in values]
        gmsh.view.addModelData(
            view_tag, step, self._model_name, 'ElementData',
            elem_tags, data, time, numComponents=1)

        if self._verbose:
            print(f"[GmshBuilder] add_scalar_element_data tag={view_tag} "
                  f"step={step} n_elems={len(elem_tags)}")

    def add_vector_element_data(self, view_tag, step, time,
                                elem_tags, values):
        """Add 3-component vector data on mesh elements.

        Parameters
        ----------
        view_tag : int
            View tag.
        step : int
            Time step index.
        time : float
            Physical time value.
        elem_tags : list of int
            GMSH element tags.
        values : list of list
            Each inner list is ``[vx, vy, vz]``.
        """
        import gmsh
        self._check_initialized()

        data = [list(v) for v in values]
        gmsh.view.addModelData(
            view_tag, step, self._model_name, 'ElementData',
            elem_tags, data, time, numComponents=3)

        if self._verbose:
            print(f"[GmshBuilder] add_vector_element_data tag={view_tag} "
                  f"step={step} n_elems={len(elem_tags)}")

    def add_list_data(self, view_tag, data_type, n_elem, data):
        """Add list-based (non-model) data to a view.

        Useful for adding data not tied to the current GMSH model
        (e.g., custom point/line/triangle data).

        Parameters
        ----------
        view_tag : int
            View tag.
        data_type : str
            List data type (e.g. ``'SP'`` for scalar point,
            ``'VP'`` for vector point, ``'ST'`` for scalar triangle).
        n_elem : int
            Number of elements in the data list.
        data : list of float
            Flat list of data values (coordinates + field values).
        """
        import gmsh
        self._check_initialized()

        gmsh.view.addListData(view_tag, data_type, n_elem, data)

        if self._verbose:
            print(f"[GmshBuilder] add_list_data tag={view_tag} "
                  f"type='{data_type}' n_elem={n_elem}")

    def add_list_data_string(self, view_tag, coord, data, style=None):
        """Add string annotation to a view.

        Parameters
        ----------
        view_tag : int
            View tag.
        coord : list of float
            Position ``[x, y, z]`` (3D) or ``[x, y]`` (2D).
        data : list of str
            String data.
        style : list of str, optional
            Style descriptors (font, size, etc.).
        """
        import gmsh
        self._check_initialized()

        if style is None:
            style = []
        gmsh.view.addListDataString(view_tag, coord, data, style)

        if self._verbose:
            print(f"[GmshBuilder] add_list_data_string tag={view_tag} "
                  f"coord={coord}")

    def add_homogeneous_data(self, view_tag, step, time, data_type,
                             n_components, data, tags):
        """Add homogeneous model data (uniform number of values per entity).

        Parameters
        ----------
        view_tag : int
            View tag.
        step : int
            Time step index.
        time : float
            Physical time value.
        data_type : str
            ``'NodeData'``, ``'ElementData'``, or ``'ElementNodeData'``.
        n_components : int
            Number of components per data point (1=scalar, 3=vector, 9=tensor).
        data : list of float
            Flat list of all data values.
        tags : list of int
            Entity tags (node or element tags) corresponding to the data.
        """
        import gmsh
        self._check_initialized()

        try:
            gmsh.view.addHomogeneousModelData(
                view_tag, step, self._model_name, data_type,
                tags, data, time, n_components)
        except Exception:
            # Older GMSH version may have different API
            if self._verbose:
                print(f"[GmshBuilder] add_homogeneous_data: API not available")
            raise

        if self._verbose:
            print(f"[GmshBuilder] add_homogeneous_data tag={view_tag} "
                  f"step={step} type='{data_type}' nc={n_components}")

    def add_field_to_view(self, view_tag, field_callable, step=0, time=0.0):
        """Evaluate a callable at all mesh nodes and store as scalar NodeData.

        The callable receives three floats ``(x, y, z)`` and must return
        a scalar value (float) or a list/array whose length determines the
        number of components added (1=scalar, 3=vector).

        Parameters
        ----------
        view_tag : int
            View tag.
        field_callable : callable
            ``f(x, y, z) -> float`` or ``f(x, y, z) -> [vx, vy, vz]``.
        step : int
            Time step index.
        time : float
            Physical time value.
        """
        import gmsh
        self._check_initialized()

        node_tags, coords, _ = gmsh.model.mesh.getNodes()
        node_tags = list(node_tags)

        if len(node_tags) == 0:
            if self._verbose:
                print(f"[GmshBuilder] add_field_to_view: no mesh nodes")
            return

        coords = np.array(coords).reshape(-1, 3)

        # Evaluate at first node to determine number of components
        sample = field_callable(coords[0, 0], coords[0, 1], coords[0, 2])
        if isinstance(sample, (int, float)):
            n_comp = 1
        else:
            n_comp = len(sample)

        data = []
        for i in range(len(node_tags)):
            val = field_callable(coords[i, 0], coords[i, 1], coords[i, 2])
            if n_comp == 1:
                data.append([float(val)])
            else:
                data.append(list(val))

        gmsh.view.addModelData(
            view_tag, step, self._model_name, 'NodeData',
            node_tags, data, time, numComponents=n_comp)

        if self._verbose:
            print(f"[GmshBuilder] add_field_to_view tag={view_tag} "
                  f"n_nodes={len(node_tags)} nc={n_comp}")

    # ------------------------------------------------------------------
    # Get / query data
    # ------------------------------------------------------------------

    def get_view_model_data(self, view_tag, step=0):
        """Retrieve model data stored in a view.

        Parameters
        ----------
        view_tag : int
            View tag.
        step : int
            Time step index.

        Returns
        -------
        dict
            ``'tags'``: list of node/element tags,
            ``'data'``: list of value lists,
            ``'time'``: time value,
            ``'n_components'``: number of components per data point.
        """
        import gmsh
        self._check_initialized()

        data_type, tags, data, time, n_comp = gmsh.view.getModelData(
            view_tag, step)

        return {
            'data_type': data_type,
            'tags': list(tags),
            'data': [list(d) for d in data],
            'time': time,
            'n_components': n_comp,
        }

    def get_view_list_data(self, view_tag):
        """Retrieve list-based data stored in a view.

        Parameters
        ----------
        view_tag : int
            View tag.

        Returns
        -------
        dict
            ``'data_type'``: list of type strings,
            ``'n_elem'``: list of element counts,
            ``'data'``: list of data arrays.
        """
        import gmsh
        self._check_initialized()

        data_types, n_elems, data = gmsh.view.getListData(view_tag)

        return {
            'data_type': list(data_types),
            'n_elem': list(n_elems),
            'data': [list(d) for d in data],
        }

    def probe_view(self, view_tag, x, y, z, step=-1):
        """Probe a view at a given spatial location.

        Parameters
        ----------
        view_tag : int
            View tag.
        x, y, z : float
            Coordinates.
        step : int
            Time step (-1 for all steps).

        Returns
        -------
        list of float
            Interpolated values at the probe location.
        """
        import gmsh
        self._check_initialized()

        values, distance = gmsh.view.probe(view_tag, x, y, z, step)
        return list(values)

    def get_view_data_range(self, view_tag):
        """Return the (min, max) value range across all model data in a view.

        Only inspects the first component of each data point.

        Parameters
        ----------
        view_tag : int
            View tag.

        Returns
        -------
        tuple
            ``(min_val, max_val)``.
        """
        import gmsh
        self._check_initialized()

        # Try to determine step count by probing step 0
        try:
            _, tags, data, _, n_comp = gmsh.view.getModelData(
                view_tag, 0)
        except Exception:
            return (0.0, 0.0)

        if len(data) == 0:
            return (0.0, 0.0)

        all_vals = []
        for d in data:
            if len(d) > 0:
                all_vals.append(d[0])

        if not all_vals:
            return (0.0, 0.0)

        return (float(min(all_vals)), float(max(all_vals)))

    # ------------------------------------------------------------------
    # View options
    # ------------------------------------------------------------------

    def set_view_option_number(self, view_tag, name, value):
        """Set a numeric option for a view.

        Parameters
        ----------
        view_tag : int
            View tag.
        name : str
            Option name (e.g. ``'RangeType'``, ``'IntervalsType'``,
            ``'NbIso'``, ``'Visible'``).
        value : float
            Option value.
        """
        import gmsh
        self._check_initialized()

        gmsh.view.option.setNumber(view_tag, name, value)

        if self._verbose:
            print(f"[GmshBuilder] set_view_option_number tag={view_tag} "
                  f"{name}={value}")

    def set_view_option_string(self, view_tag, name, value):
        """Set a string option for a view.

        Parameters
        ----------
        view_tag : int
            View tag.
        name : str
            Option name (e.g. ``'Name'``, ``'Format'``).
        value : str
            Option value.
        """
        import gmsh
        self._check_initialized()

        gmsh.view.option.setString(view_tag, name, value)

        if self._verbose:
            print(f"[GmshBuilder] set_view_option_string tag={view_tag} "
                  f"{name}='{value}'")

    def get_view_option_number(self, view_tag, name):
        """Get a numeric option from a view.

        Parameters
        ----------
        view_tag : int
            View tag.
        name : str
            Option name.

        Returns
        -------
        float
            Option value.
        """
        import gmsh
        self._check_initialized()

        return gmsh.view.option.getNumber(view_tag, name)

    def get_view_option_string(self, view_tag, name):
        """Get a string option from a view.

        Parameters
        ----------
        view_tag : int
            View tag.
        name : str
            Option name.

        Returns
        -------
        str
            Option value.
        """
        import gmsh
        self._check_initialized()

        return gmsh.view.option.getString(view_tag, name)
