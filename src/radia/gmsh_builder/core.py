"""
GmshBuilder core: context manager, state management, ID allocation.
"""


class CoreMixin:
    """Base mixin providing GMSH lifecycle and volume ID management."""

    def __init__(self, model_name='gmsh_model', verbose=True):
        self._model_name = model_name
        self._verbose = verbose
        self._volumes = {}        # {user_id: [gmsh_vol_tag, ...]}
        self._next_id = 1
        self._divisions = {}      # {gmsh_vol_tag: (nx, ny, nz)}
        self._global_divisions = None
        self._mesh_size = None
        self._initialized = False
        self._meshed = False
        self._physical_groups = {}  # {vol_id: (pg_tag, name)}
        self._materials = {}       # {vol_id: {'name': str, 'mu_r': float}}
        self._names = {}           # {vol_id: name}
        self._groups = {}          # {group_id: {'name': str, 'entities': {dim: [tags]}}}
        self._next_group_id = 1

    def __enter__(self):
        import gmsh
        gmsh.initialize()
        gmsh.model.add(self._model_name)
        gmsh.option.setNumber('General.Terminal', 1 if self._verbose else 0)
        self._initialized = True
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        import gmsh
        if self._initialized:
            gmsh.finalize()
            self._initialized = False
        return False

    def _check_initialized(self):
        if not self._initialized:
            raise RuntimeError("GmshBuilder must be used as context manager")

    def _alloc_id(self):
        uid = self._next_id
        self._next_id += 1
        return uid

    def _register_volume(self, gmsh_tags, label=None):
        """Register GMSH volume tags and return user ID."""
        uid = self._alloc_id()
        self._volumes[uid] = list(gmsh_tags) if isinstance(
            gmsh_tags, (list, tuple)) else [gmsh_tags]
        self._meshed = False
        return uid

    def _remove_volume(self, vol_id):
        """Remove a volume from tracking."""
        if vol_id in self._volumes:
            del self._volumes[vol_id]
        if vol_id in self._physical_groups:
            del self._physical_groups[vol_id]
        if vol_id in self._materials:
            del self._materials[vol_id]
        if vol_id in self._names:
            del self._names[vol_id]

    def reset(self):
        """Reset all geometry and mesh state (reset equivalent)."""
        import gmsh
        self._check_initialized()

        gmsh.clear()
        gmsh.model.add(self._model_name)
        self._volumes.clear()
        self._next_id = 1
        self._divisions.clear()
        self._global_divisions = None
        self._mesh_size = None
        self._meshed = False
        self._physical_groups.clear()
        self._materials.clear()
        self._names.clear()
        self._groups = {}
        self._next_group_id = 1

        if self._verbose:
            print("[GmshBuilder] reset")

    def compress_ids(self):
        """Re-number volume IDs to be contiguous starting from 1."""
        old_vols = dict(self._volumes)
        self._volumes.clear()
        self._next_id = 1

        id_map = {}
        for old_id in sorted(old_vols.keys()):
            new_id = self._alloc_id()
            self._volumes[new_id] = old_vols[old_id]
            id_map[old_id] = new_id

        # Remap related dicts
        new_pg = {}
        for old_id, val in self._physical_groups.items():
            if old_id in id_map:
                new_pg[id_map[old_id]] = val
        self._physical_groups = new_pg

        new_mat = {}
        for old_id, val in self._materials.items():
            if old_id in id_map:
                new_mat[id_map[old_id]] = val
        self._materials = new_mat

        new_names = {}
        for old_id, val in self._names.items():
            if old_id in id_map:
                new_names[id_map[old_id]] = val
        self._names = new_names

        if self._verbose:
            print(f"[GmshBuilder] compress_ids: {id_map}")
        return id_map

    def set_verbose(self, verbose):
        """Set verbose mode on/off.

        Parameters
        ----------
        verbose : bool
            True to enable verbose output, False to suppress.
        """
        import gmsh
        self._verbose = bool(verbose)
        if self._initialized:
            gmsh.option.setNumber('General.Terminal', 1 if self._verbose else 0)

    def get_verbose(self):
        """Return current verbose setting.

        Returns
        -------
        bool
            True if verbose mode is enabled.
        """
        return self._verbose

    def get_volume_count(self):
        """Return number of registered volumes.

        Returns
        -------
        int
            Number of volumes currently tracked.
        """
        return len(self._volumes)

    def get_all_volume_ids(self):
        """Return sorted list of all volume IDs.

        Returns
        -------
        list of int
            Sorted list of registered volume IDs.
        """
        return sorted(self._volumes.keys())

    def has_volume(self, vol_id):
        """Check if volume ID exists.

        Parameters
        ----------
        vol_id : int
            Volume ID to check.

        Returns
        -------
        bool
            True if the volume ID is registered.
        """
        return vol_id in self._volumes

    def get_gmsh_volume_tags(self, vol_id):
        """Return raw GMSH tags for a volume ID.

        Parameters
        ----------
        vol_id : int
            Volume ID to look up.

        Returns
        -------
        list of int
            List of GMSH volume tags associated with this volume ID.

        Raises
        ------
        KeyError
            If vol_id is not registered.
        """
        if vol_id not in self._volumes:
            raise KeyError(f"Volume ID {vol_id} not found")
        return list(self._volumes[vol_id])

    def set_model_name(self, name):
        """Set/rename GMSH model name.

        Parameters
        ----------
        name : str
            New model name.
        """
        import gmsh
        self._check_initialized()
        # Add new model name and make it current
        try:
            gmsh.model.setCurrent(name)
        except Exception:
            # Model doesn't exist yet, add it
            gmsh.model.add(name)
            gmsh.model.setCurrent(name)
        self._model_name = name

    def get_model_name(self):
        """Return current model name.

        Returns
        -------
        str
            The model name.
        """
        return self._model_name

    def is_initialized(self):
        """Return True if GMSH is initialized.

        Returns
        -------
        bool
            True if the GMSH context is active.
        """
        return self._initialized

    def is_meshed(self):
        """Return True if mesh has been generated.

        Returns
        -------
        bool
            True if the mesh has been generated for this model.
        """
        return self._meshed

    def snapshot(self):
        """Return a dict capturing current state.

        Returns
        -------
        dict
            Dictionary with keys: 'model_name', 'initialized', 'meshed',
            'volume_count', 'volume_ids', 'verbose'.
        """
        return {
            'model_name': self._model_name,
            'initialized': self._initialized,
            'meshed': self._meshed,
            'volume_count': len(self._volumes),
            'volume_ids': sorted(self._volumes.keys()),
            'verbose': self._verbose,
        }

    def log(self, message):
        """Print message if verbose mode is on.

        Parameters
        ----------
        message : str
            Message to print.
        """
        if self._verbose:
            print(f"[GmshBuilder] {message}")
