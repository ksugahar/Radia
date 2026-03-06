"""
Physical group / Block / Sideset / Nodeset management for GmshBuilder.

Gmsh equivalents: block, sideset, nodeset, material assignment.
"""


class PhysicalGroupMixin:
    """Mixin providing physical group management."""

    def set_physical_group(self, vol_id, name):
        """Assign physical group to a volume (Gmsh 'block add volume').

        Parameters
        ----------
        vol_id : int
            Volume ID.
        name : str
            Physical group name (e.g., 'iron_core', 'coil').
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        tags = self._volumes[vol_id]
        pg = gmsh.model.addPhysicalGroup(3, tags)
        gmsh.model.setPhysicalName(3, pg, name)
        self._physical_groups[vol_id] = (pg, name)

        if self._verbose:
            print(f"[GmshBuilder] set_physical_group vol={vol_id} name='{name}'")

    def add_block(self, vol_ids, block_id=None, name=None):
        """Add a block (Gmsh 'block add volume').

        Parameters
        ----------
        vol_ids : list of int
            Volume IDs to include.
        block_id : int, optional
            Block ID. Auto-assigned if None.
        name : str, optional
            Block name.

        Returns
        -------
        int
            Physical group tag.
        """
        import gmsh
        self._check_initialized()

        all_tags = []
        for vid in vol_ids:
            if vid not in self._volumes:
                raise KeyError(f"Volume {vid} not found")
            all_tags.extend(self._volumes[vid])

        if block_id is not None:
            pg = gmsh.model.addPhysicalGroup(3, all_tags, block_id)
        else:
            pg = gmsh.model.addPhysicalGroup(3, all_tags)

        if name:
            gmsh.model.setPhysicalName(3, pg, name)

        for vid in vol_ids:
            self._physical_groups[vid] = (pg, name or '')

        if self._verbose:
            print(f"[GmshBuilder] add_block vols={vol_ids} tag={pg} name='{name}'")
        return pg

    def add_sideset(self, surf_tags, sideset_id=None, name=None):
        """Add a sideset (Gmsh 'sideset add surface').

        Parameters
        ----------
        surf_tags : list of int
            GMSH surface tags.
        sideset_id : int, optional
            Sideset ID.
        name : str, optional
            Sideset name.

        Returns
        -------
        int
            Physical group tag.
        """
        import gmsh
        self._check_initialized()

        if sideset_id is not None:
            pg = gmsh.model.addPhysicalGroup(2, surf_tags, sideset_id)
        else:
            pg = gmsh.model.addPhysicalGroup(2, surf_tags)

        if name:
            gmsh.model.setPhysicalName(2, pg, name)

        if self._verbose:
            print(f"[GmshBuilder] add_sideset surfs={surf_tags} tag={pg} "
                  f"name='{name}'")
        return pg

    def add_nodeset(self, point_tags, nodeset_id=None, name=None):
        """Add a nodeset (Gmsh 'nodeset add vertex').

        Parameters
        ----------
        point_tags : list of int
            GMSH point tags.
        nodeset_id : int, optional
            Nodeset ID.
        name : str, optional
            Nodeset name.

        Returns
        -------
        int
            Physical group tag.
        """
        import gmsh
        self._check_initialized()

        if nodeset_id is not None:
            pg = gmsh.model.addPhysicalGroup(0, point_tags, nodeset_id)
        else:
            pg = gmsh.model.addPhysicalGroup(0, point_tags)

        if name:
            gmsh.model.setPhysicalName(0, pg, name)

        if self._verbose:
            print(f"[GmshBuilder] add_nodeset points={point_tags} tag={pg}")
        return pg

    def get_physical_groups(self, dim=None):
        """Get all physical groups (Gmsh 'get_block_id_list').

        Parameters
        ----------
        dim : int, optional
            Filter by dimension. None for all.

        Returns
        -------
        list of dict
            Physical group info.
        """
        import gmsh
        self._check_initialized()

        if dim is not None:
            groups = gmsh.model.getPhysicalGroups(dim)
        else:
            groups = gmsh.model.getPhysicalGroups()

        result = []
        for d, tag in groups:
            name = gmsh.model.getPhysicalName(d, tag)
            entities = gmsh.model.getEntitiesForPhysicalGroup(d, tag)
            result.append({
                'dim': d,
                'tag': tag,
                'name': name,
                'entities': list(entities),
            })

        return result

    def get_physical_group_name(self, tag, dim=3):
        """Get physical group name.

        Parameters
        ----------
        tag : int
            Physical group tag.
        dim : int
            Dimension.

        Returns
        -------
        str
        """
        import gmsh
        self._check_initialized()

        return gmsh.model.getPhysicalName(dim, tag)

    def get_physical_group_entities(self, tag, dim=3):
        """Get entities in a physical group.

        Parameters
        ----------
        tag : int
            Physical group tag.
        dim : int
            Dimension.

        Returns
        -------
        list of int
            Entity tags.
        """
        import gmsh
        self._check_initialized()

        return list(gmsh.model.getEntitiesForPhysicalGroup(dim, tag))

    def remove_physical_group(self, tag, dim=3):
        """Remove a physical group.

        Parameters
        ----------
        tag : int
            Physical group tag.
        dim : int
            Dimension.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.removePhysicalGroups([(dim, tag)])

    def remove_all_physical_groups(self):
        """Remove all physical groups."""
        import gmsh
        self._check_initialized()

        gmsh.model.removePhysicalGroups()

    def auto_assign_blocks(self):
        """Auto-assign each volume to its own block.

        Returns
        -------
        dict
            {vol_id: block_tag}.
        """
        import gmsh
        self._check_initialized()

        result = {}
        for vid, tags in self._volumes.items():
            pg = gmsh.model.addPhysicalGroup(3, tags)
            name = self._names.get(vid, f'volume_{vid}')
            gmsh.model.setPhysicalName(3, pg, name)
            self._physical_groups[vid] = (pg, name)
            result[vid] = pg

        if self._verbose:
            print(f"[GmshBuilder] auto_assign_blocks -> {len(result)} blocks")
        return result

    def auto_assign_sidesets(self):
        """Auto-assign each surface to its own sideset.

        Returns
        -------
        dict
            {surf_tag: sideset_tag}.
        """
        import gmsh
        self._check_initialized()

        result = {}
        surfaces = gmsh.model.getEntities(2)
        for dim, tag in surfaces:
            pg = gmsh.model.addPhysicalGroup(2, [tag])
            gmsh.model.setPhysicalName(2, pg, f'surface_{tag}')
            result[tag] = pg

        if self._verbose:
            print(f"[GmshBuilder] auto_assign_sidesets -> {len(result)} sidesets")
        return result

    def set_material(self, vol_id, name, mu_r=1.0):
        """Assign material properties to a volume.

        Parameters
        ----------
        vol_id : int
            Volume ID.
        name : str
            Material name (e.g., 'iron', 'copper').
        mu_r : float
            Relative permeability.
        """
        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        self._materials[vol_id] = {'name': name, 'mu_r': mu_r}

        if self._verbose:
            print(f"[GmshBuilder] set_material vol={vol_id} name='{name}' "
                  f"mu_r={mu_r}")

    def get_materials(self):
        """Get all material assignments.

        Returns
        -------
        dict
            {vol_id: {'name': str, 'mu_r': float}}.
        """
        return dict(self._materials)

    def assign_boundary_names(self, mapping):
        """Assign boundary condition names to surfaces.

        Parameters
        ----------
        mapping : dict
            {surf_tag: name} or {name: [surf_tags]}.
        """
        import gmsh
        self._check_initialized()

        for key, value in mapping.items():
            if isinstance(value, str):
                # {surf_tag: name}
                pg = gmsh.model.addPhysicalGroup(2, [key])
                gmsh.model.setPhysicalName(2, pg, value)
            else:
                # {name: [surf_tags]}
                pg = gmsh.model.addPhysicalGroup(2, value)
                gmsh.model.setPhysicalName(2, pg, key)

        if self._verbose:
            print(f"[GmshBuilder] assign_boundary_names: {len(mapping)} entries")

    def add_block_by_name(self, vol_ids, name):
        """Add a block with auto-assigned ID and given name.

        Convenience wrapper around add_block() that always auto-assigns
        the block ID and requires a name.

        Parameters
        ----------
        vol_ids : list of int
            Volume IDs to include in the block.
        name : str
            Block name.

        Returns
        -------
        int
            Physical group tag.
        """
        import gmsh
        self._check_initialized()

        all_tags = []
        for vid in vol_ids:
            if vid not in self._volumes:
                raise KeyError(f"Volume {vid} not found")
            all_tags.extend(self._volumes[vid])

        pg = gmsh.model.addPhysicalGroup(3, all_tags)
        gmsh.model.setPhysicalName(3, pg, name)

        for vid in vol_ids:
            self._physical_groups[vid] = (pg, name)

        if self._verbose:
            print(f"[GmshBuilder] add_block_by_name vols={vol_ids} tag={pg} "
                  f"name='{name}'")
        return pg

    def add_sideset_by_name(self, surf_tags, name):
        """Add a sideset with auto-assigned ID and given name.

        Convenience wrapper around add_sideset() that always auto-assigns
        the sideset ID and requires a name.

        Parameters
        ----------
        surf_tags : list of int
            GMSH surface tags.
        name : str
            Sideset name.

        Returns
        -------
        int
            Physical group tag.
        """
        import gmsh
        self._check_initialized()

        pg = gmsh.model.addPhysicalGroup(2, surf_tags)
        gmsh.model.setPhysicalName(2, pg, name)

        if self._verbose:
            print(f"[GmshBuilder] add_sideset_by_name surfs={surf_tags} tag={pg} "
                  f"name='{name}'")
        return pg

    def add_nodeset_by_name(self, point_tags, name):
        """Add a nodeset with auto-assigned ID and given name.

        Convenience wrapper around add_nodeset() that always auto-assigns
        the nodeset ID and requires a name.

        Parameters
        ----------
        point_tags : list of int
            GMSH point tags.
        name : str
            Nodeset name.

        Returns
        -------
        int
            Physical group tag.
        """
        import gmsh
        self._check_initialized()

        pg = gmsh.model.addPhysicalGroup(0, point_tags)
        gmsh.model.setPhysicalName(0, pg, name)

        if self._verbose:
            print(f"[GmshBuilder] add_nodeset_by_name points={point_tags} "
                  f"tag={pg} name='{name}'")
        return pg

    def rename_physical_group(self, tag, dim, new_name):
        """Rename a physical group.

        Parameters
        ----------
        tag : int
            Physical group tag.
        dim : int
            Dimension (0=nodeset, 2=sideset, 3=block).
        new_name : str
            New name for the physical group.
        """
        import gmsh
        self._check_initialized()

        # GMSH doesn't update names in-place; remove and re-add
        entities = list(gmsh.model.getEntitiesForPhysicalGroup(dim, tag))
        gmsh.model.removePhysicalGroups([(dim, tag)])
        new_tag = gmsh.model.addPhysicalGroup(dim, entities, tag)
        gmsh.model.setPhysicalName(dim, new_tag, new_name)

        # Update internal tracking if this is a 3D block
        if dim == 3:
            for vid, (pg, _old_name) in self._physical_groups.items():
                if pg == tag:
                    self._physical_groups[vid] = (pg, new_name)

        if self._verbose:
            print(f"[GmshBuilder] rename_physical_group dim={dim} tag={tag} "
                  f"new_name='{new_name}'")

    def get_block_count(self):
        """Count number of blocks (3D physical groups).

        Returns
        -------
        int
            Number of 3D physical groups.
        """
        import gmsh
        self._check_initialized()

        return len(gmsh.model.getPhysicalGroups(3))

    def get_sideset_count(self):
        """Count number of sidesets (2D physical groups).

        Returns
        -------
        int
            Number of 2D physical groups.
        """
        import gmsh
        self._check_initialized()

        return len(gmsh.model.getPhysicalGroups(2))

    def get_nodeset_count(self):
        """Count number of nodesets (0D physical groups).

        Returns
        -------
        int
            Number of 0D physical groups.
        """
        import gmsh
        self._check_initialized()

        return len(gmsh.model.getPhysicalGroups(0))

    def get_all_blocks(self):
        """Return list of all block info dicts.

        Returns
        -------
        list of dict
            Each dict has keys: 'tag', 'name', 'entities'.
        """
        import gmsh
        self._check_initialized()

        result = []
        for _dim, tag in gmsh.model.getPhysicalGroups(3):
            name = gmsh.model.getPhysicalName(3, tag)
            entities = list(gmsh.model.getEntitiesForPhysicalGroup(3, tag))
            result.append({
                'tag': tag,
                'name': name,
                'entities': entities,
            })
        return result

    def get_all_sidesets(self):
        """Return list of all sideset info dicts.

        Returns
        -------
        list of dict
            Each dict has keys: 'tag', 'name', 'entities'.
        """
        import gmsh
        self._check_initialized()

        result = []
        for _dim, tag in gmsh.model.getPhysicalGroups(2):
            name = gmsh.model.getPhysicalName(2, tag)
            entities = list(gmsh.model.getEntitiesForPhysicalGroup(2, tag))
            result.append({
                'tag': tag,
                'name': name,
                'entities': entities,
            })
        return result

    def get_all_nodesets(self):
        """Return list of all nodeset info dicts.

        Returns
        -------
        list of dict
            Each dict has keys: 'tag', 'name', 'entities'.
        """
        import gmsh
        self._check_initialized()

        result = []
        for _dim, tag in gmsh.model.getPhysicalGroups(0):
            name = gmsh.model.getPhysicalName(0, tag)
            entities = list(gmsh.model.getEntitiesForPhysicalGroup(0, tag))
            result.append({
                'tag': tag,
                'name': name,
                'entities': entities,
            })
        return result

    def get_block_name(self, block_tag):
        """Get name of a block (3D physical group).

        Parameters
        ----------
        block_tag : int
            Block tag.

        Returns
        -------
        str
            Block name.
        """
        return self.get_physical_group_name(block_tag, dim=3)

    def get_sideset_name(self, ss_tag):
        """Get name of a sideset (2D physical group).

        Parameters
        ----------
        ss_tag : int
            Sideset tag.

        Returns
        -------
        str
            Sideset name.
        """
        return self.get_physical_group_name(ss_tag, dim=2)

    def get_nodeset_name(self, ns_tag):
        """Get name of a nodeset (0D physical group).

        Parameters
        ----------
        ns_tag : int
            Nodeset tag.

        Returns
        -------
        str
            Nodeset name.
        """
        return self.get_physical_group_name(ns_tag, dim=0)

    def rename_block(self, block_tag, new_name):
        """Rename a block (3D physical group).

        Parameters
        ----------
        block_tag : int
            Block tag.
        new_name : str
            New name for the block.
        """
        self.rename_physical_group(block_tag, 3, new_name)

    def get_material(self, vol_id):
        """Get material for a specific volume.

        Parameters
        ----------
        vol_id : int
            Volume ID.

        Returns
        -------
        dict or None
            Material dict with keys 'name', 'mu_r', or None if no
            material is assigned.
        """
        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        return self._materials.get(vol_id, None)

    def set_material_by_name(self, vol_ids, name, properties):
        """Set material on multiple volumes with properties dict.

        Parameters
        ----------
        vol_ids : list of int
            Volume IDs.
        name : str
            Material name (e.g., 'iron', 'copper').
        properties : dict
            Material properties (e.g., {'mu_r': 1000.0, 'sigma': 5.8e7}).
        """
        for vid in vol_ids:
            if vid not in self._volumes:
                raise KeyError(f"Volume {vid} not found")

        for vid in vol_ids:
            mat_entry = dict(properties)
            mat_entry['name'] = name
            self._materials[vid] = mat_entry

        if self._verbose:
            print(f"[GmshBuilder] set_material_by_name vols={vol_ids} "
                  f"name='{name}' properties={properties}")

    def has_material(self, vol_id):
        """Check if volume has material assigned.

        Parameters
        ----------
        vol_id : int
            Volume ID.

        Returns
        -------
        bool
            True if the volume has a material assignment.
        """
        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        return vol_id in self._materials

    def clear_materials(self):
        """Remove all material assignments."""
        self._materials.clear()

        if self._verbose:
            print("[GmshBuilder] clear_materials: all materials removed")

    def auto_assign_all(self):
        """Auto-assign blocks and sidesets in one call.

        Calls auto_assign_blocks() and auto_assign_sidesets() sequentially.

        Returns
        -------
        dict
            {'blocks': {vol_id: block_tag},
             'sidesets': {surf_tag: sideset_tag}}.
        """
        blocks = self.auto_assign_blocks()
        sidesets = self.auto_assign_sidesets()

        if self._verbose:
            print(f"[GmshBuilder] auto_assign_all: {len(blocks)} blocks, "
                  f"{len(sidesets)} sidesets")
        return {'blocks': blocks, 'sidesets': sidesets}

    def export_material_table(self):
        """Return formatted table of all materials.

        Returns
        -------
        str
            Formatted table string with columns: Volume ID, Name,
            and all property keys.
        """
        if not self._materials:
            return "No materials assigned."

        # Collect all property keys across all materials
        all_keys = set()
        for mat in self._materials.values():
            all_keys.update(k for k in mat if k != 'name')
        all_keys = sorted(all_keys)

        # Build header
        headers = ['Vol ID', 'Name'] + all_keys
        col_widths = [max(len(h), 8) for h in headers]

        # Compute column widths from data
        rows = []
        for vid in sorted(self._materials):
            mat = self._materials[vid]
            row = [str(vid), mat.get('name', '')]
            for key in all_keys:
                val = mat.get(key, '')
                row.append(str(val))
            rows.append(row)
            for i, cell in enumerate(row):
                col_widths[i] = max(col_widths[i], len(cell))

        # Format table
        fmt = '  '.join(f'{{:<{w}}}' for w in col_widths)
        lines = [fmt.format(*headers)]
        lines.append('  '.join('-' * w for w in col_widths))
        for row in rows:
            lines.append(fmt.format(*row))

        return '\n'.join(lines)

    # ------------------------------------------------------------------
    # Block element access (methods 1-10)
    # ------------------------------------------------------------------

    def get_block_elements(self, block_tag):
        """Get all elements in a block (3D physical group).

        Retrieves entities for the physical group ``(3, block_tag)`` and
        queries the mesh elements on each entity.

        Parameters
        ----------
        block_tag : int
            Physical group tag (dimension 3).

        Returns
        -------
        list of dict
            Each dict has keys ``'type'`` (int, GMSH element type),
            ``'tags'`` (list of int, element tags), and ``'node_tags'``
            (list of int, node tags for those elements).
        """
        import gmsh
        self._check_initialized()

        vol_tags = gmsh.model.getEntitiesForPhysicalGroup(3, block_tag)
        result = []
        for vtag in vol_tags:
            elem_types, elem_tags, node_tags = gmsh.model.mesh.getElements(3, vtag)
            for etype, etags, ntags in zip(elem_types, elem_tags, node_tags):
                result.append({
                    'type': int(etype),
                    'tags': list(etags),
                    'node_tags': list(ntags),
                })

        if self._verbose:
            total = sum(len(e['tags']) for e in result)
            print(f"[GmshBuilder] get_block_elements block={block_tag} "
                  f"-> {total} elements")
        return result

    def get_block_hex_count(self, block_tag):
        """Count hexahedral elements (GMSH type 5) in a block.

        Parameters
        ----------
        block_tag : int
            Physical group tag (dimension 3).

        Returns
        -------
        int
            Number of hexahedral elements.
        """
        import gmsh
        self._check_initialized()

        vol_tags = gmsh.model.getEntitiesForPhysicalGroup(3, block_tag)
        count = 0
        for vtag in vol_tags:
            elem_types, elem_tags, _node_tags = gmsh.model.mesh.getElements(3, vtag)
            for etype, etags in zip(elem_types, elem_tags):
                if int(etype) == 5:
                    count += len(etags)

        if self._verbose:
            print(f"[GmshBuilder] get_block_hex_count block={block_tag} "
                  f"-> {count}")
        return count

    def get_block_tet_count(self, block_tag):
        """Count tetrahedral elements (GMSH type 4) in a block.

        Parameters
        ----------
        block_tag : int
            Physical group tag (dimension 3).

        Returns
        -------
        int
            Number of tetrahedral elements.
        """
        import gmsh
        self._check_initialized()

        vol_tags = gmsh.model.getEntitiesForPhysicalGroup(3, block_tag)
        count = 0
        for vtag in vol_tags:
            elem_types, elem_tags, _node_tags = gmsh.model.mesh.getElements(3, vtag)
            for etype, etags in zip(elem_types, elem_tags):
                if int(etype) == 4:
                    count += len(etags)

        if self._verbose:
            print(f"[GmshBuilder] get_block_tet_count block={block_tag} "
                  f"-> {count}")
        return count

    def get_block_node_count(self, block_tag):
        """Count nodes belonging to a block.

        Uses ``gmsh.model.mesh.getNodesForPhysicalGroup(3, block_tag)``
        to obtain the unique node tags.

        Parameters
        ----------
        block_tag : int
            Physical group tag (dimension 3).

        Returns
        -------
        int
            Number of unique nodes in the block.
        """
        import gmsh
        self._check_initialized()

        node_tags, _coords = gmsh.model.mesh.getNodesForPhysicalGroup(3, block_tag)
        count = len(node_tags)

        if self._verbose:
            print(f"[GmshBuilder] get_block_node_count block={block_tag} "
                  f"-> {count}")
        return count

    def get_block_nodes(self, block_tag):
        """Get node coordinates for a block.

        Parameters
        ----------
        block_tag : int
            Physical group tag (dimension 3).

        Returns
        -------
        numpy.ndarray
            Node coordinates with shape ``(n, 3)``.
        """
        import numpy as np
        import gmsh
        self._check_initialized()

        _node_tags, coords = gmsh.model.mesh.getNodesForPhysicalGroup(3, block_tag)
        nodes = np.array(coords, dtype=np.float64).reshape(-1, 3)

        if self._verbose:
            print(f"[GmshBuilder] get_block_nodes block={block_tag} "
                  f"-> {len(nodes)} nodes")
        return nodes

    def get_block_element_count(self, block_tag):
        """Count total elements of all types in a block.

        Parameters
        ----------
        block_tag : int
            Physical group tag (dimension 3).

        Returns
        -------
        int
            Total number of elements.
        """
        import gmsh
        self._check_initialized()

        vol_tags = gmsh.model.getEntitiesForPhysicalGroup(3, block_tag)
        count = 0
        for vtag in vol_tags:
            elem_types, elem_tags, _node_tags = gmsh.model.mesh.getElements(3, vtag)
            for etags in elem_tags:
                count += len(etags)

        if self._verbose:
            print(f"[GmshBuilder] get_block_element_count block={block_tag} "
                  f"-> {count}")
        return count

    def get_block_volumes(self, block_tag):
        """Get volume entity tags for a block.

        Parameters
        ----------
        block_tag : int
            Physical group tag (dimension 3).

        Returns
        -------
        list of int
            Volume (entity) tags belonging to this block.
        """
        import gmsh
        self._check_initialized()

        vol_tags = list(gmsh.model.getEntitiesForPhysicalGroup(3, block_tag))

        if self._verbose:
            print(f"[GmshBuilder] get_block_volumes block={block_tag} "
                  f"-> {vol_tags}")
        return vol_tags

    def get_block_surfaces(self, block_tag):
        """Get boundary surface tags of all volumes in a block.

        Parameters
        ----------
        block_tag : int
            Physical group tag (dimension 3).

        Returns
        -------
        list of int
            Surface tags forming the boundary of the block volumes.
        """
        import gmsh
        self._check_initialized()

        vol_tags = gmsh.model.getEntitiesForPhysicalGroup(3, block_tag)
        surfaces = set()
        for vtag in vol_tags:
            boundary = gmsh.model.getBoundary([(3, vtag)], oriented=False)
            for _dim, stag in boundary:
                surfaces.add(abs(stag))

        result = sorted(surfaces)

        if self._verbose:
            print(f"[GmshBuilder] get_block_surfaces block={block_tag} "
                  f"-> {len(result)} surfaces")
        return result

    def get_block_bounding_box(self, block_tag):
        """Get bounding box of all volumes in a block.

        Computes the union of bounding boxes across all volume entities
        in the physical group.

        Parameters
        ----------
        block_tag : int
            Physical group tag (dimension 3).

        Returns
        -------
        dict
            Dictionary with keys ``'xmin'``, ``'ymin'``, ``'zmin'``,
            ``'xmax'``, ``'ymax'``, ``'zmax'``.
        """
        import gmsh
        self._check_initialized()

        vol_tags = gmsh.model.getEntitiesForPhysicalGroup(3, block_tag)
        if len(vol_tags) == 0:
            raise ValueError(f"Block {block_tag} has no volume entities")

        xmin_g, ymin_g, zmin_g = float('inf'), float('inf'), float('inf')
        xmax_g, ymax_g, zmax_g = float('-inf'), float('-inf'), float('-inf')

        for vtag in vol_tags:
            xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(3, vtag)
            xmin_g = min(xmin_g, xmin)
            ymin_g = min(ymin_g, ymin)
            zmin_g = min(zmin_g, zmin)
            xmax_g = max(xmax_g, xmax)
            ymax_g = max(ymax_g, ymax)
            zmax_g = max(zmax_g, zmax)

        result = {
            'xmin': xmin_g, 'ymin': ymin_g, 'zmin': zmin_g,
            'xmax': xmax_g, 'ymax': ymax_g, 'zmax': zmax_g,
        }

        if self._verbose:
            print(f"[GmshBuilder] get_block_bounding_box block={block_tag} "
                  f"-> {result}")
        return result

    def get_block_center(self, block_tag):
        """Get center of a block as the mean of its volume centers.

        Parameters
        ----------
        block_tag : int
            Physical group tag (dimension 3).

        Returns
        -------
        list of float
            ``[x, y, z]`` center coordinates.
        """
        import gmsh
        self._check_initialized()

        vol_tags = gmsh.model.getEntitiesForPhysicalGroup(3, block_tag)
        if len(vol_tags) == 0:
            raise ValueError(f"Block {block_tag} has no volume entities")

        cx, cy, cz = 0.0, 0.0, 0.0
        total_mass = 0.0
        n = len(vol_tags)
        for vtag in vol_tags:
            m = gmsh.model.occ.getMass(3, vtag)
            x, y, z = gmsh.model.occ.getCenterOfMass(3, vtag)
            cx += m * x
            cy += m * y
            cz += m * z
            total_mass += m

        if total_mass > 0:
            center = [cx / total_mass, cy / total_mass, cz / total_mass]
        else:
            # Fallback to unweighted average if volumes have zero mass
            cx, cy, cz = 0.0, 0.0, 0.0
            for vtag in vol_tags:
                x, y, z = gmsh.model.occ.getCenterOfMass(3, vtag)
                cx += x
                cy += y
                cz += z
            center = [cx / n, cy / n, cz / n]

        if self._verbose:
            print(f"[GmshBuilder] get_block_center block={block_tag} "
                  f"-> {center}")
        return center

    # ------------------------------------------------------------------
    # Sideset element access (methods 11-18)
    # ------------------------------------------------------------------

    def get_sideset_elements(self, ss_tag):
        """Get all elements on sideset surfaces.

        Parameters
        ----------
        ss_tag : int
            Physical group tag (dimension 2).

        Returns
        -------
        list of dict
            Each dict has keys ``'type'`` (int), ``'tags'`` (list of int),
            and ``'node_tags'`` (list of int).
        """
        import gmsh
        self._check_initialized()

        surf_tags = gmsh.model.getEntitiesForPhysicalGroup(2, ss_tag)
        result = []
        for stag in surf_tags:
            elem_types, elem_tags, node_tags = gmsh.model.mesh.getElements(2, stag)
            for etype, etags, ntags in zip(elem_types, elem_tags, node_tags):
                result.append({
                    'type': int(etype),
                    'tags': list(etags),
                    'node_tags': list(ntags),
                })

        if self._verbose:
            total = sum(len(e['tags']) for e in result)
            print(f"[GmshBuilder] get_sideset_elements ss={ss_tag} "
                  f"-> {total} elements")
        return result

    def get_sideset_tri_count(self, ss_tag):
        """Count triangle elements (GMSH type 2) in a sideset.

        Parameters
        ----------
        ss_tag : int
            Physical group tag (dimension 2).

        Returns
        -------
        int
            Number of triangular elements.
        """
        import gmsh
        self._check_initialized()

        surf_tags = gmsh.model.getEntitiesForPhysicalGroup(2, ss_tag)
        count = 0
        for stag in surf_tags:
            elem_types, elem_tags, _node_tags = gmsh.model.mesh.getElements(2, stag)
            for etype, etags in zip(elem_types, elem_tags):
                if int(etype) == 2:
                    count += len(etags)

        if self._verbose:
            print(f"[GmshBuilder] get_sideset_tri_count ss={ss_tag} "
                  f"-> {count}")
        return count

    def get_sideset_quad_count(self, ss_tag):
        """Count quadrilateral elements (GMSH type 3) in a sideset.

        Parameters
        ----------
        ss_tag : int
            Physical group tag (dimension 2).

        Returns
        -------
        int
            Number of quadrilateral elements.
        """
        import gmsh
        self._check_initialized()

        surf_tags = gmsh.model.getEntitiesForPhysicalGroup(2, ss_tag)
        count = 0
        for stag in surf_tags:
            elem_types, elem_tags, _node_tags = gmsh.model.mesh.getElements(2, stag)
            for etype, etags in zip(elem_types, elem_tags):
                if int(etype) == 3:
                    count += len(etags)

        if self._verbose:
            print(f"[GmshBuilder] get_sideset_quad_count ss={ss_tag} "
                  f"-> {count}")
        return count

    def get_sideset_node_count(self, ss_tag):
        """Count unique nodes in a sideset.

        Parameters
        ----------
        ss_tag : int
            Physical group tag (dimension 2).

        Returns
        -------
        int
            Number of unique nodes.
        """
        import gmsh
        self._check_initialized()

        node_tags, _coords = gmsh.model.mesh.getNodesForPhysicalGroup(2, ss_tag)
        count = len(node_tags)

        if self._verbose:
            print(f"[GmshBuilder] get_sideset_node_count ss={ss_tag} "
                  f"-> {count}")
        return count

    def get_sideset_nodes(self, ss_tag):
        """Get node coordinates for a sideset.

        Parameters
        ----------
        ss_tag : int
            Physical group tag (dimension 2).

        Returns
        -------
        numpy.ndarray
            Node coordinates with shape ``(n, 3)``.
        """
        import numpy as np
        import gmsh
        self._check_initialized()

        _node_tags, coords = gmsh.model.mesh.getNodesForPhysicalGroup(2, ss_tag)
        nodes = np.array(coords, dtype=np.float64).reshape(-1, 3)

        if self._verbose:
            print(f"[GmshBuilder] get_sideset_nodes ss={ss_tag} "
                  f"-> {len(nodes)} nodes")
        return nodes

    def get_sideset_surfaces(self, ss_tag):
        """Get surface entity tags belonging to a sideset.

        Parameters
        ----------
        ss_tag : int
            Physical group tag (dimension 2).

        Returns
        -------
        list of int
            Surface tags.
        """
        import gmsh
        self._check_initialized()

        surf_tags = list(gmsh.model.getEntitiesForPhysicalGroup(2, ss_tag))

        if self._verbose:
            print(f"[GmshBuilder] get_sideset_surfaces ss={ss_tag} "
                  f"-> {surf_tags}")
        return surf_tags

    def get_sideset_area(self, ss_tag):
        """Compute total area of all surfaces in a sideset.

        Uses ``gmsh.model.occ.getMass`` (or ``gmsh.model.getMass``)
        to obtain each surface area.

        Parameters
        ----------
        ss_tag : int
            Physical group tag (dimension 2).

        Returns
        -------
        float
            Total area of the sideset surfaces.
        """
        import gmsh
        self._check_initialized()

        surf_tags = gmsh.model.getEntitiesForPhysicalGroup(2, ss_tag)
        total_area = 0.0
        for stag in surf_tags:
            area = gmsh.model.occ.getMass(2, stag)
            total_area += area

        if self._verbose:
            print(f"[GmshBuilder] get_sideset_area ss={ss_tag} "
                  f"-> {total_area}")
        return total_area

    def get_sideset_bounding_box(self, ss_tag):
        """Get bounding box of all surfaces in a sideset.

        Parameters
        ----------
        ss_tag : int
            Physical group tag (dimension 2).

        Returns
        -------
        dict
            Dictionary with keys ``'xmin'``, ``'ymin'``, ``'zmin'``,
            ``'xmax'``, ``'ymax'``, ``'zmax'``.
        """
        import gmsh
        self._check_initialized()

        surf_tags = gmsh.model.getEntitiesForPhysicalGroup(2, ss_tag)
        if len(surf_tags) == 0:
            raise ValueError(f"Sideset {ss_tag} has no surface entities")

        xmin_g, ymin_g, zmin_g = float('inf'), float('inf'), float('inf')
        xmax_g, ymax_g, zmax_g = float('-inf'), float('-inf'), float('-inf')

        for stag in surf_tags:
            xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(2, stag)
            xmin_g = min(xmin_g, xmin)
            ymin_g = min(ymin_g, ymin)
            zmin_g = min(zmin_g, zmin)
            xmax_g = max(xmax_g, xmax)
            ymax_g = max(ymax_g, ymax)
            zmax_g = max(zmax_g, zmax)

        result = {
            'xmin': xmin_g, 'ymin': ymin_g, 'zmin': zmin_g,
            'xmax': xmax_g, 'ymax': ymax_g, 'zmax': zmax_g,
        }

        if self._verbose:
            print(f"[GmshBuilder] get_sideset_bounding_box ss={ss_tag} "
                  f"-> {result}")
        return result

    # ------------------------------------------------------------------
    # Nodeset access (methods 19-22)
    # ------------------------------------------------------------------

    def get_nodeset_nodes(self, ns_tag):
        """Get node tags for a nodeset.

        Parameters
        ----------
        ns_tag : int
            Physical group tag (dimension 0).

        Returns
        -------
        numpy.ndarray
            1D array of node tags (int).
        """
        import numpy as np
        import gmsh
        self._check_initialized()

        point_tags = gmsh.model.getEntitiesForPhysicalGroup(0, ns_tag)
        all_node_tags = []
        for ptag in point_tags:
            node_tags, _coords, _params = gmsh.model.mesh.getNodes(0, ptag)
            all_node_tags.extend(node_tags)

        result = np.array(all_node_tags, dtype=np.int64)

        if self._verbose:
            print(f"[GmshBuilder] get_nodeset_nodes ns={ns_tag} "
                  f"-> {len(result)} nodes")
        return result

    def get_nodeset_node_count(self, ns_tag):
        """Count nodes in a nodeset.

        Parameters
        ----------
        ns_tag : int
            Physical group tag (dimension 0).

        Returns
        -------
        int
            Number of nodes.
        """
        import gmsh
        self._check_initialized()

        point_tags = gmsh.model.getEntitiesForPhysicalGroup(0, ns_tag)
        count = 0
        for ptag in point_tags:
            node_tags, _coords, _params = gmsh.model.mesh.getNodes(0, ptag)
            count += len(node_tags)

        if self._verbose:
            print(f"[GmshBuilder] get_nodeset_node_count ns={ns_tag} "
                  f"-> {count}")
        return count

    def get_nodeset_coordinates(self, ns_tag):
        """Get node coordinates for a nodeset.

        Parameters
        ----------
        ns_tag : int
            Physical group tag (dimension 0).

        Returns
        -------
        numpy.ndarray
            Node coordinates with shape ``(n, 3)``.
        """
        import numpy as np
        import gmsh
        self._check_initialized()

        point_tags = gmsh.model.getEntitiesForPhysicalGroup(0, ns_tag)
        all_coords = []
        for ptag in point_tags:
            _node_tags, coords, _params = gmsh.model.mesh.getNodes(0, ptag)
            all_coords.extend(coords)

        result = np.array(all_coords, dtype=np.float64).reshape(-1, 3)

        if self._verbose:
            print(f"[GmshBuilder] get_nodeset_coordinates ns={ns_tag} "
                  f"-> {len(result)} nodes")
        return result

    def get_nodeset_bounding_box(self, ns_tag):
        """Get bounding box of all nodes in a nodeset.

        Parameters
        ----------
        ns_tag : int
            Physical group tag (dimension 0).

        Returns
        -------
        dict
            Dictionary with keys ``'xmin'``, ``'ymin'``, ``'zmin'``,
            ``'xmax'``, ``'ymax'``, ``'zmax'``.
        """
        import numpy as np
        import gmsh
        self._check_initialized()

        coords = self.get_nodeset_coordinates(ns_tag)
        if len(coords) == 0:
            raise ValueError(f"Nodeset {ns_tag} has no nodes")

        result = {
            'xmin': float(np.min(coords[:, 0])),
            'ymin': float(np.min(coords[:, 1])),
            'zmin': float(np.min(coords[:, 2])),
            'xmax': float(np.max(coords[:, 0])),
            'ymax': float(np.max(coords[:, 1])),
            'zmax': float(np.max(coords[:, 2])),
        }

        if self._verbose:
            print(f"[GmshBuilder] get_nodeset_bounding_box ns={ns_tag} "
                  f"-> {result}")
        return result

    # ------------------------------------------------------------------
    # Material extensions (methods 23-26)
    # ------------------------------------------------------------------

    def set_material_conductivity(self, vol_id, sigma):
        """Set electrical conductivity for a volume's material.

        Stores ``sigma`` in ``self._materials[vol_id]['sigma']``.
        Creates the material entry if it does not already exist.

        Parameters
        ----------
        vol_id : int
            Volume ID.
        sigma : float
            Electrical conductivity in S/m.
        """
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        if vol_id not in self._materials:
            self._materials[vol_id] = {'name': '', 'mu_r': 1.0}

        self._materials[vol_id]['sigma'] = sigma

        if self._verbose:
            print(f"[GmshBuilder] set_material_conductivity vol={vol_id} "
                  f"sigma={sigma}")

    def set_material_density(self, vol_id, density):
        """Set mass density for a volume's material.

        Stores ``density`` in ``self._materials[vol_id]['density']``.
        Creates the material entry if it does not already exist.

        Parameters
        ----------
        vol_id : int
            Volume ID.
        density : float
            Mass density in kg/m^3.
        """
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        if vol_id not in self._materials:
            self._materials[vol_id] = {'name': '', 'mu_r': 1.0}

        self._materials[vol_id]['density'] = density

        if self._verbose:
            print(f"[GmshBuilder] set_material_density vol={vol_id} "
                  f"density={density}")

    def get_material_property(self, vol_id, key):
        """Get a specific material property for a volume.

        Parameters
        ----------
        vol_id : int
            Volume ID.
        key : str
            Property key (e.g., ``'mu_r'``, ``'sigma'``, ``'density'``).

        Returns
        -------
        value
            The property value.

        Raises
        ------
        KeyError
            If the volume has no material or the key is not found.
        """
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        if vol_id not in self._materials:
            raise KeyError(f"Volume {vol_id} has no material assigned")

        mat = self._materials[vol_id]
        if key not in mat:
            raise KeyError(f"Material for volume {vol_id} has no "
                           f"property '{key}'")

        value = mat[key]

        if self._verbose:
            print(f"[GmshBuilder] get_material_property vol={vol_id} "
                  f"key='{key}' -> {value}")
        return value

    def copy_material(self, from_vol, to_vol):
        """Copy material properties from one volume to another.

        Parameters
        ----------
        from_vol : int
            Source volume ID.
        to_vol : int
            Destination volume ID.

        Raises
        ------
        KeyError
            If either volume is not found or source has no material.
        """
        self._check_initialized()

        if from_vol not in self._volumes:
            raise KeyError(f"Source volume {from_vol} not found")
        if to_vol not in self._volumes:
            raise KeyError(f"Destination volume {to_vol} not found")
        if from_vol not in self._materials:
            raise KeyError(f"Source volume {from_vol} has no material assigned")

        self._materials[to_vol] = dict(self._materials[from_vol])

        if self._verbose:
            print(f"[GmshBuilder] copy_material from_vol={from_vol} "
                  f"to_vol={to_vol}")

    # ------------------------------------------------------------------
    # Block attributes (methods 27-30)
    # ------------------------------------------------------------------

    def set_block_attribute(self, block_tag, key, value):
        """Set a custom attribute on a block.

        Attempts to use ``gmsh.model.setAttribute(3, block_tag, key,
        [str(value)])``. If the GMSH API does not support this call,
        falls back to storing the attribute in an internal dictionary
        ``self._block_attributes``.

        Parameters
        ----------
        block_tag : int
            Physical group tag (dimension 3).
        key : str
            Attribute name.
        value
            Attribute value (will be converted to string).
        """
        import gmsh
        self._check_initialized()

        attr_name = f"block_{block_tag}_{key}"
        try:
            gmsh.model.setAttribute(attr_name, [str(value)])
        except Exception:
            pass
        if not hasattr(self, '_block_attributes'):
            self._block_attributes = {}
        if block_tag not in self._block_attributes:
            self._block_attributes[block_tag] = {}
        self._block_attributes[block_tag][key] = str(value)

        if self._verbose:
            print(f"[GmshBuilder] set_block_attribute block={block_tag} "
                  f"key='{key}' value='{value}'")

    def get_block_attribute(self, block_tag, key):
        """Get a custom attribute from a block.

        Attempts to use ``gmsh.model.getAttribute(3, block_tag, key)``.
        Falls back to the internal ``self._block_attributes`` dictionary
        if the GMSH API is unavailable.

        Parameters
        ----------
        block_tag : int
            Physical group tag (dimension 3).
        key : str
            Attribute name.

        Returns
        -------
        str
            Attribute value as a string.

        Raises
        ------
        KeyError
            If the attribute is not found.
        """
        import gmsh
        self._check_initialized()

        attr_name = f"block_{block_tag}_{key}"
        try:
            values = gmsh.model.getAttribute(attr_name)
            result = values[0] if values else ''
        except Exception:
            if not hasattr(self, '_block_attributes'):
                self._block_attributes = {}
            block_attrs = self._block_attributes.get(block_tag, {})
            if key not in block_attrs:
                raise KeyError(f"Block {block_tag} has no attribute '{key}'")
            result = block_attrs[key]

        if self._verbose:
            print(f"[GmshBuilder] get_block_attribute block={block_tag} "
                  f"key='{key}' -> '{result}'")
        return result

    def get_block_attribute_names(self, block_tag):
        """Get names of all custom attributes on a block.

        Attempts to use ``gmsh.model.getAttributeNames(3, block_tag)``.
        Falls back to the internal ``self._block_attributes`` dictionary
        if the GMSH API is unavailable.

        Parameters
        ----------
        block_tag : int
            Physical group tag (dimension 3).

        Returns
        -------
        list of str
            Attribute names.
        """
        import gmsh
        self._check_initialized()

        prefix = f"block_{block_tag}_"
        try:
            all_names = list(gmsh.model.getAttributeNames())
            names = [n[len(prefix):] for n in all_names if n.startswith(prefix)]
        except Exception:
            if not hasattr(self, '_block_attributes'):
                self._block_attributes = {}
            names = list(self._block_attributes.get(block_tag, {}).keys())

        if self._verbose:
            print(f"[GmshBuilder] get_block_attribute_names block={block_tag} "
                  f"-> {names}")
        return names

    def remove_block_attribute(self, block_tag, key):
        """Remove a custom attribute from a block.

        Attempts to use ``gmsh.model.removeAttribute(3, block_tag,
        key)``. Falls back to removing from the internal
        ``self._block_attributes`` dictionary if the GMSH API is
        unavailable.

        Parameters
        ----------
        block_tag : int
            Physical group tag (dimension 3).
        key : str
            Attribute name to remove.
        """
        import gmsh
        self._check_initialized()

        attr_name = f"block_{block_tag}_{key}"
        try:
            gmsh.model.removeAttribute(attr_name)
        except Exception:
            pass
        if not hasattr(self, '_block_attributes'):
            self._block_attributes = {}
        block_attrs = self._block_attributes.get(block_tag, {})
        block_attrs.pop(key, None)

        if self._verbose:
            print(f"[GmshBuilder] remove_block_attribute block={block_tag} "
                  f"key='{key}'")
