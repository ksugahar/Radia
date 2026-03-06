"""
Selection group management for GmshBuilder.

Gmsh equivalents: group create/delete/add/remove, set operations.
Groups are lightweight containers stored in-memory (self._groups) that
hold references to GMSH entities by dimension and tag. Curvesets are
GMSH physical groups of dimension 1 stored in the GMSH model.
"""


class GroupMixin:
    """Mixin providing selection group management and set operations."""

    # ------------------------------------------------------------------
    # Group lifecycle (10 methods)
    # ------------------------------------------------------------------

    def create_group(self, name):
        """Create a new selection group.

        Parameters
        ----------
        name : str
            Group name.

        Returns
        -------
        int
            New group ID.
        """
        self._check_initialized()

        gid = self._next_group_id
        self._groups[gid] = {
            'name': name,
            'entities': {0: [], 1: [], 2: [], 3: []},
        }
        self._next_group_id += 1

        if self._verbose:
            print(f"[GmshBuilder] create_group id={gid} name='{name}'")
        return gid

    def delete_group(self, group_id):
        """Delete a selection group.

        Parameters
        ----------
        group_id : int
            Group ID to delete.

        Raises
        ------
        KeyError
            If group_id does not exist.
        """
        self._check_initialized()

        if group_id not in self._groups:
            raise KeyError(f"Group {group_id} not found")

        del self._groups[group_id]

        if self._verbose:
            print(f"[GmshBuilder] delete_group id={group_id}")

    def rename_group(self, group_id, new_name):
        """Rename an existing group.

        Parameters
        ----------
        group_id : int
            Group ID.
        new_name : str
            New name for the group.

        Raises
        ------
        KeyError
            If group_id does not exist.
        """
        self._check_initialized()

        if group_id not in self._groups:
            raise KeyError(f"Group {group_id} not found")

        self._groups[group_id]['name'] = new_name

        if self._verbose:
            print(f"[GmshBuilder] rename_group id={group_id} "
                  f"new_name='{new_name}'")

    def get_group_name(self, group_id):
        """Return the name of a group.

        Parameters
        ----------
        group_id : int
            Group ID.

        Returns
        -------
        str
            Group name.

        Raises
        ------
        KeyError
            If group_id does not exist.
        """
        self._check_initialized()

        if group_id not in self._groups:
            raise KeyError(f"Group {group_id} not found")

        return self._groups[group_id]['name']

    def get_group_count(self):
        """Return total number of groups.

        Returns
        -------
        int
            Number of groups currently defined.
        """
        self._check_initialized()

        return len(self._groups)

    def get_all_groups(self):
        """Return summary info for all groups.

        Returns
        -------
        list of dict
            Each dict has keys 'id' (int), 'name' (str),
            'entity_count' (int).
        """
        self._check_initialized()

        result = []
        for gid in sorted(self._groups):
            grp = self._groups[gid]
            count = sum(len(tags) for tags in grp['entities'].values())
            result.append({
                'id': gid,
                'name': grp['name'],
                'entity_count': count,
            })
        return result

    def find_group_by_name(self, name):
        """Find a group by name (linear search).

        Parameters
        ----------
        name : str
            Group name to search for.

        Returns
        -------
        int or None
            Group ID if found, None otherwise.
        """
        self._check_initialized()

        for gid, grp in self._groups.items():
            if grp['name'] == name:
                return gid
        return None

    def clear_all_groups(self):
        """Remove all groups and reset the group ID counter."""
        self._check_initialized()

        self._groups.clear()
        self._next_group_id = 1

        if self._verbose:
            print("[GmshBuilder] clear_all_groups")

    def group_exists(self, group_id):
        """Check whether a group ID exists.

        Parameters
        ----------
        group_id : int
            Group ID to check.

        Returns
        -------
        bool
            True if the group exists.
        """
        self._check_initialized()

        return group_id in self._groups

    def get_group_ids(self):
        """Return sorted list of all group IDs.

        Returns
        -------
        list of int
            Sorted group IDs.
        """
        self._check_initialized()

        return sorted(self._groups.keys())

    # ------------------------------------------------------------------
    # Entity management (10 methods)
    # ------------------------------------------------------------------

    def _ensure_group(self, group_id):
        """Raise KeyError if group_id is not registered."""
        if group_id not in self._groups:
            raise KeyError(f"Group {group_id} not found")

    def group_add_volumes(self, group_id, vol_ids):
        """Add volume tags to a group.

        Parameters
        ----------
        group_id : int
            Group ID.
        vol_ids : list of int
            Volume tags to add.

        Raises
        ------
        KeyError
            If group_id does not exist.
        """
        self._check_initialized()
        self._ensure_group(group_id)

        entities = self._groups[group_id]['entities']
        if 3 not in entities:
            entities[3] = []
        existing = set(entities[3])
        for tag in vol_ids:
            if tag not in existing:
                entities[3].append(tag)
                existing.add(tag)

        if self._verbose:
            print(f"[GmshBuilder] group_add_volumes group={group_id} "
                  f"vol_ids={vol_ids}")

    def group_add_surfaces(self, group_id, surf_tags):
        """Add surface tags to a group.

        Parameters
        ----------
        group_id : int
            Group ID.
        surf_tags : list of int
            Surface tags to add.

        Raises
        ------
        KeyError
            If group_id does not exist.
        """
        self._check_initialized()
        self._ensure_group(group_id)

        entities = self._groups[group_id]['entities']
        if 2 not in entities:
            entities[2] = []
        existing = set(entities[2])
        for tag in surf_tags:
            if tag not in existing:
                entities[2].append(tag)
                existing.add(tag)

        if self._verbose:
            print(f"[GmshBuilder] group_add_surfaces group={group_id} "
                  f"surf_tags={surf_tags}")

    def group_add_curves(self, group_id, curve_tags):
        """Add curve tags to a group.

        Parameters
        ----------
        group_id : int
            Group ID.
        curve_tags : list of int
            Curve tags to add.

        Raises
        ------
        KeyError
            If group_id does not exist.
        """
        self._check_initialized()
        self._ensure_group(group_id)

        entities = self._groups[group_id]['entities']
        if 1 not in entities:
            entities[1] = []
        existing = set(entities[1])
        for tag in curve_tags:
            if tag not in existing:
                entities[1].append(tag)
                existing.add(tag)

        if self._verbose:
            print(f"[GmshBuilder] group_add_curves group={group_id} "
                  f"curve_tags={curve_tags}")

    def group_add_vertices(self, group_id, point_tags):
        """Add vertex (point) tags to a group.

        Parameters
        ----------
        group_id : int
            Group ID.
        point_tags : list of int
            Point tags to add.

        Raises
        ------
        KeyError
            If group_id does not exist.
        """
        self._check_initialized()
        self._ensure_group(group_id)

        entities = self._groups[group_id]['entities']
        if 0 not in entities:
            entities[0] = []
        existing = set(entities[0])
        for tag in point_tags:
            if tag not in existing:
                entities[0].append(tag)
                existing.add(tag)

        if self._verbose:
            print(f"[GmshBuilder] group_add_vertices group={group_id} "
                  f"point_tags={point_tags}")

    def group_remove_volumes(self, group_id, vol_ids):
        """Remove volume tags from a group.

        Parameters
        ----------
        group_id : int
            Group ID.
        vol_ids : list of int
            Volume tags to remove.

        Raises
        ------
        KeyError
            If group_id does not exist.
        """
        self._check_initialized()
        self._ensure_group(group_id)

        entities = self._groups[group_id]['entities']
        remove_set = set(vol_ids)
        entities[3] = [t for t in entities.get(3, []) if t not in remove_set]

        if self._verbose:
            print(f"[GmshBuilder] group_remove_volumes group={group_id} "
                  f"vol_ids={vol_ids}")

    def group_remove_surfaces(self, group_id, surf_tags):
        """Remove surface tags from a group.

        Parameters
        ----------
        group_id : int
            Group ID.
        surf_tags : list of int
            Surface tags to remove.

        Raises
        ------
        KeyError
            If group_id does not exist.
        """
        self._check_initialized()
        self._ensure_group(group_id)

        entities = self._groups[group_id]['entities']
        remove_set = set(surf_tags)
        entities[2] = [t for t in entities.get(2, []) if t not in remove_set]

        if self._verbose:
            print(f"[GmshBuilder] group_remove_surfaces group={group_id} "
                  f"surf_tags={surf_tags}")

    def group_remove_curves(self, group_id, curve_tags):
        """Remove curve tags from a group.

        Parameters
        ----------
        group_id : int
            Group ID.
        curve_tags : list of int
            Curve tags to remove.

        Raises
        ------
        KeyError
            If group_id does not exist.
        """
        self._check_initialized()
        self._ensure_group(group_id)

        entities = self._groups[group_id]['entities']
        remove_set = set(curve_tags)
        entities[1] = [t for t in entities.get(1, []) if t not in remove_set]

        if self._verbose:
            print(f"[GmshBuilder] group_remove_curves group={group_id} "
                  f"curve_tags={curve_tags}")

    def group_remove_vertices(self, group_id, point_tags):
        """Remove vertex (point) tags from a group.

        Parameters
        ----------
        group_id : int
            Group ID.
        point_tags : list of int
            Point tags to remove.

        Raises
        ------
        KeyError
            If group_id does not exist.
        """
        self._check_initialized()
        self._ensure_group(group_id)

        entities = self._groups[group_id]['entities']
        remove_set = set(point_tags)
        entities[0] = [t for t in entities.get(0, []) if t not in remove_set]

        if self._verbose:
            print(f"[GmshBuilder] group_remove_vertices group={group_id} "
                  f"point_tags={point_tags}")

    def group_get_entities(self, group_id, dim=None):
        """Get entities in a group.

        Parameters
        ----------
        group_id : int
            Group ID.
        dim : int, optional
            Dimension filter (0=vertex, 1=curve, 2=surface, 3=volume).
            If None, return all entities as (dim, tag) tuples.

        Returns
        -------
        list
            If dim is given, list of int (entity tags for that dimension).
            If dim is None, list of (dim, tag) tuples across all dimensions.

        Raises
        ------
        KeyError
            If group_id does not exist.
        """
        self._check_initialized()
        self._ensure_group(group_id)

        entities = self._groups[group_id]['entities']

        if dim is not None:
            return list(entities.get(dim, []))

        result = []
        for d in sorted(entities.keys()):
            for tag in entities[d]:
                result.append((d, tag))
        return result

    def group_entity_count(self, group_id, dim=None):
        """Count entities in a group.

        Parameters
        ----------
        group_id : int
            Group ID.
        dim : int, optional
            Dimension filter. If None, count all entities.

        Returns
        -------
        int
            Number of entities.

        Raises
        ------
        KeyError
            If group_id does not exist.
        """
        self._check_initialized()
        self._ensure_group(group_id)

        entities = self._groups[group_id]['entities']

        if dim is not None:
            return len(entities.get(dim, []))

        return sum(len(tags) for tags in entities.values())

    # ------------------------------------------------------------------
    # Curveset + set operations (10 methods)
    # ------------------------------------------------------------------

    def add_curveset(self, curve_tags, name=None):
        """Create a curveset (1D physical group).

        Parameters
        ----------
        curve_tags : list of int
            GMSH curve tags to include.
        name : str, optional
            Name for the curveset. If None, no name is assigned.

        Returns
        -------
        int
            GMSH physical group tag.
        """
        import gmsh
        self._check_initialized()

        pg_tag = gmsh.model.addPhysicalGroup(1, curve_tags)
        if name is not None:
            gmsh.model.setPhysicalName(1, pg_tag, name)

        if self._verbose:
            print(f"[GmshBuilder] add_curveset tag={pg_tag} "
                  f"curves={curve_tags} name='{name}'")
        return pg_tag

    def add_curveset_by_name(self, curve_tags, name):
        """Create a curveset with a required name.

        Parameters
        ----------
        curve_tags : list of int
            GMSH curve tags to include.
        name : str
            Name for the curveset.

        Returns
        -------
        int
            GMSH physical group tag.
        """
        import gmsh
        self._check_initialized()

        pg_tag = gmsh.model.addPhysicalGroup(1, curve_tags)
        gmsh.model.setPhysicalName(1, pg_tag, name)

        if self._verbose:
            print(f"[GmshBuilder] add_curveset_by_name tag={pg_tag} "
                  f"curves={curve_tags} name='{name}'")
        return pg_tag

    def get_all_curvesets(self):
        """Return info for all curvesets (1D physical groups).

        Returns
        -------
        list of dict
            Each dict has keys 'tag' (int), 'name' (str),
            'curves' (list of int).
        """
        import gmsh
        self._check_initialized()

        result = []
        for _dim, tag in gmsh.model.getPhysicalGroups(1):
            name = gmsh.model.getPhysicalName(1, tag)
            curves = list(gmsh.model.getEntitiesForPhysicalGroup(1, tag))
            result.append({
                'tag': tag,
                'name': name,
                'curves': curves,
            })
        return result

    def get_curveset_count(self):
        """Count number of curvesets (1D physical groups).

        Returns
        -------
        int
            Number of 1D physical groups.
        """
        import gmsh
        self._check_initialized()

        return len(gmsh.model.getPhysicalGroups(1))

    def get_curveset_curves(self, pg_tag):
        """Get curve tags in a curveset.

        Parameters
        ----------
        pg_tag : int
            Physical group tag of the curveset.

        Returns
        -------
        list of int
            Curve tags belonging to this curveset.
        """
        import gmsh
        self._check_initialized()

        return list(gmsh.model.getEntitiesForPhysicalGroup(1, pg_tag))

    def get_pg_for_entity(self, dim, tag):
        """Get physical group tags that contain a given entity.

        Parameters
        ----------
        dim : int
            Entity dimension (0-3).
        tag : int
            Entity tag.

        Returns
        -------
        list of int
            Physical group tags that include this entity.
        """
        import gmsh
        self._check_initialized()

        return list(gmsh.model.getPhysicalGroupsForEntity(dim, tag))

    def get_entities_by_name(self, name):
        """Find all entities belonging to physical groups with a given name.

        Iterates over all physical groups across all dimensions and
        collects entities whose physical group name matches.

        Parameters
        ----------
        name : str
            Physical group name to search for.

        Returns
        -------
        list of tuple
            List of (dim, tag) tuples for matching entities.
        """
        import gmsh
        self._check_initialized()

        result = []
        for dim, pg_tag in gmsh.model.getPhysicalGroups():
            pg_name = gmsh.model.getPhysicalName(dim, pg_tag)
            if pg_name == name:
                entities = gmsh.model.getEntitiesForPhysicalGroup(dim, pg_tag)
                for etag in entities:
                    result.append((dim, etag))
        return result

    def merge_groups(self, group_id1, group_id2):
        """Create a new group as the union of two groups.

        The original groups are kept unchanged.

        Parameters
        ----------
        group_id1 : int
            First group ID.
        group_id2 : int
            Second group ID.

        Returns
        -------
        int
            New group ID containing the union.

        Raises
        ------
        KeyError
            If either group does not exist.
        """
        self._check_initialized()
        self._ensure_group(group_id1)
        self._ensure_group(group_id2)

        name1 = self._groups[group_id1]['name']
        name2 = self._groups[group_id2]['name']
        new_gid = self.create_group(f"merge({name1}, {name2})")

        ent1 = self._groups[group_id1]['entities']
        ent2 = self._groups[group_id2]['entities']
        new_ent = self._groups[new_gid]['entities']

        all_dims = set(ent1.keys()) | set(ent2.keys())
        for d in all_dims:
            tags1 = ent1.get(d, [])
            tags2 = ent2.get(d, [])
            seen = set()
            merged = []
            for t in tags1 + tags2:
                if t not in seen:
                    merged.append(t)
                    seen.add(t)
            new_ent[d] = merged

        if self._verbose:
            print(f"[GmshBuilder] merge_groups {group_id1}+{group_id2} "
                  f"-> {new_gid}")
        return new_gid

    def subtract_groups(self, group_id1, group_id2):
        """Create a new group as the difference of two groups.

        The result contains entities in group_id1 that are NOT in
        group_id2. The original groups are kept unchanged.

        Parameters
        ----------
        group_id1 : int
            Base group ID.
        group_id2 : int
            Group ID whose entities are subtracted.

        Returns
        -------
        int
            New group ID containing the difference.

        Raises
        ------
        KeyError
            If either group does not exist.
        """
        self._check_initialized()
        self._ensure_group(group_id1)
        self._ensure_group(group_id2)

        name1 = self._groups[group_id1]['name']
        name2 = self._groups[group_id2]['name']
        new_gid = self.create_group(f"subtract({name1}, {name2})")

        ent1 = self._groups[group_id1]['entities']
        ent2 = self._groups[group_id2]['entities']
        new_ent = self._groups[new_gid]['entities']

        for d in ent1:
            remove_set = set(ent2.get(d, []))
            new_ent[d] = [t for t in ent1[d] if t not in remove_set]

        if self._verbose:
            print(f"[GmshBuilder] subtract_groups {group_id1}-{group_id2} "
                  f"-> {new_gid}")
        return new_gid

    def intersect_groups(self, group_id1, group_id2):
        """Create a new group as the intersection of two groups.

        The result contains only entities present in both groups.
        The original groups are kept unchanged.

        Parameters
        ----------
        group_id1 : int
            First group ID.
        group_id2 : int
            Second group ID.

        Returns
        -------
        int
            New group ID containing the intersection.

        Raises
        ------
        KeyError
            If either group does not exist.
        """
        self._check_initialized()
        self._ensure_group(group_id1)
        self._ensure_group(group_id2)

        name1 = self._groups[group_id1]['name']
        name2 = self._groups[group_id2]['name']
        new_gid = self.create_group(f"intersect({name1}, {name2})")

        ent1 = self._groups[group_id1]['entities']
        ent2 = self._groups[group_id2]['entities']
        new_ent = self._groups[new_gid]['entities']

        all_dims = set(ent1.keys()) | set(ent2.keys())
        for d in all_dims:
            set2 = set(ent2.get(d, []))
            new_ent[d] = [t for t in ent1.get(d, []) if t in set2]

        if self._verbose:
            print(f"[GmshBuilder] intersect_groups {group_id1}&{group_id2} "
                  f"-> {new_gid}")
        return new_gid
