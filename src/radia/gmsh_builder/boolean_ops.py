"""
Boolean operations for GmshBuilder.

Gmsh equivalents: unite, subtract, intersect, webcut, imprint, merge.
"""

from .utils import (get_bounding_box_union, create_half_space_box,
                     remove_excess_tools, dimtags_to_tags)


class BooleanMixin:
    """Mixin providing boolean and webcut operations."""

    def fuse(self, vol_ids):
        """Boolean union (Gmsh 'unite body').

        Parameters
        ----------
        vol_ids : list of int
            Volume IDs to fuse.

        Returns
        -------
        int
            New volume ID.
        """
        import gmsh
        self._check_initialized()

        if len(vol_ids) < 2:
            raise ValueError("fuse requires at least 2 volumes")

        all_tags = []
        for vid in vol_ids:
            if vid not in self._volumes:
                raise KeyError(f"Volume {vid} not found")
            all_tags.extend(self._volumes[vid])

        object_dimtags = [(3, all_tags[0])]
        tool_dimtags = [(3, t) for t in all_tags[1:]]

        out_dimtags, _ = gmsh.model.occ.fuse(object_dimtags, tool_dimtags)
        gmsh.model.occ.synchronize()

        new_tags = dimtags_to_tags(out_dimtags)
        for vid in vol_ids:
            self._remove_volume(vid)

        uid = self._register_volume(new_tags)
        if self._verbose:
            print(f"[GmshBuilder] fuse {vol_ids} -> id={uid}")
        return uid

    def cut(self, target_id, tool_id):
        """Boolean subtraction (Gmsh 'subtract').

        Parameters
        ----------
        target_id : int
            Volume to cut from.
        tool_id : int
            Volume to subtract.

        Returns
        -------
        int
            New volume ID.
        """
        import gmsh
        self._check_initialized()

        if target_id not in self._volumes:
            raise KeyError(f"Target volume {target_id} not found")
        if tool_id not in self._volumes:
            raise KeyError(f"Tool volume {tool_id} not found")

        object_dimtags = [(3, t) for t in self._volumes[target_id]]
        tool_dimtags = [(3, t) for t in self._volumes[tool_id]]

        out_dimtags, _ = gmsh.model.occ.cut(object_dimtags, tool_dimtags)
        gmsh.model.occ.synchronize()

        new_tags = dimtags_to_tags(out_dimtags)
        self._remove_volume(target_id)
        self._remove_volume(tool_id)

        uid = self._register_volume(new_tags)
        if self._verbose:
            print(f"[GmshBuilder] cut target={target_id} tool={tool_id} -> id={uid}")
        return uid

    def intersect(self, vol_id1, vol_id2):
        """Boolean intersection (Gmsh 'intersect').

        Parameters
        ----------
        vol_id1, vol_id2 : int
            Volume IDs.

        Returns
        -------
        int
            New volume ID (intersection region).
        """
        import gmsh
        self._check_initialized()

        if vol_id1 not in self._volumes:
            raise KeyError(f"Volume {vol_id1} not found")
        if vol_id2 not in self._volumes:
            raise KeyError(f"Volume {vol_id2} not found")

        object_dimtags = [(3, t) for t in self._volumes[vol_id1]]
        tool_dimtags = [(3, t) for t in self._volumes[vol_id2]]

        out_dimtags, _ = gmsh.model.occ.intersect(
            object_dimtags, tool_dimtags)
        gmsh.model.occ.synchronize()

        new_tags = dimtags_to_tags(out_dimtags)
        self._remove_volume(vol_id1)
        self._remove_volume(vol_id2)

        uid = self._register_volume(new_tags)
        if self._verbose:
            print(f"[GmshBuilder] intersect {vol_id1},{vol_id2} -> id={uid}")
        return uid

    def webcut_plane(self, vol_id, axis, position):
        """Split volume with a plane (Gmsh 'webcut with plane').

        Parameters
        ----------
        vol_id : int
            Volume ID to split.
        axis : str
            Cutting plane normal: 'x', 'y', or 'z'.
        position : float
            Position along axis.

        Returns
        -------
        list of int
            New volume IDs (typically 2).
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        axis = axis.lower()
        if axis not in ('x', 'y', 'z'):
            raise ValueError(f"axis must be 'x', 'y', or 'z', got '{axis}'")

        vol_tags = self._volumes[vol_id]
        bbox = get_bounding_box_union(vol_tags)
        span = max(bbox[3] - bbox[0], bbox[4] - bbox[1], bbox[5] - bbox[2])
        pad = span * 0.1

        tool = create_half_space_box(axis, position, bbox, pad)

        object_dimtags = [(3, t) for t in vol_tags]
        tool_dimtags = [(3, tool)]

        out_dimtags, out_map = gmsh.model.occ.fragment(
            object_dimtags, tool_dimtags)
        gmsh.model.occ.synchronize()

        target_set = set()
        for i in range(len(vol_tags)):
            for dt in out_map[i]:
                if dt[0] == 3:
                    target_set.add(dt[1])

        tool_set = set()
        for dt in out_map[len(vol_tags)]:
            if dt[0] == 3:
                tool_set.add(dt[1])

        tool_only = tool_set - target_set
        remove_excess_tools(tool_only)

        kept_tags = sorted(target_set)
        self._remove_volume(vol_id)

        new_ids = []
        for tag in kept_tags:
            uid = self._register_volume(tag)
            new_ids.append(uid)

        self._meshed = False
        if self._verbose:
            print(f"[GmshBuilder] webcut_plane vol={vol_id} axis={axis} "
                  f"pos={position} -> {new_ids}")
        return new_ids

    def webcut_cylinder(self, vol_id, center, axis, radius):
        """Split volume with a cylinder (Gmsh 'webcut with cylinder').

        Parameters
        ----------
        vol_id : int
            Volume to split.
        center : list
            Cylinder axis start point [x, y, z].
        axis : str
            Cylinder axis: 'x', 'y', or 'z'.
        radius : float
            Cylinder radius.

        Returns
        -------
        list of int
            New volume IDs.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        vol_tags = self._volumes[vol_id]
        bbox = get_bounding_box_union(vol_tags)
        span = max(bbox[3] - bbox[0], bbox[4] - bbox[1], bbox[5] - bbox[2])
        height = span * 2

        cx, cy, cz = center
        a = axis.lower()
        dx = dy = dz = 0
        if a == 'x':
            dx = height
            cx -= height / 2
        elif a == 'y':
            dy = height
            cy -= height / 2
        elif a == 'z':
            dz = height
            cz -= height / 2
        else:
            raise ValueError(f"axis must be 'x', 'y', or 'z'")

        cyl = gmsh.model.occ.addCylinder(cx, cy, cz, dx, dy, dz, radius)

        object_dimtags = [(3, t) for t in vol_tags]
        tool_dimtags = [(3, cyl)]

        out_dimtags, out_map = gmsh.model.occ.fragment(
            object_dimtags, tool_dimtags)
        gmsh.model.occ.synchronize()

        target_set = set()
        for i in range(len(vol_tags)):
            for dt in out_map[i]:
                if dt[0] == 3:
                    target_set.add(dt[1])

        tool_set = set()
        for dt in out_map[len(vol_tags)]:
            if dt[0] == 3:
                tool_set.add(dt[1])

        tool_only = tool_set - target_set
        remove_excess_tools(tool_only)

        kept_tags = sorted(target_set)
        self._remove_volume(vol_id)

        new_ids = []
        for tag in kept_tags:
            uid = self._register_volume(tag)
            new_ids.append(uid)

        self._meshed = False
        if self._verbose:
            print(f"[GmshBuilder] webcut_cylinder vol={vol_id} -> {new_ids}")
        return new_ids

    def webcut_sphere(self, vol_id, center, radius):
        """Split volume with a sphere (Gmsh 'webcut with sphere').

        Parameters
        ----------
        vol_id : int
            Volume to split.
        center : list
            Sphere center.
        radius : float
            Sphere radius.

        Returns
        -------
        list of int
            New volume IDs.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        vol_tags = self._volumes[vol_id]
        sph = gmsh.model.occ.addSphere(center[0], center[1], center[2],
                                         radius)

        object_dimtags = [(3, t) for t in vol_tags]
        tool_dimtags = [(3, sph)]

        out_dimtags, out_map = gmsh.model.occ.fragment(
            object_dimtags, tool_dimtags)
        gmsh.model.occ.synchronize()

        target_set = set()
        for i in range(len(vol_tags)):
            for dt in out_map[i]:
                if dt[0] == 3:
                    target_set.add(dt[1])

        tool_set = set()
        for dt in out_map[len(vol_tags)]:
            if dt[0] == 3:
                tool_set.add(dt[1])

        tool_only = tool_set - target_set
        remove_excess_tools(tool_only)

        kept_tags = sorted(target_set)
        self._remove_volume(vol_id)

        new_ids = []
        for tag in kept_tags:
            uid = self._register_volume(tag)
            new_ids.append(uid)

        self._meshed = False
        if self._verbose:
            print(f"[GmshBuilder] webcut_sphere vol={vol_id} -> {new_ids}")
        return new_ids

    def webcut_box(self, vol_id, center, size):
        """Split volume with a box (Gmsh 'webcut with brick').

        Parameters
        ----------
        vol_id : int
            Volume to split.
        center : list
            Box center.
        size : list
            Box dimensions [dx, dy, dz].

        Returns
        -------
        list of int
            New volume IDs.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        vol_tags = self._volumes[vol_id]
        cx, cy, cz = center
        dx, dy, dz = size
        box = gmsh.model.occ.addBox(cx - dx / 2, cy - dy / 2, cz - dz / 2,
                                      dx, dy, dz)

        object_dimtags = [(3, t) for t in vol_tags]
        tool_dimtags = [(3, box)]

        out_dimtags, out_map = gmsh.model.occ.fragment(
            object_dimtags, tool_dimtags)
        gmsh.model.occ.synchronize()

        target_set = set()
        for i in range(len(vol_tags)):
            for dt in out_map[i]:
                if dt[0] == 3:
                    target_set.add(dt[1])

        tool_set = set()
        for dt in out_map[len(vol_tags)]:
            if dt[0] == 3:
                tool_set.add(dt[1])

        tool_only = tool_set - target_set
        remove_excess_tools(tool_only)

        kept_tags = sorted(target_set)
        self._remove_volume(vol_id)

        new_ids = []
        for tag in kept_tags:
            uid = self._register_volume(tag)
            new_ids.append(uid)

        self._meshed = False
        if self._verbose:
            print(f"[GmshBuilder] webcut_box vol={vol_id} -> {new_ids}")
        return new_ids

    def imprint_and_merge(self):
        """Imprint and merge all volumes (Gmsh 'imprint all; merge all').

        Uses OCC fragment for shared interfaces.
        """
        import gmsh
        self._check_initialized()

        all_tags = []
        for tags in self._volumes.values():
            all_tags.extend(tags)

        if len(all_tags) < 2:
            return

        object_dimtags = [(3, all_tags[0])]
        tool_dimtags = [(3, t) for t in all_tags[1:]]

        out_dimtags, _ = gmsh.model.occ.fragment(
            object_dimtags, tool_dimtags)
        gmsh.model.occ.synchronize()

        new_tags = dimtags_to_tags(out_dimtags)
        self._volumes.clear()
        self._next_id = 1
        self._divisions = {}
        self._physical_groups.clear()
        self._materials.clear()
        self._names.clear()
        self._groups.clear()
        self._next_group_id = 1
        for tag in new_tags:
            self._register_volume(tag)

        self._meshed = False
        if self._verbose:
            print(f"[GmshBuilder] imprint_and_merge -> {len(new_tags)} volumes")

    def imprint_volumes(self, vol_ids):
        """Imprint specific volumes (Gmsh 'imprint vol').

        Parameters
        ----------
        vol_ids : list of int
            Volume IDs to imprint.
        """
        import gmsh
        self._check_initialized()

        all_tags = []
        for vid in vol_ids:
            if vid not in self._volumes:
                raise KeyError(f"Volume {vid} not found")
            all_tags.extend(self._volumes[vid])

        if len(all_tags) < 2:
            return

        object_dimtags = [(3, all_tags[0])]
        tool_dimtags = [(3, t) for t in all_tags[1:]]

        out_dimtags, out_map = gmsh.model.occ.fragment(
            object_dimtags, tool_dimtags)
        gmsh.model.occ.synchronize()

        # Remap volumes
        offset = 0
        for vid in vol_ids:
            tags = self._volumes[vid]
            n_tags = len(tags)
            new_tags = []
            for j in range(n_tags):
                if offset + j < len(out_map):
                    new_tags.extend([dt[1] for dt in out_map[offset + j] if dt[0] == 3])
            offset += n_tags
            self._volumes[vid] = new_tags if new_tags else tags

        self._meshed = False
        if self._verbose:
            print(f"[GmshBuilder] imprint_volumes {vol_ids}")

    def merge_volumes(self, vol_ids, tolerance=1e-6):
        """Merge volumes within tolerance (Gmsh 'merge vol tol').

        Same as imprint_volumes but with tolerance hint.
        """
        self.imprint_volumes(vol_ids)

    def section(self, vol_id, axis, position):
        """Section a volume: keep both halves (Gmsh 'section').

        Alias for webcut_plane.
        """
        return self.webcut_plane(vol_id, axis, position)

    def remove_volume(self, vol_id):
        """Delete a volume (Gmsh 'delete vol').

        Parameters
        ----------
        vol_id : int
            Volume ID to remove.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        for tag in self._volumes[vol_id]:
            try:
                gmsh.model.occ.remove([(3, tag)], recursive=True)
            except Exception:
                pass

        gmsh.model.occ.synchronize()
        self._remove_volume(vol_id)
        self._meshed = False

        if self._verbose:
            print(f"[GmshBuilder] remove_volume {vol_id}")

    def remove_all_duplicates(self):
        """Remove all duplicate entities (Gmsh 'merge all').

        Uses OCC removeAllDuplicates.
        """
        import gmsh
        self._check_initialized()

        gmsh.model.occ.removeAllDuplicates()
        gmsh.model.occ.synchronize()

        if self._verbose:
            print("[GmshBuilder] remove_all_duplicates")

    def fragment(self, vol_ids):
        """Boolean fragmentation preserving all pieces.

        Fragments all specified volumes against each other, keeping
        every resulting piece as a separate registered volume.

        Parameters
        ----------
        vol_ids : list of int
            Volume IDs to fragment.

        Returns
        -------
        list of int
            New volume IDs for all resulting pieces.
        """
        import gmsh
        self._check_initialized()

        if len(vol_ids) < 2:
            raise ValueError("fragment requires at least 2 volumes")

        all_tags = []
        for vid in vol_ids:
            if vid not in self._volumes:
                raise KeyError(f"Volume {vid} not found")
            all_tags.extend(self._volumes[vid])

        object_dimtags = [(3, all_tags[0])]
        tool_dimtags = [(3, t) for t in all_tags[1:]]

        out_dimtags, _ = gmsh.model.occ.fragment(
            object_dimtags, tool_dimtags)
        gmsh.model.occ.synchronize()

        new_tags = dimtags_to_tags(out_dimtags)
        for vid in vol_ids:
            self._remove_volume(vid)

        new_ids = []
        for tag in sorted(new_tags):
            uid = self._register_volume(tag)
            new_ids.append(uid)

        self._meshed = False
        if self._verbose:
            print(f"[GmshBuilder] fragment {vol_ids} -> {new_ids}")
        return new_ids

    def fragment_tracked(self, vol_ids):
        """Boolean fragmentation with origin tracking.

        Like fragment(), but returns a mapping from each input volume ID
        to the list of fragment volume IDs it produced.  Essential for
        volume classification after Boolean decomposition.

        Parameters
        ----------
        vol_ids : list of int
            Volume IDs to fragment.

        Returns
        -------
        dict of {int: list of int}
            Mapping from each input vol_id to resulting fragment vol_ids.
        """
        import gmsh
        self._check_initialized()

        if len(vol_ids) < 2:
            raise ValueError("fragment_tracked requires at least 2 volumes")

        # Collect all Gmsh tags and record which input index they belong to
        input_tags = []  # [(gmsh_tag, input_index), ...]
        for idx, vid in enumerate(vol_ids):
            if vid not in self._volumes:
                raise KeyError(f"Volume {vid} not found")
            for gtag in self._volumes[vid]:
                input_tags.append((gtag, idx))

        all_dimtags = [(3, gt) for gt, _ in input_tags]
        object_dimtags = all_dimtags[:1]
        tool_dimtags = all_dimtags[1:]

        _, result_map = gmsh.model.occ.fragment(object_dimtags, tool_dimtags)
        gmsh.model.occ.synchronize()

        # result_map[i] = list of (dim, tag) fragments from input i
        mapping = {vid: [] for vid in vol_ids}
        for i, (gtag, input_idx) in enumerate(input_tags):
            vid = vol_ids[input_idx]
            frag_tags = set(t[1] for t in result_map[i])
            mapping[vid].extend(frag_tags)

        # De-duplicate (a fragment can appear in multiple inputs)
        for vid in mapping:
            mapping[vid] = sorted(set(mapping[vid]))

        # Remove old volume IDs
        for vid in vol_ids:
            self._remove_volume(vid)

        # Register all unique fragments
        all_new_tags = sorted(set(
            tag for tags in mapping.values() for tag in tags))
        tag_to_uid = {}
        for tag in all_new_tags:
            uid = self._register_volume(tag)
            tag_to_uid[tag] = uid

        # Build result: input vol_id -> list of new vol_ids
        result = {}
        for vid in vol_ids:
            result[vid] = [tag_to_uid[t] for t in mapping[vid]]

        self._meshed = False
        if self._verbose:
            print(f"[GmshBuilder] fragment_tracked {vol_ids} -> {result}")
        return result

    def cut_keep_tool(self, target_id, tool_id):
        """Boolean cut keeping the tool volume.

        Subtracts tool from target but preserves the tool volume
        (removeObject=True, removeTool=False).

        Parameters
        ----------
        target_id : int
            Volume to cut from.
        tool_id : int
            Volume to subtract (kept after operation).

        Returns
        -------
        int
            New volume ID for the cut result.
        """
        import gmsh
        self._check_initialized()

        if target_id not in self._volumes:
            raise KeyError(f"Target volume {target_id} not found")
        if tool_id not in self._volumes:
            raise KeyError(f"Tool volume {tool_id} not found")

        object_dimtags = [(3, t) for t in self._volumes[target_id]]
        tool_dimtags = [(3, t) for t in self._volumes[tool_id]]

        out_dimtags, out_map = gmsh.model.occ.cut(
            object_dimtags, tool_dimtags,
            removeObject=True, removeTool=False)
        gmsh.model.occ.synchronize()

        new_tags = dimtags_to_tags(out_dimtags)
        self._remove_volume(target_id)

        # Tool is preserved by GMSH, but tags may have been renumbered
        # by OCC. Re-validate tool registration.
        remaining_tool_tags = []
        for t in self._volumes[tool_id]:
            try:
                gmsh.model.occ.getBoundingBox(3, t)
                remaining_tool_tags.append(t)
            except Exception:
                pass
        if remaining_tool_tags:
            self._volumes[tool_id] = remaining_tool_tags

        uid = self._register_volume(new_tags)
        self._meshed = False
        if self._verbose:
            print(f"[GmshBuilder] cut_keep_tool target={target_id} "
                  f"tool={tool_id} -> id={uid}")
        return uid

    def fuse_keep(self, vol_ids):
        """Fuse volumes but keep originals by copying first.

        Creates copies of all input volumes, fuses the copies,
        and returns the fused result. Original volumes are unchanged.

        Parameters
        ----------
        vol_ids : list of int
            Volume IDs to fuse (originals are preserved).

        Returns
        -------
        int
            New volume ID of the fused copy.
        """
        import gmsh
        self._check_initialized()

        if len(vol_ids) < 2:
            raise ValueError("fuse_keep requires at least 2 volumes")

        copy_tags = []
        for vid in vol_ids:
            if vid not in self._volumes:
                raise KeyError(f"Volume {vid} not found")
            for t in self._volumes[vid]:
                out = gmsh.model.occ.copy([(3, t)])
                copy_tags.extend(dimtags_to_tags(out))

        gmsh.model.occ.synchronize()

        object_dimtags = [(3, copy_tags[0])]
        tool_dimtags = [(3, t) for t in copy_tags[1:]]

        out_dimtags, _ = gmsh.model.occ.fuse(object_dimtags, tool_dimtags)
        gmsh.model.occ.synchronize()

        new_tags = dimtags_to_tags(out_dimtags)
        uid = self._register_volume(new_tags)
        self._meshed = False
        if self._verbose:
            print(f"[GmshBuilder] fuse_keep {vol_ids} -> id={uid} "
                  f"(originals preserved)")
        return uid

    def cut_with_plane(self, vol_id, point, normal):
        """Cut volume with an arbitrary plane.

        Creates a thin disk at the given point with the given normal
        direction, then uses fragment to split the volume.

        Parameters
        ----------
        vol_id : int
            Volume to cut.
        point : list
            Point on the cutting plane [x, y, z].
        normal : list
            Normal vector of the cutting plane [nx, ny, nz].

        Returns
        -------
        list of int
            New volume IDs (typically 2).
        """
        import gmsh
        import math
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        vol_tags = self._volumes[vol_id]
        bbox = get_bounding_box_union(vol_tags)
        span = max(bbox[3] - bbox[0], bbox[4] - bbox[1], bbox[5] - bbox[2])
        disk_radius = span * 2.0
        thickness = span * 0.001

        # Normalize normal vector
        nx, ny, nz = normal
        nmag = math.sqrt(nx * nx + ny * ny + nz * nz)
        if nmag < 1e-30:
            raise ValueError("Normal vector must be non-zero")
        nx, ny, nz = nx / nmag, ny / nmag, nz / nmag

        # Create a disk and extrude it thin along the normal
        px, py, pz = point
        disk = gmsh.model.occ.addDisk(px, py, pz, disk_radius, disk_radius)

        # Rotate disk so its normal aligns with the given normal.
        # Default disk normal is (0, 0, 1).
        # Compute rotation axis = cross(z, normal), angle = acos(dot(z, normal))
        dot_z = nz  # dot([0,0,1], [nx,ny,nz])
        if abs(dot_z - 1.0) > 1e-12 and abs(dot_z + 1.0) > 1e-12:
            # General case: rotate around cross product axis
            ax_x = -ny  # cross([0,0,1], [nx,ny,nz])
            ax_y = nx
            ax_z = 0.0
            angle = math.acos(max(-1.0, min(1.0, dot_z)))
            gmsh.model.occ.rotate([(2, disk)], px, py, pz,
                                  ax_x, ax_y, ax_z, angle)
        elif dot_z < 0:
            # Normal is (0, 0, -1): rotate 180 degrees around x
            gmsh.model.occ.rotate([(2, disk)], px, py, pz,
                                  1, 0, 0, math.pi)

        # Extrude disk to create a thin slab
        half_t = thickness / 2.0
        disk_copy = gmsh.model.occ.copy([(2, disk)])
        disk_copy_tag = disk_copy[0][1]
        extruded = gmsh.model.occ.extrude(
            [(2, disk)], nx * half_t, ny * half_t, nz * half_t)
        # Also extrude in the opposite direction
        extruded2 = gmsh.model.occ.extrude(
            [(2, disk_copy_tag)], -nx * half_t, -ny * half_t, -nz * half_t)

        # Collect tool volume tags from extrusion results
        tool_tags = []
        for dt in extruded + extruded2:
            if dt[0] == 3:
                tool_tags.append(dt[1])

        if not tool_tags:
            raise RuntimeError("Failed to create cutting plane tool")

        # Fuse tool volumes into one if multiple
        if len(tool_tags) > 1:
            obj_dt = [(3, tool_tags[0])]
            tl_dt = [(3, t) for t in tool_tags[1:]]
            fused_dt, _ = gmsh.model.occ.fuse(obj_dt, tl_dt)
            gmsh.model.occ.synchronize()
            tool_tags = dimtags_to_tags(fused_dt)

        object_dimtags = [(3, t) for t in vol_tags]
        tool_dimtags = [(3, t) for t in tool_tags]

        out_dimtags, out_map = gmsh.model.occ.fragment(
            object_dimtags, tool_dimtags)
        gmsh.model.occ.synchronize()

        target_set = set()
        for i in range(len(vol_tags)):
            for dt in out_map[i]:
                if dt[0] == 3:
                    target_set.add(dt[1])

        tool_set = set()
        for i in range(len(vol_tags), len(vol_tags) + len(tool_tags)):
            for dt in out_map[i]:
                if dt[0] == 3:
                    tool_set.add(dt[1])

        tool_only = tool_set - target_set
        remove_excess_tools(tool_only)

        kept_tags = sorted(target_set)
        self._remove_volume(vol_id)

        new_ids = []
        for tag in kept_tags:
            uid = self._register_volume(tag)
            new_ids.append(uid)

        self._meshed = False
        if self._verbose:
            print(f"[GmshBuilder] cut_with_plane vol={vol_id} "
                  f"point={point} normal={normal} -> {new_ids}")
        return new_ids

    def webcut_cone(self, vol_id, center, axis, r_base, r_top, height):
        """Split volume with a cone (Gmsh 'webcut with cone').

        Parameters
        ----------
        vol_id : int
            Volume to split.
        center : list
            Base center of the cone [x, y, z].
        axis : list
            Cone axis direction vector [ax, ay, az].
        r_base : float
            Radius at the base.
        r_top : float
            Radius at the top (0 for a pointed cone).
        height : float
            Height of the cone along the axis.

        Returns
        -------
        list of int
            New volume IDs.
        """
        import gmsh
        import math
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        vol_tags = self._volumes[vol_id]

        # Normalize axis
        ax, ay, az = axis
        amag = math.sqrt(ax * ax + ay * ay + az * az)
        if amag < 1e-30:
            raise ValueError("Axis vector must be non-zero")
        dx = ax / amag * height
        dy = ay / amag * height
        dz = az / amag * height

        cx, cy, cz = center
        cone = gmsh.model.occ.addCone(cx, cy, cz, dx, dy, dz,
                                       r_base, r_top)

        object_dimtags = [(3, t) for t in vol_tags]
        tool_dimtags = [(3, cone)]

        out_dimtags, out_map = gmsh.model.occ.fragment(
            object_dimtags, tool_dimtags)
        gmsh.model.occ.synchronize()

        target_set = set()
        for i in range(len(vol_tags)):
            for dt in out_map[i]:
                if dt[0] == 3:
                    target_set.add(dt[1])

        tool_set = set()
        for dt in out_map[len(vol_tags)]:
            if dt[0] == 3:
                tool_set.add(dt[1])

        tool_only = tool_set - target_set
        remove_excess_tools(tool_only)

        kept_tags = sorted(target_set)
        self._remove_volume(vol_id)

        new_ids = []
        for tag in kept_tags:
            uid = self._register_volume(tag)
            new_ids.append(uid)

        self._meshed = False
        if self._verbose:
            print(f"[GmshBuilder] webcut_cone vol={vol_id} -> {new_ids}")
        return new_ids

    def webcut_torus(self, vol_id, center, axis, R, r):
        """Split volume with a torus (Gmsh 'webcut with torus').

        Parameters
        ----------
        vol_id : int
            Volume to split.
        center : list
            Torus center [x, y, z].
        axis : list
            Torus axis direction vector [ax, ay, az].
        R : float
            Major radius (center of tube to center of torus).
        r : float
            Minor radius (tube radius).

        Returns
        -------
        list of int
            New volume IDs.
        """
        import gmsh
        import math
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        vol_tags = self._volumes[vol_id]
        cx, cy, cz = center

        # GMSH addTorus expects the torus centered at (cx,cy,cz) with
        # axis along z by default. We create it then rotate if needed.
        torus = gmsh.model.occ.addTorus(cx, cy, cz, R, r)

        # Rotate torus so its axis aligns with the given axis.
        # Default torus axis is (0, 0, 1).
        ax, ay, az = axis
        amag = math.sqrt(ax * ax + ay * ay + az * az)
        if amag < 1e-30:
            raise ValueError("Axis vector must be non-zero")
        ax, ay, az = ax / amag, ay / amag, az / amag

        dot_z = az
        if abs(dot_z - 1.0) > 1e-12 and abs(dot_z + 1.0) > 1e-12:
            rot_x = -ay
            rot_y = ax
            rot_z = 0.0
            angle = math.acos(max(-1.0, min(1.0, dot_z)))
            gmsh.model.occ.rotate([(3, torus)], cx, cy, cz,
                                  rot_x, rot_y, rot_z, angle)
        elif dot_z < 0:
            gmsh.model.occ.rotate([(3, torus)], cx, cy, cz,
                                  1, 0, 0, math.pi)

        object_dimtags = [(3, t) for t in vol_tags]
        tool_dimtags = [(3, torus)]

        out_dimtags, out_map = gmsh.model.occ.fragment(
            object_dimtags, tool_dimtags)
        gmsh.model.occ.synchronize()

        target_set = set()
        for i in range(len(vol_tags)):
            for dt in out_map[i]:
                if dt[0] == 3:
                    target_set.add(dt[1])

        tool_set = set()
        for dt in out_map[len(vol_tags)]:
            if dt[0] == 3:
                tool_set.add(dt[1])

        tool_only = tool_set - target_set
        remove_excess_tools(tool_only)

        kept_tags = sorted(target_set)
        self._remove_volume(vol_id)

        new_ids = []
        for tag in kept_tags:
            uid = self._register_volume(tag)
            new_ids.append(uid)

        self._meshed = False
        if self._verbose:
            print(f"[GmshBuilder] webcut_torus vol={vol_id} -> {new_ids}")
        return new_ids

    def subtract_list(self, target_id, tool_ids):
        """Subtract multiple tool volumes from one target.

        All tool volumes are removed after the operation.

        Parameters
        ----------
        target_id : int
            Volume to cut from.
        tool_ids : list of int
            Volumes to subtract.

        Returns
        -------
        int
            New volume ID.
        """
        import gmsh
        self._check_initialized()

        if target_id not in self._volumes:
            raise KeyError(f"Target volume {target_id} not found")

        tool_tags = []
        for tid in tool_ids:
            if tid not in self._volumes:
                raise KeyError(f"Tool volume {tid} not found")
            tool_tags.extend(self._volumes[tid])

        object_dimtags = [(3, t) for t in self._volumes[target_id]]
        tool_dimtags = [(3, t) for t in tool_tags]

        out_dimtags, _ = gmsh.model.occ.cut(object_dimtags, tool_dimtags)
        gmsh.model.occ.synchronize()

        new_tags = dimtags_to_tags(out_dimtags)
        self._remove_volume(target_id)
        for tid in tool_ids:
            self._remove_volume(tid)

        uid = self._register_volume(new_tags)
        self._meshed = False
        if self._verbose:
            print(f"[GmshBuilder] subtract_list target={target_id} "
                  f"tools={tool_ids} -> id={uid}")
        return uid

    def fuse_all(self):
        """Fuse all registered volumes into one.

        Returns
        -------
        int
            New volume ID.
        """
        import gmsh
        self._check_initialized()

        all_ids = list(self._volumes.keys())
        if len(all_ids) < 2:
            if len(all_ids) == 1:
                return all_ids[0]
            raise ValueError("No volumes to fuse")

        all_tags = []
        for vid in all_ids:
            all_tags.extend(self._volumes[vid])

        object_dimtags = [(3, all_tags[0])]
        tool_dimtags = [(3, t) for t in all_tags[1:]]

        out_dimtags, _ = gmsh.model.occ.fuse(object_dimtags, tool_dimtags)
        gmsh.model.occ.synchronize()

        new_tags = dimtags_to_tags(out_dimtags)
        for vid in all_ids:
            self._remove_volume(vid)

        uid = self._register_volume(new_tags)
        self._meshed = False
        if self._verbose:
            print(f"[GmshBuilder] fuse_all {all_ids} -> id={uid}")
        return uid

    def intersect_list(self, vol_ids):
        """Intersect multiple volumes pairwise.

        Performs sequential pairwise intersection of all given volumes,
        keeping only the region common to all.

        Parameters
        ----------
        vol_ids : list of int
            Volume IDs to intersect.

        Returns
        -------
        int
            New volume ID (common region).
        """
        import gmsh
        self._check_initialized()

        if len(vol_ids) < 2:
            raise ValueError("intersect_list requires at least 2 volumes")

        for vid in vol_ids:
            if vid not in self._volumes:
                raise KeyError(f"Volume {vid} not found")

        all_tags = []
        for vid in vol_ids:
            all_tags.extend(self._volumes[vid])

        object_dimtags = [(3, all_tags[0])]
        tool_dimtags = [(3, t) for t in all_tags[1:]]

        out_dimtags, _ = gmsh.model.occ.intersect(
            object_dimtags, tool_dimtags)
        gmsh.model.occ.synchronize()

        new_tags = dimtags_to_tags(out_dimtags)
        for vid in vol_ids:
            self._remove_volume(vid)

        uid = self._register_volume(new_tags)
        self._meshed = False
        if self._verbose:
            print(f"[GmshBuilder] intersect_list {vol_ids} -> id={uid}")
        return uid

    def remove_volumes(self, vol_ids):
        """Remove multiple volumes at once.

        Parameters
        ----------
        vol_ids : list of int
            Volume IDs to remove.
        """
        import gmsh
        self._check_initialized()

        for vid in vol_ids:
            if vid not in self._volumes:
                raise KeyError(f"Volume {vid} not found")

        for vid in vol_ids:
            for tag in self._volumes[vid]:
                try:
                    gmsh.model.occ.remove([(3, tag)], recursive=True)
                except Exception:
                    pass

        gmsh.model.occ.synchronize()
        for vid in vol_ids:
            self._remove_volume(vid)

        self._meshed = False
        if self._verbose:
            print(f"[GmshBuilder] remove_volumes {vol_ids}")

    def keep_only(self, vol_ids):
        """Remove all volumes except the specified ones.

        Parameters
        ----------
        vol_ids : list of int
            Volume IDs to keep.
        """
        import gmsh
        self._check_initialized()

        keep_set = set(vol_ids)
        for vid in keep_set:
            if vid not in self._volumes:
                raise KeyError(f"Volume {vid} not found")

        to_remove = [vid for vid in list(self._volumes.keys())
                     if vid not in keep_set]

        for vid in to_remove:
            for tag in self._volumes[vid]:
                try:
                    gmsh.model.occ.remove([(3, tag)], recursive=True)
                except Exception:
                    pass

        gmsh.model.occ.synchronize()
        for vid in to_remove:
            self._remove_volume(vid)

        self._meshed = False
        if self._verbose:
            print(f"[GmshBuilder] keep_only {vol_ids} "
                  f"(removed {len(to_remove)} volumes)")

    def split_by_volume(self, vol_id, tool_id):
        """Split vol_id by tool_id, keeping tool unchanged.

        Fragments vol_id with tool_id, then removes pieces that
        belong only to the tool (not overlapping with the target).
        The tool volume is preserved.

        Parameters
        ----------
        vol_id : int
            Volume to split.
        tool_id : int
            Splitting tool volume (preserved after operation).

        Returns
        -------
        list of int
            New volume IDs from the split target.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")
        if tool_id not in self._volumes:
            raise KeyError(f"Tool volume {tool_id} not found")

        vol_tags = self._volumes[vol_id]
        tool_tags = self._volumes[tool_id]

        # Copy tool so the original is preserved
        tool_copies = []
        for t in tool_tags:
            out = gmsh.model.occ.copy([(3, t)])
            tool_copies.extend(dimtags_to_tags(out))

        object_dimtags = [(3, t) for t in vol_tags]
        tool_dimtags = [(3, t) for t in tool_copies]

        out_dimtags, out_map = gmsh.model.occ.fragment(
            object_dimtags, tool_dimtags)
        gmsh.model.occ.synchronize()

        # Collect target pieces from out_map
        target_set = set()
        for i in range(len(vol_tags)):
            for dt in out_map[i]:
                if dt[0] == 3:
                    target_set.add(dt[1])

        # Collect tool-copy pieces from out_map
        tool_set = set()
        for i in range(len(vol_tags), len(vol_tags) + len(tool_copies)):
            for dt in out_map[i]:
                if dt[0] == 3:
                    tool_set.add(dt[1])

        # Remove pieces that belong only to the tool copy
        tool_only = tool_set - target_set
        remove_excess_tools(tool_only)

        kept_tags = sorted(target_set)
        self._remove_volume(vol_id)

        new_ids = []
        for tag in kept_tags:
            uid = self._register_volume(tag)
            new_ids.append(uid)

        self._meshed = False
        if self._verbose:
            print(f"[GmshBuilder] split_by_volume vol={vol_id} "
                  f"tool={tool_id} -> {new_ids}")
        return new_ids

    def chop(self, vol_id, tool_ids):
        """Chop a volume with multiple tools, preserving tools.

        Performs sequential cuts where each tool is preserved.
        Similar to subtract_list but tools remain registered.

        Parameters
        ----------
        vol_id : int
            Volume to chop.
        tool_ids : list of int
            Tool volumes (preserved after operation).

        Returns
        -------
        int
            New volume ID of the chopped result.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        tool_all_tags = []
        for tid in tool_ids:
            if tid not in self._volumes:
                raise KeyError(f"Tool volume {tid} not found")
            tool_all_tags.extend(self._volumes[tid])

        object_dimtags = [(3, t) for t in self._volumes[vol_id]]
        tool_dimtags = [(3, t) for t in tool_all_tags]

        out_dimtags, out_map = gmsh.model.occ.cut(
            object_dimtags, tool_dimtags,
            removeObject=True, removeTool=False)
        gmsh.model.occ.synchronize()

        new_tags = dimtags_to_tags(out_dimtags)
        self._remove_volume(vol_id)

        # Re-validate tool registrations (tags may have been renumbered)
        for tid in tool_ids:
            remaining = []
            for t in self._volumes[tid]:
                try:
                    gmsh.model.occ.getBoundingBox(3, t)
                    remaining.append(t)
                except Exception:
                    pass
            if remaining:
                self._volumes[tid] = remaining

        uid = self._register_volume(new_tags)
        self._meshed = False
        if self._verbose:
            print(f"[GmshBuilder] chop vol={vol_id} tools={tool_ids} "
                  f"-> id={uid} (tools preserved)")
        return uid

    def separate_volumes(self, vol_id):
        """Separate a multi-tag volume into individual volume IDs.

        If a registered volume ID maps to multiple GMSH volume tags,
        this method splits them into separate registered volumes.

        Parameters
        ----------
        vol_id : int
            Volume ID that may map to multiple GMSH tags.

        Returns
        -------
        list of int
            Individual volume IDs (one per GMSH tag). If vol_id
            already maps to a single tag, returns [vol_id].
        """
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        tags = self._volumes[vol_id]
        if len(tags) <= 1:
            if self._verbose:
                print(f"[GmshBuilder] separate_volumes vol={vol_id} "
                      f"-> already single tag")
            return [vol_id]

        self._remove_volume(vol_id)

        new_ids = []
        for tag in sorted(tags):
            uid = self._register_volume(tag)
            new_ids.append(uid)

        self._meshed = False
        if self._verbose:
            print(f"[GmshBuilder] separate_volumes vol={vol_id} "
                  f"-> {new_ids} ({len(tags)} tags)")
        return new_ids

    def webcut_general(self, vol_id, tool_dimtags):
        """General webcut with arbitrary GMSH dimtags as tool.

        Uses fragment + out_map + remove excess tools, following
        the same pattern as other webcut methods.

        Parameters
        ----------
        vol_id : int
            Volume to split.
        tool_dimtags : list of tuple
            GMSH (dim, tag) pairs to use as cutting tools.
            These are raw GMSH dimtags, not GmshBuilder volume IDs.

        Returns
        -------
        list of int
            New volume IDs.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        vol_tags = self._volumes[vol_id]

        object_dimtags = [(3, t) for t in vol_tags]

        out_dimtags, out_map = gmsh.model.occ.fragment(
            object_dimtags, list(tool_dimtags))
        gmsh.model.occ.synchronize()

        target_set = set()
        for i in range(len(vol_tags)):
            for dt in out_map[i]:
                if dt[0] == 3:
                    target_set.add(dt[1])

        tool_set = set()
        for i in range(len(vol_tags), len(vol_tags) + len(tool_dimtags)):
            for dt in out_map[i]:
                if dt[0] == 3:
                    tool_set.add(dt[1])

        tool_only = tool_set - target_set
        remove_excess_tools(tool_only)

        kept_tags = sorted(target_set)
        self._remove_volume(vol_id)

        new_ids = []
        for tag in kept_tags:
            uid = self._register_volume(tag)
            new_ids.append(uid)

        self._meshed = False
        if self._verbose:
            print(f"[GmshBuilder] webcut_general vol={vol_id} "
                  f"tools={tool_dimtags} -> {new_ids}")
        return new_ids
