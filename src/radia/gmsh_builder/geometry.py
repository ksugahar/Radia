"""
Geometry creation methods for GmshBuilder.

Gmsh equivalents: create brick, cylinder, sphere, cone, torus, wedge,
sweep, revolve, loft, fillet, chamfer, import step/iges/brep/stl.
"""

import os
import math


class GeometryMixin:
    """Mixin providing geometry creation methods."""

    def add_box(self, center, size):
        """Create a box (Gmsh 'create brick' equivalent).

        Parameters
        ----------
        center : list
            Center [cx, cy, cz] in meters.
        size : list
            Dimensions [Lx, Ly, Lz] in meters.

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        cx, cy, cz = center
        lx, ly, lz = size
        x0 = cx - lx / 2
        y0 = cy - ly / 2
        z0 = cz - lz / 2

        vol_tag = gmsh.model.occ.addBox(x0, y0, z0, lx, ly, lz)
        gmsh.model.occ.synchronize()

        uid = self._register_volume(vol_tag)
        if self._verbose:
            print(f"[GmshBuilder] add_box id={uid} center={center} size={size}")
        return uid

    def add_cylinder(self, center, axis, radius, height):
        """Create a cylinder (Gmsh 'create cylinder' equivalent).

        Parameters
        ----------
        center : list
            Center of base circle [x, y, z].
        axis : str
            Axis direction: 'x', 'y', or 'z'.
        radius : float
            Radius in meters.
        height : float
            Height in meters.

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        cx, cy, cz = center
        dx = dy = dz = 0
        a = axis.lower()
        if a == 'x':
            dx = height
        elif a == 'y':
            dy = height
        elif a == 'z':
            dz = height
        else:
            raise ValueError(f"axis must be 'x', 'y', or 'z', got '{axis}'")

        vol_tag = gmsh.model.occ.addCylinder(cx, cy, cz, dx, dy, dz, radius)
        gmsh.model.occ.synchronize()

        uid = self._register_volume(vol_tag)
        if self._verbose:
            print(f"[GmshBuilder] add_cylinder id={uid} axis={axis} "
                  f"r={radius} h={height}")
        return uid

    def add_sphere(self, center, radius):
        """Create a sphere (Gmsh 'create sphere' equivalent).

        Parameters
        ----------
        center : list
            Center [x, y, z].
        radius : float
            Radius in meters.

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        cx, cy, cz = center
        vol_tag = gmsh.model.occ.addSphere(cx, cy, cz, radius)
        gmsh.model.occ.synchronize()

        uid = self._register_volume(vol_tag)
        if self._verbose:
            print(f"[GmshBuilder] add_sphere id={uid} r={radius}")
        return uid

    def add_cone(self, center, axis, r_base, r_top, height):
        """Create a cone or truncated cone (Gmsh 'create cone' equivalent).

        Parameters
        ----------
        center : list
            Center of the base [x, y, z].
        axis : str
            Axis direction: 'x', 'y', or 'z'.
        r_base : float
            Base radius.
        r_top : float
            Top radius (0 for pointed cone).
        height : float
            Height.

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        cx, cy, cz = center
        dx = dy = dz = 0
        a = axis.lower()
        if a == 'x':
            dx = height
        elif a == 'y':
            dy = height
        elif a == 'z':
            dz = height
        else:
            raise ValueError(f"axis must be 'x', 'y', or 'z', got '{axis}'")

        vol_tag = gmsh.model.occ.addCone(
            cx, cy, cz, dx, dy, dz, r_base, r_top)
        gmsh.model.occ.synchronize()

        uid = self._register_volume(vol_tag)
        if self._verbose:
            print(f"[GmshBuilder] add_cone id={uid} r_base={r_base} "
                  f"r_top={r_top} h={height}")
        return uid

    def add_torus(self, center, axis, R, r):
        """Create a torus (Gmsh 'create torus' equivalent).

        Parameters
        ----------
        center : list
            Center [x, y, z].
        axis : str
            Axis direction: 'x', 'y', or 'z'.
        R : float
            Major radius (center to tube center).
        r : float
            Minor radius (tube radius).

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        cx, cy, cz = center
        a = axis.lower()

        # GMSH addTorus: z-axis is default, need rotation for x/y
        vol_tag = gmsh.model.occ.addTorus(cx, cy, cz, R, r)

        if a == 'x':
            gmsh.model.occ.rotate([(3, vol_tag)], cx, cy, cz, 0, 1, 0,
                                   math.pi / 2)
        elif a == 'y':
            gmsh.model.occ.rotate([(3, vol_tag)], cx, cy, cz, 1, 0, 0,
                                   math.pi / 2)
        elif a != 'z':
            raise ValueError(f"axis must be 'x', 'y', or 'z', got '{axis}'")

        gmsh.model.occ.synchronize()

        uid = self._register_volume(vol_tag)
        if self._verbose:
            print(f"[GmshBuilder] add_torus id={uid} R={R} r={r}")
        return uid

    def add_wedge(self, center, size, ltx=0.0):
        """Create a wedge (Gmsh 'create wedge' equivalent).

        Parameters
        ----------
        center : list
            Center [x, y, z].
        size : list
            [dx, dy, dz] full dimensions.
        ltx : float
            Top x-length (0 for triangular wedge).

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        cx, cy, cz = center
        dx, dy, dz = size
        x0 = cx - dx / 2
        y0 = cy - dy / 2
        z0 = cz - dz / 2

        vol_tag = gmsh.model.occ.addWedge(x0, y0, z0, dx, dy, dz, ltx=ltx)
        gmsh.model.occ.synchronize()

        uid = self._register_volume(vol_tag)
        if self._verbose:
            print(f"[GmshBuilder] add_wedge id={uid} size={size} ltx={ltx}")
        return uid

    def extrude(self, vol_id, dx, dy, dz, num_elements=None, heights=None,
                recombine=False):
        """Extrude surfaces of a volume (Gmsh 'sweep direction').

        Creates new volume by extruding the surfaces of the given volume.

        Parameters
        ----------
        vol_id : int
            Volume ID whose surfaces to extrude.
        dx, dy, dz : float
            Extrusion direction vector.
        num_elements : list of int or int, optional
            Number of elements per layer along the extrusion direction
            for structured meshing.  If a single int, converted to [int].
        heights : list of float, optional
            Normalized cumulative heights of each layer (values in
            [0, 1]).  Must have the same length as num_elements.
        recombine : bool, optional
            If True, recombine extruded mesh into hexahedra/quads.

        Returns
        -------
        list of int
            New volume IDs from extrusion.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        kwargs = {}
        if num_elements is not None:
            if isinstance(num_elements, int):
                num_elements = [num_elements]
            kwargs['numElements'] = num_elements
        if heights is not None:
            kwargs['heights'] = heights
        if recombine:
            kwargs['recombine'] = True

        # Get all boundary surfaces of the volume (GMSH extrude does
        # not support 3D entities directly)
        surf_dimtags = []
        for t in self._volumes[vol_id]:
            boundary = gmsh.model.getBoundary([(3, t)], oriented=False)
            surf_dimtags.extend([dt for dt in boundary if dt[0] == 2])
        out = gmsh.model.occ.extrude(surf_dimtags, dx, dy, dz, **kwargs)
        gmsh.model.occ.synchronize()

        new_tags = [dt[1] for dt in out if dt[0] == 3]
        new_ids = []
        for tag in new_tags:
            uid = self._register_volume(tag)
            new_ids.append(uid)

        if self._verbose:
            print(f"[GmshBuilder] extrude vol={vol_id} d=[{dx},{dy},{dz}] "
                  f"-> {new_ids}")
        return new_ids

    def revolve(self, vol_id, center, axis, angle, num_elements=None,
                heights=None, recombine=False):
        """Revolve surfaces of a volume (Gmsh 'sweep rotate').

        Parameters
        ----------
        vol_id : int
            Volume ID whose surfaces to revolve.
        center : list
            Rotation center [x, y, z].
        axis : list
            Rotation axis [ax, ay, az].
        angle : float
            Rotation angle in radians.
        num_elements : list of int or int, optional
            Number of elements per layer along the revolution for
            structured meshing.  If a single int, converted to [int].
        heights : list of float, optional
            Normalized cumulative heights of each layer (values in
            [0, 1]).  Must have the same length as num_elements.
        recombine : bool, optional
            If True, recombine revolved mesh into hexahedra/quads.

        Returns
        -------
        list of int
            New volume IDs.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        kwargs = {}
        if num_elements is not None:
            if isinstance(num_elements, int):
                num_elements = [num_elements]
            kwargs['numElements'] = num_elements
        if heights is not None:
            kwargs['heights'] = heights
        if recombine:
            kwargs['recombine'] = True

        # Get all boundary surfaces of the volume (GMSH revolve does
        # not support 3D entities directly)
        surf_dimtags = []
        for t in self._volumes[vol_id]:
            boundary = gmsh.model.getBoundary([(3, t)], oriented=False)
            surf_dimtags.extend([dt for dt in boundary if dt[0] == 2])
        out = gmsh.model.occ.revolve(
            surf_dimtags,
            center[0], center[1], center[2],
            axis[0], axis[1], axis[2], angle, **kwargs)
        gmsh.model.occ.synchronize()

        new_tags = [dt[1] for dt in out if dt[0] == 3]
        new_ids = []
        for tag in new_tags:
            uid = self._register_volume(tag)
            new_ids.append(uid)

        if self._verbose:
            print(f"[GmshBuilder] revolve vol={vol_id} "
                  f"angle={math.degrees(angle):.1f} deg -> {new_ids}")
        return new_ids

    def loft(self, wire_tags):
        """Create volume by lofting through sections (Gmsh 'create vol loft').

        Parameters
        ----------
        wire_tags : list of int
            GMSH wire/curve loop tags defining cross-sections.

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        out = gmsh.model.occ.addThruSections(wire_tags)
        gmsh.model.occ.synchronize()

        vol_tags = [dt[1] for dt in out if dt[0] == 3]
        if not vol_tags:
            raise RuntimeError("Loft produced no volumes")

        uid = self._register_volume(vol_tags)
        if self._verbose:
            print(f"[GmshBuilder] loft -> id={uid}")
        return uid

    # ---- Low-level geometry primitives ----

    def add_point(self, x, y, z, mesh_size=0.0):
        """Create a point (Gmsh 'create vertex').

        Returns GMSH point tag (not a volume ID).
        """
        import gmsh
        self._check_initialized()

        tag = gmsh.model.occ.addPoint(x, y, z, mesh_size)
        gmsh.model.occ.synchronize()
        return tag

    def add_line(self, p1, p2):
        """Create a line between two points. Returns GMSH curve tag."""
        import gmsh
        self._check_initialized()

        tag = gmsh.model.occ.addLine(p1, p2)
        gmsh.model.occ.synchronize()
        return tag

    def add_circle_arc(self, start, center, end):
        """Create a circle arc. Returns GMSH curve tag."""
        import gmsh
        self._check_initialized()

        tag = gmsh.model.occ.addCircleArc(start, center, end)
        gmsh.model.occ.synchronize()
        return tag

    def add_spline(self, point_tags):
        """Create a spline through points. Returns GMSH curve tag."""
        import gmsh
        self._check_initialized()

        tag = gmsh.model.occ.addSpline(point_tags)
        gmsh.model.occ.synchronize()
        return tag

    def add_bspline(self, point_tags):
        """Create a B-spline through points. Returns GMSH curve tag."""
        import gmsh
        self._check_initialized()

        tag = gmsh.model.occ.addBSpline(point_tags)
        gmsh.model.occ.synchronize()
        return tag

    def add_wire(self, curve_tags):
        """Create a wire (curve loop). Returns GMSH wire tag."""
        import gmsh
        self._check_initialized()

        tag = gmsh.model.occ.addWire(curve_tags)
        gmsh.model.occ.synchronize()
        return tag

    def add_rectangle(self, center, dx, dy):
        """Create a rectangular surface (Gmsh 'create surface rectangle').

        Returns GMSH surface tag.
        """
        import gmsh
        self._check_initialized()

        cx, cy, cz = center
        tag = gmsh.model.occ.addRectangle(cx - dx / 2, cy - dy / 2, cz,
                                            dx, dy)
        gmsh.model.occ.synchronize()
        return tag

    def add_disk(self, center, rx, ry=None):
        """Create a disk surface (Gmsh 'create surface circle').

        Returns GMSH surface tag.
        """
        import gmsh
        self._check_initialized()

        cx, cy, cz = center
        if ry is None:
            ry = rx
        tag = gmsh.model.occ.addDisk(cx, cy, cz, rx, ry)
        gmsh.model.occ.synchronize()
        return tag

    def add_plane_surface(self, wire_tags):
        """Create a plane surface from wire(s). Returns GMSH surface tag."""
        import gmsh
        self._check_initialized()

        tag = gmsh.model.occ.addPlaneSurface(wire_tags)
        gmsh.model.occ.synchronize()
        return tag

    def add_surface_filling(self, wire_tag):
        """Create a surface filling a wire boundary. Returns GMSH surface tag."""
        import gmsh
        self._check_initialized()

        tag = gmsh.model.occ.addSurfaceFilling(wire_tag)
        gmsh.model.occ.synchronize()
        return tag

    def add_surface_loop(self, surface_tags):
        """Create a surface loop (shell). Returns GMSH shell tag."""
        import gmsh
        self._check_initialized()

        tag = gmsh.model.occ.addSurfaceLoop(surface_tags)
        gmsh.model.occ.synchronize()
        return tag

    def add_volume(self, shell_tags):
        """Create a volume from shell(s). Returns volume ID."""
        import gmsh
        self._check_initialized()

        vol_tag = gmsh.model.occ.addVolume(shell_tags)
        gmsh.model.occ.synchronize()

        uid = self._register_volume(vol_tag)
        if self._verbose:
            print(f"[GmshBuilder] add_volume id={uid}")
        return uid

    def add_pipe(self, wire_tag, radius):
        """Create a pipe (tube) along a wire path.

        Creates a disk at the start of the wire and sweeps it along
        the wire path using addPipe.

        Parameters
        ----------
        wire_tag : int
            GMSH wire tag defining the pipe path.
        radius : float
            Pipe radius.

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        # Get the start point of the wire to place the disk
        bounds = gmsh.model.getParametrizationBounds(1, wire_tag)
        start_pt = gmsh.model.getValue(1, wire_tag, [bounds[0][0]])
        # Create a disk perpendicular to wire at start
        disk = gmsh.model.occ.addDisk(
            start_pt[0], start_pt[1], start_pt[2], radius, radius)
        # Sweep disk along wire
        out = gmsh.model.occ.addPipe([(2, disk)], wire_tag)
        gmsh.model.occ.synchronize()

        vol_tags = [dt[1] for dt in out if dt[0] == 3]
        if not vol_tags:
            raise RuntimeError("add_pipe produced no volumes")

        uid = self._register_volume(vol_tags)
        if self._verbose:
            print(f"[GmshBuilder] add_pipe: wire={wire_tag}, "
                  f"radius={radius} -> id={uid}")
        return uid

    # ---- Geometry modification ----

    def fillet(self, vol_id, curve_ids, radius):
        """Apply fillet to edges (Gmsh 'modify vol fillet').

        Parameters
        ----------
        vol_id : int
            Volume ID.
        curve_ids : list of int
            GMSH curve tags to fillet.
        radius : float or list of float
            Fillet radius for each curve.

        Returns
        -------
        int
            New volume ID (original is replaced).
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        vol_tags = self._volumes[vol_id]
        if isinstance(radius, (int, float)):
            radii = [radius] * len(curve_ids)
        else:
            radii = list(radius)

        out = gmsh.model.occ.fillet(vol_tags, curve_ids, radii)
        gmsh.model.occ.synchronize()

        new_tags = [dt[1] for dt in out if dt[0] == 3]
        self._remove_volume(vol_id)
        uid = self._register_volume(new_tags)

        if self._verbose:
            print(f"[GmshBuilder] fillet vol={vol_id} -> id={uid} r={radius}")
        return uid

    def chamfer(self, vol_id, curve_ids, distances):
        """Apply chamfer to edges (Gmsh 'modify vol chamfer').

        Parameters
        ----------
        vol_id : int
            Volume ID.
        curve_ids : list of int
            GMSH curve tags to chamfer.
        distances : float or list of float
            Chamfer distance(s).

        Returns
        -------
        int
            New volume ID.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        vol_tags = self._volumes[vol_id]
        if isinstance(distances, (int, float)):
            dists = [distances] * len(curve_ids)
        else:
            dists = list(distances)

        # GMSH chamfer needs surface tags adjacent to each curve
        # Get surfaces for each curve
        # getAdjacencies returns (upward, downward); for dim=1 (curve),
        # upward[0] gives adjacent surfaces
        surf_tags = []
        for ct in curve_ids:
            # Find surfaces that contain this curve
            up = gmsh.model.getAdjacencies(1, ct)
            if len(up[0]) > 0:
                surf_tags.append(up[0][0])
            else:
                surf_tags.append(0)

        out = gmsh.model.occ.chamfer(vol_tags, curve_ids, surf_tags, dists)
        gmsh.model.occ.synchronize()

        new_tags = [dt[1] for dt in out if dt[0] == 3]
        self._remove_volume(vol_id)
        uid = self._register_volume(new_tags)

        if self._verbose:
            print(f"[GmshBuilder] chamfer vol={vol_id} -> id={uid}")
        return uid

    def heal_geometry(self, vol_id=None, tolerance=1e-8):
        """Heal geometry (Gmsh 'healer autoheal').

        Parameters
        ----------
        vol_id : int, optional
            Volume to heal. If None, heals all.
        tolerance : float
            Healing tolerance.
        """
        import gmsh
        self._check_initialized()

        if vol_id is not None:
            if vol_id not in self._volumes:
                raise KeyError(f"Volume {vol_id} not found")
            dimtags = [(3, t) for t in self._volumes[vol_id]]
        else:
            dimtags = []
            for tags in self._volumes.values():
                dimtags.extend([(3, t) for t in tags])

        outDimTags = gmsh.model.occ.healShapes(
            dimtags, tolerance=tolerance)
        gmsh.model.occ.synchronize()

        # Update _volumes if healShapes changed any volume tags
        new_vol_tags = {dt[1] for dt in outDimTags if dt[0] == 3}
        if new_vol_tags:
            old_vol_tags = {dt[1] for dt in dimtags if dt[0] == 3}
            if new_vol_tags != old_vol_tags:
                for uid, tags in list(self._volumes.items()):
                    updated = []
                    changed = False
                    for t in tags:
                        if t in old_vol_tags and t not in new_vol_tags:
                            changed = True
                        else:
                            updated.append(t)
                    # Add any new tags not previously tracked
                    for nt in new_vol_tags:
                        if nt not in old_vol_tags and nt not in updated:
                            updated.append(nt)
                            changed = True
                    if changed:
                        self._volumes[uid] = updated

        self._meshed = False

        if self._verbose:
            print(f"[GmshBuilder] heal_geometry vol={vol_id}")

    def remove_small_edges(self, vol_id=None, tolerance=1e-6):
        """Remove small edges (small edge removal).

        Parameters
        ----------
        vol_id : int, optional
            Volume to clean. If None, cleans all.
        tolerance : float
            Minimum edge length to keep.
        """
        import gmsh
        self._check_initialized()

        if vol_id is not None:
            if vol_id not in self._volumes:
                raise KeyError(f"Volume {vol_id} not found")
            dimtags = [(3, t) for t in self._volumes[vol_id]]
        else:
            dimtags = []
            for tags in self._volumes.values():
                dimtags.extend([(3, t) for t in tags])

        outDimTags = gmsh.model.occ.healShapes(
            dimtags, tolerance=tolerance,
            fixSmallEdges=True, fixSmallFaces=True)
        gmsh.model.occ.synchronize()

        # Update _volumes if healShapes changed any volume tags
        new_vol_tags = {dt[1] for dt in outDimTags if dt[0] == 3}
        if new_vol_tags:
            old_vol_tags = {dt[1] for dt in dimtags if dt[0] == 3}
            if new_vol_tags != old_vol_tags:
                for uid, tags in list(self._volumes.items()):
                    updated = []
                    changed = False
                    for t in tags:
                        if t in old_vol_tags and t not in new_vol_tags:
                            changed = True
                        else:
                            updated.append(t)
                    # Add any new tags not previously tracked
                    for nt in new_vol_tags:
                        if nt not in old_vol_tags and nt not in updated:
                            updated.append(nt)
                            changed = True
                    if changed:
                        self._volumes[uid] = updated

        self._meshed = False

        if self._verbose:
            print(f"[GmshBuilder] remove_small_edges tol={tolerance}")

    # ---- Import ----

    def import_step(self, filename):
        """Import STEP file geometry (Gmsh 'import step').

        Parameters
        ----------
        filename : str
            Path to STEP file.

        Returns
        -------
        list of int
            Volume IDs for imported volumes.
        """
        import gmsh
        self._check_initialized()

        if not os.path.exists(filename):
            raise FileNotFoundError(f"STEP file not found: {filename}")

        out_dimtags = gmsh.model.occ.importShapes(filename)
        gmsh.model.occ.synchronize()

        vol_tags = [dt[1] for dt in out_dimtags if dt[0] == 3]
        new_ids = []
        for tag in vol_tags:
            uid = self._register_volume(tag)
            new_ids.append(uid)

        if self._verbose:
            print(f"[GmshBuilder] import_step '{filename}' -> "
                  f"{len(vol_tags)} vols, ids={new_ids}")
        return new_ids

    def import_iges(self, filename):
        """Import IGES file geometry (Gmsh 'import iges').

        Parameters
        ----------
        filename : str
            Path to IGES file.

        Returns
        -------
        list of int
            Volume IDs.
        """
        return self.import_step(filename)  # Same OCC import

    def import_brep(self, filename):
        """Import BREP file geometry.

        Parameters
        ----------
        filename : str
            Path to BREP file.

        Returns
        -------
        list of int
            Volume IDs.
        """
        return self.import_step(filename)  # Same OCC import

    def import_stl(self, filename):
        """Import STL file (Gmsh 'import stl').

        Parameters
        ----------
        filename : str
            Path to STL file.

        Returns
        -------
        list of int
            Volume IDs (may be empty for surface-only STL).
        """
        import gmsh
        self._check_initialized()

        if not os.path.exists(filename):
            raise FileNotFoundError(f"STL file not found: {filename}")

        # Snapshot existing volumes before merge
        existing = set(t for d, t in gmsh.model.getEntities(3))

        gmsh.merge(filename)

        # Find only newly added volumes
        all_vols = set(t for d, t in gmsh.model.getEntities(3))
        new_vols = list(all_vols - existing)
        new_ids = []
        for tag in new_vols:
            uid = self._register_volume(tag)
            new_ids.append(uid)

        if self._verbose:
            print(f"[GmshBuilder] import_stl '{filename}' -> {len(new_ids)} vols")
        return new_ids

    # ---- Extended geometry primitives ----

    def add_cylinder_vector(self, center, direction, radius, height):
        """Create a cylinder with direction specified as a vector.

        Unlike add_cylinder which takes an axis string ('x','y','z'),
        this method accepts an arbitrary direction vector.

        Parameters
        ----------
        center : list
            Base center [x, y, z] in meters.
        direction : list
            Direction vector [dx, dy, dz] (will be normalized and scaled
            by height).
        radius : float
            Radius in meters.
        height : float
            Height in meters.

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        cx, cy, cz = center
        dx, dy, dz = direction
        mag = math.sqrt(dx * dx + dy * dy + dz * dz)
        if mag == 0:
            raise ValueError("direction vector must be non-zero")
        scale = height / mag
        dx *= scale
        dy *= scale
        dz *= scale

        vol_tag = gmsh.model.occ.addCylinder(cx, cy, cz, dx, dy, dz, radius)
        gmsh.model.occ.synchronize()

        uid = self._register_volume(vol_tag)
        if self._verbose:
            print(f"[GmshBuilder] add_cylinder_vector id={uid} "
                  f"dir={direction} r={radius} h={height}")
        return uid

    def add_cone_vector(self, center, direction, r_base, r_top, height):
        """Create a cone with direction specified as a vector.

        Unlike add_cone which takes an axis string ('x','y','z'),
        this method accepts an arbitrary direction vector.

        Parameters
        ----------
        center : list
            Base center [x, y, z] in meters.
        direction : list
            Direction vector [dx, dy, dz] (will be normalized and scaled
            by height).
        r_base : float
            Base radius in meters.
        r_top : float
            Top radius in meters (0 for pointed cone).
        height : float
            Height in meters.

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        cx, cy, cz = center
        dx, dy, dz = direction
        mag = math.sqrt(dx * dx + dy * dy + dz * dz)
        if mag == 0:
            raise ValueError("direction vector must be non-zero")
        scale = height / mag
        dx *= scale
        dy *= scale
        dz *= scale

        vol_tag = gmsh.model.occ.addCone(
            cx, cy, cz, dx, dy, dz, r_base, r_top)
        gmsh.model.occ.synchronize()

        uid = self._register_volume(vol_tag)
        if self._verbose:
            print(f"[GmshBuilder] add_cone_vector id={uid} "
                  f"r_base={r_base} r_top={r_top} h={height}")
        return uid

    def add_half_space(self, axis, position):
        """Create a large box representing a half-space for boolean operations.

        The half-space extends 10 meters in each direction from the
        splitting plane. Useful as a tool body for boolean cut operations.

        Parameters
        ----------
        axis : str
            Normal axis of the splitting plane: 'x', 'y', or 'z'.
        position : float
            Position along the axis where the half-space starts.

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        half = 10.0  # 10 meters half-extent
        a = axis.lower()
        if a == 'x':
            vol_tag = gmsh.model.occ.addBox(
                position, -half, -half, half, 2 * half, 2 * half)
        elif a == 'y':
            vol_tag = gmsh.model.occ.addBox(
                -half, position, -half, 2 * half, half, 2 * half)
        elif a == 'z':
            vol_tag = gmsh.model.occ.addBox(
                -half, -half, position, 2 * half, 2 * half, half)
        else:
            raise ValueError(f"axis must be 'x', 'y', or 'z', got '{axis}'")
        gmsh.model.occ.synchronize()

        uid = self._register_volume(vol_tag)
        if self._verbose:
            print(f"[GmshBuilder] add_half_space id={uid} "
                  f"axis={axis} pos={position}")
        return uid

    def add_ellipsoid(self, center, rx, ry, rz):
        """Create an ellipsoid by scaling a sphere non-uniformly.

        Parameters
        ----------
        center : list
            Center [x, y, z] in meters.
        rx : float
            Semi-axis length in x direction.
        ry : float
            Semi-axis length in y direction.
        rz : float
            Semi-axis length in z direction.

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        cx, cy, cz = center
        # Create unit sphere at origin, then dilate and translate
        vol_tag = gmsh.model.occ.addSphere(0, 0, 0, 1.0)
        gmsh.model.occ.dilate([(3, vol_tag)], 0, 0, 0, rx, ry, rz)
        gmsh.model.occ.translate([(3, vol_tag)], cx, cy, cz)
        gmsh.model.occ.synchronize()

        uid = self._register_volume(vol_tag)
        if self._verbose:
            print(f"[GmshBuilder] add_ellipsoid id={uid} "
                  f"rx={rx} ry={ry} rz={rz}")
        return uid

    def add_circular_arc(self, center, axis, radius, start_angle, end_angle):
        """Create a circular arc curve segment.

        Parameters
        ----------
        center : list
            Center point [x, y, z] in meters.
        axis : str
            Normal axis of the arc plane: 'x', 'y', or 'z'.
        radius : float
            Arc radius in meters.
        start_angle : float
            Start angle in radians.
        end_angle : float
            End angle in radians.

        Returns
        -------
        int
            GMSH curve tag.
        """
        import gmsh
        self._check_initialized()

        cx, cy, cz = center
        a = axis.lower()

        if a == 'z':
            x1 = cx + radius * math.cos(start_angle)
            y1 = cy + radius * math.sin(start_angle)
            z1 = cz
            x2 = cx + radius * math.cos(end_angle)
            y2 = cy + radius * math.sin(end_angle)
            z2 = cz
        elif a == 'x':
            x1 = cx
            y1 = cy + radius * math.cos(start_angle)
            z1 = cz + radius * math.sin(start_angle)
            x2 = cx
            y2 = cy + radius * math.cos(end_angle)
            z2 = cz + radius * math.sin(end_angle)
        elif a == 'y':
            x1 = cx + radius * math.cos(start_angle)
            y1 = cy
            z1 = cz + radius * math.sin(start_angle)
            x2 = cx + radius * math.cos(end_angle)
            y2 = cy
            z2 = cz + radius * math.sin(end_angle)
        else:
            raise ValueError(f"axis must be 'x', 'y', or 'z', got '{axis}'")

        p_start = gmsh.model.occ.addPoint(x1, y1, z1)
        p_center = gmsh.model.occ.addPoint(cx, cy, cz)
        p_end = gmsh.model.occ.addPoint(x2, y2, z2)
        tag = gmsh.model.occ.addCircleArc(p_start, p_center, p_end)
        gmsh.model.occ.synchronize()

        if self._verbose:
            print(f"[GmshBuilder] add_circular_arc tag={tag} r={radius} "
                  f"angles=[{start_angle}, {end_angle}]")
        return tag

    def add_hollow_cylinder(self, center, axis, r_inner, r_outer, height):
        """Create a hollow cylinder (tube) by subtracting inner from outer.

        Parameters
        ----------
        center : list
            Base center [x, y, z] in meters.
        axis : str
            Axis direction: 'x', 'y', or 'z'.
        r_inner : float
            Inner radius in meters.
        r_outer : float
            Outer radius in meters.
        height : float
            Height in meters.

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        if r_inner >= r_outer:
            raise ValueError("r_inner must be less than r_outer")

        cx, cy, cz = center
        dx = dy = dz = 0
        a = axis.lower()
        if a == 'x':
            dx = height
        elif a == 'y':
            dy = height
        elif a == 'z':
            dz = height
        else:
            raise ValueError(f"axis must be 'x', 'y', or 'z', got '{axis}'")

        outer_tag = gmsh.model.occ.addCylinder(
            cx, cy, cz, dx, dy, dz, r_outer)
        inner_tag = gmsh.model.occ.addCylinder(
            cx, cy, cz, dx, dy, dz, r_inner)

        out, _ = gmsh.model.occ.cut(
            [(3, outer_tag)], [(3, inner_tag)])
        gmsh.model.occ.synchronize()

        vol_tags = [dt[1] for dt in out if dt[0] == 3]
        if not vol_tags:
            raise RuntimeError("Hollow cylinder boolean cut produced no volumes")

        uid = self._register_volume(vol_tags)
        if self._verbose:
            print(f"[GmshBuilder] add_hollow_cylinder id={uid} "
                  f"r_inner={r_inner} r_outer={r_outer} h={height}")
        return uid

    def add_hollow_sphere(self, center, r_inner, r_outer):
        """Create a hollow sphere (spherical shell).

        Parameters
        ----------
        center : list
            Center [x, y, z] in meters.
        r_inner : float
            Inner radius in meters.
        r_outer : float
            Outer radius in meters.

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        if r_inner >= r_outer:
            raise ValueError("r_inner must be less than r_outer")

        cx, cy, cz = center
        outer_tag = gmsh.model.occ.addSphere(cx, cy, cz, r_outer)
        inner_tag = gmsh.model.occ.addSphere(cx, cy, cz, r_inner)

        out, _ = gmsh.model.occ.cut(
            [(3, outer_tag)], [(3, inner_tag)])
        gmsh.model.occ.synchronize()

        vol_tags = [dt[1] for dt in out if dt[0] == 3]
        if not vol_tags:
            raise RuntimeError(
                "Hollow sphere boolean cut produced no volumes")

        uid = self._register_volume(vol_tags)
        if self._verbose:
            print(f"[GmshBuilder] add_hollow_sphere id={uid} "
                  f"r_inner={r_inner} r_outer={r_outer}")
        return uid

    def add_c_shape(self, center, axis, outer_r, inner_r, height, gap_angle):
        """Create a C-shaped magnet (hollow cylinder with angular gap).

        The gap is centered at angle=0 (positive x for z-axis).

        Parameters
        ----------
        center : list
            Base center [x, y, z] in meters.
        axis : str
            Axis direction: 'x', 'y', or 'z'.
        outer_r : float
            Outer radius in meters.
        inner_r : float
            Inner radius in meters.
        height : float
            Height in meters.
        gap_angle : float
            Gap opening angle in radians (centered at angle=0).

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        if inner_r >= outer_r:
            raise ValueError("inner_r must be less than outer_r")
        if gap_angle <= 0 or gap_angle >= 2 * math.pi:
            raise ValueError("gap_angle must be between 0 and 2*pi")

        cx, cy, cz = center
        dx = dy = dz = 0
        a = axis.lower()
        if a == 'x':
            dx = height
        elif a == 'y':
            dy = height
        elif a == 'z':
            dz = height
        else:
            raise ValueError(f"axis must be 'x', 'y', or 'z', got '{axis}'")

        # Create hollow cylinder
        outer_tag = gmsh.model.occ.addCylinder(
            cx, cy, cz, dx, dy, dz, outer_r)
        inner_tag = gmsh.model.occ.addCylinder(
            cx, cy, cz, dx, dy, dz, inner_r)
        out_hollow, _ = gmsh.model.occ.cut(
            [(3, outer_tag)], [(3, inner_tag)])

        hollow_tags = [dt[1] for dt in out_hollow if dt[0] == 3]
        if not hollow_tags:
            raise RuntimeError("Hollow cylinder creation failed")

        # Create a wedge box for the gap
        # The gap is a box extending beyond outer_r, centered at angle=0
        gap_extent = outer_r * 1.5
        half_gap = gap_angle / 2

        if a == 'z':
            # Create two cutting planes to form the gap wedge
            # Use a box rotated to cover the gap
            gap_box = gmsh.model.occ.addBox(
                cx, cy - gap_extent, cz - height * 0.1,
                gap_extent, 2 * gap_extent, height * 1.2)
            # Rotate the gap box by -half_gap and +half_gap, intersect
            # Simpler: create a wedge using two rotated boxes
            box1 = gmsh.model.occ.addBox(
                cx, cy, cz - height * 0.1,
                gap_extent, gap_extent, height * 1.2)
            gmsh.model.occ.rotate(
                [(3, box1)], cx, cy, cz, 0, 0, 1, -half_gap)
            box2 = gmsh.model.occ.addBox(
                cx, cy - gap_extent, cz - height * 0.1,
                gap_extent, gap_extent, height * 1.2)
            gmsh.model.occ.rotate(
                [(3, box2)], cx, cy, cz, 0, 0, 1, half_gap)
            # Fuse the two boxes to form the gap region
            gap_fused, _ = gmsh.model.occ.fuse(
                [(3, box1)], [(3, box2)])
            gmsh.model.occ.remove([(3, gap_box)], recursive=True)
        elif a == 'x':
            box1 = gmsh.model.occ.addBox(
                cx - height * 0.1, cy, cz,
                height * 1.2, gap_extent, gap_extent)
            gmsh.model.occ.rotate(
                [(3, box1)], cx, cy, cz, 1, 0, 0, -half_gap)
            box2 = gmsh.model.occ.addBox(
                cx - height * 0.1, cy - gap_extent, cz,
                height * 1.2, gap_extent, gap_extent)
            gmsh.model.occ.rotate(
                [(3, box2)], cx, cy, cz, 1, 0, 0, half_gap)
            gap_fused, _ = gmsh.model.occ.fuse(
                [(3, box1)], [(3, box2)])
        else:  # y
            box1 = gmsh.model.occ.addBox(
                cx, cy - height * 0.1, cz,
                gap_extent, height * 1.2, gap_extent)
            gmsh.model.occ.rotate(
                [(3, box1)], cx, cy, cz, 0, 1, 0, -half_gap)
            box2 = gmsh.model.occ.addBox(
                cx - gap_extent, cy - height * 0.1, cz,
                gap_extent, height * 1.2, gap_extent)
            gmsh.model.occ.rotate(
                [(3, box2)], cx, cy, cz, 0, 1, 0, half_gap)
            gap_fused, _ = gmsh.model.occ.fuse(
                [(3, box1)], [(3, box2)])

        # Cut the gap from the hollow cylinder
        gap_tags = [dt[1] for dt in gap_fused if dt[0] == 3]
        out_c, _ = gmsh.model.occ.cut(
            [(3, t) for t in hollow_tags],
            [(3, t) for t in gap_tags])
        gmsh.model.occ.synchronize()

        c_tags = [dt[1] for dt in out_c if dt[0] == 3]
        if not c_tags:
            raise RuntimeError("C-shape boolean cut produced no volumes")

        uid = self._register_volume(c_tags)
        if self._verbose:
            print(f"[GmshBuilder] add_c_shape id={uid} "
                  f"outer_r={outer_r} inner_r={inner_r} "
                  f"gap={math.degrees(gap_angle):.1f} deg")
        return uid

    def add_box_at(self, x0, y0, z0, dx, dy, dz):
        """Create a box from a corner point and dimensions.

        Unlike add_box which takes center and size, this method starts
        from the minimum corner (x0, y0, z0).

        Parameters
        ----------
        x0, y0, z0 : float
            Minimum corner coordinates in meters.
        dx, dy, dz : float
            Dimensions in x, y, z directions.

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        vol_tag = gmsh.model.occ.addBox(x0, y0, z0, dx, dy, dz)
        gmsh.model.occ.synchronize()

        uid = self._register_volume(vol_tag)
        if self._verbose:
            print(f"[GmshBuilder] add_box_at id={uid} "
                  f"corner=[{x0},{y0},{z0}] size=[{dx},{dy},{dz}]")
        return uid

    def add_cylinder_between(self, p1, p2, radius):
        """Create a cylinder between two arbitrary points.

        The cylinder axis runs from p1 to p2. The height is the
        distance between the two points.

        Parameters
        ----------
        p1 : list
            Start point [x, y, z] in meters.
        p2 : list
            End point [x, y, z] in meters.
        radius : float
            Cylinder radius in meters.

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        x1, y1, z1 = p1
        x2, y2, z2 = p2
        dx = x2 - x1
        dy = y2 - y1
        dz = z2 - z1
        length = math.sqrt(dx * dx + dy * dy + dz * dz)
        if length == 0:
            raise ValueError("p1 and p2 must be different points")

        vol_tag = gmsh.model.occ.addCylinder(x1, y1, z1, dx, dy, dz, radius)
        gmsh.model.occ.synchronize()

        uid = self._register_volume(vol_tag)
        if self._verbose:
            print(f"[GmshBuilder] add_cylinder_between id={uid} "
                  f"p1={p1} p2={p2} r={radius} len={length:.4f}")
        return uid

    def extrude_surface(self, surf_tag, dx, dy, dz,
                        num_elements=None, heights=None, recombine=False):
        """Extrude a surface to create a volume.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag to extrude.
        dx, dy, dz : float
            Extrusion direction and distance.
        num_elements : list of int or int, optional
            Number of elements per layer along the extrusion direction
            for structured meshing. If a single int, converted to [int].
            Example: [5] for 5 uniform layers, [3, 2] for two layers
            with 3 and 2 elements respectively.
        heights : list of float, optional
            Normalized cumulative heights of each layer (values in
            [0, 1]).  Must have the same length as num_elements.
            Example: [0.3, 1.0] for a first layer occupying 30% and a
            second layer occupying 70% of the total extrusion distance.
            If None, layers are equally spaced.
        recombine : bool, optional
            If True, recombine extruded mesh into hexahedra/quads.
            Default is False.

        Returns
        -------
        list of int
            New volume IDs from extrusion.
        """
        import gmsh
        self._check_initialized()

        kwargs = {}
        if num_elements is not None:
            if isinstance(num_elements, int):
                num_elements = [num_elements]
            kwargs['numElements'] = num_elements
        if heights is not None:
            kwargs['heights'] = heights
        if recombine:
            kwargs['recombine'] = True

        out = gmsh.model.occ.extrude([(2, surf_tag)], dx, dy, dz, **kwargs)
        gmsh.model.occ.synchronize()

        new_tags = [dt[1] for dt in out if dt[0] == 3]
        new_ids = []
        for tag in new_tags:
            uid = self._register_volume(tag)
            new_ids.append(uid)

        if self._verbose:
            ne_str = f" numElements={num_elements}" if num_elements else ""
            print(f"[GmshBuilder] extrude_surface surf={surf_tag} "
                  f"d=[{dx},{dy},{dz}]{ne_str} -> {new_ids}")
        return new_ids

    def revolve_surface(self, surf_tag, center, axis, angle,
                        num_elements=None, heights=None, recombine=False):
        """Revolve a surface around an axis to create a volume.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag to revolve.
        center : list
            Rotation center [x, y, z].
        axis : list
            Rotation axis vector [ax, ay, az].
        angle : float
            Rotation angle in radians.
        num_elements : list of int or int, optional
            Number of elements per layer along the revolution for
            structured meshing.  If a single int, converted to [int].
        heights : list of float, optional
            Normalized cumulative heights of each layer (values in
            [0, 1]).  Must have the same length as num_elements.
        recombine : bool, optional
            If True, recombine revolved mesh into hexahedra/quads.
            Default is False.

        Returns
        -------
        list of int
            New volume IDs.
        """
        import gmsh
        self._check_initialized()

        kwargs = {}
        if num_elements is not None:
            if isinstance(num_elements, int):
                num_elements = [num_elements]
            kwargs['numElements'] = num_elements
        if heights is not None:
            kwargs['heights'] = heights
        if recombine:
            kwargs['recombine'] = True

        out = gmsh.model.occ.revolve(
            [(2, surf_tag)],
            center[0], center[1], center[2],
            axis[0], axis[1], axis[2], angle, **kwargs)
        gmsh.model.occ.synchronize()

        new_tags = [dt[1] for dt in out if dt[0] == 3]
        new_ids = []
        for tag in new_tags:
            uid = self._register_volume(tag)
            new_ids.append(uid)

        if self._verbose:
            ne_str = f" numElements={num_elements}" if num_elements else ""
            print(f"[GmshBuilder] revolve_surface surf={surf_tag} "
                  f"angle={math.degrees(angle):.1f} deg{ne_str} -> {new_ids}")
        return new_ids

    def add_polygon(self, points, z=0):
        """Create a polygon surface from 2D points.

        Points are specified as 2D coordinates in the xy-plane at the
        given z height.

        Parameters
        ----------
        points : list of list
            List of [x, y] coordinates for polygon vertices.
        z : float
            Z coordinate for the polygon plane (default 0).

        Returns
        -------
        int
            GMSH surface tag.
        """
        import gmsh
        self._check_initialized()

        if len(points) < 3:
            raise ValueError("Polygon requires at least 3 points")

        # Create points
        pt_tags = []
        for pt in points:
            tag = gmsh.model.occ.addPoint(pt[0], pt[1], z)
            pt_tags.append(tag)

        # Create lines between consecutive points
        line_tags = []
        n = len(pt_tags)
        for i in range(n):
            tag = gmsh.model.occ.addLine(pt_tags[i], pt_tags[(i + 1) % n])
            line_tags.append(tag)

        # Create wire and plane surface
        wire_tag = gmsh.model.occ.addWire(line_tags)
        surf_tag = gmsh.model.occ.addPlaneSurface([wire_tag])
        gmsh.model.occ.synchronize()

        if self._verbose:
            print(f"[GmshBuilder] add_polygon surf={surf_tag} "
                  f"n_pts={len(points)} z={z}")
        return surf_tag

    def add_prism(self, base_points, height, axis='z'):
        """Create a prism by extruding a polygon.

        Parameters
        ----------
        base_points : list of list
            List of [x, y] coordinates for the base polygon vertices.
        height : float
            Extrusion height in meters.
        axis : str
            Extrusion axis: 'x', 'y', or 'z' (default 'z').

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        # Create the base polygon surface
        surf_tag = self.add_polygon(base_points, z=0)

        a = axis.lower()
        if a == 'x':
            dx, dy, dz = height, 0, 0
        elif a == 'y':
            dx, dy, dz = 0, height, 0
        elif a == 'z':
            dx, dy, dz = 0, 0, height
        else:
            raise ValueError(f"axis must be 'x', 'y', or 'z', got '{axis}'")

        out = gmsh.model.occ.extrude([(2, surf_tag)], dx, dy, dz)
        gmsh.model.occ.synchronize()

        vol_tags = [dt[1] for dt in out if dt[0] == 3]
        if not vol_tags:
            raise RuntimeError("Prism extrusion produced no volumes")

        uid = self._register_volume(vol_tags)
        if self._verbose:
            print(f"[GmshBuilder] add_prism id={uid} "
                  f"n_pts={len(base_points)} h={height} axis={axis}")
        return uid

    def copy_geometry(self, vol_id):
        """Copy a volume (convenience alias for copy in transforms).

        Parameters
        ----------
        vol_id : int
            Volume ID to copy.

        Returns
        -------
        int
            New volume ID.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        vol_tags = self._volumes[vol_id]
        dimtags = [(3, t) for t in vol_tags]
        out = gmsh.model.occ.copy(dimtags)
        gmsh.model.occ.synchronize()

        new_tags = [dt[1] for dt in out if dt[0] == 3]
        if not new_tags:
            raise RuntimeError("Copy produced no volumes")

        uid = self._register_volume(new_tags)
        if self._verbose:
            print(f"[GmshBuilder] copy_geometry vol={vol_id} -> id={uid}")
        return uid

    def add_circle(self, center, axis, radius):
        """Create a full circle curve.

        Parameters
        ----------
        center : list
            Center [x, y, z] in meters.
        axis : str
            Normal axis of the circle plane: 'x', 'y', or 'z'.
        radius : float
            Circle radius in meters.

        Returns
        -------
        int
            GMSH curve tag.
        """
        import gmsh
        self._check_initialized()

        cx, cy, cz = center
        a = axis.lower()

        # GMSH addCircle creates a circle in the xy-plane by default
        tag = gmsh.model.occ.addCircle(cx, cy, cz, radius)

        if a == 'x':
            gmsh.model.occ.rotate([(1, tag)], cx, cy, cz, 0, 1, 0,
                                   math.pi / 2)
        elif a == 'y':
            gmsh.model.occ.rotate([(1, tag)], cx, cy, cz, 1, 0, 0,
                                   math.pi / 2)
        elif a != 'z':
            raise ValueError(f"axis must be 'x', 'y', or 'z', got '{axis}'")

        gmsh.model.occ.synchronize()

        if self._verbose:
            print(f"[GmshBuilder] add_circle tag={tag} r={radius} axis={axis}")
        return tag

    def add_ellipse_arc(self, start, center, major_point, end):
        """Create an ellipse arc curve.

        Parameters
        ----------
        start : int
            GMSH point tag for the start of the arc.
        center : int
            GMSH point tag for the center of the ellipse.
        major_point : int
            GMSH point tag on the major axis.
        end : int
            GMSH point tag for the end of the arc.

        Returns
        -------
        int
            GMSH curve tag.
        """
        import gmsh
        self._check_initialized()

        tag = gmsh.model.occ.addEllipseArc(start, center, major_point, end)
        gmsh.model.occ.synchronize()

        if self._verbose:
            print(f"[GmshBuilder] add_ellipse_arc tag={tag}")
        return tag

    def add_bezier(self, point_tags):
        """Create a Bezier curve through control points.

        Parameters
        ----------
        point_tags : list of int
            GMSH point tags defining the Bezier control points.

        Returns
        -------
        int
            GMSH curve tag.
        """
        import gmsh
        self._check_initialized()

        tag = gmsh.model.occ.addBezier(point_tags)
        gmsh.model.occ.synchronize()

        if self._verbose:
            print(f"[GmshBuilder] add_bezier tag={tag} "
                  f"n_pts={len(point_tags)}")
        return tag

    def add_curve_loop(self, curve_tags):
        """Create a curve loop from curves (alias for add_wire).

        Parameters
        ----------
        curve_tags : list of int
            GMSH curve tags forming a closed loop.

        Returns
        -------
        int
            GMSH wire tag.
        """
        import gmsh
        self._check_initialized()

        tag = gmsh.model.occ.addWire(curve_tags)
        gmsh.model.occ.synchronize()

        if self._verbose:
            print(f"[GmshBuilder] add_curve_loop tag={tag} "
                  f"n_curves={len(curve_tags)}")
        return tag

    def import_mesh(self, filename):
        """Import a .msh mesh file using gmsh.merge.

        Parameters
        ----------
        filename : str
            Path to .msh file.

        Returns
        -------
        list of int
            Volume IDs for any volumes found in the mesh.
        """
        import gmsh
        self._check_initialized()

        if not os.path.exists(filename):
            raise FileNotFoundError(f"Mesh file not found: {filename}")

        # Snapshot existing volumes before merge
        existing = set(t for d, t in gmsh.model.getEntities(3))

        gmsh.merge(filename)

        # Find only newly added volumes
        all_vols = set(t for d, t in gmsh.model.getEntities(3))
        new_vols = list(all_vols - existing)
        new_ids = []
        for tag in new_vols:
            uid = self._register_volume(tag)
            new_ids.append(uid)

        if self._verbose:
            print(f"[GmshBuilder] import_mesh '{filename}' "
                  f"-> {len(new_ids)} vols")
        return new_ids

    # ---- Advanced OCC primitives ----

    def add_thick_solid(self, vol_id, exclude_face_tags, offset):
        """Create a thick solid (shell) by offsetting faces of a volume.

        Parameters
        ----------
        vol_id : int
            Volume ID of the source solid.
        exclude_face_tags : list of int
            GMSH face tags to exclude from offsetting (these faces
            become the openings in the resulting shell).
        offset : float
            Offset distance in meters (positive = outward).

        Returns
        -------
        int
            New volume ID.

        Raises
        ------
        NotImplementedError
            If the GMSH version does not support addThickSolid.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        vol_tags = self._volumes[vol_id]
        try:
            out = gmsh.model.occ.addThickSolid(
                vol_tags[0],
                exclude_face_tags,
                offset)
        except AttributeError:
            raise NotImplementedError(
                "addThickSolid is not available in this GMSH version")
        gmsh.model.occ.synchronize()

        new_tags = [dt[1] for dt in out if dt[0] == 3]
        if not new_tags:
            raise RuntimeError("addThickSolid produced no volumes")

        uid = self._register_volume(new_tags)
        if self._verbose:
            print(f"[GmshBuilder] add_thick_solid id={uid} "
                  f"source={vol_id} offset={offset}")
        return uid

    def add_thru_sections(self, wire_tags, make_solid=True,
                          make_ruled=False):
        """Create a volume by lofting through wire cross-sections.

        Parameters
        ----------
        wire_tags : list of int
            GMSH wire/curve loop tags defining cross-sections.
        make_solid : bool
            If True, create a solid volume (default True).
        make_ruled : bool
            If True, use ruled surfaces (default False).

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        out = gmsh.model.occ.addThruSections(
            wire_tags, makeSolid=make_solid, makeRuled=make_ruled)
        gmsh.model.occ.synchronize()

        vol_tags = [dt[1] for dt in out if dt[0] == 3]
        if not vol_tags:
            raise RuntimeError("addThruSections produced no volumes")

        uid = self._register_volume(vol_tags)
        if self._verbose:
            print(f"[GmshBuilder] add_thru_sections id={uid} "
                  f"n_wires={len(wire_tags)} solid={make_solid}")
        return uid

    def add_general_pipe(self, section_dimtags, wire_tag):
        """Create a pipe by sweeping section(s) along a wire path.

        Parameters
        ----------
        section_dimtags : list of tuple
            GMSH (dim, tag) pairs for the cross-section entities.
        wire_tag : int
            GMSH wire tag defining the sweep path.

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        out = gmsh.model.occ.addPipe(section_dimtags, wire_tag)
        gmsh.model.occ.synchronize()

        vol_tags = [dt[1] for dt in out if dt[0] == 3]
        if not vol_tags:
            raise RuntimeError("addPipe produced no volumes")

        uid = self._register_volume(vol_tags)
        if self._verbose:
            print(f"[GmshBuilder] add_general_pipe id={uid} "
                  f"wire={wire_tag}")
        return uid

    def add_bspline_surface(self, point_tags, n_u, n_v,
                            degree_u=3, degree_v=3):
        """Create a B-spline surface from a grid of control points.

        Parameters
        ----------
        point_tags : list of int
            GMSH point tags arranged in a n_u x n_v grid
            (row-major order).
        n_u : int
            Number of control points in the u direction.
        n_v : int
            Number of control points in the v direction.
        degree_u : int
            Degree in u direction (default 3).
        degree_v : int
            Degree in v direction (default 3).

        Returns
        -------
        int
            GMSH surface tag.
        """
        import gmsh
        self._check_initialized()

        try:
            tag = gmsh.model.occ.addBSplineSurface(
                point_tags, n_u,
                degreeU=degree_u, degreeV=degree_v)
        except AttributeError:
            raise NotImplementedError(
                "addBSplineSurface is not available in this GMSH version")
        gmsh.model.occ.synchronize()

        if self._verbose:
            print(f"[GmshBuilder] add_bspline_surface tag={tag} "
                  f"n_u={n_u} n_v={n_v}")
        return tag

    def add_bspline_filling(self, wire_tag, type_str='Coons'):
        """Create a B-spline surface filling a wire boundary.

        Parameters
        ----------
        wire_tag : int
            GMSH wire tag defining the boundary.
        type_str : str
            Filling type: 'Coons', 'Curved', etc. (default 'Coons').

        Returns
        -------
        int
            GMSH surface tag.
        """
        import gmsh
        self._check_initialized()

        try:
            tag = gmsh.model.occ.addBSplineFilling(wire_tag, type=type_str)
        except AttributeError:
            raise NotImplementedError(
                "addBSplineFilling is not available in this GMSH version")
        gmsh.model.occ.synchronize()

        if self._verbose:
            print(f"[GmshBuilder] add_bspline_filling tag={tag} "
                  f"type={type_str}")
        return tag

    def add_bezier_surface(self, point_tags, n_u, n_v):
        """Create a Bezier surface from a grid of control points.

        Parameters
        ----------
        point_tags : list of int
            GMSH point tags arranged in a n_u x n_v grid
            (row-major order).
        n_u : int
            Number of control points in the u direction.
        n_v : int
            Number of control points in the v direction.

        Returns
        -------
        int
            GMSH surface tag.

        Raises
        ------
        NotImplementedError
            If the GMSH version does not support addBezierSurface.
        """
        import gmsh
        self._check_initialized()

        try:
            tag = gmsh.model.occ.addBezierSurface(point_tags, n_u)
        except AttributeError:
            raise NotImplementedError(
                "addBezierSurface is not available in this GMSH version")
        gmsh.model.occ.synchronize()

        if self._verbose:
            print(f"[GmshBuilder] add_bezier_surface tag={tag} "
                  f"n_u={n_u} n_v={n_v}")
        return tag

    def add_bezier_filling(self, wire_tag):
        """Create a Bezier surface filling a wire boundary.

        Parameters
        ----------
        wire_tag : int
            GMSH wire tag defining the boundary.

        Returns
        -------
        int
            GMSH surface tag.

        Raises
        ------
        NotImplementedError
            If the GMSH version does not support addBezierFilling.
        """
        import gmsh
        self._check_initialized()

        try:
            tag = gmsh.model.occ.addBezierFilling(wire_tag)
        except AttributeError:
            raise NotImplementedError(
                "addBezierFilling is not available in this GMSH version")
        gmsh.model.occ.synchronize()

        if self._verbose:
            print(f"[GmshBuilder] add_bezier_filling tag={tag}")
        return tag

    def add_trimmed_surface(self, surf_tag, wire_tags):
        """Create a trimmed surface from a base surface and trimming wires.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag of the base surface.
        wire_tags : list of int
            GMSH wire tags defining the trimming boundaries.

        Returns
        -------
        int
            GMSH surface tag.
        """
        import gmsh
        self._check_initialized()

        try:
            tag = gmsh.model.occ.addTrimmedSurface(surf_tag, wire_tags)
        except AttributeError:
            raise NotImplementedError(
                "addTrimmedSurface is not available in this GMSH version")
        gmsh.model.occ.synchronize()

        if self._verbose:
            print(f"[GmshBuilder] add_trimmed_surface tag={tag} "
                  f"base={surf_tag}")
        return tag

    def add_offset_curve(self, curve_tag, offset):
        """Create a curve offset from an existing curve.

        The GMSH ``offsetCurve`` API accepts only a curve loop tag and an
        offset distance.  A direction vector is **not** supported by GMSH.

        Parameters
        ----------
        curve_tag : int
            GMSH curve loop tag of the source curve.
        offset : float
            Offset distance in meters.

        Returns
        -------
        int
            GMSH curve tag.
        """
        import gmsh
        self._check_initialized()

        try:
            tag = gmsh.model.occ.offsetCurve(curve_tag, offset)
        except AttributeError:
            raise NotImplementedError(
                "offsetCurve is not available in this GMSH version")
        gmsh.model.occ.synchronize()

        if self._verbose:
            print(f"[GmshBuilder] add_offset_curve tag={tag} "
                  f"source={curve_tag} offset={offset}")
        return tag

    def convert_to_nurbs(self, vol_id):
        """Convert volume geometry representation to NURBS.

        Parameters
        ----------
        vol_id : int
            Volume ID to convert.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        vol_tags = self._volumes[vol_id]
        try:
            for tag in vol_tags:
                gmsh.model.occ.convertToNURBS([(3, tag)])
        except AttributeError:
            raise NotImplementedError(
                "convertToNURBS is not available in this GMSH version")
        gmsh.model.occ.synchronize()

        if self._verbose:
            print(f"[GmshBuilder] convert_to_nurbs vol={vol_id}")

    def add_full_ellipse(self, center, rx, ry, axis='z'):
        """Create a full ellipse curve.

        Parameters
        ----------
        center : list
            Center [x, y, z] in meters.
        rx : float
            Semi-axis length in the first in-plane direction.
        ry : float
            Semi-axis length in the second in-plane direction.
        axis : str
            Normal axis of the ellipse plane: 'x', 'y', or 'z'
            (default 'z').

        Returns
        -------
        int
            GMSH curve tag.
        """
        import gmsh
        self._check_initialized()

        cx, cy, cz = center
        tag = gmsh.model.occ.addEllipse(cx, cy, cz, rx, ry)

        a = axis.lower()
        if a == 'x':
            gmsh.model.occ.rotate([(1, tag)], cx, cy, cz, 0, 1, 0,
                                   math.pi / 2)
        elif a == 'y':
            gmsh.model.occ.rotate([(1, tag)], cx, cy, cz, 1, 0, 0,
                                   math.pi / 2)
        elif a != 'z':
            raise ValueError(
                f"axis must be 'x', 'y', or 'z', got '{axis}'")

        gmsh.model.occ.synchronize()

        if self._verbose:
            print(f"[GmshBuilder] add_full_ellipse tag={tag} "
                  f"rx={rx} ry={ry} axis={axis}")
        return tag

    def add_fillet_2d(self, edge_tag1, edge_tag2, radius):
        """Apply a 2D fillet between two edges.

        The GMSH ``fillet2D`` API operates on a pair of edges and a
        single radius, not on a wire with multiple points.

        Parameters
        ----------
        edge_tag1 : int
            GMSH edge (curve) tag of the first edge.
        edge_tag2 : int
            GMSH edge (curve) tag of the second edge.
        radius : float
            Fillet radius in meters.

        Returns
        -------
        int
            GMSH edge tag of the fillet arc.

        Raises
        ------
        NotImplementedError
            If the GMSH version does not support fillet2D.
        """
        import gmsh
        self._check_initialized()

        try:
            out = gmsh.model.occ.fillet2D(edge_tag1, edge_tag2, radius)
        except AttributeError:
            raise NotImplementedError(
                "fillet2D is not available in this GMSH version")
        gmsh.model.occ.synchronize()

        if self._verbose:
            print(f"[GmshBuilder] add_fillet_2d edges=({edge_tag1}, "
                  f"{edge_tag2}) radius={radius}")
        return out

    # ---- Compound/assembly ----

    def add_multi_box(self, centers, sizes):
        """Create multiple boxes in batch.

        Parameters
        ----------
        centers : list of list
            List of center coordinates [[cx, cy, cz], ...].
        sizes : list of list
            List of dimensions [[Lx, Ly, Lz], ...]. If a single
            size is given, it is used for all boxes.

        Returns
        -------
        list of int
            Volume IDs for each box.
        """
        import gmsh
        self._check_initialized()

        if len(sizes) == 1:
            sizes = sizes * len(centers)
        if len(centers) != len(sizes):
            raise ValueError(
                "centers and sizes must have the same length, "
                f"got {len(centers)} and {len(sizes)}")

        vol_ids = []
        for center, size in zip(centers, sizes):
            cx, cy, cz = center
            lx, ly, lz = size
            x0 = cx - lx / 2
            y0 = cy - ly / 2
            z0 = cz - lz / 2
            vol_tag = gmsh.model.occ.addBox(x0, y0, z0, lx, ly, lz)
            uid = self._register_volume(vol_tag)
            vol_ids.append(uid)

        gmsh.model.occ.synchronize()

        if self._verbose:
            print(f"[GmshBuilder] add_multi_box n={len(vol_ids)} "
                  f"ids={vol_ids}")
        return vol_ids

    def add_coil_helix(self, center, axis, radius, pitch, n_turns,
                       wire_radius=None):
        """Create a helical coil.

        If wire_radius is given, a solid pipe (wire) is created around
        the helix. Otherwise, only the helical wire curve is created
        and registered as a volume-less entity.

        Parameters
        ----------
        center : list
            Base center [x, y, z] in meters.
        axis : str
            Helix axis: 'x', 'y', or 'z'.
        radius : float
            Helix radius in meters.
        pitch : float
            Pitch (axial advance per turn) in meters.
        n_turns : float
            Number of turns.
        wire_radius : float, optional
            Radius of the wire cross-section. If None, no solid is
            created.

        Returns
        -------
        int
            Volume ID if wire_radius is given, otherwise the GMSH
            wire tag.
        """
        import gmsh
        self._check_initialized()

        cx, cy, cz = center
        a = axis.lower()
        n_pts_per_turn = 36
        n_pts = int(n_turns * n_pts_per_turn) + 1
        total_height = pitch * n_turns

        pt_tags = []
        for i in range(n_pts):
            t = i / (n_pts - 1)
            angle = 2 * math.pi * n_turns * t
            h = total_height * t
            if a == 'z':
                px = cx + radius * math.cos(angle)
                py = cy + radius * math.sin(angle)
                pz = cz + h
            elif a == 'x':
                px = cx + h
                py = cy + radius * math.cos(angle)
                pz = cz + radius * math.sin(angle)
            elif a == 'y':
                px = cx + radius * math.cos(angle)
                py = cy + h
                pz = cz + radius * math.sin(angle)
            else:
                raise ValueError(
                    f"axis must be 'x', 'y', or 'z', got '{axis}'")
            pt_tags.append(gmsh.model.occ.addPoint(px, py, pz))

        spline_tag = gmsh.model.occ.addBSpline(pt_tags)
        wire_tag = gmsh.model.occ.addWire([spline_tag])

        if wire_radius is not None:
            # Create a circular cross-section at the start point
            if a == 'z':
                disk_tag = gmsh.model.occ.addDisk(
                    cx + radius, cy, cz, wire_radius, wire_radius)
                gmsh.model.occ.rotate(
                    [(2, disk_tag)], cx + radius, cy, cz,
                    0, 1, 0, math.pi / 2)
            elif a == 'x':
                disk_tag = gmsh.model.occ.addDisk(
                    cx, cy + radius, cz, wire_radius, wire_radius)
                gmsh.model.occ.rotate(
                    [(2, disk_tag)], cx, cy + radius, cz,
                    1, 0, 0, math.pi / 2)
            else:  # y
                disk_tag = gmsh.model.occ.addDisk(
                    cx + radius, cy, cz, wire_radius, wire_radius)

            out = gmsh.model.occ.addPipe([(2, disk_tag)], wire_tag)
            gmsh.model.occ.synchronize()

            vol_tags = [dt[1] for dt in out if dt[0] == 3]
            if not vol_tags:
                raise RuntimeError("Coil helix pipe produced no volumes")

            uid = self._register_volume(vol_tags)
            if self._verbose:
                print(f"[GmshBuilder] add_coil_helix id={uid} "
                      f"r={radius} pitch={pitch} n_turns={n_turns} "
                      f"wire_r={wire_radius}")
            return uid
        else:
            gmsh.model.occ.synchronize()
            if self._verbose:
                print(f"[GmshBuilder] add_coil_helix wire={wire_tag} "
                      f"r={radius} pitch={pitch} n_turns={n_turns}")
            return wire_tag

    def add_sector(self, center, axis, r_inner, r_outer, height, angle):
        """Create an annular sector by revolving a rectangular section.

        Parameters
        ----------
        center : list
            Base center [x, y, z] in meters.
        axis : str
            Sector axis: 'x', 'y', or 'z'.
        r_inner : float
            Inner radius in meters.
        r_outer : float
            Outer radius in meters.
        height : float
            Height in meters.
        angle : float
            Sector angle in radians.

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        cx, cy, cz = center
        a = axis.lower()

        # Create rectangular cross-section in the revolve plane
        if a == 'z':
            rect = gmsh.model.occ.addRectangle(
                cx + r_inner, cy, cz, r_outer - r_inner, height)
            ax_x, ax_y, ax_z = 0, 0, 1
        elif a == 'x':
            rect = gmsh.model.occ.addRectangle(
                cx, cy + r_inner, cz, height, r_outer - r_inner)
            # Rotate rect to align with x-axis revolve plane
            ax_x, ax_y, ax_z = 1, 0, 0
        elif a == 'y':
            rect = gmsh.model.occ.addRectangle(
                cx + r_inner, cy, cz, r_outer - r_inner, height)
            ax_x, ax_y, ax_z = 0, 1, 0
        else:
            raise ValueError(
                f"axis must be 'x', 'y', or 'z', got '{axis}'")

        out = gmsh.model.occ.revolve(
            [(2, rect)], cx, cy, cz, ax_x, ax_y, ax_z, angle)
        gmsh.model.occ.synchronize()

        vol_tags = [dt[1] for dt in out if dt[0] == 3]
        if not vol_tags:
            raise RuntimeError("Sector revolve produced no volumes")

        uid = self._register_volume(vol_tags)
        if self._verbose:
            print(f"[GmshBuilder] add_sector id={uid} "
                  f"r=[{r_inner},{r_outer}] h={height} "
                  f"angle={math.degrees(angle):.1f} deg")
        return uid

    def add_racetrack(self, center, straight_length, radius, height):
        """Create a racetrack-shaped volume (stadium extruded).

        The racetrack consists of two straight segments connected
        by semicircular ends, extruded to the given height along z.

        Parameters
        ----------
        center : list
            Center [x, y, z] in meters.
        straight_length : float
            Length of the straight portion in meters.
        radius : float
            Radius of the semicircular ends in meters.
        height : float
            Extrusion height in meters.

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        cx, cy, cz = center
        half_l = straight_length / 2

        # Create 2D racetrack (stadium shape)
        # Left semicircle center at (cx - half_l, cy)
        # Right semicircle center at (cx + half_l, cy)
        p1 = gmsh.model.occ.addPoint(cx + half_l, cy + radius, cz)
        p2 = gmsh.model.occ.addPoint(cx - half_l, cy + radius, cz)
        p3 = gmsh.model.occ.addPoint(cx - half_l, cy - radius, cz)
        p4 = gmsh.model.occ.addPoint(cx + half_l, cy - radius, cz)
        pc_left = gmsh.model.occ.addPoint(cx - half_l, cy, cz)
        pc_right = gmsh.model.occ.addPoint(cx + half_l, cy, cz)

        l_top = gmsh.model.occ.addLine(p1, p2)
        arc_left = gmsh.model.occ.addCircleArc(p2, pc_left, p3)
        l_bot = gmsh.model.occ.addLine(p3, p4)
        arc_right = gmsh.model.occ.addCircleArc(p4, pc_right, p1)

        wire = gmsh.model.occ.addWire(
            [l_top, arc_left, l_bot, arc_right])
        surf = gmsh.model.occ.addPlaneSurface([wire])

        out = gmsh.model.occ.extrude([(2, surf)], 0, 0, height)
        gmsh.model.occ.synchronize()

        vol_tags = [dt[1] for dt in out if dt[0] == 3]
        if not vol_tags:
            raise RuntimeError(
                "Racetrack extrusion produced no volumes")

        uid = self._register_volume(vol_tags)
        if self._verbose:
            print(f"[GmshBuilder] add_racetrack id={uid} "
                  f"straight={straight_length} r={radius} h={height}")
        return uid

    def add_e_shape(self, center, size, slot_depths):
        """Create an E-core shape (box with slots).

        The E-core is a box with two horizontal slots cut from one
        side, forming three prongs.

        Parameters
        ----------
        center : list
            Center [x, y, z] in meters.
        size : list
            Overall dimensions [Lx, Ly, Lz] in meters.
        slot_depths : list
            [depth, width_top, width_bot] in meters. depth is how
            far each slot extends into the body, width_top and
            width_bot are the heights of the top and bottom slots.

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        cx, cy, cz = center
        lx, ly, lz = size
        slot_depth, slot_w_top, slot_w_bot = slot_depths

        x0 = cx - lx / 2
        y0 = cy - ly / 2
        z0 = cz - lz / 2

        # Main body
        body = gmsh.model.occ.addBox(x0, y0, z0, lx, ly, lz)

        # Top slot (cut from +x side)
        slot1 = gmsh.model.occ.addBox(
            x0 + lx - slot_depth,
            y0 + ly - slot_w_top,
            z0,
            slot_depth * 1.01,
            slot_w_top,
            lz)

        # Bottom slot (cut from +x side)
        slot2 = gmsh.model.occ.addBox(
            x0 + lx - slot_depth,
            y0,
            z0,
            slot_depth * 1.01,
            slot_w_bot,
            lz)

        out, _ = gmsh.model.occ.cut(
            [(3, body)],
            [(3, slot1), (3, slot2)])
        gmsh.model.occ.synchronize()

        vol_tags = [dt[1] for dt in out if dt[0] == 3]
        if not vol_tags:
            raise RuntimeError("E-shape boolean cut produced no volumes")

        uid = self._register_volume(vol_tags)
        if self._verbose:
            print(f"[GmshBuilder] add_e_shape id={uid} "
                  f"size={size} slot_depths={slot_depths}")
        return uid

    def add_u_shape(self, center, size, slot_depth, slot_width):
        """Create a U-core shape (box with a single central slot).

        Parameters
        ----------
        center : list
            Center [x, y, z] in meters.
        size : list
            Overall dimensions [Lx, Ly, Lz] in meters.
        slot_depth : float
            Depth of the U-slot (from the top) in meters.
        slot_width : float
            Width of the U-slot in meters.

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        cx, cy, cz = center
        lx, ly, lz = size

        x0 = cx - lx / 2
        y0 = cy - ly / 2
        z0 = cz - lz / 2

        # Main body
        body = gmsh.model.occ.addBox(x0, y0, z0, lx, ly, lz)

        # Central slot (cut from +y side, centered in x)
        slot_x0 = cx - slot_width / 2
        slot = gmsh.model.occ.addBox(
            slot_x0,
            y0 + ly - slot_depth,
            z0,
            slot_width,
            slot_depth * 1.01,
            lz)

        out, _ = gmsh.model.occ.cut(
            [(3, body)], [(3, slot)])
        gmsh.model.occ.synchronize()

        vol_tags = [dt[1] for dt in out if dt[0] == 3]
        if not vol_tags:
            raise RuntimeError("U-shape boolean cut produced no volumes")

        uid = self._register_volume(vol_tags)
        if self._verbose:
            print(f"[GmshBuilder] add_u_shape id={uid} "
                  f"size={size} slot_depth={slot_depth} "
                  f"slot_width={slot_width}")
        return uid

    # ---- Import/interop ----

    def import_sat(self, filename):
        """Import ACIS SAT file geometry.

        Parameters
        ----------
        filename : str
            Path to .sat file.

        Returns
        -------
        list of int
            Volume IDs for imported volumes.
        """
        import gmsh
        self._check_initialized()

        if not os.path.exists(filename):
            raise FileNotFoundError(f"SAT file not found: {filename}")

        out_dimtags = gmsh.model.occ.importShapes(filename)
        gmsh.model.occ.synchronize()

        vol_tags = [dt[1] for dt in out_dimtags if dt[0] == 3]
        new_ids = []
        for tag in vol_tags:
            uid = self._register_volume(tag)
            new_ids.append(uid)

        if self._verbose:
            print(f"[GmshBuilder] import_sat '{filename}' -> "
                  f"{len(vol_tags)} vols, ids={new_ids}")
        return new_ids

    def import_x_t(self, filename):
        """Import Parasolid X_T file geometry.

        Parameters
        ----------
        filename : str
            Path to .x_t file.

        Returns
        -------
        list of int
            Volume IDs for imported volumes.
        """
        import gmsh
        self._check_initialized()

        if not os.path.exists(filename):
            raise FileNotFoundError(f"X_T file not found: {filename}")

        out_dimtags = gmsh.model.occ.importShapes(filename)
        gmsh.model.occ.synchronize()

        vol_tags = [dt[1] for dt in out_dimtags if dt[0] == 3]
        new_ids = []
        for tag in vol_tags:
            uid = self._register_volume(tag)
            new_ids.append(uid)

        if self._verbose:
            print(f"[GmshBuilder] import_x_t '{filename}' -> "
                  f"{len(vol_tags)} vols, ids={new_ids}")
        return new_ids

    def import_3mf(self, filename):
        """Import 3MF file geometry.

        Parameters
        ----------
        filename : str
            Path to .3mf file.

        Returns
        -------
        list of int
            Volume IDs for imported volumes.
        """
        import gmsh
        self._check_initialized()

        if not os.path.exists(filename):
            raise FileNotFoundError(f"3MF file not found: {filename}")

        out_dimtags = gmsh.model.occ.importShapes(filename)
        gmsh.model.occ.synchronize()

        vol_tags = [dt[1] for dt in out_dimtags if dt[0] == 3]
        new_ids = []
        for tag in vol_tags:
            uid = self._register_volume(tag)
            new_ids.append(uid)

        if self._verbose:
            print(f"[GmshBuilder] import_3mf '{filename}' -> "
                  f"{len(vol_tags)} vols, ids={new_ids}")
        return new_ids

    def import_from_string(self, format_str, data,
                           highestDimOnly=True):
        """Import geometry from an in-memory string representation.

        Writes the data to a temporary file with the appropriate
        extension, then imports via GMSH.

        Parameters
        ----------
        format_str : str
            Format identifier (e.g. 'step', 'iges', 'brep').
        data : str
            Geometry data as a string.
        highestDimOnly : bool
            If True, only return entities of the highest dimension
            (default True).

        Returns
        -------
        list of int
            Volume IDs for imported volumes (may be empty if the
            data contains only surfaces/curves).
        """
        import gmsh
        import tempfile
        self._check_initialized()

        ext_map = {
            'step': '.step', 'stp': '.stp',
            'iges': '.iges', 'igs': '.igs',
            'brep': '.brep',
            'stl': '.stl',
        }
        ext = ext_map.get(format_str.lower(), '.' + format_str)

        with tempfile.NamedTemporaryFile(
                mode='w', suffix=ext, delete=False) as f:
            f.write(data)
            tmp_path = f.name

        try:
            out_dimtags = gmsh.model.occ.importShapes(tmp_path,
                                                       highestDimOnly)
            gmsh.model.occ.synchronize()
        finally:
            os.remove(tmp_path)

        vol_tags = [dt[1] for dt in out_dimtags if dt[0] == 3]
        new_ids = []
        for tag in vol_tags:
            uid = self._register_volume(tag)
            new_ids.append(uid)

        if self._verbose:
            print(f"[GmshBuilder] import_from_string format={format_str} "
                  f"-> {len(new_ids)} vols")
        return new_ids

    def export_iges(self, filename):
        """Export current geometry to IGES format.

        Parameters
        ----------
        filename : str
            Output filename (should end in .iges or .igs).
        """
        import gmsh
        self._check_initialized()

        gmsh.model.occ.synchronize()
        gmsh.write(filename)

        if self._verbose:
            print(f"[GmshBuilder] export_iges '{filename}'")

    def export_x3d(self, filename):
        """Export current geometry/mesh to X3D format.

        Parameters
        ----------
        filename : str
            Output filename (should end in .x3d).
        """
        import gmsh
        self._check_initialized()

        gmsh.model.occ.synchronize()
        gmsh.write(filename)

        if self._verbose:
            print(f"[GmshBuilder] export_x3d '{filename}'")

    # ---- Sweep variants ----

    def extrude_rotate(self, surf_tags, center, axis, angle,
                       num_elements=None):
        """Revolve surfaces to create volumes (extrude + rotate).

        Parameters
        ----------
        surf_tags : list of int
            GMSH surface tags to revolve.
        center : list
            Rotation center [x, y, z].
        axis : list
            Rotation axis vector [ax, ay, az].
        angle : float
            Rotation angle in radians.
        num_elements : list of int, optional
            Number of elements along the revolution for structured
            meshing.

        Returns
        -------
        list of int
            Volume IDs.
        """
        import gmsh
        self._check_initialized()

        if num_elements is None:
            num_elements = []

        dimtags = [(2, t) for t in surf_tags]
        out = gmsh.model.occ.revolve(
            dimtags,
            center[0], center[1], center[2],
            axis[0], axis[1], axis[2], angle,
            numElements=num_elements)
        gmsh.model.occ.synchronize()

        new_tags = [dt[1] for dt in out if dt[0] == 3]
        new_ids = []
        for tag in new_tags:
            uid = self._register_volume(tag)
            new_ids.append(uid)

        if self._verbose:
            print(f"[GmshBuilder] extrude_rotate "
                  f"angle={math.degrees(angle):.1f} deg -> {new_ids}")
        return new_ids

    def extrude_along_wire(self, surf_tag, wire_tag):
        """Extrude a surface along a wire path to create a volume.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag to extrude.
        wire_tag : int
            GMSH wire tag defining the sweep path.

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        out = gmsh.model.occ.addPipe([(2, surf_tag)], wire_tag)
        gmsh.model.occ.synchronize()

        vol_tags = [dt[1] for dt in out if dt[0] == 3]
        if not vol_tags:
            raise RuntimeError(
                "extrude_along_wire produced no volumes")

        uid = self._register_volume(vol_tags)
        if self._verbose:
            print(f"[GmshBuilder] extrude_along_wire id={uid} "
                  f"surf={surf_tag} wire={wire_tag}")
        return uid

    def extrude_with_twist(self, surf_tag, dx, dy, dz, angle,
                           num_elements=None):
        """Extrude a surface with a twist rotation.

        .. note::

           GMSH does not provide a ``twist`` API in ``gmsh.model.occ``.
           Twist extrusion is not supported by the OCC kernel.  This
           method raises ``NotImplementedError``.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag to extrude.
        dx, dy, dz : float
            Extrusion direction and distance.
        angle : float
            Twist angle in radians applied over the extrusion.
        num_elements : list of int, optional
            Number of elements along the extrusion for structured
            meshing.

        Returns
        -------
        int
            Volume ID.

        Raises
        ------
        NotImplementedError
            Always raised -- GMSH OCC kernel does not support twist
            extrusion.
        """
        raise NotImplementedError(
            "Twist extrusion is not supported by the GMSH OCC kernel. "
            "gmsh.model.occ.twist does not exist. Use "
            "gmsh.model.occ.extrude for plain extrusion, or build "
            "twisted geometry via loft (add_thru_sections) with "
            "rotated cross-sections.")

    def revolve_partial(self, surf_tag, center, axis, angle):
        """Revolve a surface by a partial angle (less than 2*pi).

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag to revolve.
        center : list
            Rotation center [x, y, z].
        axis : list
            Rotation axis vector [ax, ay, az].
        angle : float
            Rotation angle in radians (should be < 2*pi).

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        out = gmsh.model.occ.revolve(
            [(2, surf_tag)],
            center[0], center[1], center[2],
            axis[0], axis[1], axis[2], angle)
        gmsh.model.occ.synchronize()

        vol_tags = [dt[1] for dt in out if dt[0] == 3]
        if not vol_tags:
            raise RuntimeError(
                "revolve_partial produced no volumes")

        uid = self._register_volume(vol_tags)
        if self._verbose:
            print(f"[GmshBuilder] revolve_partial id={uid} "
                  f"surf={surf_tag} "
                  f"angle={math.degrees(angle):.1f} deg")
        return uid

    def sweep_surface(self, surf_tag, wire_tag):
        """Sweep a surface along a wire to create a solid.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag to sweep.
        wire_tag : int
            GMSH wire tag defining the sweep path.

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        out = gmsh.model.occ.addPipe([(2, surf_tag)], wire_tag)
        gmsh.model.occ.synchronize()

        vol_tags = [dt[1] for dt in out if dt[0] == 3]
        if not vol_tags:
            raise RuntimeError("sweep_surface produced no volumes")

        uid = self._register_volume(vol_tags)
        if self._verbose:
            print(f"[GmshBuilder] sweep_surface id={uid} "
                  f"surf={surf_tag} wire={wire_tag}")
        return uid

    def extrude_to_point(self, surf_tag, point):
        """Extrude a surface towards a target point.

        The extrusion direction is computed from the center of mass
        of the surface to the target point.

        Parameters
        ----------
        surf_tag : int
            GMSH surface tag to extrude.
        point : list
            Target point [x, y, z] in meters.

        Returns
        -------
        int
            Volume ID.
        """
        import gmsh
        self._check_initialized()

        # Get center of mass of the surface
        com = gmsh.model.occ.getCenterOfMass(2, surf_tag)
        dx = point[0] - com[0]
        dy = point[1] - com[1]
        dz = point[2] - com[2]

        out = gmsh.model.occ.extrude([(2, surf_tag)], dx, dy, dz)
        gmsh.model.occ.synchronize()

        vol_tags = [dt[1] for dt in out if dt[0] == 3]
        if not vol_tags:
            raise RuntimeError(
                "extrude_to_point produced no volumes")

        uid = self._register_volume(vol_tags)
        if self._verbose:
            print(f"[GmshBuilder] extrude_to_point id={uid} "
                  f"surf={surf_tag} point={point}")
        return uid
