"""
Transform operations for GmshBuilder.

Gmsh equivalents: move, rotate, copy, reflect, scale, pattern.
"""

import math


class TransformMixin:
    """Mixin providing transform operations."""

    def translate(self, vol_id, dx, dy, dz):
        """Translate a volume (Gmsh 'move').

        Parameters
        ----------
        vol_id : int
            Volume ID.
        dx, dy, dz : float
            Translation vector in meters.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        dimtags = [(3, t) for t in self._volumes[vol_id]]
        gmsh.model.occ.translate(dimtags, dx, dy, dz)
        gmsh.model.occ.synchronize()
        self._meshed = False

        if self._verbose:
            print(f"[GmshBuilder] translate vol={vol_id} d=[{dx},{dy},{dz}]")

    def rotate(self, vol_id, center, axis, angle):
        """Rotate a volume (Gmsh 'rotate').

        Parameters
        ----------
        vol_id : int
            Volume ID.
        center : list
            Rotation center [x, y, z].
        axis : list
            Rotation axis [ax, ay, az].
        angle : float
            Angle in radians.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        dimtags = [(3, t) for t in self._volumes[vol_id]]
        gmsh.model.occ.rotate(dimtags,
                               center[0], center[1], center[2],
                               axis[0], axis[1], axis[2], angle)
        gmsh.model.occ.synchronize()
        self._meshed = False

        if self._verbose:
            print(f"[GmshBuilder] rotate vol={vol_id} "
                  f"angle={math.degrees(angle):.1f} deg")

    def copy(self, vol_id):
        """Copy a volume (Gmsh 'body copy').

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

        dimtags = [(3, t) for t in self._volumes[vol_id]]
        out = gmsh.model.occ.copy(dimtags)
        gmsh.model.occ.synchronize()

        new_tags = [dt[1] for dt in out if dt[0] == 3]
        uid = self._register_volume(new_tags)
        if self._verbose:
            print(f"[GmshBuilder] copy vol={vol_id} -> id={uid}")
        return uid

    def mirror(self, vol_id, plane_axis):
        """Mirror/reflect a volume (Gmsh 'reflect').

        Parameters
        ----------
        vol_id : int
            Volume ID.
        plane_axis : str
            Mirror plane normal: 'x', 'y', or 'z'.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        a, b, c, d = 0, 0, 0, 0
        pa = plane_axis.lower()
        if pa == 'x':
            a = 1
        elif pa == 'y':
            b = 1
        elif pa == 'z':
            c = 1
        else:
            raise ValueError(f"plane_axis must be 'x','y','z', got '{plane_axis}'")

        dimtags = [(3, t) for t in self._volumes[vol_id]]
        gmsh.model.occ.mirror(dimtags, a, b, c, d)
        gmsh.model.occ.synchronize()
        self._meshed = False

        if self._verbose:
            print(f"[GmshBuilder] mirror vol={vol_id} plane={plane_axis}")

    def scale(self, vol_id, sx, sy=None, sz=None, center=None):
        """Scale a volume (Gmsh 'scale').

        Parameters
        ----------
        vol_id : int
            Volume ID.
        sx : float
            Scale factor (or uniform factor if sy/sz not given).
        sy, sz : float, optional
            Scale factors for y and z. Defaults to sx (uniform).
        center : list, optional
            Scale center [x,y,z]. Default [0,0,0].
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        if sy is None:
            sy = sx
        if sz is None:
            sz = sx
        if center is None:
            center = [0, 0, 0]

        dimtags = [(3, t) for t in self._volumes[vol_id]]
        gmsh.model.occ.dilate(dimtags,
                               center[0], center[1], center[2],
                               sx, sy, sz)
        gmsh.model.occ.synchronize()
        self._meshed = False

        if self._verbose:
            print(f"[GmshBuilder] scale vol={vol_id} s=[{sx},{sy},{sz}]")

    def copy_translate(self, vol_id, dx, dy, dz):
        """Copy and translate in one step (Gmsh 'copy move').

        Returns
        -------
        int
            New volume ID.
        """
        new_id = self.copy(vol_id)
        self.translate(new_id, dx, dy, dz)
        return new_id

    def copy_rotate(self, vol_id, center, axis, angle):
        """Copy and rotate in one step (Gmsh 'copy rotate').

        Returns
        -------
        int
            New volume ID.
        """
        new_id = self.copy(vol_id)
        self.rotate(new_id, center, axis, angle)
        return new_id

    def copy_mirror(self, vol_id, plane_axis):
        """Copy and mirror in one step (Gmsh 'copy reflect').

        Returns
        -------
        int
            New volume ID.
        """
        new_id = self.copy(vol_id)
        self.mirror(new_id, plane_axis)
        return new_id

    def array_linear(self, vol_id, direction, n, spacing):
        """Create a linear array of copies (pattern copy).

        Parameters
        ----------
        vol_id : int
            Volume ID to replicate.
        direction : list
            Direction vector [dx, dy, dz] (normalized internally).
        n : int
            Number of copies (excluding original).
        spacing : float
            Distance between copies.

        Returns
        -------
        list of int
            New volume IDs.
        """
        import math

        length = math.sqrt(sum(d * d for d in direction))
        if length < 1e-20:
            raise ValueError("direction vector is zero")
        ux, uy, uz = [d / length for d in direction]

        new_ids = []
        for i in range(1, n + 1):
            dx = ux * spacing * i
            dy = uy * spacing * i
            dz = uz * spacing * i
            new_id = self.copy_translate(vol_id, dx, dy, dz)
            new_ids.append(new_id)

        if self._verbose:
            print(f"[GmshBuilder] array_linear vol={vol_id} n={n} -> {new_ids}")
        return new_ids

    def array_circular(self, vol_id, center, axis, n, angle=None):
        """Create a circular array of copies (pattern copy).

        Parameters
        ----------
        vol_id : int
            Volume ID to replicate.
        center : list
            Rotation center [x, y, z].
        axis : list
            Rotation axis [ax, ay, az].
        n : int
            Number of copies (excluding original).
        angle : float, optional
            Total angle in radians. Default 2*pi (full circle).

        Returns
        -------
        list of int
            New volume IDs.
        """
        if angle is None:
            angle = 2 * math.pi

        if angle >= 2 * math.pi - 1e-10:
            step = angle / (n + 1)  # full circle: avoid overlap
        else:
            step = angle / n  # partial arc: include endpoint
        new_ids = []
        for i in range(1, n + 1):
            a = step * i
            new_id = self.copy_rotate(vol_id, center, axis, a)
            new_ids.append(new_id)

        if self._verbose:
            print(f"[GmshBuilder] array_circular vol={vol_id} n={n} -> {new_ids}")
        return new_ids

    def affine_transform(self, vol_id, matrix):
        """Apply affine transformation (4x4 or flat 16-element list).

        Parameters
        ----------
        vol_id : int
            Volume ID.
        matrix : list
            Affine matrix as flat list [a11,a12,...,a44] (16 elements)
            or nested 4x4 list.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        # Flatten if nested
        if isinstance(matrix[0], (list, tuple)):
            flat = []
            for row in matrix:
                flat.extend(row)
            matrix = flat

        if len(matrix) not in (12, 16):
            raise ValueError("matrix must have 12 or 16 elements")

        dimtags = [(3, t) for t in self._volumes[vol_id]]
        gmsh.model.occ.affineTransform(dimtags, matrix)
        gmsh.model.occ.synchronize()
        self._meshed = False

        if self._verbose:
            print(f"[GmshBuilder] affine_transform vol={vol_id}")

    def translate_to(self, vol_id, target_point):
        """Move volume so its center is at target_point.

        Computes the current center via get_center() and translates
        by the difference to reach target_point.

        Parameters
        ----------
        vol_id : int
            Volume ID.
        target_point : list
            Target position [x, y, z] in meters.
        """
        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        center = self.get_center(vol_id)
        dx = target_point[0] - center[0]
        dy = target_point[1] - center[1]
        dz = target_point[2] - center[2]
        self.translate(vol_id, dx, dy, dz)

        if self._verbose:
            print(f"[GmshBuilder] translate_to vol={vol_id} "
                  f"target=[{target_point[0]},{target_point[1]},{target_point[2]}]")

    def rotate_degrees(self, vol_id, center, axis, angle_deg):
        """Rotate a volume by an angle specified in degrees.

        Convenience wrapper around rotate() that accepts degrees
        instead of radians.

        Parameters
        ----------
        vol_id : int
            Volume ID.
        center : list
            Rotation center [x, y, z].
        axis : list
            Rotation axis [ax, ay, az].
        angle_deg : float
            Angle in degrees.
        """
        angle_rad = math.radians(angle_deg)
        self.rotate(vol_id, center, axis, angle_rad)

        if self._verbose:
            print(f"[GmshBuilder] rotate_degrees vol={vol_id} "
                  f"angle={angle_deg:.1f} deg")

    def mirror_plane(self, vol_id, point, normal):
        """Mirror a volume about an arbitrary plane.

        The plane is defined by a point on the plane and its normal
        vector. Uses gmsh.model.occ.mirror(dimtags, a, b, c, d) where
        the plane equation is ax + by + cz + d = 0.

        Parameters
        ----------
        vol_id : int
            Volume ID.
        point : list
            A point on the mirror plane [x, y, z].
        normal : list
            Normal vector of the plane [nx, ny, nz].
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        a, b, c = normal[0], normal[1], normal[2]
        d = -(a * point[0] + b * point[1] + c * point[2])

        dimtags = [(3, t) for t in self._volumes[vol_id]]
        gmsh.model.occ.mirror(dimtags, a, b, c, d)
        gmsh.model.occ.synchronize()
        self._meshed = False

        if self._verbose:
            print(f"[GmshBuilder] mirror_plane vol={vol_id} "
                  f"point={point} normal={normal}")

    def copy_scale(self, vol_id, factor, center=None):
        """Copy a volume and scale the copy.

        Parameters
        ----------
        vol_id : int
            Volume ID to copy.
        factor : float
            Uniform scale factor.
        center : list, optional
            Scale center [x, y, z]. Default [0, 0, 0].

        Returns
        -------
        int
            New volume ID.
        """
        new_id = self.copy(vol_id)
        self.scale(new_id, factor, center=center)

        if self._verbose:
            print(f"[GmshBuilder] copy_scale vol={vol_id} "
                  f"factor={factor} -> id={new_id}")
        return new_id

    def array_grid(self, vol_id, nx, ny, nz, sx, sy, sz):
        """Create a 3D grid array of copies.

        Creates nx * ny * nz copies arranged in a rectangular grid
        with spacing sx, sy, sz along each axis. The original volume
        is at index (0, 0, 0) and is not duplicated.

        Parameters
        ----------
        vol_id : int
            Volume ID to replicate.
        nx, ny, nz : int
            Number of copies along x, y, z (excluding original).
        sx, sy, sz : float
            Spacing between copies in meters along x, y, z.

        Returns
        -------
        list of int
            New volume IDs.
        """
        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        new_ids = []
        for ix in range(nx + 1):
            for iy in range(ny + 1):
                for iz in range(nz + 1):
                    if ix == 0 and iy == 0 and iz == 0:
                        continue
                    dx = ix * sx
                    dy = iy * sy
                    dz = iz * sz
                    new_id = self.copy_translate(vol_id, dx, dy, dz)
                    new_ids.append(new_id)

        if self._verbose:
            print(f"[GmshBuilder] array_grid vol={vol_id} "
                  f"n=[{nx},{ny},{nz}] s=[{sx},{sy},{sz}] "
                  f"-> {len(new_ids)} copies")
        return new_ids

    def array_along_curve(self, vol_id, curve_tag, n):
        """Create copies of a volume distributed along a curve.

        Uses gmsh.model.getParametrization to find evenly spaced
        points along the curve and translates copies to those
        positions relative to the volume center.

        Parameters
        ----------
        vol_id : int
            Volume ID to replicate.
        curve_tag : int
            GMSH curve tag to distribute copies along.
        n : int
            Number of copies (excluding original).

        Returns
        -------
        list of int
            New volume IDs.
        """
        import gmsh
        self._check_initialized()

        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        bounds = gmsh.model.getParametrizationBounds(1, curve_tag)
        t_min = bounds[0][0]
        t_max = bounds[1][0]

        center = self.get_center(vol_id)

        new_ids = []
        for i in range(1, n + 1):
            t = t_min + (t_max - t_min) * i / n
            pt = gmsh.model.getValue(1, curve_tag, [t])
            dx = pt[0] - center[0]
            dy = pt[1] - center[1]
            dz = pt[2] - center[2]
            new_id = self.copy_translate(vol_id, dx, dy, dz)
            new_ids.append(new_id)

        if self._verbose:
            print(f"[GmshBuilder] array_along_curve vol={vol_id} "
                  f"curve={curve_tag} n={n} -> {new_ids}")
        return new_ids

    def align_to_axis(self, vol_id, from_axis, to_axis):
        """Rotate volume to align from_axis direction to to_axis direction.

        Computes the rotation axis via cross product and the rotation
        angle via dot product, then applies the rotation about the
        volume center.

        Parameters
        ----------
        vol_id : int
            Volume ID.
        from_axis : list
            Current direction vector [x, y, z].
        to_axis : list
            Target direction vector [x, y, z].
        """
        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        # Normalize from_axis
        f_len = math.sqrt(sum(v * v for v in from_axis))
        if f_len < 1e-20:
            raise ValueError("from_axis vector is zero")
        fx, fy, fz = from_axis[0] / f_len, from_axis[1] / f_len, from_axis[2] / f_len

        # Normalize to_axis
        t_len = math.sqrt(sum(v * v for v in to_axis))
        if t_len < 1e-20:
            raise ValueError("to_axis vector is zero")
        tx, ty, tz = to_axis[0] / t_len, to_axis[1] / t_len, to_axis[2] / t_len

        # Cross product: rotation axis
        cx = fy * tz - fz * ty
        cy = fz * tx - fx * tz
        cz = fx * ty - fy * tx
        cross_len = math.sqrt(cx * cx + cy * cy + cz * cz)

        # Dot product: rotation angle
        dot = fx * tx + fy * ty + fz * tz
        dot = max(-1.0, min(1.0, dot))
        angle = math.acos(dot)

        if cross_len < 1e-12:
            if dot > 0:
                # Vectors are parallel, no rotation needed
                if self._verbose:
                    print(f"[GmshBuilder] align_to_axis vol={vol_id} "
                          f"already aligned")
                return
            else:
                # Vectors are anti-parallel, pick arbitrary perpendicular axis
                if abs(fx) < 0.9:
                    cx, cy, cz = 0, -fz, fy
                else:
                    cx, cy, cz = -fz, 0, fx
                cross_len = math.sqrt(cx * cx + cy * cy + cz * cz)

        rot_axis = [cx / cross_len, cy / cross_len, cz / cross_len]
        center = self.get_center(vol_id)
        self.rotate(vol_id, center, rot_axis, angle)

        if self._verbose:
            print(f"[GmshBuilder] align_to_axis vol={vol_id} "
                  f"angle={math.degrees(angle):.1f} deg")

    def center_at_origin(self, vol_id):
        """Move volume center to the origin [0, 0, 0].

        Parameters
        ----------
        vol_id : int
            Volume ID.
        """
        self.translate_to(vol_id, [0.0, 0.0, 0.0])

        if self._verbose:
            print(f"[GmshBuilder] center_at_origin vol={vol_id}")

    def translate_all(self, dx, dy, dz):
        """Translate all volumes.

        Parameters
        ----------
        dx, dy, dz : float
            Translation vector in meters.
        """
        import gmsh
        self._check_initialized()

        all_dimtags = []
        for vol_id in self._volumes:
            all_dimtags.extend([(3, t) for t in self._volumes[vol_id]])

        if not all_dimtags:
            if self._verbose:
                print("[GmshBuilder] translate_all: no volumes to translate")
            return

        gmsh.model.occ.translate(all_dimtags, dx, dy, dz)
        gmsh.model.occ.synchronize()
        self._meshed = False

        if self._verbose:
            print(f"[GmshBuilder] translate_all d=[{dx},{dy},{dz}] "
                  f"({len(self._volumes)} volumes)")

    def rotate_all(self, center, axis, angle):
        """Rotate all volumes.

        Parameters
        ----------
        center : list
            Rotation center [x, y, z].
        axis : list
            Rotation axis [ax, ay, az].
        angle : float
            Angle in radians.
        """
        import gmsh
        self._check_initialized()

        all_dimtags = []
        for vol_id in self._volumes:
            all_dimtags.extend([(3, t) for t in self._volumes[vol_id]])

        if not all_dimtags:
            if self._verbose:
                print("[GmshBuilder] rotate_all: no volumes to rotate")
            return

        gmsh.model.occ.rotate(all_dimtags,
                               center[0], center[1], center[2],
                               axis[0], axis[1], axis[2], angle)
        gmsh.model.occ.synchronize()
        self._meshed = False

        if self._verbose:
            print(f"[GmshBuilder] rotate_all "
                  f"angle={math.degrees(angle):.1f} deg "
                  f"({len(self._volumes)} volumes)")

    def mirror_all(self, plane_axis):
        """Mirror all volumes about a coordinate plane.

        Parameters
        ----------
        plane_axis : str
            Mirror plane normal: 'x', 'y', or 'z'.
        """
        import gmsh
        self._check_initialized()

        a, b, c, d = 0, 0, 0, 0
        pa = plane_axis.lower()
        if pa == 'x':
            a = 1
        elif pa == 'y':
            b = 1
        elif pa == 'z':
            c = 1
        else:
            raise ValueError(f"plane_axis must be 'x','y','z', got '{plane_axis}'")

        all_dimtags = []
        for vol_id in self._volumes:
            all_dimtags.extend([(3, t) for t in self._volumes[vol_id]])

        if not all_dimtags:
            if self._verbose:
                print("[GmshBuilder] mirror_all: no volumes to mirror")
            return

        gmsh.model.occ.mirror(all_dimtags, a, b, c, d)
        gmsh.model.occ.synchronize()
        self._meshed = False

        if self._verbose:
            print(f"[GmshBuilder] mirror_all plane={plane_axis} "
                  f"({len(self._volumes)} volumes)")

    def scale_all(self, factor, center=None):
        """Scale all volumes uniformly.

        Parameters
        ----------
        factor : float
            Uniform scale factor.
        center : list, optional
            Scale center [x, y, z]. Default [0, 0, 0].
        """
        import gmsh
        self._check_initialized()

        if center is None:
            center = [0, 0, 0]

        all_dimtags = []
        for vol_id in self._volumes:
            all_dimtags.extend([(3, t) for t in self._volumes[vol_id]])

        if not all_dimtags:
            if self._verbose:
                print("[GmshBuilder] scale_all: no volumes to scale")
            return

        gmsh.model.occ.dilate(all_dimtags,
                               center[0], center[1], center[2],
                               factor, factor, factor)
        gmsh.model.occ.synchronize()
        self._meshed = False

        if self._verbose:
            print(f"[GmshBuilder] scale_all factor={factor} "
                  f"({len(self._volumes)} volumes)")

    def copy_all(self):
        """Copy all volumes.

        Returns
        -------
        dict
            Mapping of {old_vol_id: new_vol_id}.
        """
        import gmsh
        self._check_initialized()

        id_map = {}
        for vol_id in list(self._volumes.keys()):
            new_id = self.copy(vol_id)
            id_map[vol_id] = new_id

        if self._verbose:
            print(f"[GmshBuilder] copy_all: {len(id_map)} volumes copied")
        return id_map

    def swap_volumes(self, vol_id1, vol_id2):
        """Swap the GMSH tags between two volume IDs.

        After swapping, vol_id1 will reference the GMSH entities
        that vol_id2 previously referenced, and vice versa.

        Parameters
        ----------
        vol_id1 : int
            First volume ID.
        vol_id2 : int
            Second volume ID.
        """
        if vol_id1 not in self._volumes:
            raise KeyError(f"Volume {vol_id1} not found")
        if vol_id2 not in self._volumes:
            raise KeyError(f"Volume {vol_id2} not found")

        self._volumes[vol_id1], self._volumes[vol_id2] = (
            self._volumes[vol_id2], self._volumes[vol_id1])

        # Swap materials
        mat1 = self._materials.pop(vol_id1, None)
        mat2 = self._materials.pop(vol_id2, None)
        if mat2 is not None:
            self._materials[vol_id1] = mat2
        if mat1 is not None:
            self._materials[vol_id2] = mat1

        # Swap physical groups
        pg1 = self._physical_groups.pop(vol_id1, None)
        pg2 = self._physical_groups.pop(vol_id2, None)
        if pg2 is not None:
            self._physical_groups[vol_id1] = pg2
        if pg1 is not None:
            self._physical_groups[vol_id2] = pg1

        # Swap names
        name1 = self._names.pop(vol_id1, None)
        name2 = self._names.pop(vol_id2, None)
        if name2 is not None:
            self._names[vol_id1] = name2
        if name1 is not None:
            self._names[vol_id2] = name1

        if self._verbose:
            print(f"[GmshBuilder] swap_volumes {vol_id1} <-> {vol_id2}")

    def symmetrize(self, vol_id, plane_axis):
        """Create a mirror copy and fuse with original to make symmetric.

        Copies the volume, mirrors the copy about the given plane,
        then fuses both into a single symmetric volume.

        Parameters
        ----------
        vol_id : int
            Volume ID.
        plane_axis : str
            Mirror plane normal: 'x', 'y', or 'z'.

        Returns
        -------
        int
            Volume ID of the fused symmetric result.
        """
        if vol_id not in self._volumes:
            raise KeyError(f"Volume {vol_id} not found")

        new_id = self.copy(vol_id)
        self.mirror(new_id, plane_axis)
        result_id = self.fuse([vol_id, new_id])

        if self._verbose:
            print(f"[GmshBuilder] symmetrize vol={vol_id} "
                  f"plane={plane_axis} -> id={result_id}")
        return result_id
