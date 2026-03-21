"""
coil_geometry.py - OCC shape generation for Radia coil objects.

Generates netgen.occ shapes corresponding to Radia coil definitions,
enabling STEP export for Cubit import and mesh generation.

Pattern follows EMPY_COIL (K. Kameari): each coil has a geometry (.geo)
that can be exported to STEP for CAD interchange.

Usage:
    from radia.coil_geometry import CoilBlock, CoilArc, CoilRacetrack, CoilAssembly

    # Rectangular coil block
    coil1 = CoilBlock(center=[0, 0, 0], size=[0.1, 0.05, 0.2], j=1e6)

    # Arc coil (sector)
    coil2 = CoilArc(center=[0, 0, 0], r_inner=0.04, r_outer=0.06,
                     phi_start=0, phi_end=90, height=0.1, axis='z')

    # Racetrack coil
    coil3 = CoilRacetrack(center=[0, 0, 0], r_inner=0.03, r_outer=0.05,
                           straight_len=[0.2, 0.1], height=0.08, axis='z')

    # Assembly
    assembly = CoilAssembly()
    assembly.add(coil1)
    assembly.add(coil2)
    assembly.write_step("coils.step")
"""

from netgen.occ import *


class CoilBlock:
    """Rectangular coil block (corresponds to rad.ObjRecCur).

    Args:
        center: [x, y, z] center position in meters
        size: [lx, ly, lz] dimensions in meters
        j: current density in A/m^2 (for reference, not used in geometry)
        j_vec: [jx, jy, jz] current density vector (optional)
    """

    def __init__(self, center, size, j=0, j_vec=None):
        cx, cy, cz = center
        lx, ly, lz = size
        p1 = Pnt(cx - lx/2, cy - ly/2, cz - lz/2)
        p2 = Pnt(cx + lx/2, cy + ly/2, cz + lz/2)
        self.geo = Box(p1, p2)
        self.geo.mat("coil")
        self.center = center
        self.size = size
        self.j = j
        self.j_vec = j_vec or [0, 0, 0]


class CoilArc:
    """Arc coil with rectangular cross-section (corresponds to rad.ObjArcCur).

    Creates a sector of a torus by revolving a rectangle around an axis.

    Args:
        center: [x, y, z] center of rotation axis
        r_inner: inner radius in meters
        r_outer: outer radius in meters
        phi_start: start angle in degrees
        phi_end: end angle in degrees (360 for full loop)
        height: coil height in meters
        axis: rotation axis ('x', 'y', or 'z')
        j: azimuthal current density in A/m^2
    """

    def __init__(self, center, r_inner, r_outer, phi_start, phi_end,
                 height, axis='z', j=0):
        cx, cy, cz = center
        r_mid = (r_inner + r_outer) / 2
        width = r_outer - r_inner
        angle = phi_end - phi_start

        # Build cross-section rectangle in the r-z plane
        if axis == 'z':
            wp = WorkPlane(Axes(Pnt(r_inner, 0, -height/2), n=-Y, h=X))
            rot_axis = Axis(Pnt(0, 0, 0), Z)
            move_vec = Vec(cx, cy, cz)
        elif axis == 'y':
            wp = WorkPlane(Axes(Pnt(r_inner, -height/2, 0), n=-Z, h=X))
            rot_axis = Axis(Pnt(0, 0, 0), Y)
            move_vec = Vec(cx, cy, cz)
        else:  # x
            wp = WorkPlane(Axes(Pnt(-height/2, r_inner, 0), n=-Z, h=Y))
            rot_axis = Axis(Pnt(0, 0, 0), X)
            move_vec = Vec(cx, cy, cz)

        face = wp.Rectangle(width, height).Face()
        self.geo = face.Revolve(rot_axis, angle)

        if phi_start != 0:
            self.geo = self.geo.Rotate(rot_axis, phi_start)

        if cx != 0 or cy != 0 or cz != 0:
            self.geo = self.geo.Move(move_vec)

        self.geo.mat("coil")
        self.center = center
        self.r_inner = r_inner
        self.r_outer = r_outer
        self.phi_start = phi_start
        self.phi_end = phi_end
        self.height = height
        self.axis = axis
        self.j = j


class CoilLoop(CoilArc):
    """Full-loop coil (360-degree arc). Convenience wrapper for CoilArc.

    Args:
        center: [x, y, z] center position
        r_inner: inner radius
        r_outer: outer radius
        height: coil height
        axis: rotation axis ('x', 'y', or 'z')
        j: azimuthal current density
    """

    def __init__(self, center, r_inner, r_outer, height, axis='z', j=0):
        super().__init__(center, r_inner, r_outer, 0, 360, height, axis, j)


class CoilRacetrack:
    """Racetrack coil (corresponds to rad.ObjRaceTrk).

    Consists of two straight sections connected by two semicircular bends.

    Args:
        center: [x, y, z] center position
        r_inner: inner bend radius
        r_outer: outer bend radius
        straight_len: [lx, ly] straight section lengths
        height: coil height
        axis: coil normal axis ('z' means coil in XY plane)
        j: azimuthal current density
    """

    def __init__(self, center, r_inner, r_outer, straight_len, height,
                 axis='z', j=0):
        cx, cy, cz = center
        lx, ly = straight_len
        width = r_outer - r_inner
        half_h = height / 2

        # Build racetrack as union of 4 parts:
        # 2 straight blocks + 2 semicircular bends
        parts = []

        if axis == 'z':
            # Straight section along X (top, +Y side)
            if lx > 0:
                top = Box(Pnt(-lx/2, ly/2 + r_inner, -half_h),
                          Pnt(lx/2, ly/2 + r_outer, half_h))
                parts.append(top)
                # Bottom straight
                bot = Box(Pnt(-lx/2, -ly/2 - r_outer, -half_h),
                          Pnt(lx/2, -ly/2 - r_inner, half_h))
                parts.append(bot)

            # Straight section along Y (right, +X side)
            if ly > 0:
                right = Box(Pnt(lx/2 + r_inner, -ly/2, -half_h),
                            Pnt(lx/2 + r_outer, ly/2, half_h))
                parts.append(right)
                # Left straight
                left = Box(Pnt(-lx/2 - r_outer, -ly/2, -half_h),
                           Pnt(-lx/2 - r_inner, ly/2, half_h))
                parts.append(left)

            # 4 quarter-circle bends
            face = WorkPlane(Axes(Pnt(r_inner, 0, -half_h), n=-Y, h=X)) \
                .Rectangle(width, height).Face()

            # Bend at +X, +Y corner
            bend1 = face.Revolve(Axis(Pnt(0, 0, 0), Z), 90)
            bend1 = bend1.Move(Vec(lx/2, ly/2, 0))
            parts.append(bend1)

            # Bend at -X, +Y corner
            bend2 = face.Revolve(Axis(Pnt(0, 0, 0), Z), 90)
            bend2 = bend2.Rotate(Axis(Pnt(0, 0, 0), Z), 90)
            bend2 = bend2.Move(Vec(-lx/2, ly/2, 0))
            parts.append(bend2)

            # Bend at -X, -Y corner
            bend3 = face.Revolve(Axis(Pnt(0, 0, 0), Z), 90)
            bend3 = bend3.Rotate(Axis(Pnt(0, 0, 0), Z), 180)
            bend3 = bend3.Move(Vec(-lx/2, -ly/2, 0))
            parts.append(bend3)

            # Bend at +X, -Y corner
            bend4 = face.Revolve(Axis(Pnt(0, 0, 0), Z), 90)
            bend4 = bend4.Rotate(Axis(Pnt(0, 0, 0), Z), 270)
            bend4 = bend4.Move(Vec(lx/2, -ly/2, 0))
            parts.append(bend4)

        # Union all parts
        self.geo = parts[0]
        for p in parts[1:]:
            self.geo = self.geo + p

        if cx != 0 or cy != 0 or cz != 0:
            self.geo = self.geo.Move(Vec(cx, cy, cz))

        self.geo.mat("coil")
        self.center = center
        self.r_inner = r_inner
        self.r_outer = r_outer
        self.straight_len = straight_len
        self.height = height
        self.axis = axis
        self.j = j


class CoilAssembly:
    """Collection of coils with combined OCC shape.

    Provides STEP export for Cubit import.

    Usage:
        assembly = CoilAssembly()
        assembly.add(CoilBlock(...))
        assembly.add(CoilArc(...))
        assembly.write_step("coils.step")

        # Get combined shape for Netgen meshing
        geo = assembly.occ_geometry()
    """

    def __init__(self):
        self.coils = []
        self.geo = None

    def add(self, coil):
        """Add a coil to the assembly."""
        self.coils.append(coil)
        if self.geo is None:
            self.geo = coil.geo
        else:
            self.geo = Glue([self.geo, coil.geo])

    def write_step(self, filename):
        """Export assembly to STEP file for Cubit import."""
        if self.geo is None:
            raise ValueError("No coils added to assembly")
        self.geo.WriteStep(filename)
        return filename

    def occ_geometry(self):
        """Return OCCGeometry for direct Netgen meshing."""
        if self.geo is None:
            raise ValueError("No coils added to assembly")
        return OCCGeometry(self.geo)
