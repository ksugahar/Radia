#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Radia Coil Builder - Elegant fluent interface for constructing complex coil geometries.

This module provides a modern object-oriented design for defining multi-segment
coil paths with automatic state tracking and seamless Radia integration.

Example:
    >>> from radia_coil_builder import CoilBuilder
    >>>
    >>> # Create a racetrack coil
    >>> coil = (CoilBuilder(current=1000)
    ...     .set_start([0, 0, 0])
    ...     .set_cross_section(width=20, height=20)
    ...     .add_straight(100)
    ...     .add_arc(radius=50, arc_angle=180, tilt=90)
    ...     .add_straight(100)
    ...     .add_arc(radius=50, arc_angle=180, tilt=90)
    ...     .to_radia())
"""

import numpy as np
from scipy.spatial.transform import Rotation
from abc import ABC, abstractmethod


class CoilSegment(ABC):
    """
    Abstract base class for coil segments.

    All coil segments must implement end_pos and end_orientation properties
    to enable automatic state tracking in the builder pattern.
    """

    def __init__(self, current, start_pos, orientation, width, height, tilt=0):
        """
        Initialize coil segment.

        Args:
            current (float): Current in Amperes
            start_pos (array): Starting position [x, y, z] in mm
            orientation (array): 3x3 orientation matrix (row vectors)
            width (float): Cross-section width in mm
            height (float): Cross-section height in mm
            tilt (float): Tilt angle in degrees (not applied here, applied in subclass)
        """
        self.current = current
        self.start_pos = np.array(start_pos)
        self.orientation = np.array(orientation)
        self.width = width
        self.height = height

        # Extract Euler angles for Radia transformations
        rot = Rotation.from_matrix(self.orientation)
        self.euler_angles = rot.as_euler('ZXZ', degrees=True) * (-1)

    @property
    @abstractmethod
    def end_pos(self):
        """End position of the segment."""
        pass

    @property
    @abstractmethod
    def end_orientation(self):
        """End orientation matrix of the segment."""
        pass

    @property
    def center(self):
        """Geometric center point (midpoint between start and end)."""
        return (self.start_pos + self.end_pos) / 2

    @property
    def current_density(self):
        """Current density in A/mm²."""
        return self.current / (self.width * self.height)


class StraightSegment(CoilSegment):
    """
    Straight coil segment with optional tilt.

    The tilt rotates the cross-section around the Y-axis before
    extending in the local Y direction.
    """

    def __init__(self, current, start_pos, orientation, width, height, length, tilt=0):
        """
        Initialize straight segment.

        Args:
            current (float): Current in Amperes
            start_pos (array): Starting position [x, y, z] in mm
            orientation (array): 3x3 orientation matrix (row vectors)
            width (float): Cross-section width in mm
            height (float): Cross-section height in mm
            length (float): Segment length in mm
            tilt (float): Tilt angle in degrees (rotation around Y-axis)
        """
        # Apply tilt transformation to orientation
        tilt_rad = np.deg2rad(tilt)
        tilt_matrix = np.array([
            [np.cos(tilt_rad), 0, -np.sin(tilt_rad)],
            [0, 1, 0],
            [np.sin(tilt_rad), 0, np.cos(tilt_rad)]
        ])
        tilted_orientation = tilt_matrix @ orientation

        # Cross-section dimensions change with tilt
        tilted_width = abs(np.cos(tilt_rad) * width + np.sin(tilt_rad) * height)
        tilted_height = abs(-np.sin(tilt_rad) * width + np.cos(tilt_rad) * height)

        super().__init__(current, start_pos, tilted_orientation, tilted_width, tilted_height, tilt)
        self.length = length

    @property
    def end_pos(self):
        """End position: start + length * Y-direction."""
        return self.start_pos + self.length * self.orientation[1, :]

    @property
    def end_orientation(self):
        """End orientation: same as start (no rotation)."""
        return self.orientation


class ArcSegment(CoilSegment):
    """
    Arc coil segment with optional tilt.

    The arc rotates around a center point in the local XY plane.
    Tilt is applied first, then the arc rotation.
    """

    def __init__(self, current, start_pos, orientation, width, height, radius, arc_angle, tilt=0):
        """
        Initialize arc segment.

        Args:
            current (float): Current in Amperes
            start_pos (array): Starting position [x, y, z] in mm
            orientation (array): 3x3 orientation matrix (row vectors)
            width (float): Cross-section width in mm
            height (float): Cross-section height in mm
            radius (float): Arc radius in mm
            arc_angle (float): Arc angle in degrees
            tilt (float): Tilt angle in degrees (rotation around Y-axis)
        """
        # Apply tilt transformation to orientation
        tilt_rad = np.deg2rad(tilt)
        tilt_matrix = np.array([
            [np.cos(tilt_rad), 0, -np.sin(tilt_rad)],
            [0, 1, 0],
            [np.sin(tilt_rad), 0, np.cos(tilt_rad)]
        ])
        tilted_orientation = tilt_matrix @ orientation

        # Cross-section dimensions change with tilt
        tilted_width = abs(np.cos(tilt_rad) * width + np.sin(tilt_rad) * height)
        tilted_height = abs(-np.sin(tilt_rad) * width + np.cos(tilt_rad) * height)

        super().__init__(current, start_pos, tilted_orientation, tilted_width, tilted_height, tilt)
        self.radius = radius
        self.arc_angle = arc_angle

        # Arc center: start position minus radius in X-direction (row vector)
        self.arc_center = self.start_pos - self.radius * self.orientation[0, :]

    @property
    def end_pos(self):
        """End position: arc_center + radius * rotated X-direction."""
        phi_rad = np.deg2rad(self.arc_angle)
        rotation_matrix = np.array([
            [np.cos(phi_rad), np.sin(phi_rad), 0],
            [-np.sin(phi_rad), np.cos(phi_rad), 0],
            [0, 0, 1]
        ])
        end_orientation = rotation_matrix @ self.orientation
        return self.arc_center + self.radius * end_orientation[0, :]

    @property
    def end_orientation(self):
        """End orientation: rotated by arc_angle around Z-axis."""
        phi_rad = np.deg2rad(self.arc_angle)
        rotation_matrix = np.array([
            [np.cos(phi_rad), np.sin(phi_rad), 0],
            [-np.sin(phi_rad), np.cos(phi_rad), 0],
            [0, 0, 1]
        ])
        return rotation_matrix @ self.orientation


class CoilBuilder:
    """
    Fluent builder interface for creating multi-segment coil paths.

    The builder maintains current state (position, orientation, cross-section)
    and automatically updates it after each segment is added. This eliminates
    manual state tracking and reduces boilerplate code by ~75%.

    Example:
        >>> builder = CoilBuilder(current=1265)
        >>> coil_radia_objects = (builder
        ...     .set_start([218, -16.4, -81])
        ...     .set_cross_section(width=122, height=122)
        ...     .add_straight(length=32.9, tilt=0)
        ...     .add_arc(radius=121, arc_angle=64.6, tilt=90)
        ...     .add_straight(length=1018.5, tilt=90)
        ...     .to_radia())
        >>>
        >>> import radia as rad
        >>> coils = rad.ObjCnt(coil_radia_objects)
    """

    def __init__(self, current):
        """
        Initialize coil builder.

        Args:
            current (float): Current in Amperes (constant for all segments)
        """
        self.current = current
        self.segments = []

        # Initial state (identity orientation at origin)
        self._position = np.array([0.0, 0.0, 0.0])
        self._orientation = np.eye(3)
        self._width = 100.0
        self._height = 100.0

    def set_start(self, position, orientation=None):
        """
        Set starting position and orientation.

        Args:
            position (array): Starting position [x, y, z] in mm
            orientation (array, optional): 3x3 orientation matrix (row vectors).
                                          Defaults to identity (aligned with XYZ axes).

        Returns:
            self (for method chaining)
        """
        self._position = np.array(position)
        if orientation is not None:
            self._orientation = np.array(orientation)
        return self

    def set_cross_section(self, width, height):
        """
        Set cross-section dimensions for subsequent segments.

        Args:
            width (float): Width in mm
            height (float): Height in mm

        Returns:
            self (for method chaining)
        """
        self._width = width
        self._height = height
        return self

    def add_straight(self, length, tilt=0):
        """
        Add a straight segment.

        Args:
            length (float): Length in mm
            tilt (float): Tilt angle in degrees (rotation around Y-axis)

        Returns:
            self (for method chaining)
        """
        segment = StraightSegment(
            self.current,
            self._position,
            self._orientation,
            self._width,
            self._height,
            length,
            tilt
        )
        self.segments.append(segment)

        # Automatic state update
        self._position = segment.end_pos
        self._orientation = segment.end_orientation
        self._width = segment.width
        self._height = segment.height

        return self

    def add_arc(self, radius, arc_angle, tilt=0):
        """
        Add an arc segment.

        Args:
            radius (float): Arc radius in mm
            arc_angle (float): Arc angle in degrees
            tilt (float): Tilt angle in degrees (rotation around Y-axis)

        Returns:
            self (for method chaining)
        """
        segment = ArcSegment(
            self.current,
            self._position,
            self._orientation,
            self._width,
            self._height,
            radius,
            arc_angle,
            tilt
        )
        self.segments.append(segment)

        # Automatic state update
        self._position = segment.end_pos
        self._orientation = segment.end_orientation
        self._width = segment.width
        self._height = segment.height

        return self

    def to_radia(self):
        """
        Convert all segments to Radia objects.

        Returns:
            list: List of Radia object IDs (can be combined with rad.ObjCnt)
        """
        import radia as rad

        radia_objects = []
        for seg in self.segments:
            if isinstance(seg, StraightSegment):
                # Create straight current segment
                J = [0, seg.current_density, 0]  # Current density in Y-direction
                coil = rad.ObjRecCur([0, 0, 0], [seg.width, seg.length, seg.height], J)

                # Build transformation (ZXZ Euler angles + translation)
                trf = rad.TrfRot([0, 0, 0], [0, 0, 1], np.deg2rad(seg.euler_angles[2]))
                trf = rad.TrfCmbR(trf, rad.TrfRot([0, 0, 0], [1, 0, 0], np.deg2rad(seg.euler_angles[1])))
                trf = rad.TrfCmbR(trf, rad.TrfRot([0, 0, 0], [0, 0, 1], np.deg2rad(seg.euler_angles[0])))
                trf = rad.TrfCmbL(trf, rad.TrfTrsl(seg.center.tolist()))

                radia_objects.append(rad.TrfOrnt(coil, trf))

            elif isinstance(seg, ArcSegment):
                # Create arc current segment
                phi1 = np.deg2rad(seg.euler_angles[0])
                if phi1 <= 0:
                    phi1 += 2 * np.pi

                phi2 = np.deg2rad(seg.euler_angles[0] + seg.arc_angle)
                if phi1 > phi2 or phi2 <= 0:
                    phi2 += 2 * np.pi

                r_range = [seg.radius - seg.width / 2, seg.radius + seg.width / 2]
                coil = rad.ObjArcCur(
                    [0, 0, 0],
                    r_range,
                    [phi1, phi2],
                    seg.height,
                    10,  # Number of segments
                    seg.current_density,
                    "auto"
                )

                # Build transformation (ZX Euler angles + translation to arc center)
                trf = rad.TrfRot([0, 0, 0], [0, 0, 1], np.deg2rad(seg.euler_angles[2]))
                trf = rad.TrfCmbR(trf, rad.TrfRot([0, 0, 0], [1, 0, 0], np.deg2rad(seg.euler_angles[1])))
                trf = rad.TrfCmbL(trf, rad.TrfTrsl(seg.arc_center.tolist()))

                radia_objects.append(rad.TrfOrnt(coil, trf))

        return radia_objects


# Export public API
__all__ = ['CoilBuilder', 'CoilSegment', 'StraightSegment', 'ArcSegment']
