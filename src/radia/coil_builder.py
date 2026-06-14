#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Radia Coil Builder - Fluent interface for constructing complex coil geometries.

Provides a beam-optics-inspired path builder for multi-segment coils.
Each segment's end state (position, orientation) automatically becomes
the next segment's start state, enabling continuous coil path definition
without manual coordinate tracking.

Outputs:
  - to_radia(): Radia field source objects (ObjRecCur, ObjArcCur)
  - to_occ(): OCC shape for STEP export / GMSH visualization
  - write_step(): Direct STEP file export

Example:
	>>> from coil_builder import CoilBuilder
	>>>
	>>> mm = 1e-3
	>>> coil = (CoilBuilder(current=1000)
	...	 .set_start([0, 0, 0])
	...	 .set_cross_section(width=20*mm, height=20*mm)
	...	 .add_straight(100*mm)
	...	 .add_arc(radius=50*mm, arc_angle=180, tilt=90)
	...	 .add_straight(100*mm)
	...	 .add_arc(radius=50*mm, arc_angle=180, tilt=90))
	>>>
	>>> radia_objects = coil.to_radia()      # Radia field sources
	>>> occ_shape = coil.to_occ()            # OCC shape for visualization
	>>> coil.write_step("racetrack.step")    # STEP export
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

	def __init__(self, current, start_pos, orientation, width, height,
	             tilt=0, profile=None):
		"""
		Initialize coil segment.

		Args:
			current (float): Current in Amperes
			start_pos (array): Starting position [x, y, z] (in active Radia units)
			orientation (array): 3x3 orientation matrix (row vectors)
			width (float): Cross-section width (in active Radia units)
			height (float): Cross-section height (in active Radia units)
			tilt (float): Tilt angle in degrees (not applied here, applied in subclass)
			profile (Profile, optional): Explicit cross-section Profile. If
				given, replaces the implicit RectProfile(width, height) used
				for filament sampling. width/height are still stored for
				legacy consumers (to_radia, to_occ) and reporting.
		"""
		self.current = current
		self.start_pos = np.array(start_pos)
		self.orientation = np.array(orientation)
		self.width = width
		self.height = height

		# Cross-section profile: explicit if given, else RectProfile fallback
		# so StraightSegment/ArcSegment with the legacy (w, h) API continues
		# to work.
		if profile is None:
			from radia.coil_profile import RectProfile
			profile = RectProfile(width, height)
		self.profile = profile

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
		"""Current density: I / (width * height) in A/(active_unit²)."""
		return self.current / (self.width * self.height)

	@abstractmethod
	def to_occ_shape(self):
		"""Generate OCC shape for this segment."""
		pass


class StraightSegment(CoilSegment):
	"""
	Straight coil segment with optional tilt.

	The tilt rotates the cross-section around the Y-axis before
	extending in the local Y direction.
	"""

	def __init__(self, current, start_pos, orientation, width, height,
	             length, tilt=0, profile=None):
		"""
		Initialize straight segment.

		Args:
			current (float): Current in Amperes
			start_pos (array): Starting position [x, y, z] (in active Radia units)
			orientation (array): 3x3 orientation matrix (row vectors)
			width (float): Cross-section width (in active Radia units)
			height (float): Cross-section height (in active Radia units)
			length (float): Segment length (in active Radia units)
			tilt (float): Tilt angle in degrees (rotation around Y-axis)
			profile (Profile, optional): explicit cross-section profile.
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

		super().__init__(current, start_pos, tilted_orientation,
		                 tilted_width, tilted_height, tilt, profile=profile)
		self.length = length

	@property
	def end_pos(self):
		"""End position: start + length * Y-direction."""
		return self.start_pos + self.length * self.orientation[1, :]

	@property
	def end_orientation(self):
		"""End orientation: same as start (no rotation)."""
		return self.orientation

	def to_occ_shape(self, index=0):
		"""Generate OCC Box for this straight segment."""
		from netgen.occ import Box, Pnt, Axis, Vec, Z, X, Y
		# Create box at origin aligned with XYZ
		shape = Box(Pnt(-self.width/2, 0, -self.height/2),
		            Pnt(self.width/2, self.length, self.height/2))
		# Apply ZXZ Euler rotation + translation to match segment pose
		ea = self.euler_angles
		shape = shape.Rotate(Axis(Pnt(0,0,0), Z), ea[2])
		shape = shape.Rotate(Axis(Pnt(0,0,0), X), ea[1])
		shape = shape.Rotate(Axis(Pnt(0,0,0), Z), ea[0])
		shape = shape.Move(Vec(*self.start_pos))
		shape.name = "coil_straight_" + str(index)
		return shape


class LoftStraightSegment(CoilSegment):
	"""
	Straight segment with a cross-section that lofts (linear interpolation)
	from profile_start at the segment start to profile_end at the end.

	Only Stage 3a feature: profiles must be of compatible type
	(e.g. both RectProfile). Cross-type lofts (Stage 3b) are not yet
	supported.

	Filament paths through a loft segment are NOT straight lines in 3D --
	each filament's (u, v) offset from the centerline varies with path
	parameter s in [0, 1] as the profile morphs, so the filament itself
	converges/diverges with the cross-section.
	"""

	def __init__(self, current, start_pos, orientation,
	             profile_start, profile_end, length, n_sub=20, tilt=0):
		"""
		Args:
			profile_start: Profile instance at segment start (s = 0).
			profile_end: Profile instance at segment end (s = 1).
			length: straight-path length in Radia units.
			n_sub: number of sub-segments along the path for
			       filament discretization. Larger -> better
			       resolution of cross-section variation.
			tilt: Y-axis rotation (same convention as StraightSegment).
		"""
		tilt_rad = np.deg2rad(tilt)
		tilt_matrix = np.array([
			[np.cos(tilt_rad), 0, -np.sin(tilt_rad)],
			[0, 1, 0],
			[np.sin(tilt_rad), 0, np.cos(tilt_rad)],
		])
		tilted_orientation = tilt_matrix @ orientation

		# For the base class's width/height, use the END profile's bounding
		# box so CoilBuilder._width/_height follow the cross-section at the
		# next segment's start.
		w_end, h_end = profile_end.bounding_wh()
		super().__init__(current, start_pos, tilted_orientation,
		                 w_end, h_end, tilt)
		self.profile_start = profile_start
		self.profile_end = profile_end
		self.length = length
		self.n_sub = int(n_sub)

	@property
	def end_pos(self):
		return self.start_pos + self.length * self.orientation[1, :]

	@property
	def end_orientation(self):
		return self.orientation

	def profile_at(self, s):
		"""Interpolated profile at normalized path parameter s in [0, 1]."""
		return self.profile_start.interpolate(self.profile_end, s)

	def to_occ_shape(self, index=0):
		"""Loft OCC shape from profile_start at s=0 to profile_end at s=1.

		Stage 3a implementation supports only RectProfile -> RectProfile
		and emits a ThruSections loft between the two rectangle wires.
		The resulting Solid is usable for STEP export and visual
		inspection.
		"""
		from netgen.occ import (WorkPlane, Axes, Pnt, Axis, X, Y, Z, Vec,
		                         ThruSections)
		from radia.coil_profile import RectProfile
		if not (isinstance(self.profile_start, RectProfile)
		        and isinstance(self.profile_end, RectProfile)):
			raise NotImplementedError(
				"LoftStraightSegment.to_occ_shape is only implemented for "
				"RectProfile -> RectProfile (Stage 3a).")
		w0, h0 = self.profile_start.w, self.profile_start.h
		w1, h1 = self.profile_end.w, self.profile_end.h
		# ThruSections requires Wires (boundary curves), not Faces.
		wire0 = WorkPlane(Axes(Pnt(0, 0, 0), n=Y, h=X)) \
		        .Rectangle(w0, h0).Wire()
		wire1 = WorkPlane(Axes(Pnt(0, self.length, 0), n=Y, h=X)) \
		        .Rectangle(w1, h1).Wire()
		shape = ThruSections([wire0, wire1], solid=True)
		ea = self.euler_angles
		shape = shape.Rotate(Axis(Pnt(0, 0, 0), Z), ea[2])
		shape = shape.Rotate(Axis(Pnt(0, 0, 0), X), ea[1])
		shape = shape.Rotate(Axis(Pnt(0, 0, 0), Z), ea[0])
		shape = shape.Move(Vec(*self.start_pos))
		shape.name = "coil_loft_" + str(index)
		return shape


class LoftArcSegment(CoilSegment):
	"""
	Arc segment with a cross-section that lofts (linear interpolation
	in canonical (alpha, beta) parametrization) from profile_start at
	the arc entry to profile_end at the arc exit.

	Combines the geometric path of ArcSegment (rotation around arc_center
	in the local XY plane) with the cross-section variation of
	LoftStraightSegment. Each filament k at (alpha_k, beta_k) has radial
	offset u_k(s) and axial offset v_k(s) that interpolate smoothly
	between profile_start and profile_end as the arc parameter s goes
	from 0 to 1.

	Useful for IH turn geometries where the conductor cross-section
	morphs through a bend (e.g. round wire transitioning to a flat strap
	on a pancake coil turn).
	"""

	def __init__(self, current, start_pos, orientation, profile_start,
	             profile_end, radius, arc_angle, n_sub=20, tilt=0):
		"""
		Args:
			profile_start: Profile at arc entry (s = 0).
			profile_end: Profile at arc exit (s = 1).
			radius: arc center-line radius.
			arc_angle: arc angle [deg] (signed).
			n_sub: sub-segment count along the arc for filament
			       discretization (larger -> smoother cross-section
			       morph).
			tilt: Y-axis tilt [deg] (same convention as ArcSegment).
		"""
		tilt_rad = np.deg2rad(tilt)
		tilt_matrix = np.array([
			[np.cos(tilt_rad), 0, -np.sin(tilt_rad)],
			[0, 1, 0],
			[np.sin(tilt_rad), 0, np.cos(tilt_rad)],
		])
		tilted_orientation = tilt_matrix @ orientation

		# End-profile bounding (w, h) so CoilBuilder state follows.
		w_end, h_end = profile_end.bounding_wh()
		super().__init__(current, start_pos, tilted_orientation,
		                 w_end, h_end, tilt, profile=None)
		self.profile_start = profile_start
		self.profile_end = profile_end
		self.radius = radius
		self.arc_angle = arc_angle
		self.n_sub = int(n_sub)

		# Arc center: same construction as ArcSegment.
		self.arc_center = self.start_pos - self.radius * self.orientation[0, :]

	@property
	def end_pos(self):
		phi_rad = np.deg2rad(self.arc_angle)
		rotation_matrix = np.array([
			[np.cos(phi_rad), np.sin(phi_rad), 0],
			[-np.sin(phi_rad), np.cos(phi_rad), 0],
			[0, 0, 1],
		])
		end_orientation = rotation_matrix @ self.orientation
		return self.arc_center + self.radius * end_orientation[0, :]

	@property
	def end_orientation(self):
		phi_rad = np.deg2rad(self.arc_angle)
		rotation_matrix = np.array([
			[np.cos(phi_rad), np.sin(phi_rad), 0],
			[-np.sin(phi_rad), np.cos(phi_rad), 0],
			[0, 0, 1],
		])
		return rotation_matrix @ self.orientation

	def profile_at(self, s):
		return self.profile_start.interpolate(self.profile_end, s)

	def to_occ_shape(self, index=0):
		"""OCC shape for LoftArcSegment is not yet implemented. For
		visualization / STEP export, use multiple LoftStraightSegment
		pieces along the arc, or use build123d's sweep/loft directly."""
		raise NotImplementedError(
			"LoftArcSegment.to_occ_shape: not implemented yet. "
			"Use LoftStraightSegment chains for CAD export.")


class ArcSegment(CoilSegment):
	"""
	Arc coil segment with optional tilt.

	The arc rotates around a center point in the local XY plane.
	Tilt is applied first, then the arc rotation.
	"""

	def __init__(self, current, start_pos, orientation, width, height,
	             radius, arc_angle, tilt=0, profile=None):
		"""
		Initialize arc segment.

		Args:
			current (float): Current in Amperes
			start_pos (array): Starting position [x, y, z] (in active Radia units)
			orientation (array): 3x3 orientation matrix (row vectors)
			width (float): Cross-section width (in active Radia units)
			height (float): Cross-section height (in active Radia units)
			radius (float): Arc radius (in active Radia units)
			arc_angle (float): Arc angle in degrees
			tilt (float): Tilt angle in degrees (rotation around Y-axis)
			profile (Profile, optional): explicit cross-section profile.
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

		super().__init__(current, start_pos, tilted_orientation,
		                 tilted_width, tilted_height, tilt, profile=profile)
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

	def to_occ_shape(self, index=0):
		"""Generate OCC revolved shape for this arc segment."""
		from netgen.occ import WorkPlane, Axes, Pnt, Axis, Vec, Z, X, Y
		# Cross-section rectangle at (radius - width/2, 0) in local frame.
		# radius is center-line radius (distance from arc center to
		# cross-section center). Inner edge at radius - width/2.
		r_inner = self.radius - self.width / 2
		if r_inner < 1e-10:
			raise ValueError(
				f"Arc inner radius {r_inner} <= 0: radius={self.radius}, "
				f"width={self.width}. Use radius >= width/2 (center-line radius).")
		wp = WorkPlane(Axes(Pnt(r_inner, 0, -self.height/2), n=-Y, h=X))
		face = wp.Rectangle(self.width, self.height).Face()
		# Revolve around Z axis
		shape = face.Revolve(Axis(Pnt(0, 0, 0), Z), self.arc_angle)
		# Apply ZXZ Euler rotation + translation to arc center
		ea = self.euler_angles
		shape = shape.Rotate(Axis(Pnt(0, 0, 0), Z), ea[2])
		shape = shape.Rotate(Axis(Pnt(0, 0, 0), X), ea[1])
		shape = shape.Rotate(Axis(Pnt(0, 0, 0), Z), ea[0])
		shape = shape.Move(Vec(*self.arc_center))
		shape.name = "coil_arc_" + str(index)
		return shape


class CoilBuilder:
	"""
	Fluent builder interface for creating multi-segment coil paths.

	The builder maintains current state (position, orientation, cross-section)
	and automatically updates it after each segment is added. This eliminates
	manual state tracking and reduces boilerplate code by ~75%.

	Example:
		>>> mm = 1e-3  # unit conversion factor
		>>> builder = CoilBuilder(current=1265)
		>>> coil_radia_objects = (builder
		...	 .set_start([218*mm, -16.4*mm, -81*mm])
		...	 .set_cross_section(width=122*mm, height=122*mm)
		...	 .add_straight(length=32.9*mm, tilt=0)
		...	 .add_arc(radius=121*mm, arc_angle=64.6, tilt=90)
		...	 .add_straight(length=1018.5*mm, tilt=90)
		...	 .to_radia())
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
		self._width = None
		self._height = None
		# Optional explicit profile, takes precedence over (_width, _height)
		# when a segment is added. Set via set_profile().
		self._profile = None

	def set_start(self, position, orientation=None):
		"""
		Set starting position and orientation.

		Args:
			position (array): Starting position [x, y, z] (in active Radia units)
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
		Set rectangular cross-section dimensions for subsequent segments.

		Clears any explicit profile previously set via set_profile().

		Args:
			width (float): Width (in active Radia units)
			height (float): Height (in active Radia units)

		Returns:
			self (for method chaining)
		"""
		self._width = width
		self._height = height
		self._profile = None
		return self

	def set_profile(self, profile):
		"""Set an explicit cross-section Profile for subsequent segments.

		Overrides set_cross_section(). Profile bounding (w, h) still
		populates _width / _height for legacy consumers, but filament
		sampling in to_filaments uses profile.sample_at() directly.

		Args:
			profile: Profile instance (RectProfile, CircleProfile,
				AnnularProfile, etc.)

		Returns:
			self (for method chaining)
		"""
		self._profile = profile
		w, h = profile.bounding_wh()
		self._width = w
		self._height = h
		return self

	def _check_cross_section(self):
		"""Raise ValueError if cross-section has not been set."""
		if self._width is None or self._height is None:
			raise ValueError(
				"Cross-section not set. Call set_cross_section(width, height) "
				"before adding segments."
			)

	def add_straight(self, length, tilt=0, profile=None):
		"""
		Add a straight segment.

		Args:
			length (float): Length (in active Radia units)
			tilt (float): Tilt angle in degrees (rotation around Y-axis)
			profile (Profile, optional): explicit cross-section profile.
				If None, uses self._profile (from set_profile) or
				RectProfile(self._width, self._height).

		Returns:
			self (for method chaining)

		Raises:
			ValueError: If set_cross_section() has not been called.
		"""
		self._check_cross_section()
		if profile is None:
			profile = self._profile  # may still be None -> CoilSegment
			                          # builds RectProfile(w, h)
		segment = StraightSegment(
			self.current,
			self._position,
			self._orientation,
			self._width,
			self._height,
			length,
			tilt,
			profile=profile,
		)
		self.segments.append(segment)

		# Automatic state update
		self._position = segment.end_pos
		self._orientation = segment.end_orientation
		self._width = segment.width
		self._height = segment.height
		# Carry forward the explicit profile (if any) so subsequent
		# add_straight / add_arc without explicit profile= inherit it.
		if profile is not None:
			self._profile = profile

		return self

	def add_loft_straight(self, profile_end, length, n_sub=20, tilt=0,
	                     profile_start=None):
		"""Add a straight segment with a lofting cross-section.

		The cross-section interpolates linearly from profile_start at
		s=0 to profile_end at s=1. profile_start defaults to the
		current builder cross-section (a RectProfile constructed from
		the builder's _width, _height).

		Args:
			profile_end: Profile at the segment end. Same type as
				profile_start (Stage 3a: RectProfile only).
			length: straight-path length [active Radia units].
			n_sub: number of sub-segments along the path for filament
				discretization. Default 20.
			tilt: Y-axis tilt [deg].
			profile_start: explicit Profile at segment start. If None,
				uses RectProfile(self._width, self._height) matching
				the builder's current cross-section.
		Returns:
			self (for chaining).
		"""
		from radia.coil_profile import RectProfile
		if profile_start is None:
			self._check_cross_section()
			profile_start = RectProfile(self._width, self._height)

		segment = LoftStraightSegment(
			self.current, self._position, self._orientation,
			profile_start, profile_end, length, n_sub=n_sub, tilt=tilt,
		)
		self.segments.append(segment)

		self._position = segment.end_pos
		self._orientation = segment.end_orientation
		# Update builder cross-section to end profile (so the next
		# add_straight / add_arc inherits the new size).
		w_end, h_end = profile_end.bounding_wh()
		self._width = w_end
		self._height = h_end

		return self

	def add_loft_arc(self, profile_end, radius, arc_angle, n_sub=20, tilt=0,
	                  profile_start=None):
		"""Add an arc segment with a lofting cross-section.

		The cross-section interpolates linearly in canonical (alpha,
		beta) space from profile_start at s=0 (arc entry) to
		profile_end at s=1 (arc exit).

		Args:
			profile_end: Profile at arc exit. Interpolable with
				profile_start (same-type for Stage 3a rect-rect
				loft; cross-type via InterpolatedProfile otherwise).
			radius: arc center-line radius.
			arc_angle: arc angle [deg].
			n_sub: sub-segment count along the arc. Default 20.
			tilt: Y-axis tilt [deg].
			profile_start: Profile at arc entry. Defaults to current
				builder profile / rect fallback.

		Returns:
			self (for chaining).
		"""
		if profile_start is None:
			self._check_cross_section()
			if self._profile is not None:
				profile_start = self._profile
			else:
				from radia.coil_profile import RectProfile
				profile_start = RectProfile(self._width, self._height)

		segment = LoftArcSegment(
			self.current, self._position, self._orientation,
			profile_start, profile_end, radius, arc_angle,
			n_sub=n_sub, tilt=tilt,
		)
		self.segments.append(segment)

		self._position = segment.end_pos
		self._orientation = segment.end_orientation
		w_end, h_end = profile_end.bounding_wh()
		self._width = w_end
		self._height = h_end
		self._profile = profile_end

		return self

	def add_arc(self, radius, arc_angle, tilt=0, profile=None):
		"""
		Add an arc segment.

		Args:
			radius (float): Arc radius (in active Radia units)
			arc_angle (float): Arc angle in degrees
			tilt (float): Tilt angle in degrees (rotation around Y-axis)
			profile (Profile, optional): explicit cross-section profile.
				If None, uses self._profile or the rectangular fallback.

		Returns:
			self (for method chaining)

		Raises:
			ValueError: If set_cross_section() has not been called.
		"""
		self._check_cross_section()
		if profile is None:
			profile = self._profile
		segment = ArcSegment(
			self.current,
			self._position,
			self._orientation,
			self._width,
			self._height,
			radius,
			arc_angle,
			tilt,
			profile=profile,
		)
		self.segments.append(segment)

		# Automatic state update
		self._position = segment.end_pos
		self._orientation = segment.end_orientation
		self._width = segment.width
		self._height = segment.height
		if profile is not None:
			self._profile = profile

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
				phi2 = np.deg2rad(seg.euler_angles[0] + seg.arc_angle)
				j_density = seg.current_density

				# Radia requires 0 <= phi1 < phi2 <= 2*pi
				# For negative arc_angle: swap phi1/phi2 and negate current
				if phi1 > phi2:
					phi1, phi2 = phi2, phi1
					j_density = -j_density

				# Normalize to [0, 2*pi]
				while phi1 < 0:
					phi1 += 2 * np.pi
					phi2 += 2 * np.pi
				while phi1 >= 2 * np.pi:
					phi1 -= 2 * np.pi
					phi2 -= 2 * np.pi
				if phi2 <= 0:
					phi2 += 2 * np.pi

				r_range = [seg.radius - seg.width / 2, seg.radius + seg.width / 2]
				coil = rad.ObjArcCur(
					[0, 0, 0],       # center
					r_range,          # radii [r_min, r_max]
					[phi1, phi2],     # phi range
					seg.height,       # height
					10,               # nseg
					"auto",           # man_auto
					"z",              # axis (transformed by Euler angles below)
					j_density         # j (current density, sign handles direction)
				)

				# Build transformation (ZX Euler angles + translation to arc center)
				trf = rad.TrfRot([0, 0, 0], [0, 0, 1], np.deg2rad(seg.euler_angles[2]))
				trf = rad.TrfCmbR(trf, rad.TrfRot([0, 0, 0], [1, 0, 0], np.deg2rad(seg.euler_angles[1])))
				trf = rad.TrfCmbL(trf, rad.TrfTrsl(seg.arc_center.tolist()))

				radia_objects.append(rad.TrfOrnt(coil, trf))

		return radia_objects


	def to_occ(self):
		"""Convert all segments to a combined OCC shape.

		Returns:
			OCC shape (can be exported to STEP, displayed in GMSH, etc.)
		"""
		from netgen.occ import Glue
		shapes = [seg.to_occ_shape(i) for i, seg in enumerate(self.segments)]
		if len(shapes) == 0:
			raise ValueError("No segments added")
		if len(shapes) == 1:
			return shapes[0]
		return Glue(shapes)

	def write_step(self, filename):
		"""Export coil geometry as a SWEPT SOLID STEP file.

		Produces the swept rectangular cross-section glued across all
		segments.  The PEEC pipeline consumes this STEP (walker-based
		centerline + cross-section subdivision via nwinc x nhinc).

		Args:
			filename: Output .step file path

		Returns:
			filename
		"""
		shape = self.to_occ()
		shape.WriteStep(filename)
		return filename

	# ============================================================
	# Wire segment extraction (for Biot-Savart)
	# ============================================================

	def to_wire_segments(self, n_arc=20):
		"""Extract centerline wire segments for Biot-Savart field computation.

		Straight segments become a single wire segment (start -> end).
		Arc segments are discretized into n_arc straight sub-segments.

		Args:
			n_arc: Number of straight segments per arc (default 20)

		Returns:
			(segments, current) where:
				segments: list of ((x1,y1,z1), (x2,y2,z2)) endpoint pairs [m]
				current: coil current [A]
		"""
		wire_segments = []
		for seg in self.segments:
			if isinstance(seg, StraightSegment):
				wire_segments.append((tuple(seg.start_pos), tuple(seg.end_pos)))
			elif isinstance(seg, ArcSegment):
				angle_rad = np.deg2rad(seg.arc_angle)
				for i in range(n_arc):
					t1 = angle_rad * i / n_arc
					t2 = angle_rad * (i + 1) / n_arc
					p1 = seg.arc_center + seg.radius * (
						np.cos(t1) * seg.orientation[0, :] +
						np.sin(t1) * seg.orientation[1, :])
					p2 = seg.arc_center + seg.radius * (
						np.cos(t2) * seg.orientation[0, :] +
						np.sin(t2) * seg.orientation[1, :])
					wire_segments.append((tuple(p1), tuple(p2)))
		return wire_segments, self.current

	# ============================================================
	# Filament approximation for AC analysis (Tier A: analytical)
	# ============================================================

	def to_filaments(self, nw, nh, frequency=0.0, sigma=5.8e7,
	                 mu_r=1.0, n_arc=20):
		"""Generate filament paths with complex currents for AC analysis.

		Divides the rectangular cross-section into nw x nh filaments.
		Each filament follows the coil centerline offset by (u, v) in the
		local cross-section frame. The filament current is computed from
		an analytical "Tier A" model:

		    I_k = I_total * w_k / sum(w_j)
		    w_k = (1 / R_k) * exp(-(1+j) d_k / delta)

		where
		    R_k    = rho * L_k / A_cell     (DC resistance of filament k)
		    L_k    = total path length of filament k (inner arcs shorter)
		    A_cell = (width/nw) * (height/nh)  (cross-section cell area)
		    d_k    = L-inf distance from filament (u_k, v_k) to conductor
		             surface
		    delta  = sqrt(2 / (omega * mu * sigma))   (skin depth)

		Captures:
		  - "Circulating current" = path-length bias on arcs (1/R_k)
		  - A-priori skin effect via exp(-(1+j) d/delta)

		Does NOT capture:
		  - Proximity effect (coupling between filaments or to external
		    conductors/workpieces). For workpiece-less coil design this
		    approximation typically agrees well with BEM; add a workpiece
		    and the proximity effect grows, which is where this model breaks
		    down and a full PEEC / BEM solve is needed.

		Args:
		    nw (int): Number of filaments across width (>= 1)
		    nh (int): Number of filaments across height (>= 1)
		    frequency (float): Frequency [Hz]. Zero -> DC (real currents,
		                       no skin factor; only path-length bias on
		                       arc segments).
		    sigma (float): Conductivity [S/m]. Default 5.8e7 (Cu).
		    mu_r (float): Relative permeability. Default 1.0.
		    n_arc (int): Polyline sub-segments per arc (default 20).

		Returns:
		    (filament_paths, currents)
		      filament_paths : list of length N = nw*nh. Each element is a
		                       list of ((x1,y1,z1), (x2,y2,z2)) polyline
		                       segment endpoint pairs [same length unit as
		                       the builder].
		      currents       : ndarray (N,), complex128 if frequency>0 else
		                       float64. Sum over k equals self.current.

		Notes:
		    Assumes the cross-section is uniform across segments (takes
		    width/height from segments[0]). For tilted segments the local
		    orientation rotates with the segment; the filament indexing
		    (u, v) is consistent across segments because orientation is
		    continuous in CoilBuilder.
		"""
		if not self.segments:
			raise ValueError("No segments to filament-ize.")
		if nw < 1 or nh < 1:
			raise ValueError(f"nw={nw}, nh={nh}; both must be >= 1.")

		MU_0 = 4.0 * np.pi * 1e-7

		# Reference (w0, h0) from the FIRST segment's profile bounding box.
		# Used for legacy _last_filament_info reporting and the Tier A skin
		# factor fallback (approximate for non-rect profiles; PEEC is the
		# AC-truth path anyway).
		w0, h0 = self.segments[0].profile.bounding_wh()

		# Canonical (alpha, beta) in [0, 1]^2, invariant across all segments
		# and across profile types / lofts.
		i_idx, j_idx = np.meshgrid(np.arange(nw), np.arange(nh),
		                           indexing='ij')
		alphas = ((i_idx + 0.5) / nw).ravel()
		betas = ((j_idx + 0.5) / nh).ravel()
		N = alphas.size

		# Reference (u0, v0) sampled on the FIRST segment's profile, for
		# Tier A skin factor (bounding-box distance approximation) and
		# legacy reporting.
		us0, vs0 = self.segments[0].profile.sample_at(alphas, betas)

		filament_paths = [[] for _ in range(N)]
		cell_wh = [[] for _ in range(N)]
		lengths = np.zeros(N)
		R_k = np.zeros(N)

		def _check_arc_inner_radius(R_arc, us):
			r_min = R_arc + float(np.min(us))
			if r_min <= 0:
				raise ValueError(
					f"Minimum filament arc radius {r_min} <= 0: "
					f"R_arc={R_arc}, min(u)={np.min(us):g}. Reduce "
					"profile size or increase arc radius.")

		for seg in self.segments:
			x_hat = seg.orientation[0, :]
			y_hat = seg.orientation[1, :]
			z_hat = seg.orientation[2, :]

			if isinstance(seg, LoftStraightSegment):
				n_sub = seg.n_sub
				s_waypoints = np.linspace(0.0, 1.0, n_sub + 1)
				profile_at_s = [seg.profile_start.interpolate(seg.profile_end, s)
				                for s in s_waypoints]
				uv_wp = [pr.sample_at(alphas, betas) for pr in profile_at_s]
				for p_idx in range(n_sub):
					s1 = s_waypoints[p_idx]
					s2 = s_waypoints[p_idx + 1]
					u1, v1 = uv_wp[p_idx]
					u2, v2 = uv_wp[p_idx + 1]
					pr_mid = seg.profile_start.interpolate(
						seg.profile_end, 0.5 * (s1 + s2))
					dw_mid, dh_mid = pr_mid.cell_wh(nw, nh)
					A_mid = pr_mid.total_area() / (nw * nh)
					for k in range(N):
						p1 = (seg.start_pos + s1 * seg.length * y_hat
						      + u1[k] * x_hat + v1[k] * z_hat)
						p2 = (seg.start_pos + s2 * seg.length * y_hat
						      + u2[k] * x_hat + v2[k] * z_hat)
						filament_paths[k].append((tuple(p1), tuple(p2)))
						cell_wh[k].append((dw_mid, dh_mid))
						dl = float(np.linalg.norm(p2 - p1))
						lengths[k] += dl
						R_k[k] += dl / (sigma * A_mid)

			elif isinstance(seg, StraightSegment):
				us_s, vs_s = seg.profile.sample_at(alphas, betas)
				du_s, dv_s = seg.profile.cell_wh(nw, nh)
				A_cell_s = seg.profile.total_area() / (nw * nh)
				for k in range(N):
					offset = us_s[k] * x_hat + vs_s[k] * z_hat
					p1 = seg.start_pos + offset
					p2 = seg.end_pos + offset
					filament_paths[k].append((tuple(p1), tuple(p2)))
					cell_wh[k].append((du_s, dv_s))
					lengths[k] += seg.length
					R_k[k] += seg.length / (sigma * A_cell_s)

			elif isinstance(seg, ArcSegment):
				theta_total = np.deg2rad(seg.arc_angle)
				R_arc = seg.radius
				us_s, vs_s = seg.profile.sample_at(alphas, betas)
				du_s, dv_s = seg.profile.cell_wh(nw, nh)
				A_cell_s = seg.profile.total_area() / (nw * nh)
				_check_arc_inner_radius(R_arc, us_s)
				for k in range(N):
					u, v = us_s[k], vs_s[k]
					r_k = R_arc + u
					z_off = v * z_hat
					for p_idx in range(n_arc):
						t1 = theta_total * p_idx / n_arc
						t2 = theta_total * (p_idx + 1) / n_arc
						p1 = (seg.arc_center + r_k * (np.cos(t1) * x_hat
						      + np.sin(t1) * y_hat) + z_off)
						p2 = (seg.arc_center + r_k * (np.cos(t2) * x_hat
						      + np.sin(t2) * y_hat) + z_off)
						filament_paths[k].append((tuple(p1), tuple(p2)))
						cell_wh[k].append((du_s, dv_s))
					arc_len_k = abs(r_k * theta_total)
					lengths[k] += arc_len_k
					R_k[k] += arc_len_k / (sigma * A_cell_s)

			elif isinstance(seg, LoftArcSegment):
				theta_total = np.deg2rad(seg.arc_angle)
				R_arc = seg.radius
				n_sub = seg.n_sub
				s_waypoints = np.linspace(0.0, 1.0, n_sub + 1)
				profile_at_s = [seg.profile_start.interpolate(seg.profile_end, s)
				                for s in s_waypoints]
				uv_wp = [pr.sample_at(alphas, betas) for pr in profile_at_s]
				# Check inner radius across all s, all k.
				for pr_idx, (u_s, _) in enumerate(uv_wp):
					_check_arc_inner_radius(R_arc, u_s)
				for p_idx in range(n_sub):
					s1 = s_waypoints[p_idx]
					s2 = s_waypoints[p_idx + 1]
					theta1 = theta_total * s1
					theta2 = theta_total * s2
					u1, v1 = uv_wp[p_idx]
					u2, v2 = uv_wp[p_idx + 1]
					pr_mid = seg.profile_start.interpolate(
						seg.profile_end, 0.5 * (s1 + s2))
					dw_mid, dh_mid = pr_mid.cell_wh(nw, nh)
					A_mid = pr_mid.total_area() / (nw * nh)
					for k in range(N):
						r_k_1 = R_arc + u1[k]
						r_k_2 = R_arc + u2[k]
						p1 = (seg.arc_center
						      + r_k_1 * (np.cos(theta1) * x_hat
						                 + np.sin(theta1) * y_hat)
						      + v1[k] * z_hat)
						p2 = (seg.arc_center
						      + r_k_2 * (np.cos(theta2) * x_hat
						                 + np.sin(theta2) * y_hat)
						      + v2[k] * z_hat)
						filament_paths[k].append((tuple(p1), tuple(p2)))
						cell_wh[k].append((dw_mid, dh_mid))
						dl = float(np.linalg.norm(p2 - p1))
						lengths[k] += dl
						R_k[k] += dl / (sigma * A_mid)
			else:
				raise NotImplementedError(
					f"Unsupported segment type {type(seg).__name__}")

		# Tier A current weights. R_k is exact (accumulated from per-piece
		# A_cell(s) integrals) for any profile / loft. Skin factor uses a
		# bounding-box distance approximation on the FIRST segment's
		# profile -- coarse for circles (overestimates d near the corners
		# of the bounding box, so skin factor is slightly too large), but
		# Tier A is only a seed; PEEC is the AC-truth path.
		if frequency > 0:
			omega = 2.0 * np.pi * frequency
			delta = np.sqrt(2.0 / (omega * MU_0 * mu_r * sigma))
			d_k = np.minimum(w0 / 2 - np.abs(us0), h0 / 2 - np.abs(vs0))
			d_k = np.maximum(d_k, 0.0)
			skin_factor = np.exp(-(1 + 1j) * d_k / delta)
			weights = (1.0 / R_k) * skin_factor
		else:
			weights = 1.0 / R_k

		currents = self.current * weights / np.sum(weights)

		self._last_filament_info = {
			'us': us0, 'vs': vs0,
			'alphas': alphas, 'betas': betas,
			'nw': nw, 'nh': nh,
			'w': w0, 'h': h0,
			'lengths': lengths,
			'R_k': R_k,
			'cell_wh': cell_wh,
			'currents': currents,
			'frequency': frequency,
			'sigma': sigma,
			'mu_r': mu_r,
			'has_loft': any(isinstance(s, (LoftStraightSegment,
			                                 LoftArcSegment))
			                for s in self.segments),
		}
		return filament_paths, currents

	def to_filaments_peri(self, n_peri, frequency=0.0, sigma=5.8e7,
	                      mu_r=1.0, n_arc=20):
		"""Perimeter-only filament placement (thin-skin limit, IH typical).

		Places ``n_peri`` filaments equally spaced along the arc-length
		perimeter of the cross-section.  Each filament carries an equal
		share of the cross-section area (A_cell = total_area / n_peri),
		so the bundle's DC resistance matches the solid conductor.

		Appropriate regime:
		  - Thin-skin: d / delta >= 3 (skin current localised near the
		    outer surface).
		  - For DC or thick-skin (d / delta < 2) this model is inaccurate
		    because it forbids interior current paths.  Use the
		    volume-grid ``to_filaments(nw, nh, ...)`` instead.

		Args:
		    n_peri (int): Number of filaments on the perimeter (>= 1).
		    frequency, sigma, mu_r, n_arc: same as ``to_filaments``.

		Returns:
		    (filament_paths, currents) — list of length ``n_peri``,
		    complex128 if frequency>0 else float64.
		"""
		if not self.segments:
			raise ValueError("No segments to filament-ize.")
		n_peri = int(n_peri)
		if n_peri < 1:
			raise ValueError(f"n_peri={n_peri}; must be >= 1.")

		MU_0 = 4.0 * np.pi * 1e-7

		w0, h0 = self.segments[0].profile.bounding_wh()

		us0, vs0 = self.segments[0].profile.sample_perimeter(n_peri)
		N = n_peri

		filament_paths = [[] for _ in range(N)]
		cell_wh = [[] for _ in range(N)]
		lengths = np.zeros(N)
		R_k = np.zeros(N)

		def _check_arc_inner_radius(R_arc, us):
			r_min = R_arc + float(np.min(us))
			if r_min <= 0:
				raise ValueError(
					f"Minimum filament arc radius {r_min} <= 0: "
					f"R_arc={R_arc}, min(u)={np.min(us):g}. Reduce "
					"profile size or increase arc radius.")

		def _cell_side(profile):
			side = float(np.sqrt(profile.total_area() / n_peri))
			return side, side

		for seg in self.segments:
			x_hat = seg.orientation[0, :]
			y_hat = seg.orientation[1, :]
			z_hat = seg.orientation[2, :]

			if isinstance(seg, LoftStraightSegment):
				n_sub = seg.n_sub
				s_waypoints = np.linspace(0.0, 1.0, n_sub + 1)
				profile_at_s = [seg.profile_start.interpolate(seg.profile_end, s)
				                for s in s_waypoints]
				uv_wp = [pr.sample_perimeter(n_peri) for pr in profile_at_s]
				for p_idx in range(n_sub):
					s1 = s_waypoints[p_idx]
					s2 = s_waypoints[p_idx + 1]
					u1, v1 = uv_wp[p_idx]
					u2, v2 = uv_wp[p_idx + 1]
					pr_mid = seg.profile_start.interpolate(
						seg.profile_end, 0.5 * (s1 + s2))
					dw_mid, dh_mid = _cell_side(pr_mid)
					A_mid = pr_mid.total_area() / n_peri
					for k in range(N):
						p1 = (seg.start_pos + s1 * seg.length * y_hat
						      + u1[k] * x_hat + v1[k] * z_hat)
						p2 = (seg.start_pos + s2 * seg.length * y_hat
						      + u2[k] * x_hat + v2[k] * z_hat)
						filament_paths[k].append((tuple(p1), tuple(p2)))
						cell_wh[k].append((dw_mid, dh_mid))
						dl = float(np.linalg.norm(p2 - p1))
						lengths[k] += dl
						R_k[k] += dl / (sigma * A_mid)

			elif isinstance(seg, StraightSegment):
				us_s, vs_s = seg.profile.sample_perimeter(n_peri)
				du_s, dv_s = _cell_side(seg.profile)
				A_cell_s = seg.profile.total_area() / n_peri
				for k in range(N):
					offset = us_s[k] * x_hat + vs_s[k] * z_hat
					p1 = seg.start_pos + offset
					p2 = seg.end_pos + offset
					filament_paths[k].append((tuple(p1), tuple(p2)))
					cell_wh[k].append((du_s, dv_s))
					lengths[k] += seg.length
					R_k[k] += seg.length / (sigma * A_cell_s)

			elif isinstance(seg, ArcSegment):
				theta_total = np.deg2rad(seg.arc_angle)
				R_arc = seg.radius
				us_s, vs_s = seg.profile.sample_perimeter(n_peri)
				du_s, dv_s = _cell_side(seg.profile)
				A_cell_s = seg.profile.total_area() / n_peri
				_check_arc_inner_radius(R_arc, us_s)
				for k in range(N):
					u, v = us_s[k], vs_s[k]
					r_k = R_arc + u
					z_off = v * z_hat
					for p_idx in range(n_arc):
						t1 = theta_total * p_idx / n_arc
						t2 = theta_total * (p_idx + 1) / n_arc
						p1 = (seg.arc_center + r_k * (np.cos(t1) * x_hat
						      + np.sin(t1) * y_hat) + z_off)
						p2 = (seg.arc_center + r_k * (np.cos(t2) * x_hat
						      + np.sin(t2) * y_hat) + z_off)
						filament_paths[k].append((tuple(p1), tuple(p2)))
						cell_wh[k].append((du_s, dv_s))
					arc_len_k = abs(r_k * theta_total)
					lengths[k] += arc_len_k
					R_k[k] += arc_len_k / (sigma * A_cell_s)

			elif isinstance(seg, LoftArcSegment):
				theta_total = np.deg2rad(seg.arc_angle)
				R_arc = seg.radius
				n_sub = seg.n_sub
				s_waypoints = np.linspace(0.0, 1.0, n_sub + 1)
				profile_at_s = [seg.profile_start.interpolate(seg.profile_end, s)
				                for s in s_waypoints]
				uv_wp = [pr.sample_perimeter(n_peri) for pr in profile_at_s]
				for pr_idx, (u_s, _) in enumerate(uv_wp):
					_check_arc_inner_radius(R_arc, u_s)
				for p_idx in range(n_sub):
					s1 = s_waypoints[p_idx]
					s2 = s_waypoints[p_idx + 1]
					theta1 = theta_total * s1
					theta2 = theta_total * s2
					u1, v1 = uv_wp[p_idx]
					u2, v2 = uv_wp[p_idx + 1]
					pr_mid = seg.profile_start.interpolate(
						seg.profile_end, 0.5 * (s1 + s2))
					dw_mid, dh_mid = _cell_side(pr_mid)
					A_mid = pr_mid.total_area() / n_peri
					for k in range(N):
						r_k_1 = R_arc + u1[k]
						r_k_2 = R_arc + u2[k]
						p1 = (seg.arc_center
						      + r_k_1 * (np.cos(theta1) * x_hat
						                 + np.sin(theta1) * y_hat)
						      + v1[k] * z_hat)
						p2 = (seg.arc_center
						      + r_k_2 * (np.cos(theta2) * x_hat
						                 + np.sin(theta2) * y_hat)
						      + v2[k] * z_hat)
						filament_paths[k].append((tuple(p1), tuple(p2)))
						cell_wh[k].append((dw_mid, dh_mid))
						dl = float(np.linalg.norm(p2 - p1))
						lengths[k] += dl
						R_k[k] += dl / (sigma * A_mid)
			else:
				raise NotImplementedError(
					f"Unsupported segment type {type(seg).__name__}")

		# Perimeter filaments all sit on the surface: d_k = 0, so the
		# Tier A skin factor = 1 for all.  Only the path-length bias
		# (1/R_k, arcs have different r_k) drives the weight.
		if frequency > 0:
			weights = (1.0 / R_k).astype(np.complex128)
		else:
			weights = 1.0 / R_k

		currents = self.current * weights / np.sum(weights)

		self._last_filament_info = {
			'us': us0, 'vs': vs0,
			'n_peri': n_peri,
			'w': w0, 'h': h0,
			'lengths': lengths,
			'R_k': R_k,
			'cell_wh': cell_wh,
			'currents': currents,
			'frequency': frequency,
			'sigma': sigma,
			'mu_r': mu_r,
			'has_loft': any(isinstance(s, (LoftStraightSegment,
			                                 LoftArcSegment))
			                for s in self.segments),
			'placement': 'perimeter',
		}
		return filament_paths, currents

	def plot_cross_section_current(self, nw, nh, frequency=0.0,
	                               sigma=5.8e7, mu_r=1.0, fig=None):
		"""Visualize the analytical filament current distribution.

		Shows |J| (magnitude) and phase(J) on a nw x nh grid matching the
		conductor cross-section. Requires matplotlib.

		Returns:
		    matplotlib Figure (newly created if fig is None).
		"""
		import matplotlib.pyplot as plt

		# Compute filament currents (also populates _last_filament_info)
		_, currents = self.to_filaments(nw, nh, frequency=frequency,
		                                sigma=sigma, mu_r=mu_r, n_arc=1)
		info = self._last_filament_info
		w, h = info['w'], info['h']

		# Current density per cell [A / cross-section-unit^2]
		A_cell = (w / nw) * (h / nh)
		J = currents.reshape(nw, nh) / A_cell

		if fig is None:
			if frequency > 0:
				fig, axes = plt.subplots(1, 2, figsize=(10, 4))
			else:
				fig, ax0 = plt.subplots(figsize=(5, 4))
				axes = [ax0]

		extent = [-w / 2, w / 2, -h / 2, h / 2]

		# |J|
		ax = axes[0]
		im = ax.imshow(np.abs(J).T, origin='lower', extent=extent,
		               aspect='equal', cmap='viridis')
		ax.set_xlabel('u (width)')
		ax.set_ylabel('v (height)')
		ax.set_title(f'|J|  f={frequency:g} Hz,  I_total={self.current:g} A')
		fig.colorbar(im, ax=ax)

		# Phase (AC only)
		if frequency > 0 and len(axes) > 1:
			ax = axes[1]
			im = ax.imshow(np.angle(J, deg=True).T, origin='lower',
			               extent=extent, aspect='equal',
			               cmap='twilight', vmin=-180, vmax=180)
			ax.set_xlabel('u (width)')
			ax.set_ylabel('v (height)')
			omega = 2.0 * np.pi * frequency
			delta = np.sqrt(2.0 / (omega * (4e-7 * np.pi) * mu_r * sigma))
			ax.set_title(f'phase(J) [deg],  skin depth = {delta*1e3:.3g} mm')
			fig.colorbar(im, ax=ax)

		fig.tight_layout()
		return fig

	# ============================================================
	# Loop closure and symmetry
	# ============================================================

	@property
	def gap(self):
		"""Distance between end position and start position."""
		if len(self.segments) == 0:
			return 0.0
		start = self.segments[0].start_pos
		end = self._position
		return np.linalg.norm(end - start)

	@property
	def is_closed(self):
		"""Check if the coil path forms a closed loop (within tolerance)."""
		return self.gap < 1e-10

	def close(self, tolerance=1e-6):
		"""Verify loop closure or optimize arc angles to close the gap.

		If the gap is within tolerance, does nothing.
		If the gap is larger, adjusts ALL arc segment angles simultaneously
		to close the loop. This handles 3D paths (tilt, helical) where
		a single arc adjustment is insufficient.

		Args:
			tolerance: Maximum allowed gap in meters

		Returns:
			self (for method chaining)

		Raises:
			ValueError: If no arc segments exist to adjust, or
			            if optimization fails to close the gap.
		"""
		if self.gap <= tolerance:
			return self

		# Find all arc segment indices
		arc_indices = [i for i, seg in enumerate(self.segments)
		               if isinstance(seg, ArcSegment)]

		if not arc_indices:
			raise ValueError(
				f"Gap = {self.gap:.6e} m but no arc segments to adjust. "
				f"Add arc segments or close manually.")

		from scipy.optimize import minimize

		start_pos_0 = self.segments[0].start_pos
		original_angles = [self.segments[i].arc_angle for i in arc_indices]
		original_segments = [self.segments[i] for i in arc_indices]

		def _rebuild_and_gap(angles):
			"""Rebuild path with modified arc angles, return gap."""
			# Start from segment 0
			pos = self.segments[0].start_pos.copy()
			orient = self.segments[0].orientation.copy()

			for i, seg in enumerate(self.segments):
				if i in arc_indices:
					angle = angles[arc_indices.index(i)]
					test = ArcSegment(
						seg.current, pos, orient,
						seg.width, seg.height, seg.radius, angle)
				elif isinstance(seg, StraightSegment):
					test = StraightSegment(
						seg.current, pos, orient,
						seg.width, seg.height, seg.length)
				else:
					test = ArcSegment(
						seg.current, pos, orient,
						seg.width, seg.height, seg.radius, seg.arc_angle)
				pos = test.end_pos
				orient = test.end_orientation

			return np.linalg.norm(pos - start_pos_0)

		# Optimize all arc angles simultaneously
		result = minimize(
			_rebuild_and_gap,
			x0=original_angles,
			method='Nelder-Mead',
			options={'xatol': 1e-8, 'fatol': tolerance * 0.01, 'maxiter': 10000})

		if result.fun > tolerance:
			raise ValueError(
				f"Could not close loop. Residual gap = {result.fun:.6e} m. "
				f"Try adding more arc segments for more degrees of freedom.")

		# Rebuild with optimized angles
		optimized_angles = result.x
		pos = self.segments[0].start_pos.copy()
		orient = self.segments[0].orientation.copy()

		for i in range(len(self.segments)):
			seg = self.segments[i]
			if i in arc_indices:
				angle = optimized_angles[arc_indices.index(i)]
				new_seg = ArcSegment(
					seg.current, pos, orient,
					seg.width, seg.height, seg.radius, angle)
			elif isinstance(seg, StraightSegment):
				new_seg = StraightSegment(
					seg.current, pos, orient,
					seg.width, seg.height, seg.length)
			else:
				new_seg = ArcSegment(
					seg.current, pos, orient,
					seg.width, seg.height, seg.radius, seg.arc_angle)
			self.segments[i] = new_seg
			pos = new_seg.end_pos
			orient = new_seg.end_orientation

		self._position = pos
		self._orientation = orient
		return self

	def mirror(self, plane='xz'):
		"""Create a mirrored copy of the coil.

		Returns a new CoilBuilder containing mirrored segments.
		The mirror operation reverses current direction.

		Args:
			plane: Mirror plane ('xz', 'yz', or 'xy')

		Returns:
			New CoilBuilder with mirrored coil (current reversed)
		"""
		mirror_matrix = {
			'xz': np.diag([1, -1, 1]),   # mirror across XZ (flip Y)
			'yz': np.diag([-1, 1, 1]),    # mirror across YZ (flip X)
			'xy': np.diag([1, 1, -1]),    # mirror across XY (flip Z)
		}
		if plane not in mirror_matrix:
			raise ValueError(f"Unknown plane '{plane}'. Use 'xz', 'yz', or 'xy'.")

		M = mirror_matrix[plane]
		mirrored = CoilBuilder(current=-self.current)
		mirrored._width = self._width
		mirrored._height = self._height

		for seg in self.segments:
			new_start = M @ seg.start_pos
			# Mirror flips handedness. Fix by negating one row to restore
			# right-handed orientation, then negate arc_angle to compensate.
			new_orient = M @ seg.orientation
			# Restore right-handedness: ensure det > 0
			if np.linalg.det(new_orient) < 0:
				new_orient = -new_orient  # flip all axes = equivalent rotation

			if isinstance(seg, StraightSegment):
				new_seg = StraightSegment.__new__(StraightSegment)
				CoilSegment.__init__(new_seg, -seg.current, new_start,
				                     new_orient, seg.width, seg.height)
				new_seg.length = seg.length
			elif isinstance(seg, ArcSegment):
				new_seg = ArcSegment.__new__(ArcSegment)
				CoilSegment.__init__(new_seg, -seg.current, new_start,
				                     new_orient, seg.width, seg.height)
				new_seg.radius = seg.radius
				new_seg.arc_angle = -seg.arc_angle
				new_seg.arc_center = M @ seg.arc_center

			mirrored.segments.append(new_seg)

		if len(mirrored.segments) > 0:
			last = mirrored.segments[-1]
			mirrored._position = last.end_pos
			mirrored._orientation = last.end_orientation

		return mirrored

	def rotate_copies(self, axis='z', n_copies=4):
		"""Create rotational copies of the coil.

		Returns a list of CoilBuilders, each rotated by 360/n_copies degrees.
		The first element is the original (unrotated).

		Args:
			axis: Rotation axis ('x', 'y', or 'z')
			n_copies: Total number of copies (including original)

		Returns:
			List of CoilBuilder objects
		"""
		angle_step = 360.0 / n_copies
		axis_vec = {'x': np.array([1, 0, 0]),
		            'y': np.array([0, 1, 0]),
		            'z': np.array([0, 0, 1])}[axis]

		copies = [self]
		for i in range(1, n_copies):
			angle_rad = np.deg2rad(i * angle_step)
			R = Rotation.from_rotvec(angle_rad * axis_vec).as_matrix()

			rotated = CoilBuilder(current=self.current)
			rotated._width = self._width
			rotated._height = self._height

			for seg in self.segments:
				new_start = R @ seg.start_pos
				new_orient = R @ seg.orientation
				if isinstance(seg, StraightSegment):
					new_seg = StraightSegment(
						seg.current, new_start, new_orient,
						seg.width, seg.height, seg.length)
				elif isinstance(seg, ArcSegment):
					new_seg = ArcSegment(
						seg.current, new_start, new_orient,
						seg.width, seg.height, seg.radius, seg.arc_angle)
				rotated.segments.append(new_seg)

			if len(rotated.segments) > 0:
				last = rotated.segments[-1]
				rotated._position = last.end_pos
				rotated._orientation = last.end_orientation
			copies.append(rotated)

		return copies

	def combined_occ(self, others=None):
		"""Combine this coil with others into a single OCC shape.

		Args:
			others: List of CoilBuilder objects to combine with

		Returns:
			Combined OCC shape
		"""
		from netgen.occ import Glue
		shapes = [self.to_occ()]
		if others:
			for other in others:
				shapes.append(other.to_occ())
		return Glue(shapes)


# Export public API
__all__ = ['CoilBuilder', 'CoilSegment', 'StraightSegment', 'ArcSegment']
