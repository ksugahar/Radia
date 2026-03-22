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
	>>> from radia_coil_builder import CoilBuilder
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

	def __init__(self, current, start_pos, orientation, width, height, tilt=0):
		"""
		Initialize coil segment.

		Args:
			current (float): Current in Amperes
			start_pos (array): Starting position [x, y, z] (in active Radia units)
			orientation (array): 3x3 orientation matrix (row vectors)
			width (float): Cross-section width (in active Radia units)
			height (float): Cross-section height (in active Radia units)
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

	def __init__(self, current, start_pos, orientation, width, height, length, tilt=0):
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

	def to_occ_shape(self):
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
		return shape


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
			start_pos (array): Starting position [x, y, z] (in active Radia units)
			orientation (array): 3x3 orientation matrix (row vectors)
			width (float): Cross-section width (in active Radia units)
			height (float): Cross-section height (in active Radia units)
			radius (float): Arc radius (in active Radia units)
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

	def to_occ_shape(self):
		"""Generate OCC revolved shape for this arc segment."""
		from netgen.occ import WorkPlane, Axes, Pnt, Axis, Vec, Z, X, Y
		# Cross-section rectangle at (radius, 0) in local frame
		r_inner = self.radius - self.width / 2
		wp = WorkPlane(Axes(Pnt(r_inner, 0, -self.height/2), n=-Y, h=X))
		face = wp.Rectangle(self.width, self.height).Face()
		# Revolve around Z axis
		shape = face.Revolve(Axis(Pnt(0, 0, 0), Z), self.arc_angle)
		# Apply ZX Euler rotation + translation to arc center
		ea = self.euler_angles
		shape = shape.Rotate(Axis(Pnt(0, 0, 0), Z), ea[2])
		shape = shape.Rotate(Axis(Pnt(0, 0, 0), X), ea[1])
		shape = shape.Move(Vec(*self.arc_center))
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
		Set cross-section dimensions for subsequent segments.

		Args:
			width (float): Width (in active Radia units)
			height (float): Height (in active Radia units)

		Returns:
			self (for method chaining)
		"""
		self._width = width
		self._height = height
		return self

	def _check_cross_section(self):
		"""Raise ValueError if cross-section has not been set."""
		if self._width is None or self._height is None:
			raise ValueError(
				"Cross-section not set. Call set_cross_section(width, height) "
				"before adding segments."
			)

	def add_straight(self, length, tilt=0):
		"""
		Add a straight segment.

		Args:
			length (float): Length (in active Radia units)
			tilt (float): Tilt angle in degrees (rotation around Y-axis)

		Returns:
			self (for method chaining)

		Raises:
			ValueError: If set_cross_section() has not been called.
		"""
		self._check_cross_section()
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
			radius (float): Arc radius (in active Radia units)
			arc_angle (float): Arc angle in degrees
			tilt (float): Tilt angle in degrees (rotation around Y-axis)

		Returns:
			self (for method chaining)

		Raises:
			ValueError: If set_cross_section() has not been called.
		"""
		self._check_cross_section()
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
		shapes = [seg.to_occ_shape() for seg in self.segments]
		if len(shapes) == 0:
			raise ValueError("No segments added")
		if len(shapes) == 1:
			return shapes[0]
		return Glue(shapes)

	def write_step(self, filename):
		"""Export coil geometry to STEP file.

		Args:
			filename: Output .step file path

		Returns:
			filename
		"""
		shape = self.to_occ()
		shape.WriteStep(filename)
		return filename

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
