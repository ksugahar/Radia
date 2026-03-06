"""
Radia-NGSolve Field Integration

Provides arbitrary background field from NGSolve CoefficientFunction or analytical functions.

Date: 2025-10-31
"""

import radia as rad
import numpy as np


def create_cf_field_source(cf, unit='m'):
	"""
	Create Radia background field source from NGSolve CoefficientFunction.

	Parameters
	----------
	cf : ngsolve.CoefficientFunction or callable
	    NGSolve CF defining the field [Bx, By, Bz]
	    CF should be evaluable at (x, y, z) and return 3-vector
	    Can also be a Python function (x, y, z) -> [Bx, By, Bz]
	unit : str, optional
	    Unit of CF coordinates: 'm' (meters) or 'mm' (millimeters)
	    Default: 'm' (NGSolve standard)

	Returns
	-------
	int
	    Radia object key for the field source

	Examples
	--------
	>>> from ngsolve import CF, x, y, z
	>>> # Create spatially varying field
	>>> Bx = CF(0)
	>>> By = CF(0)
	>>> Bz = 1.0 + 0.01 * x  # Linear gradient in x
	>>> B_cf = CF((Bx, By, Bz))
	>>>
	>>> # Register as Radia background field
	>>> field_obj = create_cf_field_source(B_cf, unit='m')
	>>> magnet = rad.ObjRecMag([0,0,0], [5,5,5], [0,0,1])
	>>> container = rad.ObjCnt([magnet, field_obj])
	>>> rad.Solve(container, 0.0001, 10000)

	>>> # Check total field
	>>> B_total = rad.Fld(container, 'b', [0, 0, 0])

	Notes
	-----
	- Radia uses millimeters internally
	- NGSolve typically uses meters
	- This function handles unit conversion automatically
	- For best performance with NGSolve CF, evaluate on a mesh with specialcf
	"""

	# Check if cf is NGSolve CF or Python callable
	has_ngsolve = False
	try:
	    import ngsolve
	    has_ngsolve = True
	    is_ngsolve_cf = isinstance(cf, ngsolve.CoefficientFunction)
	except ImportError:
	    is_ngsolve_cf = False

	if is_ngsolve_cf:
	    # NGSolve CoefficientFunction: need special handling
	    print("Warning: Direct NGSolve CF evaluation without mesh may be limited.")
	    print("For best results, use create_analytical_field() with a function that")
	    print("evaluates the CF on a GridFunction.")

	    # Create wrapper for NGSolve CF
	    def cf_wrapper(coords):
	        """
	        Wrapper to evaluate NGSolve CF at given point.

	        Args:
	            coords: [x, y, z] in mm (Radia units)

	        Returns:
	            [Bx, By, Bz] in Tesla
	        """
	        # Convert Radia mm to NGSolve units
	        if unit == 'm':
	            x_ngs = coords[0] / 1000.0  # mm → m
	            y_ngs = coords[1] / 1000.0
	            z_ngs = coords[2] / 1000.0
	        else:  # unit == 'mm'
	            x_ngs = coords[0]
	            y_ngs = coords[1]
	            z_ngs = coords[2]

	        try:
	            # Try to evaluate CF
	            # Note: This may not work for all CFs without a mesh context
	            B_value = cf(ngsolve.mesh.MakeIntegrationPoint(x_ngs, y_ngs, z_ngs))

	            if hasattr(B_value, '__len__') and len(B_value) >= 3:
	                return [float(B_value[0]), float(B_value[1]), float(B_value[2])]
	            elif hasattr(B_value, '__len__'):
	                # Scalar or less than 3 components
	                return [0.0, 0.0, float(B_value[0]) if len(B_value) > 0 else 0.0]
	            else:
	                # Scalar CF
	                return [0.0, 0.0, float(B_value)]

	        except Exception as e:
	            print(f"Warning: CF evaluation failed at ({x_ngs}, {y_ngs}, {z_ngs}): {e}")
	            print("Returning zero field. Consider using create_analytical_field() instead.")
	            return [0.0, 0.0, 0.0]

	    return rad.ObjBckgCF(cf_wrapper)

	elif callable(cf):
	    # Python callable
	    def wrapper(coords):
	        """
	        Wrapper for Python function with unit conversion.

	        Args:
	            coords: [x, y, z] in mm (Radia units)

	        Returns:
	            [Bx, By, Bz] in Tesla
	        """
	        if unit == 'm':
	            x_m = coords[0] / 1000.0
	            y_m = coords[1] / 1000.0
	            z_m = coords[2] / 1000.0
	            return cf(x_m, y_m, z_m)
	        else:  # unit == 'mm'
	            return cf(coords[0], coords[1], coords[2])

	    return rad.ObjBckgCF(wrapper)

	else:
	    raise TypeError("cf must be NGSolve CoefficientFunction or Python callable")


def create_analytical_field(field_func, unit='mm'):
	"""
	Create background field from analytical Python function.

	Parameters
	----------
	field_func : callable
	    Python function: [Bx, By, Bz] = field_func(x, y, z)
	    Coordinates in specified units
	unit : str, optional
	    Unit of coordinates: 'm' or 'mm'
	    Default: 'mm' (Radia internal units)

	Returns
	-------
	int
	    Radia object key for the field source

	Examples
	--------
	>>> # Quadrupole field
	>>> def quadrupole_field(x, y, z):
	...     gradient = 10.0  # T/m
	...     Bx = gradient * y / 1000  # Convert mm to m for gradient
	...     By = gradient * x / 1000
	...     Bz = 0.0
	...     return [Bx, By, Bz]
	>>>
	>>> field_obj = create_analytical_field(quadrupole_field, unit='mm')
	>>>
	>>> # Use in Radia
	>>> magnet = rad.ObjRecMag([0,0,0], [10,10,10], [0,0,1])
	>>> container = rad.ObjCnt([magnet, field_obj])
	>>> rad.Solve(container, 0.0001, 10000)

	>>> # Linear gradient in z
	>>> def gradient_field(x, y, z):
	...     B0 = 1.0  # T
	...     grad = 0.01  # T/mm
	...     return [0, 0, B0 + grad * z]
	>>>
	>>> field_obj = create_analytical_field(gradient_field, unit='mm')
	"""

	def wrapper(coords):
	    """
	    Wrapper with unit conversion.

	    Args:
	        coords: [x, y, z] in mm (Radia units)

	    Returns:
	        [Bx, By, Bz] in Tesla
	    """
	    if unit == 'mm':
	        return field_func(coords[0], coords[1], coords[2])
	    else:  # unit == 'm'
	        x_m = coords[0] / 1000.0
	        y_m = coords[1] / 1000.0
	        z_m = coords[2] / 1000.0
	        return field_func(x_m, y_m, z_m)

	return rad.ObjBckgCF(wrapper)


def create_uniform_field(Bx, By, Bz):
	"""
	Create uniform background field (alternative to rad.ObjBckg).

	Parameters
	----------
	Bx, By, Bz : float
	    Uniform field components in Tesla

	Returns
	-------
	int
	    Radia object key

	Examples
	--------
	>>> # Equivalent to rad.ObjBckg([0, 0, 1.0])
	>>> field_obj = create_uniform_field(0, 0, 1.0)

	Notes
	-----
	This is provided for completeness. For uniform fields,
	rad.ObjBckg() is more efficient.
	"""

	def uniform(x, y, z):
	    return [Bx, By, Bz]

	return create_analytical_field(uniform, unit='mm')


# Predefined field templates

def dipole_field(moment, position=(0, 0, 0), unit='m'):
	"""
	Create magnetic dipole field.

	Parameters
	----------
	moment : array_like
	    Magnetic dipole moment [mx, my, mz] in A·m²
	position : array_like, optional
	    Dipole position [x0, y0, z0] in specified units
	    Default: (0, 0, 0)
	unit : str, optional
	    Unit of coordinates: 'm' or 'mm'
	    Default: 'm'

	Returns
	-------
	int
	    Radia object key

	Examples
	--------
	>>> # Dipole at origin with moment in z
	>>> field_obj = dipole_field([0, 0, 1.0], unit='m')

	Notes
	-----
	B(r) = (μ₀/4π) * (3(m·r)r/r⁵ - m/r³)
	"""

	m = np.array(moment, dtype=float)
	r0 = np.array(position, dtype=float)
	mu0_4pi = 1e-7  # T·m/A

	def dipole(x, y, z):
	    # Position relative to dipole
	    if unit == 'm':
	        r = np.array([x - r0[0], y - r0[1], z - r0[2]])
	    else:  # mm
	        r = np.array([(x - r0[0])/1000, (y - r0[1])/1000, (z - r0[2])/1000])

	    r_mag = np.linalg.norm(r)

	    if r_mag < 1e-10:
	        # Avoid singularity at dipole location
	        return [0.0, 0.0, 0.0]

	    r3 = r_mag**3
	    r5 = r_mag**5

	    m_dot_r = np.dot(m, r)

	    B = mu0_4pi * (3 * m_dot_r * r / r5 - m / r3)

	    return [float(B[0]), float(B[1]), float(B[2])]

	return create_analytical_field(dipole, unit=unit)


def solenoid_field(I, turns_per_m, radius, length, center=(0, 0, 0), unit='m'):
	"""
	Create ideal solenoid field (uniform inside, zero outside).

	Parameters
	----------
	I : float
	    Current in Amperes
	turns_per_m : float
	    Number of turns per meter
	radius : float
	    Solenoid radius in specified units
	length : float
	    Solenoid length in specified units
	center : array_like, optional
	    Solenoid center [x0, y0, z0]
	unit : str, optional
	    Unit: 'm' or 'mm'

	Returns
	-------
	int
	    Radia object key

	Examples
	--------
	>>> # 1000 A·turns/m solenoid, 0.1m radius, 1m long
	>>> field_obj = solenoid_field(10.0, 100, 0.1, 1.0, unit='m')
	"""

	mu0 = 4*np.pi*1e-7  # H/m
	B_inside = mu0 * I * turns_per_m  # T

	c = np.array(center, dtype=float)

	def solenoid(x, y, z):
	    # Check if inside solenoid
	    if unit == 'm':
	        r_perp = np.sqrt((x - c[0])**2 + (y - c[1])**2)
	        z_rel = z - c[2]
	    else:  # mm
	        r_perp = np.sqrt((x - c[0])**2 + (y - c[1])**2) / 1000
	        z_rel = (z - c[2]) / 1000
	        radius_m = radius / 1000
	        length_m = length / 1000

	    if r_perp < radius and abs(z_rel) < length/2:
	        return [0.0, 0.0, B_inside]
	    else:
	        return [0.0, 0.0, 0.0]

	return create_analytical_field(solenoid, unit=unit)


if __name__ == "__main__":
	# Test examples
	print("Radia-NGSolve Field Integration Module")
	print("=======================================")
	print()
	print("Example usage:")
	print()
	print("1. Uniform field:")
	print("   field = create_uniform_field(0, 0, 1.0)")
	print()
	print("2. Quadrupole field:")
	print("   def quad(x, y, z):")
	print("       g = 10.0  # T/m")
	print("       return [g*y, g*x, 0]")
	print("   field = create_analytical_field(quad, unit='m')")
	print()
	print("3. Dipole field:")
	print("   field = dipole_field([0, 0, 1.0], unit='m')")
	print()
	print("4. With NGSolve CF:")
	print("   from ngsolve import CF, x, y, z")
	print("   B_cf = CF((0, 0, 1 + 0.01*x))")
	print("   field = create_cf_field_source(B_cf, unit='m')")
