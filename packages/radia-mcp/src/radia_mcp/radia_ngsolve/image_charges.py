r"""Method-of-images electrostatic forces -- pure analytic closed forms.

The method of images replaces a grounded (Dirichlet) conductor by a small set of
fictitious "image" charges that enforce the zero-potential boundary, so the force on
the real charge is the bare Coulomb interaction with its images. Two classic geometries
are provided: a point charge outside a grounded conducting sphere, and an infinite line
charge parallel to a grounded conducting plane.

References:
  - J. D. Jackson, *Classical Electrodynamics*, 3rd ed., Secs. 2.1-2.3.
  - D. J. Griffiths, *Introduction to Electrodynamics*, 4th ed., Sec. 3.2.

These are tool-independent analytic helpers (no field solver involved).
"""
import math

EPS0 = 8.8541878128e-12        # vacuum permittivity [F/m]
MU0 = 4.0e-7 * math.pi          # vacuum permeability [H/m]
C0 = 299792458.0                # speed of light in vacuum [m/s]
ETA0 = 376.730313668            # vacuum wave impedance [Ohm] (= MU0*C0)


def grounded_sphere_image_force(charge_q, sphere_radius_a, distance_d, eps_r=1.0):
    r"""Force on a point charge q outside a grounded conducting sphere.

    A single image charge q' = -q*a/d sits at radius b = a^2/d from the sphere centre
    (inside the sphere), where a is the sphere radius and d the charge distance from the
    centre. The image is closer to the real charge than the centre, so the interaction is
    purely attractive. Two equivalent closed forms:

        F = (1/(4*pi*eps)) * q*q' / (d - b)^2
          = -(1/(4*pi*eps)) * q^2 * a * d / (d^2 - a^2)^2          (attractive, F < 0)

    with eps = eps_r * EPS0. Requires d > a (charge outside the sphere).
    (Jackson, *Classical Electrodynamics*, Sec. 2.3; Griffiths, Sec. 3.2.)

    Returns a dict with the signed force (attractive < 0), the image charge q' and its
    radial position b.
    """
    if sphere_radius_a <= 0.0:
        raise ValueError("sphere_radius_a must be > 0 (got %r)" % (sphere_radius_a,))
    if distance_d <= sphere_radius_a:
        raise ValueError(
            "distance_d must be > sphere_radius_a (charge must be outside the sphere); "
            "got d=%r, a=%r" % (distance_d, sphere_radius_a)
        )
    eps_r = float(eps_r)
    if eps_r <= 0.0:
        raise ValueError("eps_r must be > 0 (got %r)" % (eps_r,))

    eps = eps_r * EPS0
    a = float(sphere_radius_a)
    d = float(distance_d)
    q = float(charge_q)

    q_image = -q * a / d
    b = a * a / d
    # compact closed form (always attractive for a real charge)
    force = -(1.0 / (4.0 * math.pi * eps)) * q * q * a * d / (d * d - a * a) ** 2
    return {
        "force": force,
        "image_charge": q_image,
        "image_position": b,
    }


def line_charge_grounded_plane_image(line_charge_density, height_h, eps_r=1.0):
    r"""Force per unit length on a line charge above a grounded conducting plane.

    An infinite line charge of linear density lambda a height h above a grounded plane
    has an image line charge -lambda at depth h below the plane, a separation 2h. The
    field of a line charge falls as 1/r, giving a force per unit length

        F' = -lambda^2 / (4*pi*eps*h)                              (attractive, F' < 0)

    with eps = eps_r * EPS0. Requires h > 0.
    (Jackson, *Classical Electrodynamics*, Sec. 2.1.)

    Returns a dict with the force per unit length (attractive < 0) and the (signed)
    surface-normal field magnitude seen by the real line at its own location due to the
    image, E_image = lambda / (2*pi*eps*(2h)).
    """
    if height_h <= 0.0:
        raise ValueError("height_h must be > 0 (got %r)" % (height_h,))
    eps_r = float(eps_r)
    if eps_r <= 0.0:
        raise ValueError("eps_r must be > 0 (got %r)" % (eps_r,))

    eps = eps_r * EPS0
    lam = float(line_charge_density)
    h = float(height_h)

    force_per_length = -lam * lam / (4.0 * math.pi * eps * h)
    field_at_line = lam / (2.0 * math.pi * eps * (2.0 * h))
    return {
        "force_per_length": force_per_length,
        "field_at_line": field_at_line,
    }
