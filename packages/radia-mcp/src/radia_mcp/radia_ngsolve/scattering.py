"""High-frequency physical-optics radar cross section (RCS) of canonical PEC shapes.

The radar cross section sigma is the equivalent isotropic area that a target presents
to an incident plane wave: P_scat = sigma * S_inc / (4 pi) into 4 pi steradians.  In the
HIGH-FREQUENCY (optical) regime ka >> 1 the scattered field of a smooth conducting body
is dominated by the geometric-optics specular flash, and the physical-optics (PO)
approximation -- integrating the induced surface current Js = 2 n x H_inc over the lit
region and evaluating the radiation integral in the far field -- gives compact closed
forms for the canonical "calibration" shapes (flat plate, dihedral and trihedral corner
reflectors).  These are the textbook RCS workhorses used to size radar returns and to
calibrate measurement ranges.

This module collects those pure closed forms as analytic helpers (no field solve); they
are the wave-optics, far-field counterpart of the static/eigenvalue FE solvers elsewhere
in the package.  All quantities are SI; angles are in degrees at the public boundary and
converted internally.

References
----------
* E. F. Knott, J. F. Shaeffer, M. T. Tuley, *Radar Cross Section*, 2nd ed.,
  Artech House, 1993 -- Ch. 5 (flat plates, physical optics) and Ch. 6 (corner reflectors).
* C. A. Balanis, *Advanced Engineering Electromagnetics*, Wiley -- physical optics and
  the 2 n x H_inc surface-current approximation.
* J. D. Jackson, *Classical Electrodynamics*, 3rd ed. -- Kirchhoff / physical-optics
  diffraction by an aperture and the sinc-pattern of a finite plate.
"""
import math

EPS0 = 8.8541878128e-12        # vacuum permittivity [F/m]
MU0 = 4.0e-7 * math.pi          # vacuum permeability [H/m]
C0 = 299792458.0                # speed of light in vacuum [m/s]
ETA0 = 376.730313668            # vacuum wave impedance sqrt(MU0/EPS0) [Ohm]


def _sinc(x):
    """Unnormalised sinc sin(x)/x with the removable limit sinc(0) = 1."""
    if x == 0.0:
        return 1.0
    return math.sin(x) / x


def physical_optics_flat_plate_rcs(side_a, side_b, wavelength, theta_deg=0.0):
    """Physical-optics RCS [m^2] of a flat rectangular PEC plate a x b.

    With k = 2 pi / lambda and lit area A = a*b, the broadside (specular) maximum is

        sigma_max = 4 pi A^2 / lambda^2

    and, scanning in the plane that contains the ``a`` edge, the PO pattern is the
    aperture sinc squared modulated by the projected-area cosine:

        sigma(theta) = sigma_max * cos(theta)^2 * [ sin(k a sin theta) / (k a sin theta) ]^2 ,

    with the removable limit sinc(0) -> 1 at broadside (theta = 0), where sigma = sigma_max.
    The flat plate is the canonical PO target and the brightest of the calibration shapes.
    (Knott / Shaeffer / Tuley, *Radar Cross Section*, 2nd ed., Ch. 5.)

    Parameters
    ----------
    side_a, side_b : float
        Plate edge lengths [m] (a is the scan-plane edge), both > 0.
    wavelength : float
        Free-space wavelength lambda [m], > 0.
    theta_deg : float, default 0.0
        Aspect angle off broadside [deg], measured from the plate normal.

    Returns
    -------
    dict
        ``rcs`` -- sigma(theta) [m^2];
        ``rcs_max`` -- broadside sigma_max [m^2];
        ``rcs_dbsm`` -- 10*log10(sigma(theta)) [dBsm] (``-inf`` at a pattern null).

    Raises
    ------
    ValueError
        If any of ``side_a``, ``side_b``, ``wavelength`` is <= 0.
    """
    if side_a <= 0.0 or side_b <= 0.0:
        raise ValueError("plate edges side_a, side_b must be positive")
    if wavelength <= 0.0:
        raise ValueError("wavelength must be positive")

    k = 2.0 * math.pi / wavelength
    area = side_a * side_b
    rcs_max = 4.0 * math.pi * area * area / (wavelength * wavelength)

    theta = math.radians(theta_deg)
    u = k * side_a * math.sin(theta)
    rcs = rcs_max * math.cos(theta) ** 2 * _sinc(u) ** 2

    rcs_dbsm = 10.0 * math.log10(rcs) if rcs > 0.0 else float("-inf")
    return {"rcs": rcs, "rcs_max": rcs_max, "rcs_dbsm": rcs_dbsm}


def dihedral_corner_reflector_peak_rcs(side_a, height_h, wavelength):
    """Peak physical-optics RCS [m^2] of a 90-degree dihedral corner reflector.

    A right-angle dihedral made of two flat plates a x h returns a double-bounce flash whose
    peak (boresight) RCS is

        sigma_dihedral = 8 pi a^2 h^2 / lambda^2 ,

    exactly TWICE the prefactor of a single matched flat plate of the same lit area
    (sigma_plate,max = 4 pi (a h)^2 / lambda^2), the signature of the two-bounce retro-return.
    For completeness the square trihedral (triangular three-plate corner of edge a) peak is

        sigma_trihedral = 12 pi a^4 / lambda^2 ,

    the standard radar-range calibration reflector.  (Knott / Shaeffer / Tuley,
    *Radar Cross Section*, 2nd ed., Ch. 6, corner reflectors.)

    Parameters
    ----------
    side_a : float
        Dihedral plate width / trihedral edge length a [m], > 0.
    height_h : float
        Dihedral plate height h [m], > 0 (sets the dihedral lit area a*h).
    wavelength : float
        Free-space wavelength lambda [m], > 0.

    Returns
    -------
    dict
        ``rcs_dihedral`` -- 8 pi a^2 h^2 / lambda^2 [m^2];
        ``rcs_trihedral`` -- 12 pi a^4 / lambda^2 [m^2];
        ``rcs_dihedral_dbsm`` -- 10*log10(rcs_dihedral) [dBsm].

    Raises
    ------
    ValueError
        If any of ``side_a``, ``height_h``, ``wavelength`` is <= 0.
    """
    if side_a <= 0.0 or height_h <= 0.0:
        raise ValueError("dihedral dimensions side_a, height_h must be positive")
    if wavelength <= 0.0:
        raise ValueError("wavelength must be positive")

    lam2 = wavelength * wavelength
    rcs_dihedral = 8.0 * math.pi * side_a ** 2 * height_h ** 2 / lam2
    rcs_trihedral = 12.0 * math.pi * side_a ** 4 / lam2
    rcs_dihedral_dbsm = 10.0 * math.log10(rcs_dihedral)
    return {
        "rcs_dihedral": rcs_dihedral,
        "rcs_trihedral": rcs_trihedral,
        "rcs_dihedral_dbsm": rcs_dihedral_dbsm,
    }
