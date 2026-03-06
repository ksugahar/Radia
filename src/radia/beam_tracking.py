"""
beam_tracking.py

Beam Tracking Module with NumPy Array Support

Provides numpy-based interface for charged particle trajectory computation
in magnetic fields. Wraps the core Radia functions with numpy array output.

Designed for accelerator physics applications:
- Undulators and wigglers (synchrotron light sources)
- Permanent magnet systems (insertion devices)
- Magnetic chicanes and dipole magnets

Scope: Magnetic field only (no electric field).
Energy is conserved during tracking.

Example:
    import radia as rad
    from radia.beam_tracking import particle_trajectory, focusing_potential

    # Create undulator magnet
    magnet = rad.ObjRecMag([0, 0, 0], [100, 10, 20], [0, 0, 1.0])

    # Track electron trajectory (returns numpy array)
    trajectory = particle_trajectory(
        magnet,
        energy_gev=6.0,
        initial_conditions=[0, 0, 0, 0],  # [x0, x'0, z0, z'0]
        s_range=[-1000, 1000],
        n_points=1001
    )

    # trajectory is (1001, 5) numpy array: [[s, x, x', z, z'], ...]
    print(f"Final position: x={trajectory[-1, 1]:.6f} mm")

Future extensions (2026 roadmap):
- Hamiltonian formulation (A-field based)
- Symplectic integrators (Forest-Ruth, Yoshida)
- FMM acceleration for large field models
- GPU parallel batch tracking
"""

import numpy as np
from typing import List, Tuple, Optional, Union

# Import Radia core module
try:
    from . import radia as rad
except ImportError:
    import radia as rad

# Physical constants
ELECTRON_MASS_GEV = 0.000510998928  # Electron rest mass [GeV/c^2]
C_LIGHT = 299792458.0  # Speed of light [m/s]


def particle_trajectory(
    obj: int,
    energy_gev: float,
    initial_conditions: List[float],
    s_range: List[float],
    n_points: int
) -> np.ndarray:
    """
    Compute particle trajectory in magnetic field.

    Tracks relativistic electron using Lorentz force equation:
        dp/dt = q(v × B)

    Integration uses Runge-Kutta method (RK4 or adaptive RK5).

    Args:
        obj: Radia object key for magnetic field source
        energy_gev: Particle energy [GeV]
        initial_conditions: [x0, x'0, z0, z'0] - initial transverse position and angles
            - x0: Initial x position [length unit set by FldUnits]
            - x'0: Initial x angle dx/ds [rad]
            - z0: Initial z position [length unit]
            - z'0: Initial z angle dz/ds [rad]
        s_range: [s0, s1] - longitudinal position range [length unit]
        n_points: Number of output points along trajectory

    Returns:
        numpy array of shape (n_points, 5):
            [[s, x, x', z, z'], ...]
            where:
            - s: Longitudinal position [length unit]
            - x, z: Transverse positions [length unit]
            - x', z': Transverse angles [rad]

    Example:
        >>> import radia as rad
        >>> from radia.beam_tracking import particle_trajectory
        >>> rad.FldUnits('mm')
        >>> magnet = rad.ObjRecMag([0, 0, 0], [100, 10, 20], [0, 0, 1.0])
        >>> traj = particle_trajectory(magnet, 6.0, [0, 0, 0, 0], [-1000, 1000], 101)
        >>> print(traj.shape)
        (101, 5)
    """
    if n_points < 2:
        raise ValueError("n_points must be >= 2")

    if len(initial_conditions) != 4:
        raise ValueError("initial_conditions must have 4 elements: [x0, x'0, z0, z'0]")

    if len(s_range) != 2:
        raise ValueError("s_range must have 2 elements: [s0, s1]")

    # Call Radia's FldPtcTrj function
    result = rad.FldPtcTrj(obj, energy_gev, initial_conditions, s_range, n_points)

    # Convert to numpy array
    if isinstance(result, (int, float)):
        # Single value returned (error or edge case)
        return np.array([[result, 0, 0, 0, 0]])

    # Result is list of lists: [[s, x, x', z, z'], ...]
    return np.array(result)


def focusing_potential(
    obj: int,
    p1: List[float],
    p2: List[float],
    n_points: int
) -> float:
    """
    Compute focusing potential along a line.

    Integrates transverse magnetic field components for beam optics analysis.
    Used for calculating focusing strength in periodic magnetic structures.

    Args:
        obj: Radia object key for magnetic field source
        p1: Start point [x1, y1, z1] [length unit]
        p2: End point [x2, y2, z2] [length unit]
        n_points: Number of integration points

    Returns:
        Focusing potential value

    Example:
        >>> import radia as rad
        >>> from radia.beam_tracking import focusing_potential
        >>> magnet = rad.ObjRecMag([0, 0, 0], [100, 10, 20], [0, 0, 1.0])
        >>> pot = focusing_potential(magnet, [-500, 0, 0], [500, 0, 0], 100)
    """
    if len(p1) != 3 or len(p2) != 3:
        raise ValueError("p1 and p2 must be 3D points [x, y, z]")

    return rad.FldFocPot(obj, p1, p2, n_points)


def shim_signature(
    obj: int,
    field_id: str,
    displacement: List[float],
    p1: List[float],
    p2: List[float],
    n_points: int,
    integration_dir: Optional[List[float]] = None
) -> np.ndarray:
    """
    Compute shim signature.

    Computes field variation from displacing a magnetic element.
    Useful for field correction analysis in undulators/wigglers.

    Args:
        obj: Radia object key
        field_id: Field component to compute:
            - 'bx', 'by', 'bz': Magnetic field components
            - 'ibx', 'iby', 'ibz': Field integrals
        displacement: Displacement vector [dx, dy, dz]
        p1: Start point [x, y, z]
        p2: End point [x, y, z]
        n_points: Number of evaluation points
        integration_dir: Integration direction [vx, vy, vz] for field integrals
            Default: [0, 0, 0] (no integration)

    Returns:
        numpy array of field variation values at each evaluation point

    Example:
        >>> import radia as rad
        >>> from radia.beam_tracking import shim_signature
        >>> magnet = rad.ObjRecMag([0, 0, 0], [100, 10, 20], [0, 0, 1.0])
        >>> sig = shim_signature(magnet, 'by', [0, 0.1, 0], [-500, 0, 0], [500, 0, 0], 101)
    """
    if len(displacement) != 3:
        raise ValueError("displacement must be 3D vector [dx, dy, dz]")

    if len(p1) != 3 or len(p2) != 3:
        raise ValueError("p1 and p2 must be 3D points [x, y, z]")

    if integration_dir is None:
        integration_dir = [0, 0, 0]
    elif len(integration_dir) != 3:
        raise ValueError("integration_dir must be 3D vector [vx, vy, vz]")

    result = rad.FldShimSig(obj, field_id, displacement, p1, p2, n_points, integration_dir)

    return np.array(result)


def focusing_kick_periodic(
    obj: int,
    start_point: List[float],
    longitudinal_dir: List[float],
    period: float,
    n_periods: int,
    transverse_dir: List[float],
    range1: float,
    n_points1: int,
    range2: float,
    n_points2: int,
    n_harmonics: int = 21,
    n_long_points: int = 100,
    d1: float = 0.0,
    d2: float = 0.0,
    energy_gev: float = 0.0,
    units: str = 'T2m2',
    output_format: str = 'fix'
) -> dict:
    """
    Compute second order kick matrices for periodic magnetic field.

    Computes 2nd order kick matrices for particle tracking codes.
    Output format compatible with BETA, TRACY, etc.

    Args:
        obj: Radia object key for magnetic field source
        start_point: Start point [x, y, z] [length unit]
        longitudinal_dir: Longitudinal direction vector [nx, ny, nz]
        period: Period length [length unit]
        n_periods: Number of periods
        transverse_dir: First transverse direction vector [nx, ny, nz]
        range1: Range in first transverse direction [length unit]
        n_points1: Number of points in first transverse direction
        range2: Range in second transverse direction [length unit]
        n_points2: Number of points in second transverse direction
        n_harmonics: Maximum number of harmonics to treat (default 21)
        n_long_points: Number of longitudinal points per period (default 100)
        d1: Derivative step in first transverse direction (0 = auto)
        d2: Derivative step in second transverse direction (0 = auto)
        energy_gev: Particle energy [GeV] (required if units='rad' or 'microrad')
        units: Output units - 'T2m2', 'rad', or 'microrad' (default 'T2m2')
        output_format: Output format - 'fix' (fixed-width) or 'tab' (tab-delimited)

    Returns:
        dict with keys:
            - 'kick1': 2D numpy array of kicks in first transverse direction
            - 'kick2': 2D numpy array of kicks in second transverse direction
            - 'int_bt_e2': 2D numpy array of integrated B^2
            - 'coord1': 1D numpy array of first transverse coordinates
            - 'coord2': 1D numpy array of second transverse coordinates
            - 'text': Formatted text output string

    Example:
        >>> import radia as rad
        >>> from radia.beam_tracking import focusing_kick_periodic
        >>> # Create simple undulator period
        >>> und = rad.ObjRecMag([0, 0, 5], [40, 20, 10], [0, 0, 1.2])
        >>> result = focusing_kick_periodic(
        ...     und,
        ...     start_point=[0, -100, 0],
        ...     longitudinal_dir=[0, 1, 0],
        ...     period=50.0,
        ...     n_periods=1,
        ...     transverse_dir=[1, 0, 0],
        ...     range1=10.0, n_points1=5,
        ...     range2=10.0, n_points2=5
        ... )
    """
    if len(start_point) != 3:
        raise ValueError("start_point must be 3D point [x, y, z]")

    if len(longitudinal_dir) != 3 or len(transverse_dir) != 3:
        raise ValueError("Direction vectors must be 3D [nx, ny, nz]")

    result = rad.FldFocKickPer(
        obj,
        start_point,
        longitudinal_dir,
        period,
        n_periods,
        transverse_dir,
        range1, n_points1,
        range2, n_points2,
        '',  # comment
        n_harmonics,
        n_long_points,
        d1, d2,
        units,
        energy_gev,
        output_format
    )

    # Parse result: [kick1_matrix, kick2_matrix, int_bt_e2_matrix, coord1_list, coord2_list, text_string]
    return {
        'kick1': np.array(result[0]),
        'kick2': np.array(result[1]),
        'int_bt_e2': np.array(result[2]),
        'coord1': np.array(result[3]),
        'coord2': np.array(result[4]),
        'text': result[5] if len(result) > 5 else ''
    }


def batch_trajectory(
    obj: int,
    energy_gev: float,
    initial_conditions_list: List[List[float]],
    s_range: List[float],
    n_points: int
) -> List[np.ndarray]:
    """
    Compute multiple particle trajectories (batch processing).

    Useful for phase space analysis or error studies.

    Args:
        obj: Radia object key for magnetic field source
        energy_gev: Particle energy [GeV]
        initial_conditions_list: List of initial conditions
            Each element is [x0, x'0, z0, z'0]
        s_range: [s0, s1] - longitudinal position range [length unit]
        n_points: Number of output points per trajectory

    Returns:
        List of numpy arrays, one per particle.
        Each array has shape (n_points, 5): [[s, x, x', z, z'], ...]

    Example:
        >>> import radia as rad
        >>> from radia.beam_tracking import batch_trajectory
        >>> import numpy as np
        >>> magnet = rad.ObjRecMag([0, 0, 0], [100, 10, 20], [0, 0, 1.0])
        >>> # Track 5 particles with different initial x positions
        >>> init_conds = [[x, 0, 0, 0] for x in np.linspace(-5, 5, 5)]
        >>> trajectories = batch_trajectory(magnet, 6.0, init_conds, [-1000, 1000], 101)
        >>> print(f"Tracked {len(trajectories)} particles")
    """
    trajectories = []
    for ic in initial_conditions_list:
        traj = particle_trajectory(obj, energy_gev, ic, s_range, n_points)
        trajectories.append(traj)
    return trajectories


# Export public API
__all__ = [
    'particle_trajectory',
    'focusing_potential',
    'shim_signature',
    'focusing_kick_periodic',
    'batch_trajectory',
    'ELECTRON_MASS_GEV',
    'C_LIGHT'
]
