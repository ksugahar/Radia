"""
ESIM Workpiece Module for Induction Heating Analysis

This module provides the ESIMWorkpiece class for modeling conductive ferromagnetic
workpieces in induction heating simulations using the Effective Surface Impedance
Method (ESIM).

The ESIM workpiece uses pre-computed ESI tables (from esim_cell_problem.py) to
efficiently handle nonlinear material behavior without full 3D eddy current
computation.

Reference:
    K. Hollaus, M. Kaltenbacher, J. Schoberl, "A Nonlinear Effective Surface
    Impedance in a Magnetic Scalar Potential Formulation," IEEE Trans. Magnetics,
    2025, DOI: 10.1109/TMAG.2025.3613932

Author: Radia Development Team
Date: 2026-01-08
"""

import numpy as np
from scipy.constants import mu_0

try:
    from .esim_cell_problem import ESITable, generate_esi_table_from_bh_curve
except ImportError:
    from esim_cell_problem import ESITable, generate_esi_table_from_bh_curve


class SurfacePanel:
    """
    Represents a surface panel on the ESIM workpiece.

    Each panel has a center, normal vector, area, and vertices.
    The tangential magnetic field H_t and surface impedance Z are
    computed during the coupled solve.
    """

    def __init__(self, panel_id, center, normal, area, vertices):
        """
        Initialize a surface panel.

        Parameters:
            panel_id: Unique identifier for the panel
            center: [x, y, z] center coordinates [m]
            normal: [nx, ny, nz] outward normal vector (unit)
            area: Panel area [m^2]
            vertices: List of vertex coordinates [[x1,y1,z1], ...]
        """
        self.panel_id = panel_id
        self.center = np.array(center)
        self.normal = np.array(normal)
        self.normal = self.normal / np.linalg.norm(self.normal)  # Normalize
        self.area = area
        self.vertices = np.array(vertices)

        # Solution quantities (set during solve)
        self.H_tangential = 0.0 + 0.0j  # Complex tangential field [A/m]
        self.Z_surface = 0.0 + 0.0j     # Complex surface impedance [Ohm]
        self.P_loss = 0.0               # Active power loss [W]
        self.Q_loss = 0.0               # Reactive power [var]


class ESIMWorkpiece:
    """
    ESIM-based workpiece for induction heating analysis.

    Uses pre-computed ESI table for nonlinear ferromagnetic material behavior.
    The workpiece is discretized into surface panels for field computation.
    """

    def __init__(self, esi_table, geometry='block', **kwargs):
        """
        Initialize an ESIM workpiece.

        Parameters:
            esi_table: ESITable object with pre-computed Z(H0) data
            geometry: Geometry type ('block', 'cylinder', 'custom')
            **kwargs: Geometry-specific parameters

        For 'block' geometry:
            center: [x, y, z] center coordinates [m]
            dimensions: [Lx, Ly, Lz] block dimensions [m]
            panels_per_side: Number of panels per side (default: 5)

        For 'cylinder' geometry:
            center: [x, y, z] center of top face [m]
            radius: Cylinder radius [m]
            height: Cylinder height [m]
            panels_radial: Panels in radial direction
            panels_axial: Panels in axial direction

        For 'custom' geometry:
            panels: List of SurfacePanel objects
        """
        self.esi_table = esi_table
        self.geometry_type = geometry
        self.panels = []

        if geometry == 'block':
            self._create_block_panels(**kwargs)
        elif geometry == 'cylinder':
            self._create_cylinder_panels(**kwargs)
        elif geometry == 'custom':
            self.panels = kwargs.get('panels', [])
        else:
            raise ValueError(f"Unknown geometry type: {geometry}")

    def _create_block_panels(self, center, dimensions, panels_per_side=5):
        """
        Create surface panels for a rectangular block.

        Parameters:
            center: [x, y, z] center coordinates
            dimensions: [Lx, Ly, Lz] block dimensions
            panels_per_side: Number of panels per side
        """
        cx, cy, cz = center
        Lx, Ly, Lz = dimensions
        n = panels_per_side

        panel_id = 0

        # Top face (z = cz + Lz/2, normal = +z)
        z_top = cz + Lz / 2
        dx = Lx / n
        dy = Ly / n
        for i in range(n):
            for j in range(n):
                x = cx - Lx/2 + (i + 0.5) * dx
                y = cy - Ly/2 + (j + 0.5) * dy
                vertices = [
                    [x - dx/2, y - dy/2, z_top],
                    [x + dx/2, y - dy/2, z_top],
                    [x + dx/2, y + dy/2, z_top],
                    [x - dx/2, y + dy/2, z_top],
                ]
                panel = SurfacePanel(
                    panel_id=panel_id,
                    center=[x, y, z_top],
                    normal=[0, 0, 1],
                    area=dx * dy,
                    vertices=vertices
                )
                self.panels.append(panel)
                panel_id += 1

        # Bottom face (z = cz - Lz/2, normal = -z)
        z_bot = cz - Lz / 2
        for i in range(n):
            for j in range(n):
                x = cx - Lx/2 + (i + 0.5) * dx
                y = cy - Ly/2 + (j + 0.5) * dy
                vertices = [
                    [x - dx/2, y - dy/2, z_bot],
                    [x + dx/2, y - dy/2, z_bot],
                    [x + dx/2, y + dy/2, z_bot],
                    [x - dx/2, y + dy/2, z_bot],
                ]
                panel = SurfacePanel(
                    panel_id=panel_id,
                    center=[x, y, z_bot],
                    normal=[0, 0, -1],
                    area=dx * dy,
                    vertices=vertices
                )
                self.panels.append(panel)
                panel_id += 1

        # Front face (y = cy + Ly/2, normal = +y)
        y_front = cy + Ly / 2
        dx = Lx / n
        dz = Lz / n
        for i in range(n):
            for k in range(n):
                x = cx - Lx/2 + (i + 0.5) * dx
                z = cz - Lz/2 + (k + 0.5) * dz
                vertices = [
                    [x - dx/2, y_front, z - dz/2],
                    [x + dx/2, y_front, z - dz/2],
                    [x + dx/2, y_front, z + dz/2],
                    [x - dx/2, y_front, z + dz/2],
                ]
                panel = SurfacePanel(
                    panel_id=panel_id,
                    center=[x, y_front, z],
                    normal=[0, 1, 0],
                    area=dx * dz,
                    vertices=vertices
                )
                self.panels.append(panel)
                panel_id += 1

        # Back face (y = cy - Ly/2, normal = -y)
        y_back = cy - Ly / 2
        for i in range(n):
            for k in range(n):
                x = cx - Lx/2 + (i + 0.5) * dx
                z = cz - Lz/2 + (k + 0.5) * dz
                vertices = [
                    [x - dx/2, y_back, z - dz/2],
                    [x + dx/2, y_back, z - dz/2],
                    [x + dx/2, y_back, z + dz/2],
                    [x - dx/2, y_back, z + dz/2],
                ]
                panel = SurfacePanel(
                    panel_id=panel_id,
                    center=[x, y_back, z],
                    normal=[0, -1, 0],
                    area=dx * dz,
                    vertices=vertices
                )
                self.panels.append(panel)
                panel_id += 1

        # Right face (x = cx + Lx/2, normal = +x)
        x_right = cx + Lx / 2
        dy = Ly / n
        dz = Lz / n
        for j in range(n):
            for k in range(n):
                y = cy - Ly/2 + (j + 0.5) * dy
                z = cz - Lz/2 + (k + 0.5) * dz
                vertices = [
                    [x_right, y - dy/2, z - dz/2],
                    [x_right, y + dy/2, z - dz/2],
                    [x_right, y + dy/2, z + dz/2],
                    [x_right, y - dy/2, z + dz/2],
                ]
                panel = SurfacePanel(
                    panel_id=panel_id,
                    center=[x_right, y, z],
                    normal=[1, 0, 0],
                    area=dy * dz,
                    vertices=vertices
                )
                self.panels.append(panel)
                panel_id += 1

        # Left face (x = cx - Lx/2, normal = -x)
        x_left = cx - Lx / 2
        for j in range(n):
            for k in range(n):
                y = cy - Ly/2 + (j + 0.5) * dy
                z = cz - Lz/2 + (k + 0.5) * dz
                vertices = [
                    [x_left, y - dy/2, z - dz/2],
                    [x_left, y + dy/2, z - dz/2],
                    [x_left, y + dy/2, z + dz/2],
                    [x_left, y - dy/2, z + dz/2],
                ]
                panel = SurfacePanel(
                    panel_id=panel_id,
                    center=[x_left, y, z],
                    normal=[-1, 0, 0],
                    area=dy * dz,
                    vertices=vertices
                )
                self.panels.append(panel)
                panel_id += 1

        # Store geometry parameters
        self.center = np.array(center)
        self.dimensions = np.array(dimensions)

    def _create_cylinder_panels(self, center, radius, height,
                                 panels_radial=8, panels_axial=5):
        """
        Create surface panels for a cylinder.

        Parameters:
            center: [x, y, z] center of top face
            radius: Cylinder radius
            height: Cylinder height (extends in -z direction)
            panels_radial: Number of panels around circumference
            panels_axial: Number of panels along axis
        """
        cx, cy, cz = center
        R = radius
        H = height
        n_r = panels_radial
        n_a = panels_axial

        panel_id = 0

        # Top face (z = cz, normal = +z)
        # Divide into radial segments (pizza slices)
        dtheta = 2 * np.pi / n_r
        for i in range(n_r):
            theta1 = i * dtheta
            theta2 = (i + 1) * dtheta
            theta_mid = (theta1 + theta2) / 2

            # Approximate as triangle from center to arc
            x_mid = cx + R * 0.67 * np.cos(theta_mid)  # Centroid of sector
            y_mid = cy + R * 0.67 * np.sin(theta_mid)

            vertices = [
                [cx, cy, cz],
                [cx + R * np.cos(theta1), cy + R * np.sin(theta1), cz],
                [cx + R * np.cos(theta2), cy + R * np.sin(theta2), cz],
            ]
            area = 0.5 * R * R * dtheta

            panel = SurfacePanel(
                panel_id=panel_id,
                center=[x_mid, y_mid, cz],
                normal=[0, 0, 1],
                area=area,
                vertices=vertices
            )
            self.panels.append(panel)
            panel_id += 1

        # Bottom face (z = cz - H, normal = -z)
        z_bot = cz - H
        for i in range(n_r):
            theta1 = i * dtheta
            theta2 = (i + 1) * dtheta
            theta_mid = (theta1 + theta2) / 2

            x_mid = cx + R * 0.67 * np.cos(theta_mid)
            y_mid = cy + R * 0.67 * np.sin(theta_mid)

            vertices = [
                [cx, cy, z_bot],
                [cx + R * np.cos(theta1), cy + R * np.sin(theta1), z_bot],
                [cx + R * np.cos(theta2), cy + R * np.sin(theta2), z_bot],
            ]
            area = 0.5 * R * R * dtheta

            panel = SurfacePanel(
                panel_id=panel_id,
                center=[x_mid, y_mid, z_bot],
                normal=[0, 0, -1],
                area=area,
                vertices=vertices
            )
            self.panels.append(panel)
            panel_id += 1

        # Side surface
        dz = H / n_a
        for i in range(n_r):
            theta1 = i * dtheta
            theta2 = (i + 1) * dtheta
            theta_mid = (theta1 + theta2) / 2

            for k in range(n_a):
                z1 = cz - k * dz
                z2 = cz - (k + 1) * dz
                z_mid = (z1 + z2) / 2

                x_mid = cx + R * np.cos(theta_mid)
                y_mid = cy + R * np.sin(theta_mid)

                # Normal points radially outward
                nx = np.cos(theta_mid)
                ny = np.sin(theta_mid)

                vertices = [
                    [cx + R * np.cos(theta1), cy + R * np.sin(theta1), z1],
                    [cx + R * np.cos(theta2), cy + R * np.sin(theta2), z1],
                    [cx + R * np.cos(theta2), cy + R * np.sin(theta2), z2],
                    [cx + R * np.cos(theta1), cy + R * np.sin(theta1), z2],
                ]
                area = R * dtheta * dz

                panel = SurfacePanel(
                    panel_id=panel_id,
                    center=[x_mid, y_mid, z_mid],
                    normal=[nx, ny, 0],
                    area=area,
                    vertices=vertices
                )
                self.panels.append(panel)
                panel_id += 1

        # Store geometry parameters
        self.center = np.array(center)
        self.radius = radius
        self.height = height

    @property
    def num_panels(self):
        """Return the number of surface panels."""
        return len(self.panels)

    @property
    def total_surface_area(self):
        """Return the total surface area [m^2]."""
        return sum(p.area for p in self.panels)

    def get_panel_centers(self):
        """Return array of panel center coordinates."""
        return np.array([p.center for p in self.panels])

    def get_panel_normals(self):
        """Return array of panel normal vectors."""
        return np.array([p.normal for p in self.panels])

    def set_tangential_field(self, panel_id, H_tangential):
        """
        Set the tangential magnetic field for a panel.

        Parameters:
            panel_id: Panel identifier
            H_tangential: Complex tangential field [A/m]
        """
        self.panels[panel_id].H_tangential = H_tangential

        # Update surface impedance from ESI table
        H_abs = abs(H_tangential)
        self.panels[panel_id].Z_surface = self.esi_table.get_impedance(H_abs)

    def update_all_impedances(self):
        """Update surface impedances for all panels based on current H field."""
        for panel in self.panels:
            H_abs = abs(panel.H_tangential)
            panel.Z_surface = self.esi_table.get_impedance(H_abs)

    def compute_power_losses(self):
        """
        Compute power losses for all panels.

        Returns:
            P_total: Total active power loss [W]
            Q_total: Total reactive power [var]
        """
        P_total = 0.0
        Q_total = 0.0

        for panel in self.panels:
            H_abs = abs(panel.H_tangential)
            P_prime, Q_prime = self.esi_table.get_power_loss(H_abs)

            panel.P_loss = P_prime * panel.area
            panel.Q_loss = Q_prime * panel.area

            P_total += panel.P_loss
            Q_total += panel.Q_loss

        return P_total, Q_total

    def get_power_distribution(self):
        """
        Get the power loss distribution over the workpiece surface.

        Returns:
            power_data: List of dicts with panel power information
        """
        power_data = []
        for panel in self.panels:
            power_data.append({
                'panel_id': panel.panel_id,
                'center': panel.center.tolist(),
                'area': panel.area,
                'P_loss': panel.P_loss,
                'Q_loss': panel.Q_loss,
                'P_density': panel.P_loss / panel.area if panel.area > 0 else 0,
            })
        return power_data

    def get_summary(self):
        """
        Get a summary of the workpiece state.

        Returns:
            summary: Dict with workpiece information
        """
        P_total, Q_total = self.compute_power_losses()

        # Find max power density
        max_P_density = 0.0
        for panel in self.panels:
            P_density = panel.P_loss / panel.area if panel.area > 0 else 0
            max_P_density = max(max_P_density, P_density)

        return {
            'geometry_type': self.geometry_type,
            'num_panels': self.num_panels,
            'total_surface_area': self.total_surface_area,
            'P_total': P_total,
            'Q_total': Q_total,
            'max_P_density': max_P_density,
        }


def create_esim_block(center, dimensions, bh_curve, sigma, frequency,
                      panels_per_side=5, esi_n_points=30):
    """
    Convenience function to create an ESIM block workpiece.

    Parameters:
        center: [x, y, z] center coordinates [m]
        dimensions: [Lx, Ly, Lz] block dimensions [m]
        bh_curve: BH curve data [[H1, B1], [H2, B2], ...]
        sigma: Conductivity [S/m]
        frequency: Operating frequency [Hz]
        panels_per_side: Number of panels per side
        esi_n_points: Number of points in ESI table

    Returns:
        workpiece: ESIMWorkpiece object
    """
    # Generate ESI table
    esi_table = generate_esi_table_from_bh_curve(
        bh_curve, sigma, frequency, n_points=esi_n_points
    )

    # Create workpiece
    workpiece = ESIMWorkpiece(
        esi_table=esi_table,
        geometry='block',
        center=center,
        dimensions=dimensions,
        panels_per_side=panels_per_side
    )

    return workpiece


def create_esim_cylinder(center, radius, height, bh_curve, sigma, frequency,
                         panels_radial=8, panels_axial=5, esi_n_points=30):
    """
    Convenience function to create an ESIM cylinder workpiece.

    Parameters:
        center: [x, y, z] center of top face [m]
        radius: Cylinder radius [m]
        height: Cylinder height [m]
        bh_curve: BH curve data [[H1, B1], [H2, B2], ...]
        sigma: Conductivity [S/m]
        frequency: Operating frequency [Hz]
        panels_radial: Panels around circumference
        panels_axial: Panels along axis
        esi_n_points: Number of points in ESI table

    Returns:
        workpiece: ESIMWorkpiece object
    """
    # Generate ESI table
    esi_table = generate_esi_table_from_bh_curve(
        bh_curve, sigma, frequency, n_points=esi_n_points
    )

    # Create workpiece
    workpiece = ESIMWorkpiece(
        esi_table=esi_table,
        geometry='cylinder',
        center=center,
        radius=radius,
        height=height,
        panels_radial=panels_radial,
        panels_axial=panels_axial
    )

    return workpiece


# Example usage and test
if __name__ == "__main__":
    print("ESIM Workpiece Test")
    print("=" * 60)

    # BH curve for steel
    bh_curve_steel = [
        [0, 0],
        [100, 0.2],
        [250, 0.5],
        [500, 0.9],
        [1000, 1.3],
        [2500, 1.6],
        [5000, 1.8],
        [10000, 1.95],
        [50000, 2.1],
    ]

    sigma_steel = 2e6  # S/m
    freq = 50000  # 50 kHz

    print("Creating ESIM block workpiece...")
    workpiece = create_esim_block(
        center=[0, 0, -0.01],
        dimensions=[0.1, 0.1, 0.02],  # 100mm x 100mm x 20mm
        bh_curve=bh_curve_steel,
        sigma=sigma_steel,
        frequency=freq,
        panels_per_side=5
    )

    print(f"Number of panels: {workpiece.num_panels}")
    print(f"Total surface area: {workpiece.total_surface_area*1e4:.2f} cm^2")
    print()

    # Simulate uniform tangential field
    print("Simulating uniform H_t = 1000 A/m on all panels...")
    H_t_uniform = 1000.0  # A/m
    for panel in workpiece.panels:
        workpiece.set_tangential_field(panel.panel_id, H_t_uniform)

    # Compute power losses
    P_total, Q_total = workpiece.compute_power_losses()
    print(f"Total active power: {P_total:.2f} W")
    print(f"Total reactive power: {Q_total:.2f} var")
    print()

    # Get summary
    summary = workpiece.get_summary()
    print("Summary:")
    for key, value in summary.items():
        if isinstance(value, float):
            print(f"  {key}: {value:.4g}")
        else:
            print(f"  {key}: {value}")

    print()
    print("Test completed!")
