"""
ESIM VTK Export Utilities

Class-based VTK export following NGSolve's VTKOutput pattern.
Exports ESIM workpiece power distribution and field data to VTK format
for visualization in ParaView and other VTK-compatible tools.

Reference:
    NGSolve VTKOutput: https://docu.ngsolve.org/latest/i-tutorials/appendix-vtk/vtk.html

Author: Radia Development Team
Date: 2026-01-08
"""

import numpy as np


class ESIMVTKOutput:
    """
    VTK output class for ESIM analysis results.

    Follows NGSolve's VTKOutput pattern with:
    - Constructor sets up the output configuration
    - Do() method writes the actual file
    - Support for time series output

    Example:
        vtk = ESIMVTKOutput(
            workpiece=workpiece,
            coefs=['PowerDensity', 'H_tangential'],
            filename='esim_result',
            subdivision=0
        )
        vtk.Do()  # Writes esim_result.vtk

        # Time series
        for t in times:
            solver.solve()
            vtk.Do(time=t)  # Writes esim_result_0.vtk, esim_result_1.vtk, ...
    """

    # Available coefficient names for workpiece panels
    AVAILABLE_COEFS = [
        'PowerDensity',      # [W/m^2] Power loss per unit area
        'PowerLoss',         # [W] Total power loss per panel
        'H_tangential',      # [A/m] Tangential magnetic field magnitude
        'Z_magnitude',       # [Ohm] Surface impedance magnitude
        'Z_real',            # [Ohm] Surface impedance real part
        'Z_imag',            # [Ohm] Surface impedance imaginary part
        'ReactivePowerDensity',  # [var/m^2] Reactive power per unit area
        'ReactivePower',     # [var] Total reactive power per panel
        'SkinDepth',         # [m] Skin depth at panel
        'PanelNormal',       # [-] Panel normal vector (3 components)
    ]

    def __init__(self, workpiece=None, coil=None, coefs=None, names=None,
                 filename='esim_output', subdivision=0, floatsize='double',
                 legacy=False, grid_params=None):
        """
        Initialize ESIM VTK output.

        Parameters:
            workpiece: ESIMWorkpiece object (optional)
            coil: InductionHeatingCoil object (optional, for field export)
            coefs: List of coefficient names to export (default: all available)
                   For workpiece: 'PowerDensity', 'H_tangential', 'Z_magnitude', etc.
                   For coil field: 'B_field', 'B_magnitude'
            names: List of custom names for coefficients (default: same as coefs)
            filename: Base filename without extension
            subdivision: Subdivision level (0=no subdivision, n>0 subdivides panels)
            floatsize: 'single' or 'double' precision (default: 'double')
            legacy: Use legacy VTK format (default: False, uses XML VTU)
            grid_params: Dict for 3D field grid (optional):
                'x_range': [x_min, x_max]
                'y_range': [y_min, y_max]
                'z_range': [z_min, z_max]
                'nx', 'ny', 'nz': Grid points per direction
        """
        self.workpiece = workpiece
        self.coil = coil
        self.filename = filename
        self.subdivision = subdivision
        self.floatsize = floatsize
        self.legacy = legacy
        self.grid_params = grid_params

        # Default coefficients
        if coefs is None:
            if workpiece is not None:
                coefs = ['PowerDensity', 'H_tangential', 'Z_magnitude', 'PanelNormal']
            else:
                coefs = []
        self.coefs = coefs

        # Default names (same as coefs)
        if names is None:
            names = coefs
        self.names = names

        # Output counter for time series
        self._output_index = 0
        self._time_values = []

    def Do(self, time=None, vb='VOL'):
        """
        Write VTK output file.

        Parameters:
            time: Optional time stamp for time series
            vb: 'VOL' for volume, 'BND' for boundary (default: 'VOL')

        Returns:
            str: The filename that was written
        """
        # Determine filename with index for time series
        if time is not None:
            output_filename = f"{self.filename}_{self._output_index}"
            self._time_values.append(time)
            self._output_index += 1
        else:
            output_filename = self.filename

        files_written = []

        # Export workpiece if available
        if self.workpiece is not None:
            if self.legacy:
                fname = self._write_workpiece_legacy(output_filename)
            else:
                fname = self._write_workpiece_vtu(output_filename)
            files_written.append(fname)

        # Export coil field if grid_params provided
        if self.coil is not None and self.grid_params is not None:
            fname = self._write_field_grid(output_filename)
            files_written.append(fname)

        # Write PVD collection file for time series
        if len(self._time_values) > 1:
            self._write_pvd_collection()

        return files_written[0] if len(files_written) == 1 else files_written

    def _get_float_format(self):
        """Get format string for float output."""
        if self.floatsize == 'single':
            return '{:.6e}'
        else:
            return '{:.15e}'

    def _get_panel_data(self, panel, coef_name):
        """Extract coefficient value from panel."""
        if coef_name == 'PowerDensity':
            return panel.P_loss / panel.area if panel.area > 0 else 0.0
        elif coef_name == 'PowerLoss':
            return panel.P_loss
        elif coef_name == 'H_tangential':
            return abs(panel.H_tangential)
        elif coef_name == 'Z_magnitude':
            return abs(panel.Z_surface)
        elif coef_name == 'Z_real':
            return panel.Z_surface.real
        elif coef_name == 'Z_imag':
            return panel.Z_surface.imag
        elif coef_name == 'ReactivePowerDensity':
            return panel.Q_loss / panel.area if panel.area > 0 else 0.0
        elif coef_name == 'ReactivePower':
            return panel.Q_loss
        elif coef_name == 'SkinDepth':
            return panel.skin_depth if hasattr(panel, 'skin_depth') else 0.0
        elif coef_name == 'PanelNormal':
            return panel.normal  # Returns 3-component vector
        else:
            raise ValueError(f"Unknown coefficient: {coef_name}")

    def _write_workpiece_legacy(self, filename):
        """Write workpiece in legacy VTK format (.vtk)."""
        panels = self.workpiece.panels
        fmt = self._get_float_format()

        # Collect vertices and connectivity
        all_vertices = []
        cell_connectivity = []
        vertex_offset = 0

        for panel in panels:
            vertices = panel.vertices
            if self.subdivision > 0:
                vertices = self._subdivide_polygon(vertices, self.subdivision)
            n_verts = len(vertices)
            for v in vertices:
                all_vertices.append(v)
            cell_connectivity.append([n_verts] + list(range(vertex_offset, vertex_offset + n_verts)))
            vertex_offset += n_verts

        n_points = len(all_vertices)
        n_cells = len(panels)

        output_file = f"{filename}.vtk"
        with open(output_file, 'w') as f:
            # Header
            f.write("# vtk DataFile Version 3.0\n")
            f.write("ESIM Workpiece - Radia\n")
            f.write("ASCII\n")
            f.write("DATASET POLYDATA\n\n")

            # Points
            f.write(f"POINTS {n_points} double\n")
            for v in all_vertices:
                f.write(f"{fmt.format(v[0])} {fmt.format(v[1])} {fmt.format(v[2])}\n")

            # Polygons
            total_ints = sum(len(c) for c in cell_connectivity)
            f.write(f"\nPOLYGONS {n_cells} {total_ints}\n")
            for conn in cell_connectivity:
                f.write(" ".join(str(c) for c in conn) + "\n")

            # Cell data
            f.write(f"\nCELL_DATA {n_cells}\n")

            for coef, name in zip(self.coefs, self.names):
                if coef == 'PanelNormal':
                    # Vector data
                    f.write(f"\nNORMALS {name} double\n")
                    for panel in panels:
                        n = self._get_panel_data(panel, coef)
                        f.write(f"{fmt.format(n[0])} {fmt.format(n[1])} {fmt.format(n[2])}\n")
                else:
                    # Scalar data
                    f.write(f"\nSCALARS {name} double 1\n")
                    f.write("LOOKUP_TABLE default\n")
                    for panel in panels:
                        val = self._get_panel_data(panel, coef)
                        f.write(f"{fmt.format(val)}\n")

        print(f"Exported: {output_file} ({n_cells} panels, {n_points} vertices)")
        return output_file

    def _write_workpiece_vtu(self, filename):
        """Write workpiece in VTK XML format (.vtu)."""
        panels = self.workpiece.panels
        fmt = self._get_float_format()

        # Collect vertices and connectivity
        all_vertices = []
        cell_connectivity = []
        cell_offsets = []
        vertex_offset = 0
        current_offset = 0

        for panel in panels:
            vertices = panel.vertices
            if self.subdivision > 0:
                vertices = self._subdivide_polygon(vertices, self.subdivision)
            n_verts = len(vertices)
            for v in vertices:
                all_vertices.append(v)
            for i in range(n_verts):
                cell_connectivity.append(vertex_offset + i)
            current_offset += n_verts
            cell_offsets.append(current_offset)
            vertex_offset += n_verts

        n_points = len(all_vertices)
        n_cells = len(panels)

        output_file = f"{filename}.vtu"
        with open(output_file, 'w') as f:
            # XML header
            f.write('<?xml version="1.0"?>\n')
            f.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
            f.write('  <UnstructuredGrid>\n')
            f.write(f'    <Piece NumberOfPoints="{n_points}" NumberOfCells="{n_cells}">\n')

            # Points
            f.write('      <Points>\n')
            f.write('        <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
            for v in all_vertices:
                f.write(f"          {fmt.format(v[0])} {fmt.format(v[1])} {fmt.format(v[2])}\n")
            f.write('        </DataArray>\n')
            f.write('      </Points>\n')

            # Cells
            f.write('      <Cells>\n')
            # Connectivity
            f.write('        <DataArray type="Int32" Name="connectivity" format="ascii">\n')
            f.write('          ' + ' '.join(str(c) for c in cell_connectivity) + '\n')
            f.write('        </DataArray>\n')
            # Offsets
            f.write('        <DataArray type="Int32" Name="offsets" format="ascii">\n')
            f.write('          ' + ' '.join(str(o) for o in cell_offsets) + '\n')
            f.write('        </DataArray>\n')
            # Cell types (7 = VTK_POLYGON)
            f.write('        <DataArray type="UInt8" Name="types" format="ascii">\n')
            f.write('          ' + ' '.join(['7'] * n_cells) + '\n')
            f.write('        </DataArray>\n')
            f.write('      </Cells>\n')

            # Cell data
            f.write('      <CellData>\n')
            for coef, name in zip(self.coefs, self.names):
                if coef == 'PanelNormal':
                    f.write(f'        <DataArray type="Float64" Name="{name}" NumberOfComponents="3" format="ascii">\n')
                    for panel in panels:
                        n = self._get_panel_data(panel, coef)
                        f.write(f"          {fmt.format(n[0])} {fmt.format(n[1])} {fmt.format(n[2])}\n")
                    f.write('        </DataArray>\n')
                else:
                    f.write(f'        <DataArray type="Float64" Name="{name}" format="ascii">\n')
                    for panel in panels:
                        val = self._get_panel_data(panel, coef)
                        f.write(f"          {fmt.format(val)}\n")
                    f.write('        </DataArray>\n')
            f.write('      </CellData>\n')

            # Close tags
            f.write('    </Piece>\n')
            f.write('  </UnstructuredGrid>\n')
            f.write('</VTKFile>\n')

        print(f"Exported: {output_file} ({n_cells} panels, {n_points} vertices)")
        return output_file

    def _write_field_grid(self, filename):
        """Write 3D field grid in VTK format."""
        params = self.grid_params
        x_min, x_max = params['x_range']
        y_min, y_max = params['y_range']
        z_min, z_max = params['z_range']
        nx = params.get('nx', 21)
        ny = params.get('ny', 21)
        nz = params.get('nz', 21)

        x = np.linspace(x_min, x_max, nx)
        y = np.linspace(y_min, y_max, ny)
        z = np.linspace(z_min, z_max, nz)

        # Compute B field
        print(f"Computing B field on {nx}x{ny}x{nz} grid...")
        B_data = np.zeros((nx, ny, nz, 3))
        total = nx * ny * nz
        count = 0

        for i, xi in enumerate(x):
            for j, yj in enumerate(y):
                for k, zk in enumerate(z):
                    B = self.coil.compute_field_at_point([xi, yj, zk])
                    B_data[i, j, k, :] = np.real(B)
                    count += 1
                    if count % 1000 == 0:
                        print(f"  Progress: {count}/{total} ({100*count/total:.1f}%)")

        fmt = self._get_float_format()
        output_file = f"{filename}_field.vts"

        with open(output_file, 'w') as f:
            f.write('<?xml version="1.0"?>\n')
            f.write('<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">\n')
            f.write(f'  <StructuredGrid WholeExtent="0 {nx-1} 0 {ny-1} 0 {nz-1}">\n')
            f.write(f'    <Piece Extent="0 {nx-1} 0 {ny-1} 0 {nz-1}">\n')

            # Points
            f.write('      <Points>\n')
            f.write('        <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        f.write(f"          {fmt.format(x[i])} {fmt.format(y[j])} {fmt.format(z[k])}\n")
            f.write('        </DataArray>\n')
            f.write('      </Points>\n')

            # Point data
            f.write('      <PointData>\n')
            # B field vector
            f.write('        <DataArray type="Float64" Name="B_field" NumberOfComponents="3" format="ascii">\n')
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        Bx, By, Bz = B_data[i, j, k, :]
                        f.write(f"          {fmt.format(Bx)} {fmt.format(By)} {fmt.format(Bz)}\n")
            f.write('        </DataArray>\n')
            # B magnitude
            f.write('        <DataArray type="Float64" Name="B_magnitude" format="ascii">\n')
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        B_mag = np.linalg.norm(B_data[i, j, k, :])
                        f.write(f"          {fmt.format(B_mag)}\n")
            f.write('        </DataArray>\n')
            f.write('      </PointData>\n')

            f.write('    </Piece>\n')
            f.write('  </StructuredGrid>\n')
            f.write('</VTKFile>\n')

        print(f"Exported: {output_file} ({nx}x{ny}x{nz} = {total} points)")
        return output_file

    def _write_pvd_collection(self):
        """Write PVD collection file for time series."""
        pvd_file = f"{self.filename}.pvd"
        with open(pvd_file, 'w') as f:
            f.write('<?xml version="1.0"?>\n')
            f.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n')
            f.write('  <Collection>\n')
            for i, t in enumerate(self._time_values):
                ext = '.vtk' if self.legacy else '.vtu'
                f.write(f'    <DataSet timestep="{t}" file="{self.filename}_{i}{ext}"/>\n')
            f.write('  </Collection>\n')
            f.write('</VTKFile>\n')
        print(f"Exported time series collection: {pvd_file}")

    def _subdivide_polygon(self, vertices, level):
        """Subdivide polygon for higher resolution output."""
        if level <= 0:
            return vertices

        # Simple subdivision: add midpoints
        new_vertices = []
        n = len(vertices)
        for i in range(n):
            v1 = np.array(vertices[i])
            v2 = np.array(vertices[(i + 1) % n])
            new_vertices.append(vertices[i])
            new_vertices.append(((v1 + v2) / 2).tolist())

        if level > 1:
            return self._subdivide_polygon(new_vertices, level - 1)
        return new_vertices


# Convenience functions (backward compatible with old API)

def export_esim_workpiece_vtk(workpiece, filename='esim_workpiece', include_fields=True):
    """
    Export ESIM workpiece to VTK format (legacy function).

    Parameters:
        workpiece: ESIMWorkpiece object
        filename: Output filename without extension
        include_fields: Include field/power data (default: True)

    Returns:
        str: Filename written
    """
    if include_fields:
        coefs = ['PowerDensity', 'PowerLoss', 'H_tangential', 'Z_magnitude',
                 'ReactivePowerDensity', 'PanelNormal']
    else:
        coefs = []

    vtk = ESIMVTKOutput(
        workpiece=workpiece,
        coefs=coefs,
        filename=filename,
        legacy=True
    )
    return vtk.Do()


def export_esim_coil_field_vtk(coil, grid_params, filename='esim_coil_field'):
    """
    Export coil B field on 3D grid to VTK format (legacy function).

    Parameters:
        coil: InductionHeatingCoil object
        grid_params: Dict with x_range, y_range, z_range, nx, ny, nz
        filename: Output filename without extension

    Returns:
        str: Filename written
    """
    vtk = ESIMVTKOutput(
        coil=coil,
        grid_params=grid_params,
        filename=filename
    )
    return vtk.Do()


def export_esim_combined_vtk(solver, grid_params=None, base_filename='esim_analysis'):
    """
    Export complete ESIM analysis to VTK files (legacy function).

    Parameters:
        solver: ESIMCoupledSolver object
        grid_params: Optional dict for 3D field grid
        base_filename: Base filename for outputs

    Returns:
        list: Filenames written
    """
    vtk = ESIMVTKOutput(
        workpiece=solver.workpiece,
        coil=solver.coil if grid_params else None,
        grid_params=grid_params,
        filename=base_filename
    )
    return vtk.Do()


# Example usage
if __name__ == "__main__":
    import sys
    import os
    sys.path.insert(0, os.path.dirname(__file__))

    from esim_coupled_solver import InductionHeatingCoil, ESIMCoupledSolver
    from esim_workpiece import create_esim_block

    print("ESIM VTK Export Test")
    print("=" * 60)
    print()

    # Steel BH curve
    bh_curve = [
        [0, 0], [100, 0.2], [500, 0.9], [1000, 1.3],
        [5000, 1.8], [50000, 2.1],
    ]

    sigma = 2e6
    freq = 50000

    # Create coil
    coil = InductionHeatingCoil(
        coil_type='spiral',
        center=[0, 0, 0.02],
        inner_radius=0.03,
        outer_radius=0.05,
        pitch=0.005,
        num_turns=3,
        axis=[0, 0, 1],
    )
    coil.set_current(100)

    # Create workpiece
    workpiece = create_esim_block(
        center=[0, 0, -0.01],
        dimensions=[0.08, 0.08, 0.02],
        bh_curve=bh_curve,
        sigma=sigma,
        frequency=freq,
        panels_per_side=8
    )

    # Solve
    print("Solving ESIM coupled problem...")
    solver = ESIMCoupledSolver(coil, workpiece, freq)
    result = solver.solve(tol=1e-4, max_iter=10, verbose=False)
    print(f"  Power: P = {result['P_total']:.1f} W, Q = {result['Q_total']:.1f} var")
    print()

    # Test 1: Legacy VTK format
    print("Test 1: Legacy VTK format (.vtk)")
    vtk1 = ESIMVTKOutput(
        workpiece=workpiece,
        coefs=['PowerDensity', 'H_tangential', 'Z_magnitude'],
        filename='test_legacy',
        legacy=True
    )
    vtk1.Do()
    print()

    # Test 2: VTK XML format (.vtu)
    print("Test 2: VTK XML format (.vtu)")
    vtk2 = ESIMVTKOutput(
        workpiece=workpiece,
        coefs=['PowerDensity', 'H_tangential', 'Z_magnitude', 'PanelNormal'],
        names=['P_density [W/m2]', 'H_t [A/m]', '|Z| [Ohm]', 'Normal'],
        filename='test_xml',
        legacy=False
    )
    vtk2.Do()
    print()

    # Test 3: Time series
    print("Test 3: Time series output")
    vtk3 = ESIMVTKOutput(
        workpiece=workpiece,
        coefs=['PowerDensity'],
        filename='test_timeseries',
        legacy=False
    )
    for t in [0.0, 0.5, 1.0]:
        # In real use, you would update the solution here
        vtk3.Do(time=t)
    print()

    # Test 4: Backward-compatible function
    print("Test 4: Legacy function (backward compatible)")
    export_esim_workpiece_vtk(workpiece, filename='test_compat')
    print()

    print("=" * 60)
    print("All tests completed!")
    print()
    print("Open .vtk or .vtu files in ParaView for visualization.")
