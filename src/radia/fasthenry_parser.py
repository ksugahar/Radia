"""
fasthenry_parser.py

FastHenry .inp file parser for Radia PEEC

Parses FastHenry input files and converts to PEECBuilder topology.
Supports magnetic material blocks for coupled PEEC+MMM analysis.

Supported directives:
  .Units  - Length unit (m, mm, um, cm, etc.)
  N<name> - Node definition: N1 x=0 y=0 z=0
  E<name> - Segment definition: E1 N1 N2 w=1 h=1 [nwinc=3] [nhinc=3] [sigma=5.8e7]
  .external - Port definition: .external N1 N2
  .freq    - Frequency sweep: .freq fmin=1e3 fmax=1e9 ndec=10
  .default - Default segment parameters
  .equiv   - Node merge (equivalence)
  .magnetic / .endmagnetic - Magnetic material block (box/hexahedron/mesh)

Usage:
    from fasthenry_parser import FastHenryParser

    parser = FastHenryParser()
    parser.parse_file('inductor.inp')

    # Get Radia PEEC builder
    builder = parser.to_peec_builder()
    topo = builder.build_topology()

    # Get frequencies
    freqs = parser.get_frequencies()

    # Solve (auto-detects magnetic blocks and uses coupled solver)
    result = parser.solve(freqs)

Part of Radia project
"""

import re
import os
import sys
import numpy as np


class FastHenryParser:
    """
    Parser for FastHenry .inp input files.

    Converts FastHenry node-segment models to Radia PEECBuilder topology.
    """

    # Unit scale factors (to meters)
    UNIT_SCALES = {
        'm': 1.0,
        'meters': 1.0,
        'cm': 1e-2,
        'mm': 1e-3,
        'um': 1e-6,
        'in': 0.0254,
        'mils': 0.0254e-3,
    }

    def __init__(self):
        self.reset()

    def reset(self):
        """Reset parser state."""
        self.unit_scale = 1.0  # Default: meters
        self.unit_name = 'm'

        # Nodes: name -> {'x': float, 'y': float, 'z': float}
        self.nodes = {}

        # Segments: name -> {'n1': str, 'n2': str, 'w': float, 'h': float,
        #                     'sigma': float, 'nwinc': int, 'nhinc': int, ...}
        self.segments = {}

        # Ports: list of (node_pos_name, node_neg_name)
        self.ports = []

        # Frequency specification
        self.freq_spec = None  # {'fmin': float, 'fmax': float, 'ndec': int}

        # Default segment parameters
        self.defaults = {
            'w': 1e-3,
            'h': 1e-3,
            'sigma': 5.8e7,
            'nwinc': 1,
            'nhinc': 1,
            'rho': None,  # Resistivity (alternative to sigma)
            'rw': None,  # Wire radius (for circular cross-section)
        }

        # Node equivalences
        self.equiv_groups = []  # list of [node_name, ...]

        # Segment ordering (preserve insertion order)
        self.segment_order = []

        # Magnetic material blocks
        self.magnetic_blocks = []  # list of dicts with block specs

        # Comments / title
        self.title = ''

    def parse_file(self, filename):
        """
        Parse a FastHenry .inp file.

        Args:
            filename: Path to .inp file

        Raises:
            FileNotFoundError: If file does not exist
            ValueError: If file contains syntax errors
        """
        if not os.path.exists(filename):
            raise FileNotFoundError(f"FastHenry file not found: {filename}")

        with open(filename, 'r') as f:
            lines = f.readlines()

        self.reset()
        self._current_file_dir = os.path.dirname(os.path.abspath(filename))
        self._parse_lines(lines)

    def parse_string(self, text):
        """
        Parse FastHenry format from a string.

        Args:
            text: FastHenry format text
        """
        lines = text.splitlines(keepends=True)
        self.reset()
        self._current_file_dir = os.getcwd()
        self._parse_lines(lines)

    def _parse_lines(self, lines):
        """Parse list of lines."""
        # Join continuation lines (lines ending with '+')
        joined = []
        current = ''
        for line in lines:
            stripped = line.strip()
            if stripped.endswith('+'):
                current += stripped[:-1] + ' '
            else:
                current += stripped
                joined.append(current)
                current = ''
        if current:
            joined.append(current)

        # Parse with multi-line block support
        idx = 0
        while idx < len(joined):
            line = joined[idx]
            stripped = line.strip()

            # Check for .magnetic block start
            if stripped and stripped.lower().startswith('.magnetic'):
                # Collect lines until .endmagnetic
                block_lines = []
                idx += 1
                while idx < len(joined):
                    bline = joined[idx].strip()
                    if bline.lower().startswith('.endmagnetic'):
                        break
                    block_lines.append(bline)
                    idx += 1
                self._parse_magnetic_block(block_lines)
            else:
                self._parse_line(line)
            idx += 1

        # Apply node equivalences
        self._apply_equivalences()

    def _parse_line(self, line):
        """Parse a single line."""
        # Strip comments (everything after '*')
        if '*' in line:
            line = line[:line.index('*')]

        line = line.strip()
        if not line:
            return

        # Title line (first non-empty, non-directive line)
        lower = line.lower()

        # Directives start with '.'
        if lower.startswith('.units'):
            self._parse_units(line)
        elif lower.startswith('.default'):
            self._parse_default(line)
        elif lower.startswith('.external'):
            self._parse_external(line)
        elif lower.startswith('.freq'):
            self._parse_freq(line)
        elif lower.startswith('.equiv'):
            self._parse_equiv(line)
        elif lower.startswith('.end'):
            return  # End of file
        elif lower.startswith('.'):
            # Unknown directive - skip
            return
        elif line[0].upper() == 'N' and not line[0:2].upper().startswith('NW'):
            # Node definition: N<name> x=... y=... z=...
            self._parse_node(line)
        elif line[0].upper() == 'E':
            # Segment definition: E<name> N1 N2 w=... h=... ...
            self._parse_segment(line)
        # else: title or comment - skip

    def _parse_units(self, line):
        """Parse .Units directive."""
        # .Units mm
        parts = line.split()
        if len(parts) >= 2:
            unit = parts[1].lower()
            if unit in self.UNIT_SCALES:
                self.unit_scale = self.UNIT_SCALES[unit]
                self.unit_name = unit
            else:
                raise ValueError(f"Unknown unit: {parts[1]}")

    def _parse_default(self, line):
        """Parse .default directive for default segment parameters."""
        # .default w=1 h=1 sigma=5.8e7 nwinc=3 nhinc=3
        params = self._parse_key_value_pairs(line)
        for key, val in params.items():
            key_lower = key.lower()
            if key_lower == 'w':
                self.defaults['w'] = float(val) * self.unit_scale
            elif key_lower == 'h':
                self.defaults['h'] = float(val) * self.unit_scale
            elif key_lower == 'sigma':
                self.defaults['sigma'] = float(val)
            elif key_lower == 'rho':
                self.defaults['rho'] = float(val)
            elif key_lower == 'nwinc':
                self.defaults['nwinc'] = int(float(val))
            elif key_lower == 'nhinc':
                self.defaults['nhinc'] = int(float(val))
            elif key_lower == 'rw':
                self.defaults['rw'] = float(val) * self.unit_scale

    def _parse_node(self, line):
        """Parse node definition: N<name> x=... y=... z=..."""
        parts = line.split()
        name = parts[0]  # N<name>

        # Parse key=value pairs
        params = self._parse_key_value_pairs(line, start_idx=1)

        x = float(params.get('x', 0)) * self.unit_scale
        y = float(params.get('y', 0)) * self.unit_scale
        z = float(params.get('z', 0)) * self.unit_scale

        self.nodes[name] = {'x': x, 'y': y, 'z': z}

    def _parse_segment(self, line):
        """Parse segment definition: E<name> N1 N2 w=... h=... [nwinc=...] [nhinc=...] [sigma=...]"""
        parts = line.split()
        name = parts[0]  # E<name>

        # N1 and N2 are the next two tokens (not key=value)
        n1 = None
        n2 = None
        remaining_parts = []

        idx = 1
        while idx < len(parts):
            if '=' in parts[idx]:
                remaining_parts.append(parts[idx])
            else:
                if n1 is None:
                    n1 = parts[idx]
                elif n2 is None:
                    n2 = parts[idx]
                else:
                    remaining_parts.append(parts[idx])
            idx += 1

        if n1 is None or n2 is None:
            raise ValueError(f"Segment {name} missing node references: {line}")

        # Parse key=value pairs from remaining
        params = {}
        for part in remaining_parts:
            if '=' in part:
                key, val = part.split('=', 1)
                params[key.strip().lower()] = val.strip()

        # Also parse from the original line
        line_params = self._parse_key_value_pairs(line, start_idx=1)
        for k, v in line_params.items():
            if k.lower() not in ('x', 'y', 'z'):
                params[k.lower()] = v

        # Apply defaults
        w = float(params.get('w', self.defaults['w'] / self.unit_scale)) * self.unit_scale
        h = float(params.get('h', self.defaults['h'] / self.unit_scale)) * self.unit_scale

        sigma = self.defaults['sigma']
        if 'sigma' in params:
            sigma = float(params['sigma'])
        elif 'rho' in params:
            sigma = 1.0 / float(params['rho'])
        elif self.defaults['rho'] is not None:
            sigma = 1.0 / self.defaults['rho']

        nwinc = int(float(params.get('nwinc', self.defaults['nwinc'])))
        nhinc = int(float(params.get('nhinc', self.defaults['nhinc'])))

        # Check for circular wire (rw parameter)
        rw = None
        if 'rw' in params:
            rw = float(params['rw']) * self.unit_scale
        elif self.defaults.get('rw') is not None:
            rw = self.defaults['rw']

        seg = {
            'n1': n1,
            'n2': n2,
            'w': w,
            'h': h,
            'sigma': sigma,
            'nwinc': nwinc,
            'nhinc': nhinc,
        }
        if rw is not None:
            seg['rw'] = rw
            seg['cross_section'] = 'circular'
        else:
            seg['cross_section'] = 'rectangular'

        self.segments[name] = seg
        self.segment_order.append(name)

    def _parse_external(self, line):
        """Parse .external directive: .external N1 N2 [portname]"""
        parts = line.split()
        if len(parts) >= 3:
            n1 = parts[1]
            n2 = parts[2]
            self.ports.append((n1, n2))

    def _parse_freq(self, line):
        """Parse .freq directive: .freq fmin=1e3 fmax=1e9 ndec=10"""
        params = self._parse_key_value_pairs(line)

        fmin = float(params.get('fmin', 1e3))
        fmax = float(params.get('fmax', 1e9))
        ndec = int(float(params.get('ndec', 10)))

        self.freq_spec = {'fmin': fmin, 'fmax': fmax, 'ndec': ndec}

    def _parse_equiv(self, line):
        """Parse .equiv directive: .equiv N1 N2 N3 ..."""
        parts = line.split()
        if len(parts) >= 3:
            # All nodes after .equiv are equivalent
            equiv_nodes = parts[1:]
            self.equiv_groups.append(equiv_nodes)

    def _apply_equivalences(self):
        """Apply node equivalences: merge nodes in equiv groups."""
        for group in self.equiv_groups:
            if len(group) < 2:
                continue
            # Keep first node, replace all others
            master = group[0]
            for slave in group[1:]:
                # Replace slave with master in all segments
                for seg_name in self.segment_order:
                    seg = self.segments[seg_name]
                    if seg['n1'] == slave:
                        seg['n1'] = master
                    if seg['n2'] == slave:
                        seg['n2'] = master
                # Replace in ports
                new_ports = []
                for p1, p2 in self.ports:
                    if p1 == slave:
                        p1 = master
                    if p2 == slave:
                        p2 = master
                    new_ports.append((p1, p2))
                self.ports = new_ports

    def _parse_magnetic_block(self, lines):
        """
        Parse a .magnetic / .endmagnetic block.

        Supported formats:

        Box form (structured hex mesh):
            .magnetic
              type=box
              center=0,0,0
              size=0.03,0.03,0.03
              divisions=3,3,3
              mu_r=1000
              mu_r_imag=0
            .endmagnetic

        Hexahedron form (explicit 8 vertices):
            .magnetic
              type=hexahedron
              vertices=x1,y1,z1, x2,y2,z2, ..., x8,y8,z8
              mu_r=1000
              mu_r_imag=0
            .endmagnetic

        Mesh file form (GMSH .msh import):
            .magnetic
              type=mesh
              file=ferrite_core.msh
              mu_r=1000
              mu_r_imag=0
              physical_group=1    (optional: GMSH physical group filter)
            .endmagnetic

        Args:
            lines: List of lines between .magnetic and .endmagnetic
        """
        params = {}
        for line in lines:
            # Strip comments
            if '*' in line:
                line = line[:line.index('*')]
            line = line.strip()
            if not line:
                continue
            # Parse key=value
            if '=' in line:
                key, val = line.split('=', 1)
                params[key.strip().lower()] = val.strip()

        block_type = params.get('type', 'box').lower()
        mu_r = float(params.get('mu_r', 1000.0))
        mu_r_imag = float(params.get('mu_r_imag', 0.0))

        block = {
            'type': block_type,
            'mu_r': mu_r,
            'mu_r_imag': mu_r_imag,
        }

        if block_type == 'box':
            # Parse center, size, divisions
            center = [float(v) * self.unit_scale
                      for v in params.get('center', '0,0,0').split(',')]
            size = [float(v) * self.unit_scale
                    for v in params.get('size', '0.01,0.01,0.01').split(',')]
            divisions = [int(v)
                         for v in params.get('divisions', '1,1,1').split(',')]

            block['center'] = center
            block['size'] = size
            block['divisions'] = divisions

        elif block_type == 'hexahedron':
            # Parse 24 floats (8 vertices x 3 coords)
            vert_str = params.get('vertices', '')
            vals = [float(v) * self.unit_scale for v in vert_str.split(',')]
            if len(vals) != 24:
                raise ValueError(
                    f".magnetic hexahedron requires 24 vertex values "
                    f"(8 vertices x 3 coords), got {len(vals)}")
            vertices = [vals[i*3:(i+1)*3] for i in range(8)]
            block['vertices'] = vertices

        elif block_type == 'mesh':
            # GMSH .msh file import
            mesh_file = params.get('file', '')
            if not mesh_file:
                raise ValueError(
                    ".magnetic type=mesh requires 'file' parameter "
                    "(path to GMSH .msh file)")
            # Resolve relative path against input file directory
            if hasattr(self, '_current_file_dir') and self._current_file_dir:
                mesh_file = os.path.join(self._current_file_dir, mesh_file)
            block['file'] = mesh_file
            # Optional physical group filter
            pg = params.get('physical_group', None)
            block['physical_group'] = int(pg) if pg is not None else None

        else:
            raise ValueError(f"Unknown .magnetic type: {block_type}")

        self.magnetic_blocks.append(block)

    def _parse_key_value_pairs(self, line, start_idx=0):
        """
        Extract key=value pairs from a line.

        Returns dict of {key: value_string}.
        """
        params = {}
        # Find all key=value patterns
        pattern = r'(\w+)\s*=\s*([^\s,]+)'
        matches = re.findall(pattern, line)
        for key, val in matches:
            params[key.lower()] = val
        return params

    def get_frequencies(self):
        """
        Get frequency array from .freq directive.

        Returns:
            numpy array of frequencies [Hz], or None if no .freq directive
        """
        if self.freq_spec is None:
            return None

        fmin = self.freq_spec['fmin']
        fmax = self.freq_spec['fmax']
        ndec = self.freq_spec['ndec']

        if fmin <= 0 or fmax <= 0:
            return np.array([fmax])

        if fmin == fmax:
            return np.array([fmin])

        # Log-spaced: ndec points per decade
        n_decades = np.log10(fmax / fmin)
        n_points = int(n_decades * ndec) + 1
        return np.logspace(np.log10(fmin), np.log10(fmax), n_points)

    def get_summary(self):
        """
        Get a summary of the parsed model.

        Returns:
            dict with model statistics
        """
        return {
            'units': self.unit_name,
            'n_nodes': len(self.nodes),
            'n_segments': len(self.segments),
            'n_ports': len(self.ports),
            'n_equiv_groups': len(self.equiv_groups),
            'n_magnetic_blocks': len(self.magnetic_blocks),
            'freq_spec': self.freq_spec,
            'defaults': {k: v for k, v in self.defaults.items() if v is not None},
        }

    def to_peec_builder(self):
        """
        Convert parsed FastHenry model to a Radia PEECBuilder.

        Returns:
            PyPEECBuilder instance with topology set up

        Raises:
            ImportError: If peec_matrices module not available
            ValueError: If model has undefined node references
        """
        # Import PEECBuilder
        try:
            from peec_matrices import PEECBuilder as PyPEECBuilder
        except ImportError:
            # Try adding src/radia to path
            src_radia = os.path.join(os.path.dirname(__file__))
            if src_radia not in sys.path:
                sys.path.insert(0, src_radia)
            from peec_matrices import PEECBuilder as PyPEECBuilder

        builder = PyPEECBuilder()

        # Map node names to builder node IDs
        node_id_map = {}

        # Add all nodes
        for name, coords in self.nodes.items():
            nid = builder.add_node_at(coords['x'], coords['y'], coords['z'])
            node_id_map[name] = nid

        # Add all segments
        for seg_name in self.segment_order:
            seg = self.segments[seg_name]
            n1_name = seg['n1']
            n2_name = seg['n2']

            if n1_name not in node_id_map:
                raise ValueError(f"Segment {seg_name}: undefined node {n1_name}")
            if n2_name not in node_id_map:
                raise ValueError(f"Segment {seg_name}: undefined node {n2_name}")

            n1_id = node_id_map[n1_name]
            n2_id = node_id_map[n2_name]

            cross_type = 0  # rectangular
            if seg.get('cross_section') == 'circular':
                cross_type = 1

            builder.add_connected_segment(
                n1_id, n2_id,
                seg['w'], seg['h'],
                sigma=seg['sigma'],
                cross_section_type=cross_type,
                nwinc=seg['nwinc'],
                nhinc=seg['nhinc']
            )

        # Add ports
        for p_pos, p_neg in self.ports:
            if p_pos not in node_id_map:
                raise ValueError(f"Port: undefined node {p_pos}")
            if p_neg not in node_id_map:
                raise ValueError(f"Port: undefined node {p_neg}")
            builder.add_port(node_id_map[p_pos], node_id_map[p_neg])

        return builder

    def build_magnetic_objects(self):
        """
        Create Radia magnetic objects from parsed .magnetic blocks.

        Returns:
            list of Radia object handles, or empty list if no magnetic blocks

        Raises:
            ImportError: If radia module not available
        """
        if not self.magnetic_blocks:
            return []

        import sys
        import os
        try:
            sys.path.insert(0, os.path.dirname(__file__))
            from _radia_pybind import (ObjHexahedron, ObjCnt, MatLin, MatApl,
                                       FldUnits)
        except ImportError:
            import radia as rad
            ObjHexahedron = rad.ObjHexahedron
            ObjCnt = rad.ObjCnt
            MatLin = rad.MatLin
            MatApl = rad.MatApl

        objects = []

        for block in self.magnetic_blocks:
            mu_r = block['mu_r']

            if block['type'] == 'box':
                # Create structured hex mesh
                center = np.array(block['center'])
                size = np.array(block['size'])
                divs = block['divisions']

                # Element sizes
                dx = size[0] / divs[0]
                dy = size[1] / divs[1]
                dz = size[2] / divs[2]

                # Starting corner
                corner = center - size / 2.0

                sub_objs = []
                for iz in range(divs[2]):
                    for iy in range(divs[1]):
                        for ix in range(divs[0]):
                            x0 = corner[0] + ix * dx
                            y0 = corner[1] + iy * dy
                            z0 = corner[2] + iz * dz

                            verts = [
                                [x0, y0, z0],
                                [x0 + dx, y0, z0],
                                [x0 + dx, y0 + dy, z0],
                                [x0, y0 + dy, z0],
                                [x0, y0, z0 + dz],
                                [x0 + dx, y0, z0 + dz],
                                [x0 + dx, y0 + dy, z0 + dz],
                                [x0, y0 + dy, z0 + dz],
                            ]
                            obj = ObjHexahedron(verts, [0, 0, 0])
                            mat = MatLin(mu_r)
                            MatApl(obj, mat)
                            sub_objs.append(obj)

                if len(sub_objs) == 1:
                    objects.append(sub_objs[0])
                else:
                    objects.append(ObjCnt(sub_objs))

            elif block['type'] == 'hexahedron':
                obj = ObjHexahedron(block['vertices'], [0, 0, 0])
                mat = MatLin(mu_r)
                MatApl(obj, mat)
                objects.append(obj)

            elif block['type'] == 'mesh':
                # Import from GMSH .msh file
                from gmsh_mesh_import import gmsh_to_radia
                obj = gmsh_to_radia(
                    block['file'],
                    mu_r=mu_r,
                    physical_group=block.get('physical_group'),
                    unit_scale=1.0,  # Already in parser units
                )
                objects.append(obj)

        return objects

    def solve(self, freqs=None, Zs_func=None, solver_method=0,
              solver_prec=0.0001, solver_maxiter=1000):
        """
        Parse, build, and solve in one step.

        Automatically uses CoupledPEECSolver when .magnetic blocks are present.

        Args:
            freqs: Frequency array [Hz], or None to use .freq directive
            Zs_func: Optional surface impedance function(freq) -> Zs array
            solver_method: Radia solver method for coupling (0=LU, 1=BiCGSTAB)
            solver_prec: Radia solver precision
            solver_maxiter: Radia solver max iterations

        Returns:
            dict with:
                'freqs': frequency array [Hz]
                'Z_port': complex impedance array
                'R': real part (resistance)
                'L': imaginary part / omega (inductance)
                'topology': build_topology() result dict
                'Delta_L': coupling matrix (if magnetic blocks present)
        """
        builder = self.to_peec_builder()
        topo = builder.build_topology()

        if freqs is None:
            freqs = self.get_frequencies()
        if freqs is None:
            freqs = np.array([0.0])  # DC only
        freqs = np.asarray(freqs)

        Delta_L = None

        if self.magnetic_blocks:
            # Coupled PEEC + MMM solve
            from peec_coupled import CoupledPEECSolver

            # Extract mu_r_imag and mu_r_real from magnetic blocks
            # Use first block's values (all blocks typically share same material)
            mu_r_imag = self.magnetic_blocks[0].get('mu_r_imag', 0.0)
            mu_r_real = self.magnetic_blocks[0].get('mu_r', 1000.0)

            mag_objects = self.build_magnetic_objects()
            solver = CoupledPEECSolver(topo, mag_objects,
                                       mu_r_imag=mu_r_imag)
            solver.compute_coupling_matrix(
                solver_method=solver_method,
                solver_prec=solver_prec,
                solver_maxiter=solver_maxiter,
                mu_r_real=mu_r_real)
            Delta_L = solver.Delta_L
            Z_port = solver.frequency_sweep(freqs, Zs_func)
        else:
            # Standard PEEC solve (no magnetic coupling)
            from peec_topology import PEECCircuitSolver
            solver = PEECCircuitSolver(topo)
            Z_port = solver.frequency_sweep(freqs, Zs_func)

        # Extract R and L
        R = np.real(Z_port)
        omega = 2.0 * np.pi * freqs
        L = np.zeros_like(freqs)
        mask = omega > 0
        L[mask] = np.imag(Z_port[mask]) / omega[mask]

        result = {
            'freqs': freqs,
            'Z_port': Z_port,
            'R': R,
            'L': L,
            'topology': topo,
        }
        if Delta_L is not None:
            result['Delta_L'] = Delta_L

        return result

    def __repr__(self):
        summary = self.get_summary()
        parts = [
            f"nodes={summary['n_nodes']}",
            f"segments={summary['n_segments']}",
            f"ports={summary['n_ports']}",
            f"units={summary['units']}",
        ]
        if summary['n_magnetic_blocks'] > 0:
            parts.append(f"magnetic_blocks={summary['n_magnetic_blocks']}")
        return f"FastHenryParser({', '.join(parts)})"
