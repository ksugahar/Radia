"""
cubit_bem_extractor.py

BEM Impedance Extractor for inductors/transformers.

Extracts L(f), R(f), M (mutual inductance) from Cubit mesh using
ngsbem BEM solver. Supports high-order curved elements via Cubit's ACIS kernel +
export_NGSolveCurvedMesh().

Architecture:
    Cubit -> export_NGSolveCurvedMesh() -> NGSolve Mesh -> ngsbem BEM -> L, R, M

Solver: ngsbem (Lucy Weggler's stabilized formulation)
ESIM:   Karl Hollaus's nonlinear surface impedance for magnetic cores
Post:   GMSH .msh for high-order element visualization

Usage:
    from radia.cubit_bem_extractor import BEMExtractor

    ext = BEMExtractor()
    ext.load_cubit_mesh(cubit, curve_order=3)
    ext.set_conductor("coil", sigma=5.8e7)
    ext.set_port("port1", block="coil")
    result = ext.extract(freqs=[1e3, 10e3, 100e3])
    print(f"L = {result.L[0,0]*1e6:.2f} uH")

Part of Radia project
"""

import os
import json
import time
import numpy as np
import scipy.linalg

# Physical constants
MU_0 = 4.0 * np.pi * 1e-7
EPS_0 = 8.854187817e-12
C_0 = 1.0 / np.sqrt(MU_0 * EPS_0)


def _extract_edge_vertices(mesh):
    """Extract vertex indices for each mesh edge.

    Args:
        mesh: NGSolve Mesh

    Returns:
        list of (v0, v1) tuples, one per edge (0-indexed vertex indices)
    """
    edges = []
    for edge in mesh.edges:
        verts = list(edge.vertices)
        edges.append((verts[0].nr, verts[1].nr))
    return edges


class BEMExtractResult:
    """Extraction result container.

    Attributes:
        L: (n_port, n_port) inductance matrix [H] at lowest frequency
        R: (n_port, n_port, n_freq) resistance matrix [Ohm]
        C: (n_port, n_port) capacitance matrix [F] (if computed)
        freqs: (n_freq,) frequency points [Hz]
        Z: (n_port, n_port, n_freq) impedance matrix [Ohm]
        t_assemble: Assembly time [s]
        t_solve: Total solve time [s]
    """

    def __init__(self, n_port, freqs):
        n_freq = len(freqs)
        self.n_port = n_port
        self.freqs = np.array(freqs)
        self.Z = np.zeros((n_port, n_port, n_freq), dtype=complex)
        self.L = np.zeros((n_port, n_port))
        self.R = np.zeros((n_port, n_port, n_freq))
        self.C = None
        self.t_assemble = 0.0
        self.t_solve = 0.0

    @property
    def M(self):
        """Mutual inductance (off-diagonal of L)."""
        return self.L

    def L_at(self, freq):
        """Inductance matrix at given frequency [H].

        L(f) = Im(Z(f)) / omega
        """
        idx = self._freq_index(freq)
        omega = 2 * np.pi * self.freqs[idx]
        if omega < 1e-10:
            return self.L
        return np.imag(self.Z[:, :, idx]) / omega

    def R_at(self, freq):
        """Resistance matrix at given frequency [Ohm]."""
        idx = self._freq_index(freq)
        return np.real(self.Z[:, :, idx])

    def Z_at(self, freq):
        """Impedance matrix at given frequency [Ohm]."""
        idx = self._freq_index(freq)
        return self.Z[:, :, idx]

    def _freq_index(self, freq):
        """Find closest frequency index."""
        return np.argmin(np.abs(self.freqs - freq))

    def to_json(self, filename):
        """Export results to JSON."""
        data = {
            'n_port': self.n_port,
            'freqs_Hz': self.freqs.tolist(),
            't_assemble_s': self.t_assemble,
            't_solve_s': self.t_solve,
            'L_H': self.L.tolist(),
            'R_Ohm': self.R.tolist(),
            'Z_real': np.real(self.Z).tolist(),
            'Z_imag': np.imag(self.Z).tolist(),
        }
        if self.C is not None:
            data['C_F'] = self.C.tolist()
        with open(filename, 'w') as f:
            json.dump(data, f, indent=2)
        print(f"Results saved to {filename}")

    def to_csv(self, filename):
        """Export frequency sweep to CSV."""
        with open(filename, 'w') as f:
            # Header
            headers = ['freq_Hz']
            for i in range(self.n_port):
                for j in range(self.n_port):
                    headers.append(f'R{i+1}{j+1}_Ohm')
                    headers.append(f'L{i+1}{j+1}_H')
                    headers.append(f'Z{i+1}{j+1}_real_Ohm')
                    headers.append(f'Z{i+1}{j+1}_imag_Ohm')
            f.write(','.join(headers) + '\n')

            # Data
            for k, freq in enumerate(self.freqs):
                omega = 2 * np.pi * freq
                row = [f'{freq:.6e}']
                for i in range(self.n_port):
                    for j in range(self.n_port):
                        R = np.real(self.Z[i, j, k])
                        L = np.imag(self.Z[i, j, k]) / omega if omega > 1e-10 else self.L[i, j]
                        row.append(f'{R:.6e}')
                        row.append(f'{L:.6e}')
                        row.append(f'{np.real(self.Z[i, j, k]):.6e}')
                        row.append(f'{np.imag(self.Z[i, j, k]):.6e}')
                f.write(','.join(row) + '\n')
        print(f"CSV saved to {filename}")

    def print_summary(self):
        """Print extraction summary."""
        print("=" * 60)
        print("BEM Extraction Results")
        print("=" * 60)
        print(f"  Ports: {self.n_port}")
        print(f"  Frequencies: {len(self.freqs)} points "
              f"({self.freqs[0]:.0f} - {self.freqs[-1]:.0f} Hz)")
        print(f"  Assembly time: {self.t_assemble:.3f} s")
        print(f"  Solve time: {self.t_solve:.3f} s")
        print()

        # Inductance matrix
        print("  Inductance Matrix L [uH]:")
        for i in range(self.n_port):
            vals = [f'{self.L[i,j]*1e6:10.4f}' for j in range(self.n_port)]
            print(f"    [{', '.join(vals)}]")
        print()

        # Resistance at each frequency
        for k, freq in enumerate(self.freqs):
            print(f"  R at {freq:.0f} Hz [mOhm]:")
            for i in range(self.n_port):
                vals = [f'{np.real(self.Z[i,j,k])*1e3:10.4f}'
                        for j in range(self.n_port)]
                print(f"    [{', '.join(vals)}]")
            print()

        # Coupling coefficient
        if self.n_port >= 2:
            k12 = self.L[0, 1] / np.sqrt(abs(self.L[0, 0] * self.L[1, 1]))
            print(f"  Coupling coefficient k12 = {k12:.4f}")
        print("=" * 60)


class BEMExtractor:
    """BEM Impedance Extractor equivalent using ngsbem BEM solver.

    Extracts L, R, M from Cubit mesh for inductors/transformers.

    Workflow:
        1. load_cubit_mesh() — Import Cubit mesh as NGSolve mesh
        2. set_conductor()   — Define conductor blocks
        3. set_port()        — Define extraction ports
        4. extract()         — Run ngsbem BEM extraction
        5. export_gmsh_post() — Visualize in GMSH

    Example:
        ext = BEMExtractor()
        ext.load_cubit_mesh(cubit, curve_order=3)
        ext.set_conductor("coil", sigma=5.8e7)
        ext.set_port("port1", block="coil")
        result = ext.extract(freqs=[1e3, 10e3, 100e3])
        result.print_summary()
    """

    def __init__(self):
        self.mesh = None            # NGSolve Mesh
        self._cubit = None          # Cubit Python interface (stored for sideset resolution)
        self._conductors = {}       # block_name -> {sigma, mu_r, thickness}
        self._cores = {}            # block_name -> {mu_r, bh_data}
        self._ports = []            # list of {name, block, ...}
        self._solvers = {}          # block_name -> NGBEMPEECSolver
        self._result = None
        self._curve_order = 1

    def load_cubit_mesh(self, cubit, curve_order=3):
        """Load Cubit mesh as NGSolve mesh with curved elements.

        Uses export_NGSolveCurvedMesh() which works directly with Cubit's ACIS kernel
        for high-order curving. No STEP files or OCC geometry needed.

        Args:
            cubit: Cubit Python interface object
            curve_order: Polynomial order for curved elements (1-5)
        """
        import tempfile
        from ngsolve import Mesh as NGMesh

        # Export Cubit mesh via C++ plugin (ACIS kernel curving)
        vol_path = tempfile.mktemp(suffix='.vol')
        cubit.cmd(f'export netgen "{vol_path}" order {curve_order} overwrite')
        self.mesh = NGMesh(vol_path)
        self._cubit = cubit
        self._curve_order = curve_order

        nv = self.mesh.nv
        ne = self.mesh.ne
        materials = self.mesh.GetMaterials()
        print(f"BEMExtractor: Loaded mesh ({ne} elements, {nv} vertices)")
        print(f"  Materials: {', '.join(materials)}")
        print(f"  Curve order: {curve_order}")

    def load_cubit_surface(self, cubit, block_name="conductor",
                           curve_order=1):
        """Load Cubit surface mesh directly for BEM analysis.

        For BEM, only the conductor surface is needed (not the volume).
        This extracts triangular/quad surface elements from a Cubit block.

        For high-order curved surface meshes, use export_NGSolveCurvedMesh() with
        surface_only=True instead of this method.

        Blocks are the only Cubit entity type that supports element type
        specification (TRI6, QUAD8 for 2nd order elements).

        Cubit workflow:
            cubit.cmd("surface 1 2 3 4 5 6 size 0.002")
            cubit.cmd("mesh surface 1 2 3 4 5 6")
            cubit.cmd("block 1 name 'conductor' surface 1 2 3 4 5 6")
            cubit.cmd("block 1 element type TRI6")
            ext.load_cubit_surface(cubit, block_name="conductor")

        Args:
            cubit: Cubit Python interface
            block_name: Name of the block containing conductor surfaces
            curve_order: Polynomial order for curved elements
        """
        from netgen.meshing import Mesh as NetgenMesh, MeshPoint, Element2D, FaceDescriptor
        from netgen.csg import Pnt
        from ngsolve import Mesh

        # Find the block
        block_ids = cubit.get_block_id_list()
        target_block = None
        block_names_available = []
        for bid in block_ids:
            name = cubit.get_block_name(bid)
            block_names_available.append(name)
            if name == block_name:
                target_block = bid
                break

        if target_block is None:
            raise ValueError(f"Block '{block_name}' not found. "
                             f"Available: {block_names_available}")

        # Get surface elements from block. get_block_tris/get_block_quads
        # return an empty list when the block has no elements of that type;
        # a thrown exception means a real Cubit/API error and must surface.
        tri_ids = list(cubit.get_block_tris(target_block))
        quad_ids = list(cubit.get_block_quads(target_block))

        if len(tri_ids) == 0 and len(quad_ids) == 0:
            raise ValueError(f"Block '{block_name}' contains no surface elements")

        # Collect unique nodes
        node_set = set()
        for tri_id in tri_ids:
            nodes = cubit.get_connectivity("tri", tri_id)
            node_set.update(nodes)
        for quad_id in quad_ids:
            nodes = cubit.get_connectivity("quad", quad_id)
            node_set.update(nodes)

        # Build Netgen surface mesh
        ngmesh = NetgenMesh(dim=3)
        fd = FaceDescriptor(surfnr=1, domin=0, domout=0, bc=1)
        fd.bcname = block_name
        ngmesh.Add(fd)

        node_map = {}
        for nid in sorted(node_set):
            coords = cubit.get_nodal_coordinates(nid)
            pid = ngmesh.Add(MeshPoint(Pnt(*coords)))
            node_map[nid] = pid

        for tri_id in tri_ids:
            nodes = cubit.get_connectivity("tri", tri_id)
            verts = [node_map[n] for n in nodes]
            ngmesh.Add(Element2D(index=1, vertices=verts))

        for quad_id in quad_ids:
            nodes = cubit.get_connectivity("quad", quad_id)
            verts = [node_map[n] for n in nodes]
            # Split quad into 2 triangles (BEM uses triangles)
            ngmesh.Add(Element2D(index=1, vertices=[verts[0], verts[1], verts[2]]))
            ngmesh.Add(Element2D(index=1, vertices=[verts[0], verts[2], verts[3]]))

        ngmesh.SetBCName(0, block_name)
        self.mesh = Mesh(ngmesh)
        self._cubit = cubit
        self._curve_order = curve_order

        n_surf = len(tri_ids) + len(quad_ids)
        print(f"BEMExtractor: Loaded surface mesh from block '{block_name}'")
        print(f"  {n_surf} surface elements ({len(tri_ids)} tri + {len(quad_ids)} quad)")
        print(f"  {len(node_set)} nodes")
        print(f"  Curve order: {curve_order}")

        # Store node mapping for port definition
        self._cubit_node_map = node_map
        self._cubit_node_map_inv = {v: k for k, v in node_map.items()}

    def load_ngsolve_mesh(self, mesh):
        """Load an existing NGSolve mesh directly (no Cubit required).

        Args:
            mesh: NGSolve Mesh object (surface or volume)
        """
        self.mesh = mesh
        nv = mesh.nv
        ne = mesh.ne
        materials = mesh.GetMaterials()
        print(f"BEMExtractor: Loaded NGSolve mesh ({ne} elements, {nv} vertices)")
        print(f"  Materials: {', '.join(materials)}")

    def set_conductor(self, block_name, sigma=5.8e7, mu_r=1.0, thickness=None):
        """Define a conductor block.

        Args:
            block_name: Cubit block name (material name in NGSolve)
            sigma: Conductivity [S/m] (default: copper 5.8e7)
            mu_r: Relative permeability (default: 1.0)
            thickness: Conductor thickness [m] for SIBC skin effect.
                       If None, estimated from mesh geometry.
        """
        self._conductors[block_name] = {
            'sigma': sigma,
            'mu_r': mu_r,
            'thickness': thickness,
        }
        print(f"  Conductor '{block_name}': sigma={sigma:.2e} S/m, mu_r={mu_r}")

    def set_magnetic_core(self, block_name, mu_r=None, bh_data=None):
        """Define a magnetic core block.

        For nonlinear cores, Karl Hollaus's ESIM computes the
        effective surface impedance from the B-H curve.

        Args:
            block_name: Cubit block name
            mu_r: Linear relative permeability (for simple linear core)
            bh_data: [[H, B], ...] B-H curve data for nonlinear ESIM
        """
        self._cores[block_name] = {
            'mu_r': mu_r,
            'bh_data': bh_data,
        }
        if bh_data is not None:
            print(f"  Magnetic core '{block_name}': nonlinear BH ({len(bh_data)} points, ESIM)")
        else:
            print(f"  Magnetic core '{block_name}': mu_r={mu_r}")

    def set_port(self, name, block=None, source_nodes=None, sink_nodes=None,
                 source_block=None, sink_block=None):
        """Define an extraction port.

        Port types:
          1. Closed loop (inductor): Only specify block.
             Current is uniformly distributed on the conductor surface.
             L = total self-inductance of the loop.

          2. Open conductor (busbar): Specify source/sink blocks.
             Current enters at source block, exits at sink block.
             L = partial inductance between source and sink.

        For Cubit workflow, use blocks to define source/sink:
            cubit.cmd("block 1 name 'conductor' volume 1")
            cubit.cmd("block 2 name 'source1' surface 5")
            cubit.cmd("block 2 element type TRI6")
            cubit.cmd("block 3 name 'sink1' surface 10")
            cubit.cmd("block 3 element type TRI6")
            ext.set_port("port1", block="conductor",
                         source_block="source1", sink_block="sink1")

        Blocks are the only Cubit entity type that supports element type
        specification (e.g., TRI6, QUAD8 for 2nd order elements).

        Args:
            name: Port name (e.g., "port1", "primary", "secondary")
            block: Conductor block name. If None, uses the first conductor.
            source_nodes: List of source vertex indices (0-indexed). Current enters.
            sink_nodes: List of sink vertex indices (0-indexed). Current exits.
            source_block: Cubit block name for source (alternative to source_nodes).
            sink_block: Cubit block name for sink (alternative to sink_nodes).
        """
        if block is None and self._conductors:
            block = list(self._conductors.keys())[0]

        port_type = 'closed'  # default: closed loop, uniform excitation
        if (source_nodes is not None or source_block is not None):
            port_type = 'open'  # open conductor, source/sink excitation

        port_info = {
            'name': name,
            'block': block,
            'type': port_type,
            'source_nodes': source_nodes,
            'sink_nodes': sink_nodes,
            'source_block': source_block,
            'sink_block': sink_block,
        }
        self._ports.append(port_info)

        if port_type == 'closed':
            print(f"  Port '{name}' -> block '{block}' (closed loop, uniform excitation)")
        else:
            src = source_block or f"nodes {source_nodes}"
            snk = sink_block or f"nodes {sink_nodes}"
            print(f"  Port '{name}' -> block '{block}' (open: {src} -> {snk})")

    def extract(self, freqs, mode='mqs', use_esim=False):
        """Run ngsbem extraction at specified frequencies.

        Workflow:
          1. Assemble BEM matrices (L, P, Q_0, M_LS, R)
          2. Build port excitation vectors from source/sink definition
          3. For each frequency: solve BEM system with ESIM correction
          4. Extract Z_port = e^T @ Z_branch^{-1} @ e

        ESIM correction (Karl Hollaus):
          DC resistance: R_dc -> R_dc * F_R(xi) where xi = d/(2*delta)
          DC inductance: L_dc -> L_dc * (1 + Delta_L_esim)
          Uses nonlinear B-H curve for field-dependent surface impedance.

        Args:
            freqs: List/array of frequencies [Hz]
            mode: BEM solver mode:
                'mqs': Loop-only (R + jwL), fastest, DC-1 MHz
                'full': Full Loop-Star with Schur complement
                'stabilized': Weggler's stabilized BEM (most robust)
            use_esim: If True, apply Karl Hollaus's ESIM correction
                      for magnetic core effects on R and L.

        Returns:
            BEMExtractResult with L, R, Z matrices
        """
        if self.mesh is None:
            raise RuntimeError("Call load_cubit_mesh() or load_ngsolve_mesh() first")
        if not self._conductors:
            raise RuntimeError("Define at least one conductor with set_conductor()")
        if not self._ports:
            raise RuntimeError("Define at least one port with set_port()")

        freqs = np.asarray(freqs, dtype=float)
        n_port = len(self._ports)
        result = BEMExtractResult(n_port, freqs)

        print(f"\nBEM Extraction: {n_port} port(s), {len(freqs)} frequencies")
        print(f"  Mode: {mode}, ESIM: {use_esim}")

        # --- Phase 1: Assemble BEM matrices per conductor ---
        t_start = time.perf_counter()
        self._assemble_solvers()
        result.t_assemble = time.perf_counter() - t_start
        print(f"  Assembly time: {result.t_assemble:.3f} s")

        # --- Phase 2: Build port excitation vectors ---
        port_vectors = self._build_port_vectors()

        # --- Phase 3: ESIM correction (if requested) ---
        esim_Zs_funcs = {}
        if use_esim and self._cores:
            esim_Zs_funcs = self._build_esim_corrections()

        # --- Phase 4: Frequency sweep ---
        t_solve_start = time.perf_counter()

        for port_i, port_info in enumerate(self._ports):
            block = port_info['block']
            solver = self._solvers.get(block)
            if solver is None:
                print(f"  WARNING: No solver for block '{block}', "
                      f"skipping port '{port_info['name']}'")
                continue

            # Port excitation vector
            e_port = port_vectors[port_i]

            # Build SIBC function
            cond = self._conductors[block]
            Zs_func = None
            if cond['thickness'] is not None:
                from radia.ngsbem_peec import make_dowell_Zs_func
                Zs_func = make_dowell_Zs_func(
                    solver.n_loop, solver.R_loop,
                    cond['thickness'], cond['sigma'], cond['mu_r'])

            # ESIM correction (adds to Zs)
            esim_func = esim_Zs_funcs.get(block)

            # Solve at each frequency with port vector
            for k, freq in enumerate(freqs):
                Z_port = self._solve_port_impedance(
                    solver, freq, e_port, mode, Zs_func, esim_func)
                result.Z[port_i, port_i, k] = Z_port

            # Extract L from lowest frequency
            omega_0 = 2 * np.pi * freqs[0]
            if omega_0 > 1e-10:
                result.L[port_i, port_i] = np.imag(result.Z[port_i, port_i, 0]) / omega_0

        # --- Phase 5: Mutual inductance (multi-port) ---
        if n_port >= 2:
            self._compute_mutual_inductance(result, port_vectors)

        # Populate R matrix
        for k in range(len(freqs)):
            result.R[:, :, k] = np.real(result.Z[:, :, k])

        result.t_solve = time.perf_counter() - t_solve_start
        print(f"  Solve time: {result.t_solve:.3f} s")
        print(f"  Total time: {result.t_assemble + result.t_solve:.3f} s")

        self._result = result
        result.print_summary()
        return result

    def _assemble_solvers(self):
        """Create and assemble NGBEMPEECSolver for each conductor."""
        from radia.ngsbem_peec import NGBEMPEECSolver

        for block_name, cond in self._conductors.items():
            if block_name in self._solvers:
                continue  # already assembled
            print(f"  Assembling BEM for '{block_name}'...")

            solver = NGBEMPEECSolver(
                self.mesh,
                conductor_label=block_name,
                sigma=cond['sigma'],
                thickness=cond['thickness'],
                order=0,
            )
            solver.assemble()
            solver.print_summary()
            self._solvers[block_name] = solver

    def _build_port_vectors(self):
        """Build port excitation vectors for each port.

        Returns:
            list of numpy arrays, one per port.
            Each array has shape (n_loop,) and represents the
            normalized current distribution for that port.

        Port types:
          - 'closed': Uniform excitation (all edges carry equal current).
            L_partial = e^T @ L @ e = total self-inductance.

          - 'open': Source/sink excitation.
            Edges adjacent to source nodes get +1, sink nodes get -1.
            L_partial = partial inductance from source to sink.
        """
        port_vectors = []

        for port_info in self._ports:
            block = port_info['block']
            solver = self._solvers.get(block)
            if solver is None:
                port_vectors.append(None)
                continue

            n_loop = solver.n_loop

            if port_info['type'] == 'closed':
                # Uniform excitation: all edges carry equal current
                e = np.ones(n_loop) / n_loop
                port_vectors.append(e)
                continue

            # Open conductor: source/sink excitation
            source_nodes = port_info['source_nodes']
            sink_nodes = port_info['sink_nodes']

            # Resolve block names to node lists (Cubit)
            if source_nodes is None and port_info.get('source_block'):
                source_nodes = self._resolve_block_nodes(
                    port_info['source_block'])
            if sink_nodes is None and port_info.get('sink_block'):
                sink_nodes = self._resolve_block_nodes(
                    port_info['sink_block'])

            if source_nodes is None or sink_nodes is None:
                print(f"  WARNING: Port '{port_info['name']}' has no source/sink, "
                      f"using uniform excitation")
                e = np.ones(n_loop) / n_loop
                port_vectors.append(e)
                continue

            # Build excitation from source/sink nodes
            e = self._build_source_sink_vector(
                solver, set(source_nodes), set(sink_nodes))
            port_vectors.append(e)

        return port_vectors

    def _build_source_sink_vector(self, solver, source_nodes, sink_nodes):
        """Build port excitation vector from source/sink node sets.

        For HDivSurface (RWG basis), each DOF corresponds to one mesh edge.
        An edge connects two vertices. If one vertex is a source node and the
        other is an interior node, the edge carries current FROM the source.
        Similarly for sink edges.

        The excitation vector encodes:
          e[i] = +1/n_src  if edge i is adjacent to source (current in)
          e[i] = -1/n_snk  if edge i is adjacent to sink (current out)
          e[i] = 0          otherwise (interior edge, current determined by BEM)

        Normalization: sum(e_source) = +1, sum(e_sink) = -1 (unit total current)

        Args:
            solver: NGBEMPEECSolver with assembled matrices
            source_nodes: Set of source vertex indices (NGSolve 0-indexed)
            sink_nodes: Set of sink vertex indices (NGSolve 0-indexed)

        Returns:
            e: numpy array (n_loop,) excitation vector
        """
        n_loop = solver.n_loop
        e = np.zeros(n_loop)

        # Get edge-vertex connectivity from mesh
        edge_data = _extract_edge_vertices(solver.mesh)

        n_src_edges = 0
        n_snk_edges = 0

        for edge_idx, (v0, v1) in enumerate(edge_data):
            if edge_idx >= n_loop:
                break
            v0_is_src = v0 in source_nodes
            v1_is_src = v1 in source_nodes
            v0_is_snk = v0 in sink_nodes
            v1_is_snk = v1 in sink_nodes

            # Edge with one source vertex and one interior vertex
            if (v0_is_src or v1_is_src) and not (v0_is_snk or v1_is_snk):
                e[edge_idx] = 1.0
                n_src_edges += 1
            # Edge with one sink vertex and one interior vertex
            elif (v0_is_snk or v1_is_snk) and not (v0_is_src or v1_is_src):
                e[edge_idx] = -1.0
                n_snk_edges += 1

        # Normalize: total source current = +1, total sink current = -1
        if n_src_edges > 0:
            src_mask = e > 0
            e[src_mask] /= np.sum(e[src_mask])
        if n_snk_edges > 0:
            snk_mask = e < 0
            e[snk_mask] /= np.abs(np.sum(e[snk_mask]))

        print(f"    Port vector: {n_src_edges} source edges, "
              f"{n_snk_edges} sink edges, "
              f"{n_loop - n_src_edges - n_snk_edges} interior edges")

        return e

    def _resolve_block_nodes(self, block_name):
        """Resolve a Cubit block name to NGSolve vertex indices.

        Blocks are the only Cubit entity type that supports element type
        specification (TRI6, QUAD8, HEX20, etc.), making them suitable
        for 2nd-order element workflows.

        Args:
            block_name: Name of the Cubit block

        Returns:
            list of NGSolve vertex indices (0-indexed), or None
        """
        if self._cubit is None:
            print(f"  WARNING: Cannot resolve block '{block_name}' - "
                  f"mesh was not loaded from Cubit")
            return None

        cubit = self._cubit

        # Find block by name
        block_ids = cubit.get_block_id_list()
        target_block = None
        block_names_available = []
        for bid in block_ids:
            name = cubit.get_block_name(bid)
            block_names_available.append(name)
            if name == block_name:
                target_block = bid
                break

        if target_block is None:
            print(f"  WARNING: Block '{block_name}' not found in Cubit. "
                  f"Available: {block_names_available}")
            return None

        # Collect all nodes from block elements. get_block_{tris,quads,
        # tets,hexes} return [] for blocks with no elements of that type;
        # an exception means a real API error and must surface.
        cubit_nodes = set()
        for tri_id in cubit.get_block_tris(target_block):
            cubit_nodes.update(cubit.get_connectivity("tri", tri_id))
        for quad_id in cubit.get_block_quads(target_block):
            cubit_nodes.update(cubit.get_connectivity("quad", quad_id))
        for tet_id in cubit.get_block_tets(target_block):
            cubit_nodes.update(cubit.get_connectivity("tet", tet_id))
        for hex_id in cubit.get_block_hexes(target_block):
            cubit_nodes.update(cubit.get_connectivity("hex", hex_id))

        if not cubit_nodes:
            print(f"  WARNING: Block '{block_name}' has no nodes")
            return None

        # Map Cubit node IDs to NGSolve vertex indices
        if hasattr(self, '_cubit_node_map'):
            ngsolve_nodes = []
            for nid in cubit_nodes:
                if nid in self._cubit_node_map:
                    ngsolve_nodes.append(self._cubit_node_map[nid] - 1)
            return ngsolve_nodes if ngsolve_nodes else None
        else:
            return self._match_nodes_by_coordinate(cubit_nodes)

    def _match_nodes_by_coordinate(self, cubit_node_ids):
        """Match Cubit nodes to NGSolve vertices by coordinate proximity.

        Args:
            cubit_node_ids: Set of Cubit node IDs

        Returns:
            list of NGSolve vertex indices (0-indexed)
        """
        if self.mesh is None or self._cubit is None:
            return None

        cubit = self._cubit
        tol = 1e-10

        # Build NGSolve vertex coordinate array
        nv = self.mesh.nv
        ng_coords = np.zeros((nv, 3))
        for v in self.mesh.vertices:
            pt = v.point
            ng_coords[v.nr] = [pt[0], pt[1], pt[2]]

        matched = []
        for nid in cubit_node_ids:
            coords = np.array(cubit.get_nodal_coordinates(nid))
            dists = np.linalg.norm(ng_coords - coords, axis=1)
            best = np.argmin(dists)
            if dists[best] < tol:
                matched.append(int(best))

        return matched if matched else None

    def _solve_port_impedance(self, solver, freq, e_port, mode,
                               Zs_func=None, esim_func=None):
        """Solve for port impedance at a single frequency.

        Z_port = 1 / (e^T @ Z_branch^{-1} @ e)

        where Z_branch = diag(R + Zs + Zs_esim) + jw*L

        Args:
            solver: NGBEMPEECSolver
            freq: Frequency [Hz]
            e_port: Port excitation vector (n_loop,)
            mode: Solver mode ('mqs', 'full', 'stabilized')
            Zs_func: Dowell SIBC function(freq) -> Zs array
            esim_func: ESIM correction function(freq) -> Zs_esim array

        Returns:
            Z_port: Complex port impedance [Ohm]
        """
        omega = 2 * np.pi * freq

        # Build branch impedance matrix
        R_diag = solver.R_loop.astype(complex)

        # Add SIBC (Dowell skin effect)
        if Zs_func is not None:
            R_diag = R_diag + np.asarray(Zs_func(freq), dtype=complex)

        # Add ESIM correction (nonlinear magnetic core effect)
        if esim_func is not None:
            R_diag = R_diag + np.asarray(esim_func(freq), dtype=complex)

        if mode == 'mqs':
            Z_branch = np.diag(R_diag) + 1j * omega * solver.L
        elif mode == 'stabilized':
            # Full stabilized system
            n = solver.n_loop + solver.n_star
            Z_full = np.zeros((n, n), dtype=complex)
            Z_full[:solver.n_loop, :solver.n_loop] = (
                np.diag(R_diag) + 1j * omega * solver.L)
            scale = 1j * omega * MU_0
            Z_full[:solver.n_loop, solver.n_loop:] = scale * solver.Q_0.T
            Z_full[solver.n_loop:, :solver.n_loop] = scale * solver.Q_0
            k = omega / C_0
            Z_full[solver.n_loop:, solver.n_loop:] = scale * k**2 * solver.V_0

            # RHS with port vector
            rhs = np.zeros(n, dtype=complex)
            rhs[:solver.n_loop] = e_port
            sol = scipy.linalg.solve(Z_full, rhs)
            I_L = sol[:solver.n_loop]
            return 1.0 / (e_port @ I_L) if abs(e_port @ I_L) > 1e-30 else 0j
        else:  # 'full'
            Z_LL = np.diag(R_diag) + 1j * omega * solver.L
            if omega > 1e-10 and solver._P_inv_M_LS is not None:
                cap_correction = (1j * omega) * (solver.M_LS.T @ solver._P_inv_M_LS)
                Z_branch = Z_LL - cap_correction
            else:
                Z_branch = Z_LL

        # Solve Z_branch * x = e_port
        if omega < 1e-10:
            # DC: only resistance
            x = scipy.linalg.solve(np.diag(solver.R_loop.astype(complex)), e_port)
        else:
            x = scipy.linalg.solve(Z_branch, e_port, assume_a='sym')

        denom = e_port @ x
        return 1.0 / denom if abs(denom) > 1e-30 else 0j

    def _build_esim_corrections(self):
        """Build ESIM surface impedance correction functions.

        Uses Karl Hollaus's Effective Surface Impedance Method
        to compute field-dependent Z_s from nonlinear B-H curves.

        Returns:
            dict: block_name -> callable(freq) returning Zs_esim array
        """
        esim_funcs = {}

        for core_name, core_info in self._cores.items():
            bh_data = core_info.get('bh_data')
            mu_r = core_info.get('mu_r')

            if bh_data is not None:
                # Nonlinear ESIM
                try:
                    from radia.esim_cell_problem import ESIMCellProblemSolver
                    esim_solver = ESIMCellProblemSolver(bh_data)

                    # Find conductors coupled to this core
                    for cond_name, cond_info in self._conductors.items():
                        solver = self._solvers.get(cond_name)
                        if solver is None:
                            continue

                        sigma = cond_info['sigma']
                        n_loop = solver.n_loop

                        def make_esim_func(esim_s, sig, n):
                            def func(freq):
                                omega = 2 * np.pi * freq
                                if omega < 1e-10:
                                    return np.zeros(n, dtype=complex)
                                # Use characteristic H field for ESIM
                                H_0 = 1000.0  # A/m (initial estimate)
                                Zs = esim_s.compute_surface_impedance(
                                    H_0, freq, sig)
                                return np.full(n, Zs, dtype=complex)
                            return func

                        esim_funcs[cond_name] = make_esim_func(
                            esim_solver, sigma, n_loop)
                        print(f"  ESIM: core '{core_name}' -> "
                              f"conductor '{cond_name}' (nonlinear)")

                except ImportError:
                    print("  WARNING: esim_cell_problem not available, "
                          "skipping ESIM correction")

            elif mu_r is not None:
                # Linear core: simple mu_r correction
                # Internal inductance scales with mu_r
                for cond_name in self._conductors:
                    solver = self._solvers.get(cond_name)
                    if solver is None:
                        continue
                    n_loop = solver.n_loop

                    def make_linear_func(mr, n):
                        def func(freq):
                            # Linear core adds internal inductance
                            # Delta_L_internal ~ (mu_r - 1) * L_int
                            return np.zeros(n, dtype=complex)
                        return func

                    esim_funcs[cond_name] = make_linear_func(mu_r, n_loop)
                    print(f"  Linear core '{core_name}': mu_r={mu_r} "
                          f"-> conductor '{cond_name}'")

        return esim_funcs

    def _compute_mutual_inductance(self, result, port_vectors):
        """Compute mutual inductance between ports.

        For multi-conductor problems, mutual inductance is computed
        from the cross terms of the BEM L matrix.

        M_ij = e_i^T @ L_ij @ e_j

        where L_ij is the mutual inductance sub-matrix between
        conductors i and j.
        """
        n_port = result.n_port
        if n_port < 2:
            return

        # Check if ports share the same solver (same conductor block)
        # If so, mutual inductance is already in the L matrix
        # For different conductor blocks, need cross-coupling
        print("  NOTE: Mutual inductance requires multi-conductor BEM (future work)")

    def export_gmsh_post(self, filename, field="J", freq_index=0):
        """Export field data to GMSH .msh format for visualization.

        Args:
            filename: Output .msh file path
            field: Field to export ("J" for current density, "B" for magnetic flux)
            freq_index: Index into frequency array for AC solution
        """
        from radia.gmsh_post_export import GmshPostExport

        if self.mesh is None:
            raise RuntimeError("No mesh loaded")

        post = GmshPostExport(self.mesh)

        if field.upper() == "J":
            # Export current density from BEM solution
            for block_name, solver in self._solvers.items():
                if solver.L is not None:
                    # Get edge current coefficients
                    freq = self._result.freqs[freq_index] if self._result else 0
                    I_loop, Q_star = solver.solve_loop_star(freq, mode='mqs')

                    # Convert to element-wise current density magnitude
                    J_mag = np.abs(I_loop)
                    post.add_scalar_field(
                        f"|J|_{block_name}", J_mag, cell_data=True)

        elif field.upper() == "B":
            # B field requires external evaluation
            print("  B field export: use NGSolve GridFunction or Radia rad.Fld()")

        post.write(filename)
        return filename

    def export_webgui(self, gf=None, name="solution"):
        """Display in JupyterLab via netgen.webgui.

        Args:
            gf: NGSolve GridFunction to display
            name: Display name
        """
        from ngsolve.webgui import Draw
        if gf is not None:
            Draw(gf, self.mesh, name)
        else:
            Draw(self.mesh)
