"""
Export methods for GmshBuilder.

Gmsh equivalents: export mesh, export step, export stl.
"""

import os
import tempfile


class ExportMixin:
    """Mixin providing mesh/geometry export."""

    def export(self, filename, version=2.2):
        """Export mesh to GMSH .msh file (Gmsh 'export mesh').

        Parameters
        ----------
        filename : str
            Output file path.
        version : float
            MSH format version (2.2 recommended for Radia/NGSolve).
        """
        import gmsh
        self._check_initialized()

        if not self._meshed:
            self.generate()

        gmsh.option.setNumber('Mesh.MshFileVersion', version)
        gmsh.write(filename)

        if self._verbose:
            print(f"[GmshBuilder] export -> {filename}")

    def merge_file(self, filename):
        """Merge a file into the current model.

        Supports STEP, BREP, MSH, STL, and other formats supported by GMSH.

        Parameters
        ----------
        filename : str
            File path to merge.
        """
        import gmsh
        self._check_initialized()

        gmsh.merge(filename)
        self._meshed = False

        if self._verbose:
            print(f"[GmshBuilder] merge_file <- {filename}")

    def to_radia(self, mu_r=1000, magnetization=None):
        """Generate mesh and convert to Radia objects.

        Parameters
        ----------
        mu_r : float or None
            Relative permeability. None for permanent magnets.
        magnetization : list, optional
            Initial magnetization [Mx, My, Mz] in A/m.

        Returns
        -------
        int
            Radia container object ID.
        """
        import gmsh
        from radia.gmsh_mesh_import import gmsh_to_radia

        self._check_initialized()

        if not self._meshed:
            self.generate()

        tmp_path = os.path.join(tempfile.gettempdir(), '_cubitmesh_temp.msh')
        gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
        gmsh.write(tmp_path)

        try:
            container = gmsh_to_radia(tmp_path, mu_r=mu_r,
                                       magnetization=magnetization)
        finally:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)

        if self._verbose:
            print(f"[GmshBuilder] to_radia -> container={container}")

        return container

    def to_radia_per_volume(self, materials=None):
        """Convert to Radia with per-volume material assignment.

        Parameters
        ----------
        materials : dict, optional
            {vol_id: {'mu_r': float, 'magnetization': [Mx,My,Mz]}}.
            Falls back to self._materials if not provided.

        Returns
        -------
        int
            Radia container object ID.
        """
        import gmsh
        from radia.gmsh_mesh_import import gmsh_to_radia

        self._check_initialized()

        if not self._meshed:
            self.generate()

        if materials is None:
            materials = {}
            for vid, mat in self._materials.items():
                materials[vid] = {'mu_r': mat.get('mu_r', 1000)}

        # If no per-volume materials defined, use default
        if not materials:
            return self.to_radia()

        # Auto-assign physical groups for material separation
        self.auto_assign_blocks()

        tmp_path = os.path.join(tempfile.gettempdir(),
                                '_cubitmesh_pervol_temp.msh')
        gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
        gmsh.write(tmp_path)

        try:
            import sys
            sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
            radia_objects = []

            for vid, mat_info in materials.items():
                if vid not in self._physical_groups:
                    continue
                pg_tag = self._physical_groups[vid][0]
                mu_r = mat_info.get('mu_r', 1000)
                mag = mat_info.get('magnetization', None)

                obj = gmsh_to_radia(tmp_path, mu_r=mu_r,
                                     magnetization=mag,
                                     physical_group=pg_tag)
                radia_objects.append(obj)

            if not radia_objects:
                container = gmsh_to_radia(tmp_path)
            elif len(radia_objects) == 1:
                container = radia_objects[0]
            else:
                try:
                    from radia._radia_pybind import ObjCnt
                except ImportError:
                    import radia as rad
                    ObjCnt = rad.ObjCnt
                container = ObjCnt(radia_objects)
        finally:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)

        if self._verbose:
            print(f"[GmshBuilder] to_radia_per_volume -> container={container}")

        return container

    def to_ngsolve_surface(self, label="conductor", mesh_size=None):
        """Generate surface mesh and convert to NGSolve Mesh for ngbem.

        Parameters
        ----------
        label : str
            BC name for conductor surface.
        mesh_size : float, optional
            Target element size.

        Returns
        -------
        ngsolve.Mesh
        """
        import gmsh
        from radia.gmsh_mesh_import import gmsh_surface_to_ngsolve

        self._check_initialized()

        if not self._meshed:
            self.generate_surface(mesh_size=mesh_size)

        tmp_path = os.path.join(tempfile.gettempdir(),
                                '_cubitmesh_surface_temp.msh')
        gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
        gmsh.write(tmp_path)

        try:
            mesh = gmsh_surface_to_ngsolve(tmp_path, label=label)
        finally:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)

        if self._verbose:
            print(f"[GmshBuilder] to_ngsolve_surface -> label='{label}'")
        return mesh

    def to_ngsolve_volume(self):
        """Export mesh and load as NGSolve volume mesh.

        Uses ReadGmsh to preserve physical group names as material labels.

        Returns
        -------
        ngsolve.Mesh
        """
        import gmsh
        from netgen.read_gmsh import ReadGmsh
        from ngsolve import Mesh

        self._check_initialized()

        if not self._meshed:
            self.generate()

        tmp_path = os.path.join(tempfile.gettempdir(),
                                '_gmshbuilder_ngsolve_temp.msh')
        gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
        gmsh.write(tmp_path)

        try:
            ngmesh = ReadGmsh(tmp_path)
            mesh = Mesh(ngmesh)
        finally:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)

        if self._verbose:
            print("[GmshBuilder] to_ngsolve_volume")
        return mesh

    def export_step(self, filename):
        """Export geometry to STEP file (Gmsh 'export step').

        Parameters
        ----------
        filename : str
            Output STEP file path.
        """
        import gmsh
        self._check_initialized()

        gmsh.write(filename)
        if self._verbose:
            print(f"[GmshBuilder] export_step -> {filename}")

    def export_brep(self, filename):
        """Export geometry to BREP file.

        Parameters
        ----------
        filename : str
            Output BREP file path.
        """
        import gmsh
        self._check_initialized()

        gmsh.write(filename)
        if self._verbose:
            print(f"[GmshBuilder] export_brep -> {filename}")

    def export_stl(self, filename):
        """Export mesh to STL file (Gmsh 'export stl').

        Parameters
        ----------
        filename : str
            Output STL file path.
        """
        import gmsh
        self._check_initialized()

        if not self._meshed:
            self.generate_surface()

        gmsh.write(filename)
        if self._verbose:
            print(f"[GmshBuilder] export_stl -> {filename}")

    def export_vtk(self, filename):
        """Export mesh to VTK file (Gmsh 'export vtk').

        Parameters
        ----------
        filename : str
            Output VTK file path.
        """
        import gmsh
        self._check_initialized()

        if not self._meshed:
            self.generate()

        gmsh.write(filename)
        if self._verbose:
            print(f"[GmshBuilder] export_vtk -> {filename}")

    def export_vtu(self, filename):
        """Export mesh to VTU file.

        Parameters
        ----------
        filename : str
            Output VTU file path.
        """
        import gmsh
        self._check_initialized()

        if not self._meshed:
            self.generate()

        gmsh.write(filename)
        if self._verbose:
            print(f"[GmshBuilder] export_vtu -> {filename}")

    def to_numpy(self):
        """Export mesh as numpy arrays.

        Returns
        -------
        dict
            'nodes': (n_nodes, 3), connectivity arrays by element type.
        """
        return self.get_elements_as_numpy()

    def to_pyvista(self):
        """Export mesh as PyVista UnstructuredGrid.

        Returns
        -------
        pyvista.UnstructuredGrid
        """
        import pyvista as pv
        import numpy as np

        self._check_initialized()

        if not self._meshed:
            self.generate()

        data = self.get_elements_as_numpy()
        points = data['node_coords']

        cells = []
        cell_types = []

        # VTK cell types
        VTK_HEXAHEDRON = 12
        VTK_TETRA = 10
        VTK_WEDGE = 13

        if 'hex' in data:
            for row in data['hex']:
                cells.append(np.concatenate([[8], row]))
                cell_types.append(VTK_HEXAHEDRON)

        if 'tet' in data:
            for row in data['tet']:
                cells.append(np.concatenate([[4], row]))
                cell_types.append(VTK_TETRA)

        if 'wedge' in data:
            for row in data['wedge']:
                cells.append(np.concatenate([[6], row]))
                cell_types.append(VTK_WEDGE)

        if cells:
            cells = np.concatenate(cells)
            cell_types = np.array(cell_types, dtype=np.uint8)
            grid = pv.UnstructuredGrid(cells, cell_types, points)
        else:
            grid = pv.UnstructuredGrid()

        return grid

    def show(self, dim=-1):
        """Show model in GMSH GUI (Gmsh 'display').

        Parameters
        ----------
        dim : int
            -1 for all, 2 for surface, 3 for volume.
        """
        import gmsh
        self._check_initialized()

        gmsh.fltk.run()

    def save_session(self, filename):
        """Save geometry as BREP (Gmsh 'save').

        Parameters
        ----------
        filename : str
            Output .brep file path.
        """
        self.export_brep(filename)

    # ------------------------------------------------------------------
    # Additional export methods
    # ------------------------------------------------------------------

    def export_msh(self, filename, version=4.1):
        """Export mesh to MSH file with specified format version.

        Unlike ``export()`` which defaults to MSH 2.2 for Radia/NGSolve
        compatibility, this method defaults to MSH 4.1 (latest).

        Parameters
        ----------
        filename : str
            Output .msh file path.
        version : float
            MSH format version (default 4.1).
        """
        import gmsh
        self._check_initialized()

        if not self._meshed:
            self.generate()

        gmsh.option.setNumber('Mesh.MshFileVersion', version)
        gmsh.write(filename)

        if self._verbose:
            print(f"[GmshBuilder] export_msh -> {filename} (v{version})")

    def export_mesh2d(self, filename, version=2.2):
        """Export surface (2D) mesh only.

        Generates a surface mesh if no mesh exists, then writes it out.

        Parameters
        ----------
        filename : str
            Output .msh file path.
        version : float
            MSH format version (default 2.2).
        """
        import gmsh
        self._check_initialized()

        if not self._meshed:
            self.generate_surface()

        gmsh.option.setNumber('Mesh.MshFileVersion', version)
        gmsh.write(filename)

        if self._verbose:
            print(f"[GmshBuilder] export_mesh2d -> {filename}")

    def export_ply(self, filename):
        """Export mesh to PLY file.

        Parameters
        ----------
        filename : str
            Output .ply file path.
        """
        import gmsh
        self._check_initialized()

        if not self._meshed:
            self.generate()

        gmsh.write(filename)

        if self._verbose:
            print(f"[GmshBuilder] export_ply -> {filename}")

    def export_pos(self, filename):
        """Export mesh to Gmsh POS file for visualization.

        Parameters
        ----------
        filename : str
            Output .pos file path.
        """
        import gmsh
        self._check_initialized()

        if not self._meshed:
            self.generate()

        gmsh.write(filename)

        if self._verbose:
            print(f"[GmshBuilder] export_pos -> {filename}")

    def export_geo(self, filename):
        """Export geometry as GEO file (Gmsh geometry script).

        Parameters
        ----------
        filename : str
            Output .geo or .geo_unrolled file path.
        """
        import gmsh
        self._check_initialized()

        gmsh.write(filename)

        if self._verbose:
            print(f"[GmshBuilder] export_geo -> {filename}")

    def export_unv(self, filename):
        """Export mesh to IDEAS Universal (.unv) format.

        Parameters
        ----------
        filename : str
            Output .unv file path.
        """
        import gmsh
        self._check_initialized()

        if not self._meshed:
            self.generate()

        gmsh.write(filename)

        if self._verbose:
            print(f"[GmshBuilder] export_unv -> {filename}")

    def to_dict(self):
        """Export mesh as a plain Python dictionary.

        Returns
        -------
        dict
            ``'nodes'`` : list of [x, y, z] coordinate lists.
            ``'elements'`` : dict mapping element type name to list of
            node-index lists (0-based).
            ``'volumes'`` : dict mapping volume id to its GMSH tags.
        """
        import gmsh
        self._check_initialized()

        if not self._meshed:
            self.generate()

        # Nodes
        tags, coords, _ = gmsh.model.mesh.getNodes()
        n_nodes = len(tags)
        coords_arr = list(zip(coords[0::3], coords[1::3], coords[2::3]))
        tag_to_idx = {int(t): i for i, t in enumerate(tags)}

        nodes = [[c[0], c[1], c[2]] for c in coords_arr]

        # Elements (3D)
        elements = {}
        elem_types, elem_tags, elem_nodes = gmsh.model.mesh.getElements(3)
        for i, et in enumerate(elem_types):
            name, _, _, n_per_elem, _, _ = \
                gmsh.model.mesh.getElementProperties(et)
            raw_nodes = elem_nodes[i]
            n_elems = len(elem_tags[i])

            conns = []
            for e in range(n_elems):
                offset = e * n_per_elem
                conn = [tag_to_idx.get(int(raw_nodes[offset + k]),
                                       int(raw_nodes[offset + k]))
                        for k in range(n_per_elem)]
                conns.append(conn)

            # Friendly name
            if 'Hexahedron' in name:
                key = 'hex'
            elif 'Tetrahedron' in name:
                key = 'tet'
            elif 'Prism' in name or 'Wedge' in name:
                key = 'wedge'
            else:
                key = name
            elements[key] = conns

        # Also collect 2D elements if no 3D found
        if not elements:
            elem_types2, elem_tags2, elem_nodes2 = \
                gmsh.model.mesh.getElements(2)
            for i, et in enumerate(elem_types2):
                name, _, _, n_per_elem, _, _ = \
                    gmsh.model.mesh.getElementProperties(et)
                raw_nodes = elem_nodes2[i]
                n_elems = len(elem_tags2[i])
                conns = []
                for e in range(n_elems):
                    offset = e * n_per_elem
                    conn = [tag_to_idx.get(int(raw_nodes[offset + k]),
                                           int(raw_nodes[offset + k]))
                            for k in range(n_per_elem)]
                    conns.append(conn)
                if 'Triangle' in name:
                    key = 'tri'
                elif 'Quadrilateral' in name:
                    key = 'quad'
                else:
                    key = name
                elements[key] = conns

        # Volume info
        volumes = {}
        if hasattr(self, '_volumes'):
            for vid, vtags in self._volumes.items():
                volumes[vid] = list(vtags)

        return {
            'nodes': nodes,
            'elements': elements,
            'volumes': volumes,
        }

    def to_meshio(self):
        """Export mesh via meshio.

        Writes a temporary MSH file, reads it with meshio, and returns
        the resulting ``meshio.Mesh`` object.

        Returns
        -------
        meshio.Mesh
        """
        import gmsh
        import meshio

        self._check_initialized()

        if not self._meshed:
            self.generate()

        tmp_path = os.path.join(tempfile.gettempdir(),
                                '_cubitmesh_meshio_temp.msh')
        gmsh.option.setNumber('Mesh.MshFileVersion', 4.1)
        gmsh.write(tmp_path)

        try:
            mesh = meshio.read(tmp_path)
        finally:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)

        if self._verbose:
            print("[GmshBuilder] to_meshio -> meshio.Mesh")

        return mesh

    def to_radia_with_material(self, vol_id, mu_r=1000, magnetization=None):
        """Export a single volume to Radia with material properties.

        Parameters
        ----------
        vol_id : int
            Volume identifier (as returned by ``add_box`` etc.).
        mu_r : float
            Relative permeability for linear soft-iron.
        magnetization : list, optional
            Magnetization [Mx, My, Mz] in A/m for permanent magnets.

        Returns
        -------
        int
            Radia container object ID.
        """
        import gmsh
        from radia.gmsh_mesh_import import gmsh_to_radia

        self._check_initialized()

        if not self._meshed:
            self.generate()

        # Ensure the volume has a physical group
        if not hasattr(self, '_physical_groups') or \
                vol_id not in self._physical_groups:
            self.auto_assign_blocks()

        if vol_id not in self._physical_groups:
            raise ValueError(
                f"Volume {vol_id} not found in physical groups")

        pg_tag = self._physical_groups[vol_id][0]

        tmp_path = os.path.join(tempfile.gettempdir(),
                                '_cubitmesh_radia_vol_temp.msh')
        gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
        gmsh.write(tmp_path)

        try:
            container = gmsh_to_radia(tmp_path, mu_r=mu_r,
                                       magnetization=magnetization,
                                       physical_group=pg_tag)
        finally:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)

        if self._verbose:
            print(f"[GmshBuilder] to_radia_with_material vol={vol_id} "
                  f"mu_r={mu_r} -> container={container}")

        return container

    def export_abaqus(self, filename):
        """Export mesh to Abaqus INP file.

        Parameters
        ----------
        filename : str
            Output .inp file path.
        """
        import gmsh
        self._check_initialized()

        if not self._meshed:
            self.generate()

        gmsh.write(filename)

        if self._verbose:
            print(f"[GmshBuilder] export_abaqus -> {filename}")

    def export_cgns(self, filename):
        """Export mesh to CGNS file.

        Parameters
        ----------
        filename : str
            Output .cgns file path.
        """
        import gmsh
        self._check_initialized()

        if not self._meshed:
            self.generate()

        gmsh.write(filename)

        if self._verbose:
            print(f"[GmshBuilder] export_cgns -> {filename}")

    def export_med(self, filename):
        """Export mesh to MED file.

        Parameters
        ----------
        filename : str
            Output .med file path.
        """
        import gmsh
        self._check_initialized()

        if not self._meshed:
            self.generate()

        gmsh.write(filename)

        if self._verbose:
            print(f"[GmshBuilder] export_med -> {filename}")

    def get_export_formats(self):
        """Return list of supported export format strings.

        Returns
        -------
        list of str
            Supported file extensions (without leading dot).
        """
        return [
            'msh', 'step', 'brep', 'stl', 'vtk', 'vtu',
            'ply', 'inp', 'unv', 'med', 'cgns', 'geo',
        ]

    def export_auto(self, filename):
        """Auto-detect export format from file extension and export.

        Parameters
        ----------
        filename : str
            Output file path. Extension determines format.

        Raises
        ------
        ValueError
            If the file extension is not recognized.
        """
        ext = os.path.splitext(filename)[1].lower().lstrip('.')

        dispatch = {
            'msh': self.export_msh,
            'step': self.export_step,
            'stp': self.export_step,
            'brep': self.export_brep,
            'stl': self.export_stl,
            'vtk': self.export_vtk,
            'vtu': self.export_vtu,
            'ply': self.export_ply,
            'pos': self.export_pos,
            'geo': self.export_geo,
            'geo_unrolled': self.export_geo,
            'unv': self.export_unv,
            'inp': self.export_abaqus,
            'cgns': self.export_cgns,
            'med': self.export_med,
        }

        handler = dispatch.get(ext)
        if handler is None:
            supported = ', '.join(sorted(dispatch.keys()))
            raise ValueError(
                f"Unsupported extension '.{ext}'. "
                f"Supported: {supported}")

        handler(filename)

    def to_scipy_sparse(self):
        """Export mesh connectivity as a scipy sparse adjacency matrix.

        Builds a symmetric node-to-node adjacency matrix from element
        connectivity (edges shared by elements). Diagonal entries are
        zero.

        Returns
        -------
        scipy.sparse.csr_matrix
            Sparse adjacency matrix of shape (n_nodes, n_nodes).
        """
        import numpy as np
        from scipy.sparse import csr_matrix

        self._check_initialized()

        if not self._meshed:
            self.generate()

        data = self.get_elements_as_numpy()
        n_nodes = len(data['node_coords'])

        rows = []
        cols = []

        for key in ('hex', 'tet', 'wedge', 'tri', 'quad'):
            if key not in data:
                continue
            conn = data[key]
            n_elems, n_per_elem = conn.shape
            for e in range(n_elems):
                nodes = conn[e]
                # Add edges between all pairs of nodes in element
                for a in range(n_per_elem):
                    for b in range(a + 1, n_per_elem):
                        ni, nj = int(nodes[a]), int(nodes[b])
                        rows.append(ni)
                        cols.append(nj)
                        rows.append(nj)
                        cols.append(ni)

        if rows:
            rows = np.array(rows, dtype=np.int32)
            cols = np.array(cols, dtype=np.int32)
            vals = np.ones(len(rows), dtype=np.float64)
            adj = csr_matrix((vals, (rows, cols)),
                             shape=(n_nodes, n_nodes))
            # Deduplicate (sum -> clip to 1)
            adj.data[:] = 1.0
        else:
            adj = csr_matrix((n_nodes, n_nodes), dtype=np.float64)

        if self._verbose:
            print(f"[GmshBuilder] to_scipy_sparse -> ({n_nodes}, {n_nodes}), "
                  f"nnz={adj.nnz}")

        return adj
