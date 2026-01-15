"""
mmm_ngsolve.py - NGSolve CoefficientFunction wrapper for mmm_core

This module provides NGSolve CoefficientFunction integration for the
standalone MMM (Magnetic Moment Method) solver.

Features:
- MMMFieldCF: CoefficientFunction for B/H field evaluation
- Direct NGSolve mesh integration via MMMBuilder
- Efficient batch field evaluation
- Nonlinear material support with B-H curve
- Permanent magnet support with fixed remanent magnetization

Usage:
    from ngsolve import *
    from mmm_ngsolve import MMMFieldCF, MMMSystem

    # Create mesh
    mesh = Mesh(...)

    # Linear material:
    system = MMMSystem.from_ngsolve_mesh(mesh, mu_r=1000)
    system.solve([0, 0, 1e5])

    # Nonlinear material with B-H curve:
    bh_curve = [[0, 0], [100, 0.5], [1000, 1.2], [10000, 1.8]]
    system = MMMSystem.from_ngsolve_mesh(mesh, bh_curve=bh_curve)
    result = system.solve_nonlinear([0, 0, 1e5])

    # Permanent magnet (no solve needed):
    pm_system = MMMSystem.from_ngsolve_mesh(pm_mesh, Mr=[0, 0, 954930])
    B_cf = pm_system.get_b_field_cf()  # Direct field from fixed magnetization

    # Mixed PM + soft iron:
    system = MMMSystem()
    system.add_pm_from_mesh(pm_mesh, Mr=[0, 0, 954930])
    system.add_soft_from_mesh(iron_mesh, mu_r=1000)
    system.build()
    system.solve()  # Solves for soft iron magnetization only

    # Get CoefficientFunction for field
    B_cf = system.get_b_field_cf()

Part of Radia project - pybind11 module: mmm_core
"""

import numpy as np

# Import mmm_core (C++ module)
try:
    from . import mmm_core
except ImportError:
    import mmm_core


class MMMSystem:
    """
    MMM (Magnetic Moment Method) system for magnetostatic analysis.

    Provides a high-level interface for:
    - Building MMM interaction matrix from NGSolve mesh
    - Solving the MMM linear system
    - Computing B/H fields at arbitrary points
    - NGSolve CoefficientFunction integration
    - Permanent magnet support (fixed remanent magnetization)

    Example:
        from ngsolve import Mesh
        from netgen.occ import Box, OCCGeometry
        from mmm_ngsolve import MMMSystem

        # Create mesh
        geo = OCCGeometry(Box((0,0,0), (1,1,1)))
        mesh = Mesh(geo.GenerateMesh(maxh=0.2))

        # Build MMM system with soft iron
        system = MMMSystem.from_ngsolve_mesh(mesh, mu_r=1000)
        system.solve([0, 0, 1e5])  # H_ext = 100 kA/m in z

        # Get B field CoefficientFunction
        B_cf = system.get_b_field_cf()

        # Permanent magnet example (no solve needed)
        pm_system = MMMSystem.from_ngsolve_mesh(pm_mesh, Mr=[0, 0, 954930])
        B_cf = pm_system.get_b_field_cf()

        # Mixed PM + soft iron
        mixed = MMMSystem()
        mixed.add_pm_from_mesh(pm_mesh, Mr=[0, 0, 954930])
        mixed.add_soft_from_mesh(iron_mesh, mu_r=1000)
        mixed.build()
        mixed.solve()
    """

    def __init__(self):
        """Create empty MMM system. Use from_* class methods to construct."""
        self.builder = mmm_core.MMMBuilder()
        self.solver = mmm_core.MMMSolver()
        self.field_computer = mmm_core.MMMFieldComputer()
        self.nonlinear_solver = None  # Created when bh_curve is set
        self.hacapk_solver = None     # Created when HACApK method is used

        self.N = None           # Interaction matrix
        self.dof_offset = None  # DOF offset array
        self.mu_r = None        # Relative permeability (scalar or array)
        self.bh_curve = None    # B-H curve for nonlinear material
        self.M = None           # Solved magnetization
        self.chi = None         # Solved chi per element (nonlinear)

        self._is_built = False
        self._is_solved = False
        self._is_nonlinear = False
        self._is_pm_only = False  # True if system has only PM elements (no solve needed)
        self._has_pm = False      # True if system has any PM elements
        self._has_soft = False    # True if system has any soft magnetic elements

    @classmethod
    def from_ngsolve_mesh(cls, mesh, mu_r=1.0, bh_curve=None, material_filter=None, Mr=None):
        """
        Create MMM system from NGSolve mesh.

        Args:
            mesh: NGSolve Mesh object
            mu_r: Relative permeability (scalar or dict {material_name: mu_r})
                  Ignored if bh_curve or Mr is provided.
            bh_curve: B-H curve as [[H1, B1], [H2, B2], ...] for nonlinear material
                      H in A/m, B in Tesla. If provided, enables nonlinear solver.
            material_filter: Optional material name filter (str or list)
            Mr: Remanent magnetization [Mx, My, Mz] in A/m for permanent magnet.
                If provided, creates PM elements (no solve needed for field computation).

        Returns:
            MMMSystem instance ready for solving (soft iron) or field computation (PM)

        Note:
            For permanent magnets (Mr provided), no solve() is needed - the field
            can be computed directly using get_b_field_cf() or compute_b_field().
        """
        system = cls()

        # Extract mesh data
        vertices, elements, element_materials = cls._extract_mesh(mesh, material_filter)

        # Determine if this is a permanent magnet system
        if Mr is not None:
            # Permanent magnet with fixed remanent magnetization
            Mr = np.asarray(Mr, dtype=np.float64)
            if Mr.shape != (3,):
                raise ValueError("Mr must be [Mx, My, Mz] shape (3,)")

            # Add PM elements to builder
            if elements.shape[1] == 4:
                system.builder.add_pm_tetrahedra_from_mesh(vertices, elements, Mr)
            elif elements.shape[1] == 8:
                system.builder.add_pm_hexahedra_from_mesh(vertices, elements, Mr)
            else:
                raise ValueError(f"Unsupported element type with {elements.shape[1]} vertices")

            system._has_pm = True
            system._is_pm_only = True  # Only PM, no soft elements

        else:
            # Soft magnetic material
            if elements.shape[1] == 4:
                system.builder.add_tetrahedra_from_mesh(vertices, elements)
            elif elements.shape[1] == 8:
                system.builder.add_hexahedra_from_mesh(vertices, elements)
            else:
                raise ValueError(f"Unsupported element type with {elements.shape[1]} vertices")

            system._has_soft = True

            # Store mu_r or bh_curve
            if bh_curve is not None:
                system.bh_curve = np.asarray(bh_curve, dtype=np.float64)
                system._is_nonlinear = True
                system.mu_r = None
            elif isinstance(mu_r, dict):
                n_elem = len(elements)
                system.mu_r = np.ones(n_elem)
                for i, mat_name in enumerate(element_materials):
                    if mat_name in mu_r:
                        system.mu_r[i] = mu_r[mat_name]
            else:
                system.mu_r = float(mu_r)

        # Build interaction matrix
        system.N, system.dof_offset = system.builder.build()
        system.field_computer.set_elements_from_builder(system.builder)

        # For PM-only system, mark as solved (no linear solve needed)
        if system._is_pm_only:
            # Create M array with fixed magnetization from PM elements
            # PM elements have dof=0, but field computation uses Mr directly
            system.M = np.zeros(system.builder.total_dof, dtype=np.float64)
            system._is_solved = True  # Field computation ready without solve
        else:
            # Setup solver for soft iron
            system.solver.set_matrix(system.N, system.dof_offset)

            if system._is_nonlinear:
                system.nonlinear_solver = mmm_core.MMMNonlinearSolver()
                system.nonlinear_solver.set_matrix(system.N, system.dof_offset)
                system.nonlinear_solver.set_bh_curve(system.bh_curve)

        system._is_built = True
        return system

    @classmethod
    def from_arrays(cls, vertices, elements, mu_r=1.0, bh_curve=None, Mr=None):
        """
        Create MMM system from numpy arrays.

        Args:
            vertices: (n_verts, 3) vertex coordinates [m]
            elements: (n_elem, 4) or (n_elem, 8) vertex indices (0-based)
            mu_r: Relative permeability (scalar or array per element)
                  Ignored if bh_curve or Mr is provided.
            bh_curve: B-H curve as [[H1, B1], [H2, B2], ...] for nonlinear material
            Mr: Remanent magnetization [Mx, My, Mz] in A/m for permanent magnet.

        Returns:
            MMMSystem instance ready for solving (soft iron) or field computation (PM)
        """
        system = cls()

        vertices = np.asarray(vertices, dtype=np.float64)
        elements = np.asarray(elements, dtype=np.int32)

        if Mr is not None:
            # Permanent magnet
            Mr = np.asarray(Mr, dtype=np.float64)
            if Mr.shape != (3,):
                raise ValueError("Mr must be [Mx, My, Mz] shape (3,)")

            if elements.shape[1] == 4:
                system.builder.add_pm_tetrahedra_from_mesh(vertices, elements, Mr)
            elif elements.shape[1] == 8:
                system.builder.add_pm_hexahedra_from_mesh(vertices, elements, Mr)
            else:
                raise ValueError(f"Unsupported element type with {elements.shape[1]} vertices")

            system._has_pm = True
            system._is_pm_only = True
        else:
            # Soft magnetic material
            if elements.shape[1] == 4:
                system.builder.add_tetrahedra_from_mesh(vertices, elements)
            elif elements.shape[1] == 8:
                system.builder.add_hexahedra_from_mesh(vertices, elements)
            else:
                raise ValueError(f"Unsupported element type with {elements.shape[1]} vertices")

            system._has_soft = True

            if bh_curve is not None:
                system.bh_curve = np.asarray(bh_curve, dtype=np.float64)
                system._is_nonlinear = True
                system.mu_r = None
            else:
                system.mu_r = mu_r

        # Build
        system.N, system.dof_offset = system.builder.build()
        system.field_computer.set_elements_from_builder(system.builder)

        if system._is_pm_only:
            system.M = np.zeros(system.builder.total_dof, dtype=np.float64)
            system._is_solved = True
        else:
            system.solver.set_matrix(system.N, system.dof_offset)

            if system._is_nonlinear:
                system.nonlinear_solver = mmm_core.MMMNonlinearSolver()
                system.nonlinear_solver.set_matrix(system.N, system.dof_offset)
                system.nonlinear_solver.set_bh_curve(system.bh_curve)

        system._is_built = True
        return system

    @staticmethod
    def _extract_mesh(mesh, material_filter=None):
        """
        Extract vertex and element data from NGSolve mesh.

        Args:
            mesh: NGSolve Mesh object
            material_filter: Optional material name filter

        Returns:
            vertices: (n_verts, 3) array
            elements: (n_elem, n_verts_per_elem) array
            element_materials: list of material names
        """
        from ngsolve import VOL, ET

        # Get vertices
        nv = mesh.nv
        vertices = np.zeros((nv, 3), dtype=np.float64)
        for i, v in enumerate(mesh.vertices):
            pt = v.point
            vertices[i] = [pt[0], pt[1], pt[2]]

        # Get 3D elements using VOL region
        elements_list = []
        element_materials = []

        for el in mesh.Elements(VOL):
            mat_name = el.mat

            # Apply material filter if specified
            if material_filter is not None:
                if isinstance(material_filter, str):
                    if mat_name != material_filter:
                        continue
                elif isinstance(material_filter, (list, tuple)):
                    if mat_name not in material_filter:
                        continue

            # Check element type
            if el.type == ET.TET:
                n_verts = 4
            elif el.type == ET.HEX:
                n_verts = 8
            else:
                # Skip unsupported element types
                continue

            # Extract vertex indices (0-based)
            verts = [v.nr for v in el.vertices]
            if len(verts) != n_verts:
                raise ValueError(f"Element has unexpected vertex count: {len(verts)} != {n_verts}")
            elements_list.append(verts)
            element_materials.append(mat_name)

        if len(elements_list) == 0:
            raise ValueError("No elements found (check material_filter)")

        elements = np.array(elements_list, dtype=np.int32)

        return vertices, elements, element_materials

    def add_pm_from_mesh(self, mesh, Mr, material_filter=None):
        """
        Add permanent magnet elements from NGSolve mesh.

        Use this method to build mixed PM + soft iron systems.

        Args:
            mesh: NGSolve Mesh object
            Mr: Remanent magnetization [Mx, My, Mz] in A/m
            material_filter: Optional material name filter (str or list)

        Example:
            system = MMMSystem()
            system.add_pm_from_mesh(pm_mesh, Mr=[0, 0, 954930])
            system.add_soft_from_mesh(iron_mesh, mu_r=1000)
            system.build()
            system.solve()
        """
        if self._is_built:
            raise RuntimeError("System already built. Create a new MMMSystem.")

        vertices, elements, _ = self._extract_mesh(mesh, material_filter)
        Mr = np.asarray(Mr, dtype=np.float64)
        if Mr.shape != (3,):
            raise ValueError("Mr must be [Mx, My, Mz] shape (3,)")

        if elements.shape[1] == 4:
            self.builder.add_pm_tetrahedra_from_mesh(vertices, elements, Mr)
        elif elements.shape[1] == 8:
            self.builder.add_pm_hexahedra_from_mesh(vertices, elements, Mr)
        else:
            raise ValueError(f"Unsupported element type with {elements.shape[1]} vertices")

        self._has_pm = True

    def add_soft_from_mesh(self, mesh, mu_r=1000, material_filter=None):
        """
        Add soft magnetic elements from NGSolve mesh.

        Use this method to build mixed PM + soft iron systems.

        Args:
            mesh: NGSolve Mesh object
            mu_r: Relative permeability (scalar)
            material_filter: Optional material name filter (str or list)

        Example:
            system = MMMSystem()
            system.add_pm_from_mesh(pm_mesh, Mr=[0, 0, 954930])
            system.add_soft_from_mesh(iron_mesh, mu_r=1000)
            system.build()
            system.solve()
        """
        if self._is_built:
            raise RuntimeError("System already built. Create a new MMMSystem.")

        vertices, elements, _ = self._extract_mesh(mesh, material_filter)

        if elements.shape[1] == 4:
            self.builder.add_tetrahedra_from_mesh(vertices, elements)
        elif elements.shape[1] == 8:
            self.builder.add_hexahedra_from_mesh(vertices, elements)
        else:
            raise ValueError(f"Unsupported element type with {elements.shape[1]} vertices")

        self._has_soft = True
        if self.mu_r is None:
            self.mu_r = float(mu_r)
        # Note: For mixed materials, mu_r should be per-element (future enhancement)

    def add_pm_from_arrays(self, vertices, elements, Mr):
        """
        Add permanent magnet elements from numpy arrays.

        Args:
            vertices: (n_verts, 3) vertex coordinates [m]
            elements: (n_elem, 4) or (n_elem, 8) vertex indices (0-based)
            Mr: Remanent magnetization [Mx, My, Mz] in A/m
        """
        if self._is_built:
            raise RuntimeError("System already built. Create a new MMMSystem.")

        vertices = np.asarray(vertices, dtype=np.float64)
        elements = np.asarray(elements, dtype=np.int32)
        Mr = np.asarray(Mr, dtype=np.float64)
        if Mr.shape != (3,):
            raise ValueError("Mr must be [Mx, My, Mz] shape (3,)")

        if elements.shape[1] == 4:
            self.builder.add_pm_tetrahedra_from_mesh(vertices, elements, Mr)
        elif elements.shape[1] == 8:
            self.builder.add_pm_hexahedra_from_mesh(vertices, elements, Mr)
        else:
            raise ValueError(f"Unsupported element type with {elements.shape[1]} vertices")

        self._has_pm = True

    def add_soft_from_arrays(self, vertices, elements, mu_r=1000):
        """
        Add soft magnetic elements from numpy arrays.

        Args:
            vertices: (n_verts, 3) vertex coordinates [m]
            elements: (n_elem, 4) or (n_elem, 8) vertex indices (0-based)
            mu_r: Relative permeability (scalar)
        """
        if self._is_built:
            raise RuntimeError("System already built. Create a new MMMSystem.")

        vertices = np.asarray(vertices, dtype=np.float64)
        elements = np.asarray(elements, dtype=np.int32)

        if elements.shape[1] == 4:
            self.builder.add_tetrahedra_from_mesh(vertices, elements)
        elif elements.shape[1] == 8:
            self.builder.add_hexahedra_from_mesh(vertices, elements)
        else:
            raise ValueError(f"Unsupported element type with {elements.shape[1]} vertices")

        self._has_soft = True
        if self.mu_r is None:
            self.mu_r = float(mu_r)

    def build(self):
        """
        Build interaction matrix after adding elements.

        Call this after using add_pm_from_mesh/add_soft_from_mesh.

        Example:
            system = MMMSystem()
            system.add_pm_from_mesh(pm_mesh, Mr=[0, 0, 954930])
            system.add_soft_from_mesh(iron_mesh, mu_r=1000)
            system.build()
            system.solve()
        """
        if self._is_built:
            raise RuntimeError("System already built.")

        if self.builder.num_elements == 0:
            raise RuntimeError("No elements added. Use add_pm_from_mesh or add_soft_from_mesh first.")

        # Build interaction matrix
        self.N, self.dof_offset = self.builder.build()
        self.field_computer.set_elements_from_builder(self.builder)

        # Determine system type
        self._is_pm_only = self._has_pm and not self._has_soft

        if self._is_pm_only:
            # PM-only: no solve needed
            self.M = np.zeros(self.builder.total_dof, dtype=np.float64)
            self._is_solved = True
        else:
            # Has soft elements: setup solver
            self.solver.set_matrix(self.N, self.dof_offset)

        self._is_built = True

    def solve(self, H_ext=None, method='lu', tol=1e-8, max_iter=1000,
              hacapk_params=None):
        """
        Solve MMM system with external field.

        For PM-only systems, no solve is needed and this method returns immediately.
        For mixed PM + soft iron systems, the soft iron magnetization is computed
        considering both H_ext and the field from PM elements.

        Args:
            H_ext: External field [Hx, Hy, Hz] in A/m, or (n_elem, 3) array.
                   For PM-only systems, this can be None.
                   For mixed systems, if None, only PM field acts on soft iron.
            method: 'lu', 'bicgstab', or 'hacapk'
                    - 'lu': Dense LU decomposition (LAPACK dgesv)
                    - 'bicgstab': Iterative BiCGSTAB solver
                    - 'hacapk': H-matrix accelerated BiCGSTAB (O(N log N))
            tol: Convergence tolerance for iterative methods
            max_iter: Maximum iterations for iterative methods
            hacapk_params: HACApK parameters dict (optional)
                           - aca_eps: ACA+ tolerance (default: 1e-4)
                           - leaf_size: Minimum cluster size (default: 32)
                           - eta: Admissibility parameter (default: 2.0)
                           - max_rank: Maximum ACA rank (default: 200)
                           - print_level: 0=silent, 1=standard (default: 0)

        Returns:
            M: Magnetization vector (total_dof,) in A/m
        """
        if not self._is_built:
            raise RuntimeError("System not built. Use from_ngsolve_mesh or from_arrays first.")

        # PM-only system: already solved (field from fixed Mr)
        if self._is_pm_only:
            self._is_solved = True
            return self.M

        n_elem = self.builder.num_elements
        total_dof = self.builder.total_dof

        # Compute 1/chi from mu_r
        # chi = mu_r - 1, so 1/chi = 1/(mu_r - 1)
        if self.mu_r is None:
            raise RuntimeError("mu_r not set. For soft iron, provide mu_r.")

        if np.isscalar(self.mu_r):
            chi = self.mu_r - 1.0
            if chi <= 0:
                raise ValueError("mu_r must be > 1 for magnetic materials")
            inv_chi = np.full(n_elem, 1.0 / chi, dtype=np.float64)
        else:
            chi = np.array(self.mu_r) - 1.0
            if np.any(chi <= 0):
                raise ValueError("mu_r must be > 1 for all elements")
            inv_chi = 1.0 / chi

        # Handle H_ext (can be None for mixed systems where PM provides field)
        if H_ext is None:
            H_ext = np.zeros(3, dtype=np.float64)

        # Expand H_ext to total_dof
        H_ext = np.asarray(H_ext, dtype=np.float64)
        if H_ext.ndim == 1 and len(H_ext) == 3:
            # Uniform field: broadcast to all DOF
            H_ext_full = np.zeros(total_dof, dtype=np.float64)
            for i in range(n_elem):
                dof_start = self.dof_offset[i]
                dof_end = self.dof_offset[i + 1]
                n_dof = dof_end - dof_start
                if n_dof == 3:
                    # Tetrahedron: direct copy
                    H_ext_full[dof_start:dof_end] = H_ext
                elif n_dof == 0:
                    # PM element: no DOF
                    pass
                else:
                    # Hexahedron (6 DOF): need to compute normal components
                    # For now, simplified: use uniform value
                    H_ext_full[dof_start:dof_end] = np.tile(H_ext[:2], 2)[:n_dof]
            H_ext = H_ext_full
        elif H_ext.shape[0] != total_dof:
            raise ValueError(f"H_ext must have shape (3,) or ({total_dof},)")

        # Solve
        if method == 'lu':
            self.M = self.solver.solve_lu(inv_chi, H_ext, chi_per_element=True)
        elif method == 'bicgstab':
            self.M, iterations = self.solver.solve_bicgstab(
                inv_chi, H_ext, tol, max_iter, chi_per_element=True
            )
        elif method == 'hacapk':
            # HACApK H-matrix solver is not yet available in mmm_core
            # Fall back to BiCGSTAB for now
            raise NotImplementedError(
                "HACApK solver is not yet available in mmm_core. "
                "Use method='lu' or method='bicgstab' instead. "
                "HACApK integration is planned for a future release."
            )
        else:
            raise ValueError(f"Unknown method: {method}. Use 'lu' or 'bicgstab'.")

        self._is_solved = True
        return self.M

    def _solve_hacapk(self, inv_chi, H_ext, tol, max_iter, hacapk_params):
        """
        Internal method to solve using HACApK H-matrix acceleration.

        NOTE: HACApK solver is not yet available in mmm_core.
        This method is a placeholder for future implementation.
        """
        raise NotImplementedError(
            "HACApK solver is not yet available in mmm_core. "
            "Use method='lu' or method='bicgstab' instead."
        )

    def get_hacapk_stats(self):
        """
        Get H-matrix statistics (after solving with method='hacapk').

        NOTE: HACApK solver is not yet available in mmm_core.
        """
        raise NotImplementedError(
            "HACApK solver is not yet available in mmm_core. "
            "H-matrix statistics are not available."
        )

    def solve_nonlinear(self, H_ext, tol=1e-4, max_iter=100, relax=1.0,
                        linear_solver='lu'):
        """
        Solve MMM system with nonlinear material (B-H curve).

        Uses Newton-Raphson iteration to find self-consistent solution where
        chi depends on local H field magnitude.

        Args:
            H_ext: External field [Hx, Hy, Hz] in A/m
            tol: Convergence tolerance for B-field relative change
            max_iter: Maximum Newton-Raphson iterations
            relax: Relaxation factor (0 < relax <= 1). Use < 1 for stability.
            linear_solver: 'lu' (default) or 'bicgstab' for inner linear solve

        Returns:
            dict with keys:
                'M': Magnetization vector (total_dof,) in A/m
                'chi': Chi per element (n_elem,)
                'iterations': Number of iterations
                'converged': True if converged
                'error': Final convergence error

        Raises:
            RuntimeError: If system not built or not configured for nonlinear
        """
        if not self._is_built:
            raise RuntimeError("System not built. Use from_ngsolve_mesh or from_arrays first.")

        if not self._is_nonlinear:
            raise RuntimeError("System not configured for nonlinear solve. "
                               "Provide bh_curve in from_ngsolve_mesh or from_arrays.")

        # Set solver parameters
        solver_code = 0 if linear_solver == 'lu' else 1
        self.nonlinear_solver.set_parameters(tol, max_iter, relax, solver_code)

        # Convert H_ext to numpy array
        H_ext = np.asarray(H_ext, dtype=np.float64)
        if H_ext.ndim != 1 or len(H_ext) != 3:
            raise ValueError("H_ext must be [Hx, Hy, Hz] for nonlinear solve")

        # Solve
        result = self.nonlinear_solver.solve(H_ext)

        # Store results
        self.M = result['M']
        self.chi = result['chi']
        self._is_solved = True

        return result

    def get_local_h(self):
        """
        Get local H field magnitude per element (after nonlinear solve).

        Returns:
            H_local: (n_elements,) local H magnitude in A/m
        """
        if not self._is_nonlinear or self.nonlinear_solver is None:
            raise RuntimeError("Only available after solve_nonlinear()")
        return self.nonlinear_solver.get_local_h()

    def get_local_b(self):
        """
        Get local B field magnitude per element (after nonlinear solve).

        Returns:
            B_local: (n_elements,) local B magnitude in Tesla
        """
        if not self._is_nonlinear or self.nonlinear_solver is None:
            raise RuntimeError("Only available after solve_nonlinear()")
        return self.nonlinear_solver.get_local_b()

    def compute_b_field(self, obs_points):
        """
        Compute B field at observation points.

        Args:
            obs_points: (n_points, 3) observation points [m]

        Returns:
            B: (n_points, 3) magnetic flux density [T]
        """
        if not self._is_solved:
            raise RuntimeError("System not solved. Call solve() first.")

        obs_points = np.asarray(obs_points, dtype=np.float64)
        if obs_points.ndim == 1:
            obs_points = obs_points.reshape(1, 3)

        return self.field_computer.compute_b_field(self.M, obs_points)

    def compute_h_field(self, obs_points):
        """
        Compute H field at observation points.

        Args:
            obs_points: (n_points, 3) observation points [m]

        Returns:
            H: (n_points, 3) magnetic field intensity [A/m]
        """
        if not self._is_solved:
            raise RuntimeError("System not solved. Call solve() first.")

        obs_points = np.asarray(obs_points, dtype=np.float64)
        if obs_points.ndim == 1:
            obs_points = obs_points.reshape(1, 3)

        return self.field_computer.compute_h_field(self.M, obs_points)

    def get_b_field_cf(self):
        """
        Get NGSolve CoefficientFunction for B field.

        Returns:
            MMMFieldCF: CoefficientFunction for B field
        """
        if not self._is_solved:
            raise RuntimeError("System not solved. Call solve() first.")

        return MMMFieldCF(self, field_type='b')

    def get_h_field_cf(self):
        """
        Get NGSolve CoefficientFunction for H field.

        Returns:
            MMMFieldCF: CoefficientFunction for H field
        """
        if not self._is_solved:
            raise RuntimeError("System not solved. Call solve() first.")

        return MMMFieldCF(self, field_type='h')

    @property
    def num_elements(self):
        """Number of elements"""
        return self.builder.num_elements

    @property
    def total_dof(self):
        """Total degrees of freedom"""
        return self.builder.total_dof

    @property
    def has_pm(self):
        """True if system contains permanent magnet elements"""
        return self._has_pm

    @property
    def has_soft(self):
        """True if system contains soft magnetic elements"""
        return self._has_soft

    @property
    def is_pm_only(self):
        """True if system contains only PM elements (no solve needed)"""
        return self._is_pm_only

    def get_element_centers(self):
        """Get element center coordinates (n_elements, 3) in meters"""
        return self.field_computer.get_element_centers()

    def get_element_volumes(self):
        """Get element volumes (n_elements,) in m^3"""
        return self.field_computer.get_element_volumes()

    def info(self):
        """
        Print system information.

        Returns:
            str: System information summary
        """
        lines = [
            f"MMMSystem:",
            f"  Elements: {self.num_elements}",
            f"  Total DOF: {self.total_dof}",
            f"  Has PM: {self._has_pm}",
            f"  Has Soft: {self._has_soft}",
            f"  PM Only: {self._is_pm_only}",
            f"  Built: {self._is_built}",
            f"  Solved: {self._is_solved}",
        ]
        if self._is_nonlinear:
            lines.append(f"  Nonlinear: True (B-H curve)")
        if self.mu_r is not None:
            if np.isscalar(self.mu_r):
                lines.append(f"  mu_r: {self.mu_r}")
            else:
                lines.append(f"  mu_r: array ({len(self.mu_r)} elements)")
        info_str = "\n".join(lines)
        print(info_str)
        return info_str


class MMMFieldCF:
    """
    NGSolve CoefficientFunction wrapper for MMM field computation.

    This class implements a CoefficientFunction-like interface that can be
    used with NGSolve's Set() and Integrate() methods.

    Note: This is a pure Python implementation. For maximum performance,
    consider using the field values directly with numpy arrays.
    """

    def __init__(self, system, field_type='b'):
        """
        Create field CoefficientFunction.

        Args:
            system: MMMSystem instance (must be solved)
            field_type: 'b' for B field (Tesla), 'h' for H field (A/m)
        """
        self.system = system
        self.field_type = field_type.lower()

        if self.field_type not in ('b', 'h'):
            raise ValueError("field_type must be 'b' or 'h'")

        # For NGSolve compatibility
        self.dim = 3

    def __call__(self, mesh_point):
        """
        Evaluate field at a point.

        This method is called by NGSolve when evaluating the CoefficientFunction.

        Args:
            mesh_point: MappedIntegrationPoint or (x, y, z) coordinates

        Returns:
            Field value at the point
        """
        # Extract coordinates
        if hasattr(mesh_point, 'point'):
            # NGSolve MappedIntegrationPoint
            pt = mesh_point.point
            x, y, z = pt[0], pt[1], pt[2]
        elif hasattr(mesh_point, '__len__') and len(mesh_point) >= 3:
            x, y, z = mesh_point[0], mesh_point[1], mesh_point[2]
        else:
            raise ValueError("Cannot extract coordinates from input")

        obs_point = np.array([[x, y, z]], dtype=np.float64)

        if self.field_type == 'b':
            field = self.system.compute_b_field(obs_point)
        else:
            field = self.system.compute_h_field(obs_point)

        return tuple(field[0])

    def Evaluate(self, obs_points):
        """
        Batch evaluate field at multiple points.

        Args:
            obs_points: (n_points, 3) array of observation points

        Returns:
            (n_points, 3) array of field values
        """
        obs_points = np.asarray(obs_points, dtype=np.float64)
        if obs_points.ndim == 1:
            obs_points = obs_points.reshape(1, 3)

        if self.field_type == 'b':
            return self.system.compute_b_field(obs_points)
        else:
            return self.system.compute_h_field(obs_points)

    def __repr__(self):
        field_name = "B" if self.field_type == 'b' else "H"
        return f"MMMFieldCF({field_name} field, {self.system.num_elements} elements)"


# Convenience function for quick setup
def create_mmm_from_mesh(mesh, mu_r=1000, H_ext=None, material_filter=None):
    """
    Convenience function to create and solve MMM system from NGSolve mesh.

    Args:
        mesh: NGSolve Mesh object
        mu_r: Relative permeability
        H_ext: External field [Hx, Hy, Hz] in A/m (optional)
        material_filter: Optional material name filter

    Returns:
        MMMSystem: Solved system (if H_ext provided) or ready system
    """
    system = MMMSystem.from_ngsolve_mesh(mesh, mu_r=mu_r, material_filter=material_filter)

    if H_ext is not None:
        system.solve(H_ext)

    return system
