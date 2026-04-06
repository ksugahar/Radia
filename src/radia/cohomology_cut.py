"""Total scalar potential solver with Gmsh cohomology cuts.

Implements the total scalar potential formulation for magnetostatics:

    H = -grad(phi) + sum_k I_k h_k

where phi is the scalar potential (H1), h_k are cohomology basis
functions (HCurl, curl-free), and I_k = N_k * I_coil are the
ampere-turns for each coil.

The cohomology basis functions h_k are computed automatically by Gmsh's
built-in homology/cohomology solver.  Each h_k satisfies:

    curl(h_k) = 0          (irrotational)
    oint h_k . dl = delta_jk   (unit circulation around loop j)

This eliminates the need for:
    - Biot-Savart source field computation
    - VoxelCoefficient interpolation (no Gibbs artifacts at coil surfaces)
    - Explicit cut surface definition (automatic via cohomology)

Weak form (I_k known):
    integral(mu * grad(phi) . grad(v)) = sum_k I_k * integral(mu * h_k . grad(v))

This is a standard H1 Poisson problem with a topological source term.

Usage:
    import gmsh
    from radia.cohomology_cut import CohomologyCutSolver

    solver = CohomologyCutSolver()
    solver.setup_from_gmsh()      # after user defines Gmsh geometry
    solver.solve([100.0], {'iron': 1000.0})  # NI=100, mu_r=1000
    B_cf = solver.get_B()

Reference:
    Pellikka et al., "Homology and Cohomology Computation in Finite
    Element Modeling," SIAM J. Sci. Comput., 2013.
"""

import numpy as np

MU_0 = 4 * np.pi * 1e-7


class CohomologyCutSolver:
    """Total scalar potential solver with automatic cohomology cuts.

    Workflow:
        1. User creates Gmsh geometry (with coil holes, iron regions)
        2. Call setup_from_gmsh() to mesh, compute cohomology, transfer to NGSolve
        3. Call solve() with ampere-turns and material properties
        4. Access results via get_H(), get_B(), get_phi()
    """

    def __init__(self):
        self._mesh = None
        self._h_basis = []
        self._phi_gf = None
        self._H_cf = None
        self._B_cf = None
        self._mu_cf = None
        self._gmsh_to_ngs_vertex = None

    def setup_from_gmsh(self, domain_physical_name='domain',
                        boundary_physical_name='outer'):
        """Mesh, compute cohomology, and transfer to NGSolve.

        Call this after defining Gmsh geometry and physical groups.
        The Gmsh model must have:
          - A 3D physical group for the domain (air + iron, excluding coil)
          - 2D physical groups for boundary surfaces

        Parameters
        ----------
        domain_physical_name : str
            Name of the 3D physical group defining the computation domain.
        boundary_physical_name : str
            Name of the 2D physical group for outer boundary (Dirichlet BC).

        Returns
        -------
        int
            Number of cohomology generators (= number of independent coils).
        """
        import gmsh

        # Find domain physical group tag
        domain_tag = self._find_physical_group(3, domain_physical_name)

        # Request cohomology (H^1 generators)
        gmsh.model.mesh.addHomologyRequest('Cohomology', [domain_tag], [], [1])

        # Mesh
        gmsh.model.mesh.generate(3)

        # Compute cohomology
        dimTags = gmsh.model.mesh.computeHomology()

        if not dimTags:
            raise RuntimeError(
                "No cohomology generators found. "
                "The domain may be simply connected (no coil holes).")

        # Extract cohomology edge chains from Gmsh
        cohom_chains = self._extract_cohomology_chains(dimTags)

        # Transfer mesh to NGSolve
        self._transfer_mesh_to_ngsolve()

        # Map cohomology edges to NGSolve HCurl DOFs
        self._build_hcurl_basis(cohom_chains)

        n_coils = len(self._h_basis)
        return n_coils

    def solve(self, NI_list, mu_r_dict=None, order=2,
              dirichlet='outer', kelvin_region=None, kelvin_radius=None,
              kelvin_center=None):
        """Solve the total scalar potential problem.

        Parameters
        ----------
        NI_list : list of float
            Ampere-turns for each coil [NI_1, NI_2, ...].
            Length must match number of cohomology generators.
        mu_r_dict : dict or None
            Material name -> relative permeability, e.g. {'iron': 1000}.
            Regions not listed default to mu_r=1 (air).
        order : int
            H1 polynomial order (default: 2).
        dirichlet : str
            Boundary label for Dirichlet BC (phi=0).
        kelvin_region : str or None
            Material name for the Kelvin transform shell region.
            If set, applies weight (R/r)^2 to the bilinear form there.
        kelvin_radius : float or None
            Kelvin sphere radius R (in meters). Required if kelvin_region is set.
        kelvin_center : tuple or None
            Center of the Kelvin sphere (x0, y0, z0). Default: origin.
        """
        from ngsolve import (H1, GridFunction, BilinearForm, LinearForm,
                             grad, dx, x, y, z, sqrt)

        if self._mesh is None:
            raise RuntimeError("Call setup_from_gmsh() first")

        if len(NI_list) != len(self._h_basis):
            raise ValueError(
                f"NI_list has {len(NI_list)} entries but "
                f"{len(self._h_basis)} cohomology generators exist")

        if kelvin_region and not kelvin_radius:
            raise ValueError("kelvin_radius required when kelvin_region is set")

        # Build mu CoefficientFunction
        self._mu_cf = self._build_mu_cf(mu_r_dict)

        # FEM spaces
        fes = H1(self._mesh, order=order, dirichlet=dirichlet)
        phi = fes.TrialFunction()
        v = fes.TestFunction()

        # Choose solver strategy based on problem size
        ndofs = fes.ndof
        use_iterative = ndofs > 200000

        # Bilinear form: a(phi, v) = integral(mu * grad(phi) . grad(v))
        a = BilinearForm(fes)
        if kelvin_region:
            # Physical region: standard bilinear form
            all_mats = list(set(self._mesh.GetMaterials()))
            phys_mats = [m for m in all_mats if m != kelvin_region]
            for mat in phys_mats:
                a += self._mu_cf * grad(phi) * grad(v) * dx(mat)
            # Kelvin shell: weight = (R/r)^2 (3D Kelvin transform)
            cx, cy, cz = kelvin_center if kelvin_center else (0, 0, 0)
            r_cf = sqrt((x - cx)**2 + (y - cy)**2 + (z - cz)**2)
            R = kelvin_radius
            kelvin_weight = (R * R) / (r_cf * r_cf + 1e-30)
            a += MU_0 * kelvin_weight * grad(phi) * grad(v) * dx(kelvin_region)
        else:
            a += self._mu_cf * grad(phi) * grad(v) * dx

        # BDDC preconditioner must be registered before Assemble
        pre = None
        if use_iterative:
            from ngsolve import Preconditioner
            pre = Preconditioner(a, 'bddc')
        a.Assemble()

        # Source: f(v) = sum_k NI_k * integral(mu * h_k . grad(v))
        # Source only in physical region (not in Kelvin shell)
        f = LinearForm(fes)
        if kelvin_region:
            for k, NI_k in enumerate(NI_list):
                if abs(NI_k) > 0:
                    for mat in phys_mats:
                        f += NI_k * self._mu_cf * self._h_basis[k] * grad(v) * dx(mat)
        else:
            for k, NI_k in enumerate(NI_list):
                if abs(NI_k) > 0:
                    f += NI_k * self._mu_cf * self._h_basis[k] * grad(v) * dx
        f.Assemble()

        # Solve
        self._phi_gf = GridFunction(fes, name='phi_total')
        if use_iterative:
            from ngsolve.krylovspace import CGSolver
            print(f"   Using BDDC+CG ({ndofs} DOFs)...", flush=True)
            inv = CGSolver(mat=a.mat, pre=pre.mat, maxiter=2000,
                           printrates=False, tol=1e-10)
            self._phi_gf.vec.data = inv * f.vec
            print(f"   CG converged ({ndofs} DOFs)", flush=True)
        else:
            self._phi_gf.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

        # Build result fields: H = -grad(phi) + sum NI_k * h_k
        from ngsolve import CF
        H_source = CF((0, 0, 0))
        for k, NI_k in enumerate(NI_list):
            if abs(NI_k) > 0:
                H_source = H_source + NI_k * self._h_basis[k]
        self._H_cf = -grad(self._phi_gf) + H_source
        self._B_cf = self._mu_cf * self._H_cf

        return self._phi_gf

    def solve_nonlinear(self, NI_list, bh_data, order=2,
                        dirichlet='outer', iron_domain='iron',
                        tol=1e-4, maxiter=50, relax=1.0):
        """Solve with nonlinear B-H curve (Newton-Raphson).

        Parameters
        ----------
        NI_list : list of float
            Ampere-turns for each coil.
        bh_data : list of [H, B] pairs
            B-H curve data [[H1, B1], [H2, B2], ...].
        order : int
            H1 polynomial order.
        dirichlet : str
            Boundary label for Dirichlet BC.
        iron_domain : str
            Material name for iron region.
        tol : float
            Convergence tolerance.
        maxiter : int
            Maximum iterations.
        relax : float
            Under-relaxation factor (0 < relax <= 1).
        """
        from ngsolve import (H1, GridFunction, BilinearForm, LinearForm,
                             grad, dx, Integrate, InnerProduct, CF)

        if self._mesh is None:
            raise RuntimeError("Call setup_from_gmsh() first")

        bh = np.array(bh_data)
        H_data, B_data = bh[:, 0], bh[:, 1]

        fes = H1(self._mesh, order=order, dirichlet=dirichlet)

        self._phi_gf = GridFunction(fes, name='phi_total')
        phi_new = GridFunction(fes)

        mu_r_current = 1000.0  # initial guess

        for iteration in range(maxiter):
            mu_r_dict = {iron_domain: mu_r_current}
            self._mu_cf = self._build_mu_cf(mu_r_dict)

            phi = fes.TrialFunction()
            v = fes.TestFunction()

            a = BilinearForm(fes)
            a += self._mu_cf * grad(phi) * grad(v) * dx
            a.Assemble()

            f = LinearForm(fes)
            for k, NI_k in enumerate(NI_list):
                if abs(NI_k) > 0:
                    f += NI_k * self._mu_cf * self._h_basis[k] * grad(v) * dx
            f.Assemble()

            phi_new.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

            # Under-relaxation
            diff_vec = phi_new.vec.CreateVector()
            diff_vec.data = phi_new.vec - self._phi_gf.vec
            self._phi_gf.vec.data += relax * diff_vec

            # Update H and compute mu_r from B-H curve
            H_source = CF((0, 0, 0))
            for k, NI_k in enumerate(NI_list):
                if abs(NI_k) > 0:
                    H_source = H_source + NI_k * self._h_basis[k]
            H_cf = -grad(self._phi_gf) + H_source
            H_mag_sq = Integrate(InnerProduct(H_cf, H_cf) * dx(iron_domain),
                                 self._mesh)
            vol_iron = Integrate(CF(1) * dx(iron_domain), self._mesh)
            H_avg = np.sqrt(H_mag_sq / max(vol_iron, 1e-30))

            B_interp = np.interp(H_avg, H_data, B_data)
            mu_r_new = B_interp / (MU_0 * H_avg) if H_avg > 1 else 1000.0
            mu_r_new = max(1.0, mu_r_new)

            change = abs(mu_r_new - mu_r_current) / max(mu_r_current, 1.0)
            mu_r_current = mu_r_new

            if change < tol:
                break

        self._mu_cf = self._build_mu_cf({iron_domain: mu_r_current})
        self._H_cf = -grad(self._phi_gf) + H_source
        self._B_cf = self._mu_cf * self._H_cf

        return self._phi_gf

    # ------------------------------------------------------------------
    # Result access
    # ------------------------------------------------------------------

    def get_H(self):
        """Get H field CoefficientFunction (A/m, dim=3)."""
        if self._H_cf is None:
            raise RuntimeError("Call solve() first")
        return self._H_cf

    def get_B(self):
        """Get B field CoefficientFunction (Tesla, dim=3)."""
        if self._B_cf is None:
            raise RuntimeError("Call solve() first")
        return self._B_cf

    def get_phi(self):
        """Get scalar potential GridFunction."""
        if self._phi_gf is None:
            raise RuntimeError("Call solve() first")
        return self._phi_gf

    def get_mesh(self):
        """Get the NGSolve mesh."""
        return self._mesh

    def get_cohomology_basis(self):
        """Get list of HCurl GridFunctions (cohomology basis)."""
        return self._h_basis

    def project_to_hdiv(self, order=None):
        """Project B field onto HDiv GridFunction (ensures div(B)=0).

        Parameters
        ----------
        order : int, optional
            HDiv polynomial order. Default: 2.

        Returns
        -------
        ngsolve.GridFunction
            B field in HDiv space.
        """
        from ngsolve import HDiv, GridFunction

        if self._B_cf is None:
            raise RuntimeError("Call solve() first")
        if order is None:
            order = 2
        fes_hdiv = HDiv(self._mesh, order=order)
        B_gf = GridFunction(fes_hdiv, name='B')
        B_gf.Set(self._B_cf)
        return B_gf

    # ------------------------------------------------------------------
    # Internal methods
    # ------------------------------------------------------------------

    def _find_physical_group(self, dim, name):
        """Find Gmsh physical group tag by name."""
        import gmsh
        pgs = gmsh.model.getPhysicalGroups(dim)
        for d, tag in pgs:
            pg_name = gmsh.model.getPhysicalName(d, tag)
            if pg_name == name:
                return tag
        raise ValueError(f"Physical group '{name}' (dim={dim}) not found")

    def _extract_cohomology_chains(self, dimTags):
        """Extract cohomology edge chains from Gmsh.

        Returns list of lists of (node1, node2) tuples.
        """
        import gmsh

        chains = []
        for dim, tag in dimTags:
            if dim != 1:
                continue
            entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
            edges = []
            for ent in entities:
                elem_types, _, node_tags_list = gmsh.model.mesh.getElements(1, ent)
                for i, _ in enumerate(elem_types):
                    nt = node_tags_list[i]
                    n_elems = len(nt) // 2
                    for j in range(n_elems):
                        edges.append((int(nt[2*j]), int(nt[2*j+1])))
            chains.append(edges)
        return chains

    def _transfer_mesh_to_ngsolve(self):
        """Export Gmsh mesh and import into NGSolve via .vol file."""
        import gmsh
        import tempfile
        import os

        tmpfile = os.path.join(tempfile.gettempdir(), '_cohom_mesh.vol')
        gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
        # Write .msh, convert to .vol via Netgen, then load
        msh_tmp = os.path.join(tempfile.gettempdir(), '_cohom_mesh.msh')
        gmsh.write(msh_tmp)

        try:
            from netgen.meshing import Mesh as NetgenMesh
            from ngsolve import Mesh
            ngmesh = NetgenMesh()
            ngmesh.Load(msh_tmp)
            ngmesh.Save(tmpfile)
            self._mesh = Mesh(tmpfile)
        finally:
            for f in [msh_tmp, tmpfile]:
                if os.path.exists(f):
                    os.remove(f)

        # Build Gmsh -> NGSolve vertex mapping
        self._build_vertex_map()

    def _build_vertex_map(self):
        """Build Gmsh node ID -> NGSolve vertex index mapping."""
        import gmsh
        from scipy.spatial import cKDTree

        # Gmsh node coordinates
        node_ids, coords, _ = gmsh.model.mesh.getNodes()
        gmsh_coords = {}
        for i, nid in enumerate(node_ids):
            gmsh_coords[int(nid)] = np.array(
                [coords[3*i], coords[3*i+1], coords[3*i+2]])

        # NGSolve vertex coordinates
        ngs_coords = np.zeros((self._mesh.nv, 3))
        for i, v in enumerate(self._mesh.vertices):
            pt = v.point
            ngs_coords[i] = [pt[0], pt[1], pt[2]]

        # KDTree matching
        tree = cKDTree(ngs_coords)
        self._gmsh_to_ngs_vertex = {}
        for gmsh_id, coord in gmsh_coords.items():
            dist, idx = tree.query(coord)
            if dist < 1e-8:
                self._gmsh_to_ngs_vertex[gmsh_id] = idx

        n_mapped = len(self._gmsh_to_ngs_vertex)
        n_total = len(gmsh_coords)
        if n_mapped < n_total:
            raise RuntimeError(
                f"Only {n_mapped}/{n_total} Gmsh nodes mapped to NGSolve. "
                "Mesh transfer may have failed.")

    def _build_hcurl_basis(self, cohom_chains):
        """Map Gmsh cohomology edge chains to NGSolve HCurl GridFunctions."""
        from ngsolve import HCurl, GridFunction

        # Build NGSolve edge lookup: (min_v, max_v) -> edge_nr
        ngs_edge_map = {}
        for edge in self._mesh.edges:
            verts = edge.vertices
            v0, v1 = int(verts[0].nr), int(verts[1].nr)
            key = (min(v0, v1), max(v0, v1))
            ngs_edge_map[key] = int(edge.nr)

        fes_hcurl = HCurl(self._mesh, order=0)

        self._h_basis = []
        for chain in cohom_chains:
            h_gf = GridFunction(fes_hcurl)
            h_gf.vec[:] = 0

            mapped = 0
            for gn1, gn2 in chain:
                if gn1 in self._gmsh_to_ngs_vertex and gn2 in self._gmsh_to_ngs_vertex:
                    nv1 = self._gmsh_to_ngs_vertex[gn1]
                    nv2 = self._gmsh_to_ngs_vertex[gn2]
                    key = (min(nv1, nv2), max(nv1, nv2))
                    if key in ngs_edge_map:
                        edge_nr = ngs_edge_map[key]
                        sign = 1.0 if nv1 < nv2 else -1.0
                        h_gf.vec[edge_nr] = sign
                        mapped += 1

            if mapped < len(chain):
                raise RuntimeError(
                    f"Only {mapped}/{len(chain)} cohomology edges mapped. "
                    "Gmsh-NGSolve edge matching failed.")

            self._h_basis.append(h_gf)

    def _build_mu_cf(self, mu_r_dict):
        """Build domain-wise mu CoefficientFunction."""
        if mu_r_dict is None:
            from ngsolve import CF
            return CF(MU_0)

        iron_dict = {name: MU_0 * mu_r for name, mu_r in mu_r_dict.items()}
        return self._mesh.MaterialCF(iron_dict, default=MU_0)
