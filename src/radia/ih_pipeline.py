"""ih_pipeline.py

High-level CAD-first IH pipeline. Caches the workpiece BEM-SIBC solver
once, then evaluates many coil designs against it cheaply.

The expensive operators (DL, SL, M, K) depend only on the workpiece
mesh and ngsolve H1 order -- NOT on the coil geometry or operating
frequency. For parametric / LLM-driven coil sweeps, building the
context once and calling .evaluate(coil) per design drops the per-
design cost from ~14 s to ~3 s on the 166-node ih_bem_sample
workpiece.

Usage:
    from radia.ih_pipeline import IHWorkpieceContext

    ctx = IHWorkpieceContext(
        wp_mesh, frequency=7e3, sigma=2e6, mu_r=100.0)

    for coil in candidate_coils:
        r = ctx.evaluate(coil, nw=5, nh=8)
        print(r['H_t_rms'], r['P_total'])
"""

import time
import numpy as np

from radia.bem_sibc_solver import (ScalarBIESIBCSolver,
                                    compute_phi_inc_from_filaments)
from radia.peec_bundle import build_bundle_solver, filament_currents

MU_0 = 4.0 * np.pi * 1e-7
DEFAULT_CU_SIGMA = 5.8e7


class IHWorkpieceContext:
    """Cached workpiece state for cheap per-coil IH evaluation.

    The BEM-SIBC solver (O(N_vert^2) assembly of DL, SL, M, K) is
    built once in __init__. Each .evaluate(coil) call runs PEEC on
    the coil filaments, projects phi_inc on the workpiece surface
    and does one linear SIBC solve -- no BEM rebuild.

    Attributes:
        mesh: NGSolve workpiece mesh.
        wp_nodes: (N, 3) surface node coordinates for phi_inc.
        frequency, omega: operating frequency [Hz] and [rad/s].
        sigma, mu_r: workpiece material.
        delta: skin depth in the workpiece [m].
        Z_s: complex surface impedance used in the SIBC term.
        solver: cached ScalarBIESIBCSolver (the expensive object).
        t_build: BEM assembly time [s] for reporting.
    """

    def __init__(self, mesh, frequency, sigma, mu_r, order=1,
                 wp_nodes=None):
        """Initialize context and assemble BEM operators.

        Args:
            mesh: NGSolve Mesh (workpiece surface mesh).
            frequency: operating frequency [Hz].
            sigma: workpiece conductivity [S/m].
            mu_r: workpiece relative permeability.
            order: H1 polynomial order (default 1).
            wp_nodes: optional (N, 3) observation nodes for phi_inc.
                If None, extracted from the mesh vertices.
        """
        self.mesh = mesh
        self.frequency = float(frequency)
        self.omega = 2.0 * np.pi * self.frequency
        self.sigma = float(sigma)
        self.mu_r = float(mu_r)
        self.order = order

        if wp_nodes is None:
            wp_nodes = np.array(
                [[mesh.vertices[i].point[j] for j in range(3)]
                 for i in range(mesh.nv)])
        self.wp_nodes = np.asarray(wp_nodes, dtype=float)

        mu = MU_0 * self.mu_r
        rho = 1.0 / self.sigma
        delta = np.sqrt(2.0 * rho / (self.omega * mu))
        self.delta = delta
        self.Z_s = complex(1.0, 1.0) * rho / delta

        t0 = time.perf_counter()
        self.solver = ScalarBIESIBCSolver(mesh, order=order)
        self.t_build = time.perf_counter() - t0

    def evaluate(self, coil, nw, nh, sigma_coil=DEFAULT_CU_SIGMA,
                 I_total=None, n_quad=20):
        """Evaluate one coil design against the cached workpiece.

        Args:
            coil: CoilBuilder instance (already assembled).
            nw, nh: filament grid density on the coil cross-section.
            sigma_coil: coil conductor conductivity [S/m].
            I_total: total port current [A]. Default: coil.current.
            n_quad: Gauss-Legendre order for phi_inc path integrals.

        Returns:
            dict with keys
                H_t_rms, P_total, Q_total -- physics metrics
                area                       -- workpiece surface area
                paths, I_peec              -- filament geometry + AC I
                n_fil, n_loop              -- solver size
                timings                    -- per-step ms breakdown
        """
        timings = {}
        if I_total is None:
            I_total = coil.current

        t0 = time.perf_counter()
        paths, _ = coil.to_filaments(nw=nw, nh=nh,
                                     frequency=self.frequency,
                                     sigma=sigma_coil)
        cell_wh = coil._last_filament_info['cell_wh']
        timings['to_filaments_ms'] = (time.perf_counter() - t0) * 1e3

        t0 = time.perf_counter()
        bundle, seg_of_fil, *_ = build_bundle_solver(
            paths, dw=None, dh=None, sigma=sigma_coil, cell_wh=cell_wh)
        timings['build_bundle_ms'] = (time.perf_counter() - t0) * 1e3

        t0 = time.perf_counter()
        I_branch = bundle.compute_branch_currents(self.frequency,
                                                  [complex(I_total)])
        I_peec = filament_currents(I_branch, seg_of_fil)
        timings['branch_currents_ms'] = (time.perf_counter() - t0) * 1e3

        t0 = time.perf_counter()
        phi = compute_phi_inc_from_filaments(
            self.wp_nodes, paths, I_peec, n_quad=n_quad)
        timings['phi_inc_ms'] = (time.perf_counter() - t0) * 1e3

        t0 = time.perf_counter()
        result = self.solver.solve(phi, Z_s=self.Z_s, omega=self.omega)
        timings['sibc_solve_ms'] = (time.perf_counter() - t0) * 1e3

        area = result['area']
        H_rms = result['H_t_rms']
        P_tot = result['P_density'] * area
        Q_tot = 0.5 * self.Z_s.imag * H_rms ** 2 * area
        timings['total_ms'] = sum(timings.values())

        return {
            'H_t_rms': H_rms,
            'P_total': P_tot,
            'Q_total': Q_tot,
            'area': area,
            'paths': paths,
            'I_peec': I_peec,
            'n_fil': len(paths),
            'n_loop': bundle.n_loop,
            'timings': timings,
        }
