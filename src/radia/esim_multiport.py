"""
esim_multiport.py — Multi-port ESIM surface impedance tensor

Extends the scalar ESIM (Z_s = E_t / H_t, single direction) to a
two-port framework with a 2×2 surface impedance tensor:

  [E_tx]   [Z_xx  Z_xy] [H_tx]
  [E_ty] = [Z_yx  Z_yy] [H_ty]

Port 1: x-tangential direction (rolling direction for electrical steel)
Port 2: y-tangential direction (transverse direction)

Dual (admittance) formulation:
  [H_tx]   [Y_xx  Y_xy] [E_tx]
  [H_ty] = [Y_yx  Y_yy] [E_ty]

  Y_s = Z_s^{-1}

  Use dual formulation when |Z_s| >> 1 (thick conductor, high frequency)
  for better numerical conditioning.

Circuit representation (two-port network):
  [V1]   [Z11  Z12] [I1]    V_i = E_t_i (electric tangential = "voltage")
  [V2] = [Z21  Z22] [I2]    I_i = H_t_i (magnetic tangential = "current")

  Z_ij = E_i / H_j  with all other ports open (H_k=0 for k≠j)

NGSolve Robin BC:
  Scalar:  a += (jω / Z_s) * u.Trace() * v.Trace() * ds(bnd)
  Tensor:  a += (jω / Z_xx) * u.Trace()[0] * v.Trace()[0] * ds(bnd)
         + (jω / Z_yy) * u.Trace()[1] * v.Trace()[1] * ds(bnd)
  (Valid for flat workpiece with normal in z-direction)

Theory:
  Kotiuga (1987), Bossavit (1998): thick cuts for multi-port T-Omega
  Hollaus (TU Wien): T-Omega dual for laminates
  Ren et al. (CEFC 2026): systematic multi-port dual formulation

TODO (requires Ren 2026 paper):
  - Off-diagonal Z_xy, Z_yx for general anisotropy
  - Topological multi-port: ring/torus workpiece via cohomology cuts
  - Cauer ladder expansion of Z_s(omega) tensor (CLN dual)
"""

import numpy as np
from .esim_cell_problem import ESIMCellProblemSolver, ESIMFiniteSlabSolver


class ESIMPort:
    """
    One port of the multi-port ESIM system.

    A port corresponds to one tangential direction on the workpiece surface,
    characterised by its material parameters (sigma, mu) and physical direction.

    For anisotropic materials (e.g. grain-oriented electrical steel),
    each principal axis has its own port with different mu.
    """

    def __init__(self, name, sigma, direction=(1.0, 0.0, 0.0),
                 bh_curve=None, mu_r=None, complex_mu=None):
        """
        Parameters
        ----------
        name : str
            Port label, e.g. 'rolling', 'transverse', 'x', 'y'.
        sigma : float
            Electrical conductivity [S/m].
        direction : tuple of 3 floats
            Tangential direction unit vector in 3D.
        bh_curve : list of [H, B] pairs, optional
            B-H curve for nonlinear real permeability. H in A/m, B in Tesla.
        mu_r : float, optional
            Constant relative permeability (used if bh_curve is None).
        complex_mu : tuple or list, optional
            Complex permeability specification (passed to ESIMCellProblemSolver).
        """
        if bh_curve is None and mu_r is None and complex_mu is None:
            raise ValueError("Specify at least one of: bh_curve, mu_r, complex_mu")
        self.name = name
        self.sigma = float(sigma)
        self.direction = np.asarray(direction, dtype=float)
        self.bh_curve = bh_curve
        self.mu_r = mu_r
        self.complex_mu = complex_mu

    def _bh_curve_for_solver(self):
        """Return bh_curve list compatible with ESIMCellProblemSolver."""
        if self.bh_curve is not None:
            return self.bh_curve
        if self.mu_r is not None:
            # Convert constant mu_r to a two-point B-H curve
            mu0 = 4e-7 * np.pi
            H_pts = [1.0, 1e6]
            B_pts = [mu0 * self.mu_r * H for H in H_pts]
            return list(zip(H_pts, B_pts))
        return None  # complex_mu case: no bh_curve


class ESIMMultiportSolver:
    """
    Two-port surface impedance solver for anisotropic conductors.

    Computes the 2x2 complex tensor Z_s and its inverse Y_s = Z_s^{-1}.

    For isotropic materials:  Z_s = Z_scalar * I  (diagonal, Z_xx = Z_yy)
    For orthotropic (principal axes aligned): Z_xy = Z_yx = 0, Z_xx != Z_yy
    For general anisotropy: off-diagonal terms non-zero (TODO: requires Ren 2026)

    Usage
    -----
    # Grain-oriented electrical steel (IEC 60404-2 grade M100-23P)
    solver = ESIMMultiportSolver.rolling_steel(
        freq=10e3, sigma=2e6,
        bh_rolling=BH_ROLLING, bh_transverse=BH_TRANSVERSE
    )
    Z_s, Y_s = solver.solve([H0_x, H0_y])  # H0 in A/m

    # Isotropic workpiece (soft iron)
    solver = ESIMMultiportSolver.isotropic(freq=10e3, sigma=1e7, bh_curve=BH)
    Z_s, Y_s = solver.solve([H0, H0])
    """

    def __init__(self, freq, ports):
        """
        Parameters
        ----------
        freq : float
            Operating frequency [Hz].
        ports : list of ESIMPort, length 2
            [port_1 (x-direction), port_2 (y-direction)].
        """
        if len(ports) != 2:
            raise ValueError(
                f"ESIMMultiportSolver requires exactly 2 ports, got {len(ports)}"
            )
        self.freq = float(freq)
        self.ports = ports

    # ------------------------------------------------------------------
    # Main interface
    # ------------------------------------------------------------------

    def solve(self, H0_list, **solver_kwargs):
        """
        Compute Z_s tensor for given surface field operating points.

        Parameters
        ----------
        H0_list : sequence of 2 floats
            [H0_x, H0_y] — tangential field amplitudes at each port [A/m].
        **solver_kwargs
            Passed to ESIMCellProblemSolver.solve() (tol, max_iter, relaxation).

        Returns
        -------
        Z_s : np.ndarray, shape (2,2), dtype complex
            Surface impedance tensor [Ohm].
            [E_tx; E_ty] = Z_s @ [H_tx; H_ty]
        Y_s : np.ndarray, shape (2,2), dtype complex
            Surface admittance tensor [S] = Z_s^{-1}.
            [H_tx; H_ty] = Y_s @ [E_tx; E_ty]

        Notes
        -----
        Off-diagonal Z_xy, Z_yx are set to zero (orthotropic assumption).
        For general anisotropy with cross-coupling permeability, the full
        2D cell problem is required (Ren CEFC 2026 -- not yet implemented).
        """
        if len(H0_list) != 2:
            raise ValueError("H0_list must have length 2 (one per port)")

        Z_diag = np.zeros(2, dtype=complex)
        results = []

        solvers = getattr(self, '_solvers', None)
        for j, (port, H0) in enumerate(zip(self.ports, H0_list)):
            cell_solver = solvers[j] if solvers is not None else self._make_cell_solver(port)
            res = cell_solver.solve(float(H0), **solver_kwargs)
            Z_diag[j] = res['Z']
            results.append(res)

        Z_s = np.diag(Z_diag)

        # Off-diagonal: zero for orthotropic materials.
        # TODO: Ren (2026) multi-port dual formulation for Z_xy, Z_yx.

        det = Z_s[0, 0] * Z_s[1, 1]  # diagonal case
        if abs(det) < 1e-300:
            raise RuntimeError("Z_s is singular — check material parameters")
        Y_s = np.linalg.inv(Z_s)

        return Z_s, Y_s

    def solve_primal(self, H0_list, **kw):
        """Return Z_s (primal / impedance formulation)."""
        Z_s, _ = self.solve(H0_list, **kw)
        return Z_s

    def solve_dual(self, H0_list, **kw):
        """
        Return Y_s = Z_s^{-1} (dual / admittance formulation).

        Preferred for high-frequency or low-conductivity cases where
        |Z_s| >> 1 causes poor conditioning in the Robin BC coefficient
        (jomega / Z_s -> 0, bilinear form loses the surface term).
        """
        _, Y_s = self.solve(H0_list, **kw)
        return Y_s

    def robin_coeffs(self, H0_list, omega, formulation='primal', **kw):
        """
        Return NGSolve Robin BC coefficients for the tensor SIBC.

        For a flat workpiece (normal in z-direction):
          a += coeff_xx * u.Trace()[0] * v.Trace()[0] * ds(bnd)
          a += coeff_yy * u.Trace()[1] * v.Trace()[1] * ds(bnd)

        Parameters
        ----------
        H0_list : sequence of 2 floats
            Operating-point field amplitudes [A/m].
        omega : float
            Angular frequency 2*pi*f [rad/s].
        formulation : 'primal' or 'dual'
            'primal': coeff = j*omega / Z_s  (standard NGSolve convention)
            'dual':   coeff = j*omega * Y_s  (equivalent, better for large |Z_s|)

        Returns
        -------
        coeff_xx, coeff_yy : complex
            Robin BC coefficients for x and y tangential components.
        """
        Z_s, Y_s = self.solve(H0_list, **kw)
        if formulation == 'primal':
            coeff_xx = 1j * omega / Z_s[0, 0]
            coeff_yy = 1j * omega / Z_s[1, 1]
        elif formulation == 'dual':
            coeff_xx = 1j * omega * Y_s[0, 0]
            coeff_yy = 1j * omega * Y_s[1, 1]
        else:
            raise ValueError(f"formulation must be 'primal' or 'dual', got {formulation!r}")
        return coeff_xx, coeff_yy

    # ------------------------------------------------------------------
    # Convenience constructors
    # ------------------------------------------------------------------

    @classmethod
    def isotropic(cls, freq, sigma, bh_curve=None, mu_r=None, complex_mu=None):
        """
        Isotropic material: Z_xx = Z_yy, Z_xy = Z_yx = 0.

        Parameters
        ----------
        freq : float
        sigma : float  [S/m]
        bh_curve / mu_r / complex_mu : material permeability (one required)
        """
        port_x = ESIMPort('x', sigma, (1, 0, 0),
                          bh_curve=bh_curve, mu_r=mu_r, complex_mu=complex_mu)
        port_y = ESIMPort('y', sigma, (0, 1, 0),
                          bh_curve=bh_curve, mu_r=mu_r, complex_mu=complex_mu)
        return cls(freq, [port_x, port_y])

    @classmethod
    def rolling_steel(cls, freq, sigma,
                      bh_rolling, bh_transverse,
                      mu_r_rolling=None, mu_r_transverse=None):
        """
        Grain-oriented electrical steel (anisotropic).

        Rolling direction (x, easy axis): higher mu, lower Z_s.
        Transverse direction (y, hard axis): lower mu, higher Z_s.

        Parameters
        ----------
        freq : float
        sigma : float  [S/m]
        bh_rolling : list of [H, B]   -- B-H curve for rolling direction
        bh_transverse : list of [H, B] -- B-H curve for transverse direction
        mu_r_rolling, mu_r_transverse : float, optional
            Constant mu_r fallback if bh curves are None.
        """
        port_rolling = ESIMPort(
            'rolling', sigma, (1, 0, 0),
            bh_curve=bh_rolling, mu_r=mu_r_rolling
        )
        port_transverse = ESIMPort(
            'transverse', sigma, (0, 1, 0),
            bh_curve=bh_transverse, mu_r=mu_r_transverse
        )
        return cls(freq, [port_rolling, port_transverse])

    @classmethod
    def from_z_scalars(cls, freq, Z_xx, Z_yy):
        """
        Build a solver from pre-computed scalar impedances.

        Useful for testing or when Z_s values come from measurement.
        """

        class _ConstantSolver:
            def __init__(self, Z_val):
                self._Z = Z_val
            def solve(self, H0, **kw):
                return {'Z': self._Z, 'converged': True, 'iterations': 0}

        obj = cls.__new__(cls)
        obj.freq = float(freq)
        obj._solvers = [_ConstantSolver(Z_xx), _ConstantSolver(Z_yy)]
        obj.ports = [
            ESIMPort('x', 1.0, (1, 0, 0), mu_r=1.0),
            ESIMPort('y', 1.0, (0, 1, 0), mu_r=1.0),
        ]
        return obj

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _make_cell_solver(self, port):
        """Create ESIMCellProblemSolver for one port."""
        return ESIMCellProblemSolver(
            bh_curve=port._bh_curve_for_solver(),
            sigma=port.sigma,
            frequency=self.freq,
            complex_mu=port.complex_mu,
        )


# ------------------------------------------------------------------
# NGSolve integration helper
# ------------------------------------------------------------------

def make_tensor_robin_cf(Z_s, omega, mesh, bnd_name,
                         normal_direction='z', formulation='primal'):
    """
    Build NGSolve CoefficientFunctions for a tensor Robin BC on a flat surface.

    For a flat workpiece with outward normal in 'normal_direction':
      a += coeff_xx * u.Trace()[idx_x] * v.Trace()[idx_x] * ds(bnd_name)
      a += coeff_yy * u.Trace()[idx_y] * v.Trace()[idx_y] * ds(bnd_name)

    Parameters
    ----------
    Z_s : np.ndarray, shape (2,2), complex
        Surface impedance tensor from ESIMMultiportSolver.solve().
    omega : float
        Angular frequency [rad/s].
    mesh : ngsolve.Mesh
        The NGSolve mesh (used only to verify bnd_name exists).
    bnd_name : str
        Boundary label for the SIBC surface.
    normal_direction : 'x', 'y', or 'z'
        Outward normal direction of the flat workpiece surface.
    formulation : 'primal' or 'dual'

    Returns
    -------
    coeff_xx, coeff_yy : complex
        Robin BC coefficients for the two tangential components.
    idx_x, idx_y : int
        Indices into u.Trace() for the two tangential components.
    """
    _normal_to_tangential_indices = {'x': (1, 2), 'y': (0, 2), 'z': (0, 1)}
    if normal_direction not in _normal_to_tangential_indices:
        raise ValueError(f"normal_direction must be 'x', 'y', or 'z'")

    if bnd_name not in mesh.GetBoundaries():
        raise ValueError(
            f"Boundary '{bnd_name}' not found. Available: {list(mesh.GetBoundaries())}"
        )

    idx_x, idx_y = _normal_to_tangential_indices[normal_direction]

    if formulation == 'primal':
        coeff_xx = 1j * omega / Z_s[0, 0]
        coeff_yy = 1j * omega / Z_s[1, 1]
    elif formulation == 'dual':
        Y_s = np.linalg.inv(Z_s)
        coeff_xx = 1j * omega * Y_s[0, 0]
        coeff_yy = 1j * omega * Y_s[1, 1]
    else:
        raise ValueError(f"formulation must be 'primal' or 'dual'")

    return coeff_xx, coeff_yy, idx_x, idx_y
