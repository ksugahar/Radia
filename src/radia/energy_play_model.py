"""
Energy-Based B-input Play Hysteresis Model

Decomposes the standard B-input Play model into reversible and irreversible
components, enabling:
- Convex energy functional (thermodynamically consistent)
- Fast inverse (Picard 2-3 iterations vs Newton ~100)
- Direct integration with Hantila polarization method
- Compatible with 2-scalar (Omega-reduced Omega) FEM formulation

Reference: docs/research/ENERGY_BASED_B_INPUT_DERIVATION.md

Usage:
    from radia.hysteresis_io import load_hys, build_shape_functions
    from radia.energy_play_model import EnergyBasedPlayModel

    loops = load_hys("material.hys")
    eta, f_k_tables, _ = build_shape_functions(loops)
    model = EnergyBasedPlayModel(eta, f_k_tables)

    # Forward: B -> H (same speed as standard Play)
    H = model.forward(B)

    # Inverse: H -> B (Picard, 2-3 iter instead of Newton ~100)
    B = model.inverse(H)

    # For Hantila solver:
    nu_rev = model.nu_rev          # constant reluctivity (LU once)
    H_irr = model.irreversible(B)  # nonlinear residual (iteration)
"""

import numpy as np
from scipy.interpolate import interp1d


class EnergyBasedPlayModel:
    """Energy-based B-input Play model with reversible/irreversible separation.

    Decomposes H(B) = nu_rev * B + H_irr(B, history)
    where nu_rev is constant and H_irr uses Play operators with g_k >= 0.

    Args:
        eta: ndarray shape (K,), play thresholds in Tesla
        f_k_tables: list of (r_array, f_array) tuples, original shape functions
        nu_rev: float, optional. Reversible reluctivity. If None, computed
                automatically as min slope of f_0.
    """

    def __init__(self, eta, f_k_tables, nu_rev=None):
        self.K = len(f_k_tables)
        self.eta = np.array(eta, dtype=float)
        self.f_k_tables = f_k_tables

        # Compute nu_rev from f_0 if not provided
        if nu_rev is None:
            self.nu_rev = self._compute_nu_rev()
        else:
            self.nu_rev = float(nu_rev)

        # Build irreversible shape functions g_k
        self.g_k_tables = self._build_irreversible_tables()
        self.g_k_interp = [interp1d(r, g, kind='linear', fill_value='extrapolate')
                           for r, g in self.g_k_tables]

        # Original f_k interpolators (for standard forward)
        self.f_k_interp = [interp1d(r, f, kind='linear', fill_value='extrapolate')
                           for r, f in self.f_k_tables]

        # Play operator states (per-element, scalar for now)
        self._p = np.zeros(self.K)

    def _compute_nu_rev(self):
        """Compute reversible reluctivity for convex decomposition.

        nu_rev must be large enough that ALL g_k have non-negative slopes.
        This means nu_rev >= max total slope of H(B) at any B.

        For Picard convergence: contraction ratio = (nu_total - nu_rev)/nu_rev.
        We choose nu_rev = max slope of the total H(B) curve, which gives
        contraction ratio = 0 at the steepest point.
        """
        # Compute total slope dH/dB = sum of all f_k slopes
        # nu_rev should equal the maximum of sum_k f_k'
        total_max_slope = 0.0
        for k in range(self.K):
            rk, fk = self.f_k_tables[k]
            if len(rk) >= 2:
                slopes = np.diff(fk) / np.diff(rk)
                total_max_slope += float(np.max(np.abs(slopes)))

        return total_max_slope

    def _build_irreversible_tables(self):
        """Build g_k tables: g_0 = f_0 - nu_rev*r, g_k = f_k for k>=1."""
        g_tables = []
        for k, (r, f) in enumerate(self.f_k_tables):
            if k == 0:
                g = f - self.nu_rev * r
            else:
                g = f.copy()
            g_tables.append((r.copy(), g))
        return g_tables

    def reset_state(self):
        """Reset play operator states to zero."""
        self._p = np.zeros(self.K)

    def _play_operator(self, B, k):
        """Evaluate play operator p_k(B) with threshold eta_k."""
        if k == 0:
            return B  # eta_0 = 0
        eta_k = self.eta[k]
        self._p[k] = np.clip(B, self._p[k] - eta_k, self._p[k] + eta_k)
        return self._p[k]

    def forward(self, B):
        """Evaluate H(B) = nu_rev * B + H_irr(B).

        O(K) direct evaluation (same speed as standard Play).

        Args:
            B: float, magnetic flux density

        Returns:
            H: float, magnetic field intensity
        """
        return self.nu_rev * B + self.irreversible(B)

    def irreversible(self, B):
        """Evaluate irreversible component H_irr(B).

        Args:
            B: float, magnetic flux density

        Returns:
            H_irr: float, irreversible field
        """
        H_irr = 0.0
        for k in range(self.K):
            pk = self._play_operator(B, k)
            H_irr += float(self.g_k_interp[k](pk))
        return H_irr

    def inverse(self, H, B_init=None, max_iter=50, tol=1e-10, method='newton'):
        """Solve H -> B using Newton or Picard iteration.

        Newton: B^{n+1} = B^n - F(B^n) / F'(B^n)
          where F(B) = nu_rev*B + H_irr(B) - H
                F'(B) = nu_rev + dH_irr/dB

        Picard: B^{n+1} = (H - H_irr(B^n)) / nu_rev

        Args:
            H: float, target magnetic field intensity
            B_init: float, initial guess. If None, uses H/nu_rev.
            max_iter: maximum iterations
            tol: convergence tolerance
            method: 'newton' (default, quadratic convergence) or 'picard'

        Returns:
            B: float, magnetic flux density
        """
        if self.nu_rev < 1e-30:
            raise ValueError("nu_rev too small for inverse")

        p_save = self._p.copy()

        if B_init is None:
            B = H / self.nu_rev if self.nu_rev > 0 else 0.0
        else:
            B = B_init

        for it in range(max_iter):
            self._p = p_save.copy()

            H_irr = self.irreversible(B)
            F = self.nu_rev * B + H_irr - H

            if abs(F) < tol:
                break

            if method == 'newton':
                # Analytical Jacobian: dH_irr/dB = sum of active g_k slopes
                dH_irr_dB = self._jacobian(B)
                dF_dB = self.nu_rev + dH_irr_dB
                if abs(dF_dB) < 1e-30:
                    # Fallback to Picard
                    B = (H - H_irr) / self.nu_rev
                else:
                    B = B - F / dF_dB
            else:
                # Picard
                B = (H - H_irr) / self.nu_rev

        # Final state update with converged B
        self._p = p_save.copy()
        self.irreversible(B)
        return B

    def _jacobian(self, B):
        """Compute dH_irr/dB (analytical Jacobian for Newton).

        dH_irr/dB = sum_k g_k'(p_k) * dp_k/dB

        dp_k/dB = 1 if play operator is active (not at limit)
        dp_k/dB = 0 if play operator is stuck at limit
        """
        dH_dB = 0.0
        for k in range(self.K):
            pk = self._p[k]
            eta_k = self.eta[k] if k > 0 else 0.0

            # dp_k/dB: active if B is within [p_k - eta_k, p_k + eta_k]
            if eta_k == 0 or (pk - eta_k < B < pk + eta_k):
                dp_dB = 1.0
            else:
                dp_dB = 0.0

            if dp_dB > 0:
                # g_k'(p_k) via finite difference on interpolator
                rk, gk = self.g_k_tables[k]
                eps = max(1e-10, abs(pk) * 1e-8)
                g_plus = float(self.g_k_interp[k](pk + eps))
                g_minus = float(self.g_k_interp[k](pk - eps))
                dg_dp = (g_plus - g_minus) / (2 * eps)
                dH_dB += dg_dp * dp_dB

        return dH_dB

    def commit_state(self):
        """Save current play operator states (call after converged step)."""
        self._committed_p = self._p.copy()

    def restore_state(self):
        """Restore play operator states to last committed."""
        if hasattr(self, '_committed_p'):
            self._p = self._committed_p.copy()

    @property
    def contraction_ratio(self):
        """Estimated contraction ratio for Picard iteration.

        rho = max|dH_irr/dB| / nu_rev
        Should be < 1 for convergence, < 0.5 for fast convergence.
        """
        max_slope = 0.0
        for k in range(self.K):
            rk, gk = self.g_k_tables[k]
            if len(rk) >= 2:
                slopes = np.diff(gk) / np.diff(rk)
                max_slope += float(np.max(np.abs(slopes)))
        return max_slope / self.nu_rev if self.nu_rev > 0 else float('inf')

    def energy(self, B):
        """Compute total free energy W(B) = W_rev + W_irr.

        Args:
            B: float

        Returns:
            W: float, energy density (J/m^3)
        """
        W_rev = 0.5 * self.nu_rev * B**2
        W_irr = 0.0
        for k in range(self.K):
            pk = self._play_operator(B, k)
            rk, gk = self.g_k_tables[k]
            # Integrate g_k from 0 to pk
            mask = rk <= pk
            if np.any(mask):
                _trapz = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
                W_irr += float(_trapz(gk[mask], rk[mask]))
        return W_rev + W_irr


class VectorEnergyPlayModel:
    """Vectorial energy-based Play model for 3D problems.

    Uses component-wise decomposition:
    H_i = nu_rev * B_i + H_irr_i(B_i)  for i = x, y, z

    Each component has independent play operator states.
    """

    def __init__(self, eta, f_k_tables, nu_rev=None):
        self.models = [
            EnergyBasedPlayModel(eta, f_k_tables, nu_rev)
            for _ in range(3)  # x, y, z
        ]
        self.nu_rev = self.models[0].nu_rev

    def forward(self, Bx, By, Bz):
        """Evaluate H(B) component-wise."""
        Hx = self.models[0].forward(Bx)
        Hy = self.models[1].forward(By)
        Hz = self.models[2].forward(Bz)
        return Hx, Hy, Hz

    def irreversible(self, Bx, By, Bz):
        """Evaluate H_irr(B) component-wise."""
        return (self.models[0].irreversible(Bx),
                self.models[1].irreversible(By),
                self.models[2].irreversible(Bz))

    def inverse(self, Hx, Hy, Hz):
        """Solve H -> B component-wise using Picard."""
        Bx = self.models[0].inverse(Hx)
        By = self.models[1].inverse(Hy)
        Bz = self.models[2].inverse(Hz)
        return Bx, By, Bz

    def reset_state(self):
        for m in self.models:
            m.reset_state()

    def commit_state(self):
        for m in self.models:
            m.commit_state()

    def restore_state(self):
        for m in self.models:
            m.restore_state()
