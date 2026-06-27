"""Bergqvist Quadruplet Model (BQM) - Table-based implementation.

H-input vector hysteresis model with table-based anhysteretic function.
Reference: Prigozhin/Bergqvist energy-based variational model.
           Alexander Sauseng (TU Graz), hystereseoperator.m

Step 1: H-input with table-based M_an (replacing analytical atan)
Step 2: B-input version (inverting H->B relationship)
Step 3: Play model shape function -> energy model table conversion
         (Potter-Schmulian data -> B-input Play -> Egger energy tables)
"""

import numpy as np
from scipy.interpolate import interp1d

MU_0 = 4.0 * np.pi * 1e-7


class BQMModel:
    """Bergqvist Quadruplet Model with table-based anhysteretic function.

    Parameters
    ----------
    M_an : callable
        Anhysteretic magnetization function: |u| [A/m] -> |m| [A/m]
    N : int
        Number of pseudo-particles.
    k_l : array_like, shape (N,)
        Dissipation constants [A/m] for each particle.
    omega : array_like, shape (N,)
        Volumetric weights for each particle (sum should be 1).
    ndim : int
        Spatial dimension (2 or 3).
    """

    def __init__(self, M_an, N, k_l, omega, ndim=2):
        self.M_an = M_an
        self.N = N
        self.k_l = np.asarray(k_l, dtype=float)
        self.omega = np.asarray(omega, dtype=float)
        self.ndim = ndim

        # State variables (per particle)
        self.hrev = [np.zeros(ndim) for _ in range(N)]
        self.m = [np.zeros(ndim) for _ in range(N)]
        self.h_old = np.zeros(ndim)

    def reset(self):
        """Reset all internal state to zero."""
        for i in range(self.N):
            self.hrev[i] = np.zeros(self.ndim)
            self.m[i] = np.zeros(self.ndim)
        self.h_old = np.zeros(self.ndim)

    def _dMan(self, u_abs, du=1e-4):
        """Numerical derivative of M_an at u_abs."""
        if u_abs < du:
            return (self.M_an(du) - self.M_an(0.0)) / du
        return (self.M_an(u_abs + du) - self.M_an(u_abs - du)) / (2 * du)

    def evaluate(self, h_act, flag=1):
        """Evaluate BQM operator for given H field.

        Parameters
        ----------
        h_act : array_like, shape (ndim,)
            Applied H field [A/m].
        flag : int
            1 = energy minimization (Newton), 0 = no energy min.

        Returns
        -------
        m_tot : ndarray, shape (ndim,)
            Total magnetization [A/m].
        """
        h_act = np.asarray(h_act, dtype=float)
        m_tot = np.zeros(self.ndim)

        for j in range(self.N):
            k = self.k_l[j]
            hrev_old = self.hrev[j]
            m_old = self.m[j]

            diff = h_act - hrev_old

            if k <= 0:
                # k=0: purely anhysteretic, hrev always tracks h_act
                hrev_new = h_act.copy()
            else:
                normalized_sq = (diff[0] / k) ** 2 + (diff[1] / k) ** 2
                if self.ndim == 3:
                    normalized_sq += (diff[2] / k) ** 2

                if normalized_sq <= 1.0:
                    # Reversible: hrev unchanged
                    hrev_new = hrev_old.copy()
                else:
                    # Irreversible: solve for direction
                    hrev_new = self._solve_irreversible(
                        h_act, k, hrev_old, m_old, flag
                    )

            # Compute magnetization from hrev
            hrev_norm = np.linalg.norm(hrev_new)
            if hrev_norm > 1e-30:
                m_j = self.M_an(hrev_norm) * hrev_new / hrev_norm
            else:
                m_j = np.zeros(self.ndim)

            m_tot += self.omega[j] * m_j

            # Update state
            self.hrev[j] = hrev_new
            self.m[j] = m_j

        self.h_old = h_act.copy()
        return m_tot

    def _solve_irreversible(self, h_act, k, hrev_old, m_old, flag):
        """Solve for hrev when dissipation is reached (2D Newton)."""
        diff = h_act - hrev_old
        diff_norm = np.linalg.norm(diff)

        if k > 0 and diff_norm > 0:
            i_phi_init = (hrev_old - h_act) / (k * diff_norm / k)
        else:
            i_phi_init = np.zeros(self.ndim)

        if self.ndim == 2:
            phi = np.arctan2(i_phi_init[1], i_phi_init[0])
        else:
            # 3D: use direction vector directly
            phi = np.arctan2(i_phi_init[1], i_phi_init[0])

        if flag == 1 and k > 0:
            for _ in range(10):
                cos_phi = np.cos(phi)
                sin_phi = np.sin(phi)

                u = h_act + k * np.array([cos_phi, sin_phi])
                if self.ndim == 3:
                    u = np.append(u, h_act[2])  # Simplified 3D

                u_abs = np.linalg.norm(u)
                u_strich = k * np.array([-sin_phi, cos_phi])
                u_strich_strich = k * np.array([-cos_phi, -sin_phi])

                Man = self.M_an(u_abs)
                Man_strich = self._dMan(u_abs)

                # g' = (M_an/|u| * u - m_old)^T * u'
                m_dir = Man / u_abs * u[:2] - m_old[:2]
                g_strich = np.dot(m_dir, u_strich)

                # g'' (second derivative of objective)
                g_strich_strich = (
                    (u_abs * Man_strich - Man) / u_abs ** 3
                    * np.dot(u[:2], u_strich) ** 2
                    + Man / u_abs * np.dot(u_strich, u_strich)
                    + np.dot(m_dir, u_strich_strich)
                )

                if abs(g_strich_strich) < 1e-8:
                    break
                phi_new = phi - g_strich / g_strich_strich
                if abs(phi_new - phi) < 1e-4:
                    break
                phi = phi_new

        i_phi = np.array([np.cos(phi), np.sin(phi)])
        if self.ndim == 3:
            i_phi = np.append(i_phi, 0.0)

        hirr = k * (-i_phi)
        hrev = h_act - hirr
        return hrev

    def _save_state(self):
        """Save internal state for later restoration."""
        return ([h.copy() for h in self.hrev],
                [m.copy() for m in self.m],
                self.h_old.copy())

    def _restore_state(self, state):
        """Restore internal state from saved snapshot."""
        hrev, m, h_old = state
        self.hrev = [h.copy() for h in hrev]
        self.m = [m.copy() for m in m]
        self.h_old = h_old.copy()

    def run_trajectory(self, H_trajectory, flag=1):
        """H-input: run BQM over a time series of H values.

        Parameters
        ----------
        H_trajectory : ndarray, shape (n_steps, ndim)
            H field values at each time step [A/m].
        flag : int
            1 = energy minimization.

        Returns
        -------
        M_trajectory : ndarray, shape (n_steps, ndim)
            Magnetization at each time step [A/m].
        B_trajectory : ndarray, shape (n_steps, ndim)
            B field at each time step [T].
        """
        n_steps = len(H_trajectory)
        M = np.zeros((n_steps, self.ndim))
        B = np.zeros((n_steps, self.ndim))

        for i in range(n_steps):
            M[i] = self.evaluate(H_trajectory[i], flag=flag)
            B[i] = MU_0 * (H_trajectory[i] + M[i])

        return M, B

    def run_trajectory_binput(self, B_trajectory, flag=1):
        """B-input: same algorithm, B replaces H as input.

        The SAME BQM algorithm (dissipation check, Newton angle,
        anhysteretic function) operates identically. Only the
        interpretation changes:

        H-input BQM:  input=H, output=M, B=mu_0*(H+M)
        B-input BQM:  input=B, output=H, M=B/mu_0-H

        For B-input, M_an should map |brev| (Tesla) -> H contribution (A/m),
        and k_l should be in Tesla (B-space dissipation thresholds).

        Parameters
        ----------
        B_trajectory : ndarray, shape (n_steps, ndim)
            B field values at each time step [T].
        flag : int
            1 = energy minimization.

        Returns
        -------
        H_trajectory : ndarray, shape (n_steps, ndim)
            H field at each time step [A/m].
        M_trajectory : ndarray, shape (n_steps, ndim)
            Magnetization at each time step [A/m].
        """
        n_steps = len(B_trajectory)
        H = np.zeros((n_steps, self.ndim))
        M = np.zeros((n_steps, self.ndim))

        for i in range(n_steps):
            # evaluate() is generic: takes input vector, returns output vector
            # In B-input mode: input=B, output is interpreted as H
            H[i] = self.evaluate(B_trajectory[i], flag=flag)
            M[i] = B_trajectory[i] / MU_0 - H[i]

        return H, M

    def run_trajectory_inverse(self, H_trajectory, flag=1, tol=1e-10):
        """Inverse of B-input BQM: given H trajectory, compute B.

        Uses Newton iteration (fsolve-style) to find B such that
        BQM_forward(B) = H_target at each time step.
        Mirrors prototype4.m approach.

        Parameters
        ----------
        H_trajectory : ndarray, shape (n_steps, ndim)
            Target H field at each time step [A/m].
        flag : int
            1 = energy minimization (Newton angle).
        tol : float
            Convergence tolerance for Newton residual.

        Returns
        -------
        B_trajectory : ndarray, shape (n_steps, ndim)
            B field at each time step [T].
        M_trajectory : ndarray, shape (n_steps, ndim)
            Magnetization at each time step [A/m].
        """
        n_steps = len(H_trajectory)
        B_out = np.zeros((n_steps, self.ndim))
        M_out = np.zeros((n_steps, self.ndim))

        # Initial guess for B
        b = np.zeros(self.ndim)

        for i in range(n_steps):
            H_target = H_trajectory[i]

            # Save state before Newton iterations
            state_saved = self._save_state()

            # Newton iteration to solve: BQM(B) - H_target = 0
            for n_iter in range(50):
                # Evaluate BQM at current B guess
                self._restore_state(state_saved)
                H_computed = self.evaluate(b, flag=flag)
                residual = H_computed - H_target

                res_norm = np.linalg.norm(residual)
                if res_norm < tol:
                    break

                # Finite difference Jacobian dH/dB (2x2 for 2D)
                eps_fd = max(1e-8, 1e-6 * np.linalg.norm(b))
                J = np.zeros((self.ndim, self.ndim))
                for d in range(self.ndim):
                    b_plus = b.copy()
                    b_plus[d] += eps_fd
                    self._restore_state(state_saved)
                    H_plus = self.evaluate(b_plus, flag=flag)
                    J[:, d] = (H_plus - H_computed) / eps_fd

                # Solve J * delta_b = -residual
                try:
                    delta_b = np.linalg.solve(J, -residual)
                except np.linalg.LinAlgError:
                    break

                # Armijo backtracking
                tau = 1.0
                for ls in range(20):
                    b_trial = b + tau * delta_b
                    self._restore_state(state_saved)
                    H_trial = self.evaluate(b_trial, flag=flag)
                    res_trial = np.linalg.norm(H_trial - H_target)
                    if res_trial < res_norm:
                        break
                    tau *= 0.5

                b = b + tau * delta_b

            # Final evaluation at converged B to update state
            self._restore_state(state_saved)
            H_final = self.evaluate(b, flag=flag)

            B_out[i] = b.copy()
            M_out[i] = b / MU_0 - H_final

        return B_out, M_out


# =====================================================================
# Egger Energy-Based Inverse: H -> B
# =====================================================================

class EggerInverse:
    """Egger energy-based Inverse operator: H -> B.

    Each operator k independently minimizes:
      phi_k(J_k) = U_k(|J_k|) - H * J_k + chi_k * |J_k - J_k_prev|

    Then B = mu_0 * H + sum_k J_k.

    This is O(K) per step (each J_k solved independently, no Schur).
    Follows rad_material_impl.cpp SolveInverseK().

    Parameters
    ----------
    K : int
        Number of operators.
    chi : array_like, shape (K,)
        Pinning strengths for each operator.
    f_k_tables : list of (r_array, f_array) tuples
        Shape function tables. f_k(r) = U_k'(r).
        r_array: |J| grid [0, r_max], f_array: f_k values.
    ndim : int
        Spatial dimension (2 or 3).
    """

    def __init__(self, K, chi, f_k_tables, ndim=2):
        self.K = K
        self.chi = np.asarray(chi, dtype=float)
        self.ndim = ndim
        self.eps = 1e-8  # smoothing for |x|

        # Precompute tables: f (shape), U (integral of f), df (derivative of f)
        self.tables = []
        for k in range(K):
            r, f = f_k_tables[k]
            r = np.asarray(r, dtype=float)
            f = np.asarray(f, dtype=float)
            n = len(r)
            r_max = r[-1]

            # U = integral_0^r f(s) ds (trapezoidal cumulative)
            dr = np.diff(r)
            f_avg = 0.5 * (f[:-1] + f[1:])
            U = np.zeros(n)
            U[1:] = np.cumsum(f_avg * dr)

            # df = f'(r) (central differences, forward/backward at ends)
            df = np.zeros(n)
            if n >= 2:
                df[0] = (f[1] - f[0]) / (r[1] - r[0]) if r[1] > r[0] else 0.0
                df[-1] = (f[-1] - f[-2]) / (r[-1] - r[-2]) if r[-1] > r[-2] else 0.0
            for i in range(1, n - 1):
                df[i] = (f[i+1] - f[i-1]) / (r[i+1] - r[i-1])

            # Store raw arrays for np.interp (faster than scipy.interp1d)
            self.tables.append({
                'r': r, 'f': f, 'U': U, 'df': df,
                'n': n, 'r_max': r_max,
            })

        # State: J_k_prev per operator
        self.J_prev = [np.zeros(ndim) for _ in range(K)]

    def reset(self):
        """Reset all internal state."""
        for k in range(self.K):
            self.J_prev[k] = np.zeros(self.ndim)

    def _Uk(self, k, r):
        t = self.tables[k]
        return float(np.interp(r, t['r'], t['U']))

    def _fk(self, k, r):
        t = self.tables[k]
        return float(np.interp(r, t['r'], t['f']))

    def _dfk(self, k, r):
        t = self.tables[k]
        return float(np.interp(r, t['r'], t['df']))

    def _norm_eps(self, x):
        """Smoothed norm |x|_eps = sqrt(x*x + eps^2)."""
        return np.sqrt(np.dot(x, x) + self.eps**2)

    def _grad_norm_eps(self, x):
        """Gradient of |x|_eps."""
        return x / self._norm_eps(x)

    def _hess_norm_eps(self, x):
        """Hessian of |x|_eps."""
        ne = self._norm_eps(x)
        n = len(x)
        I = np.eye(n)
        return (I - np.outer(x, x) / ne**2) / ne

    def _grad_Uk(self, k, J):
        """Gradient of U_k(|J|) w.r.t. J."""
        J_norm = np.linalg.norm(J)
        if J_norm < 1e-30:
            return np.zeros(self.ndim)
        return self._fk(k, J_norm) * J / J_norm

    def _hess_Uk(self, k, J):
        """Hessian of U_k(|J|) w.r.t. J."""
        J_norm = np.linalg.norm(J)
        n = self.ndim
        if J_norm < 1e-30:
            return self._dfk(k, 0.0) * np.eye(n)
        fk = self._fk(k, J_norm)
        dfk = self._dfk(k, J_norm)
        e = J / J_norm
        # d^2 U / dJ_i dJ_j = df/dr * e_i*e_j + f/r * (I - e_i*e_j)
        return dfk * np.outer(e, e) + (fk / J_norm) * (np.eye(n) - np.outer(e, e))

    def solve_k(self, k, H):
        """Solve for J_k given H (Newton with Armijo backtracking).

        Minimizes phi_k(J_k) = U_k(|J_k|) - H*J_k + chi_k*|J_k - J_k_prev|_eps

        Follows SolveInverseK() in rad_material_impl.cpp.
        """
        J_k = self.J_prev[k].copy()
        chi_k = self.chi[k]
        r_max = self.tables[k]['r_max']

        for it in range(30):
            diff = J_k - self.J_prev[k]
            # Gradient: grad U_k - H + chi_k * grad |diff|_eps
            grad = self._grad_Uk(k, J_k) - H + chi_k * self._grad_norm_eps(diff)
            grad_norm = np.linalg.norm(grad)
            if grad_norm < 1e-12:
                break

            # Hessian: Hess U_k + chi_k * Hess |diff|_eps
            hess = self._hess_Uk(k, J_k) + chi_k * self._hess_norm_eps(diff)
            try:
                dJ = np.linalg.solve(hess, -grad)
            except np.linalg.LinAlgError:
                break

            # Armijo backtracking
            tau = 1.0
            J_norm_k = np.linalg.norm(J_k)
            F_curr = (self._Uk(k, J_norm_k) - np.dot(H, J_k)
                      + chi_k * self._norm_eps(diff))
            dir_deriv = np.dot(grad, dJ)
            for ls in range(20):
                J_trial = J_k + tau * dJ
                diff_trial = J_trial - self.J_prev[k]
                F_new = (self._Uk(k, np.linalg.norm(J_trial))
                         - np.dot(H, J_trial)
                         + chi_k * self._norm_eps(diff_trial))
                if F_new <= F_curr + 0.1 * tau * dir_deriv:
                    break
                tau *= 0.5

            J_k = J_k + tau * dJ

            # Clamp to saturation
            Jmax = 0.9999 * r_max
            J_k = np.clip(J_k, -Jmax, Jmax)

        return J_k

    def evaluate(self, H):
        """Compute B from H using Egger Inverse.

        Each J_k is solved independently, then B = mu_0*H + sum J_k.

        Parameters
        ----------
        H : array_like, shape (ndim,)
            Applied H field [A/m].

        Returns
        -------
        B : ndarray, shape (ndim,)
            Magnetic flux density [T].
        """
        H = np.asarray(H, dtype=float)
        J_total = np.zeros(self.ndim)

        for k in range(self.K):
            J_k = self.solve_k(k, H)
            self.J_prev[k] = J_k.copy()
            J_total += J_k

        return MU_0 * H + J_total

    def run_trajectory(self, H_trajectory):
        """Run Egger Inverse over a trajectory.

        Parameters
        ----------
        H_trajectory : ndarray, shape (n_steps, ndim)
            H field values at each time step [A/m].

        Returns
        -------
        B_trajectory : ndarray, shape (n_steps, ndim)
            B field at each time step [T].
        M_trajectory : ndarray, shape (n_steps, ndim)
            Magnetization at each time step [A/m].
        """
        n_steps = len(H_trajectory)
        B = np.zeros((n_steps, self.ndim))
        M = np.zeros((n_steps, self.ndim))

        for i in range(n_steps):
            B[i] = self.evaluate(H_trajectory[i])
            M[i] = B[i] / MU_0 - H_trajectory[i]

        return B, M


def make_egger_inverse(M_an=None, K=20, k_max=140.0, ms=1.23e6, A=38.0,
                       H_max=10000.0, n_table=1001, ndim=2):
    """Create Egger energy-based Inverse model (H -> B).

    Derives energy model shape functions from H-input anhysteretic M_an.

    The key derivation:
      J(H) = mu_0 * M_an(H)   (polarization on anhysteretic curve)
      J_k(H) = J(H) / K       (equal partition among K operators)
      f_k(r) = J_k^{-1}(r)    (shape function = inverse of J_k)
      r_max = mu_0 * ms / K   (maximum |J_k| at saturation)
      chi_k = k_l[k]          (pinning strength in A/m)

    Parameters
    ----------
    M_an : callable, optional
        H-input anhysteretic: |H| -> |M| [A/m]. Default: atan(ms, A).
    K : int
        Number of operators.
    k_max : float
        Maximum pinning strength [A/m] (H-input dissipation).
    ms : float
        Saturation magnetization [A/m].
    A : float
        Curvature parameter [A/m].
    H_max : float
        Maximum H for table construction [A/m].
    n_table : int
        Number of table points.
    ndim : int
        Spatial dimension.
    """
    if M_an is None:
        M_an = make_atan_Man(ms, A)

    # r_max per operator = mu_0 * ms / K (saturation bound)
    r_max = MU_0 * ms / K

    # Build f_k on UNIFORM r grid (critical for interpolation accuracy)
    # f_k(r) = M_an^{-1}(K*r/mu_0): maps |J_k| -> H
    # We invert M_an numerically: build H(M) table, then evaluate at K*r/mu_0
    H_dense = np.linspace(0, H_max, 10 * n_table)
    M_dense = np.array([M_an(h) for h in H_dense])
    M_to_H = interp1d(M_dense, H_dense, kind='linear',
                       bounds_error=False,
                       fill_value=(H_dense[0], H_dense[-1]))

    r_grid = np.linspace(0, r_max * 0.9999, n_table)
    f_grid = np.array([float(M_to_H(K * r / MU_0)) for r in r_grid])

    # Build shape function tables (same for all K operators)
    f_k_tables = []
    for k in range(K):
        f_k_tables.append((r_grid.copy(), f_grid.copy()))

    # Pinning: k_l[k] in A/m (same as H-input BQM dissipation)
    chi = np.array([k_max * k / (K - 1) if K > 1 else 0.0
                    for k in range(K)])

    return EggerInverse(K, chi, f_k_tables, ndim=ndim)


# =====================================================================
# Anhysteretic function constructors
# =====================================================================

def make_atan_Man(ms, A):
    """Create analytical atan anhysteretic function (BQM standard).

    M_an(|u|) = (2*ms/pi) * atan(|u|/A)

    Parameters
    ----------
    ms : float
        Saturation magnetization [A/m].
    A : float
        Curvature parameter [A/m].
    """
    def M_an(u_abs):
        return (2.0 * ms / np.pi) * np.arctan(u_abs / A)
    return M_an


def make_table_Man(u_values, m_values):
    """Create table-based anhysteretic function via linear interpolation.

    Parameters
    ----------
    u_values : array_like
        |u| grid points [A/m], monotonically increasing.
    m_values : array_like
        |m| values at grid points [A/m].

    Returns
    -------
    M_an : callable
        |u| -> |m| function.
    """
    u_values = np.asarray(u_values, dtype=float)
    m_values = np.asarray(m_values, dtype=float)
    interp = interp1d(u_values, m_values, kind='linear',
                      bounds_error=False,
                      fill_value=(m_values[0], m_values[-1]))

    def M_an(u_abs):
        return float(interp(u_abs))
    return M_an


def atan_to_table(ms, A, n_points=1001, u_max=None):
    """Discretize atan anhysteretic function into a table.

    Parameters
    ----------
    ms : float
        Saturation magnetization [A/m].
    A : float
        Curvature parameter [A/m].
    n_points : int
        Number of table points.
    u_max : float, optional
        Maximum |u| value. Default: 100*A.

    Returns
    -------
    u_values : ndarray
        |u| grid.
    m_values : ndarray
        |m| grid.
    """
    if u_max is None:
        u_max = 100.0 * A
    u_values = np.linspace(0, u_max, n_points)
    m_values = (2.0 * ms / np.pi) * np.arctan(u_values / A)
    return u_values, m_values


# =====================================================================
# Default material (Sauseng parameters)
# =====================================================================

def invert_Man_to_Han(M_an, ms, H_max=10000.0, n_points=2001):
    """Compute B-input shape function H_an from H-input shape function M_an.

    The H-input anhysteretic curve: B(H) = mu_0*(H + M_an(|H|))
    The B-input shape function is H_an(|B|) = B^{-1}(|B|), the inverse.

    Parameters
    ----------
    M_an : callable
        H-input anhysteretic function: |H| (A/m) -> |M| (A/m).
    ms : float
        Saturation magnetization [A/m] (for setting B range).
    H_max : float
        Maximum H value for computing the curve [A/m].
    n_points : int
        Number of table points.

    Returns
    -------
    H_an : callable
        B-input shape function: |B| (T) -> |H| (A/m).
    B_sat : float
        Approximate B saturation value [T].
    """
    H_vals = np.linspace(0, H_max, n_points)
    B_vals = np.array([MU_0 * (h + M_an(h)) for h in H_vals])

    B_sat = B_vals[-1]

    # Invert: swap axes (B -> H)
    interp = interp1d(B_vals, H_vals, kind='linear',
                      bounds_error=False,
                      fill_value=(H_vals[0], H_vals[-1]))

    def H_an(b_abs):
        return float(interp(b_abs))
    return H_an, B_sat


def make_binput_material(H_an=None, N=20, k_max_tesla=None,
                         ms=1.23e6, A=38.0, k_max_hinput=140.0,
                         ndim=2):
    """Create B-input BQM model.

    Same BQM algorithm, but input is B (Tesla) and output is H (A/m).
    Shape function H_an maps |brev| (T) -> H contribution (A/m).
    Dissipation k_l in Tesla.

    Parameters
    ----------
    H_an : callable, optional
        B-input shape function. If None, derived from atan M_an.
    N : int
        Number of pseudo-particles.
    k_max_tesla : float, optional
        Maximum dissipation in Tesla. If None, derived from H-input k_max.
    ms : float
        Saturation magnetization [A/m] (for deriving H_an).
    A : float
        Curvature parameter [A/m] (for deriving H_an).
    k_max_hinput : float
        H-input max dissipation [A/m] (for scaling to Tesla).
    ndim : int
        Spatial dimension.
    """
    if H_an is None:
        M_an_hinput = make_atan_Man(ms, A)
        H_an, B_sat = invert_Man_to_Han(M_an_hinput, ms)

    if k_max_tesla is None:
        # Scale k_l from H-space (A/m) to B-space (T)
        # At low field: dB/dH ~ mu_0*(1 + chi_initial)
        # chi_initial = dM_an/dH at H=0 = (2*ms)/(pi*A)
        chi_0 = (2.0 * ms) / (np.pi * A)
        k_max_tesla = MU_0 * (1.0 + chi_0) * k_max_hinput

    k_l = np.array([k_max_tesla * i / (N - 1) for i in range(N)])
    omega = np.ones(N) / N

    return BQMModel(H_an, N, k_l, omega, ndim=ndim)


def make_default_material(M_an=None, N=20, k_max=140.0, ndim=2):
    """Create BQM model with default Sauseng parameters.

    Parameters
    ----------
    M_an : callable, optional
        Anhysteretic function. Default: atan with ms=1.23e6, A=38.
    N : int
        Number of pseudo-particles.
    k_max : float
        Maximum dissipation constant [A/m].
    ndim : int
        Spatial dimension.
    """
    if M_an is None:
        M_an = make_atan_Man(ms=1.23e6, A=38.0)

    k_l = np.array([k_max * i / (N - 1) for i in range(N)])
    omega = np.ones(N) / N

    return BQMModel(M_an, N, k_l, omega, ndim=ndim)


# =====================================================================
# Step 3: B-input Play Model -> Energy Model Conversion
# =====================================================================

def compute_shape_functions(HdataLUT, Bplay):
    """Compute play model shape functions from descending branch data.

    Python translation of ShapeFunction.m (JMAG identification method).
    Decomposes descending H(B) branches into play operator contributions
    using nested finite differences.

    Parameters
    ----------
    HdataLUT : ndarray, shape (nB_max, K)
        H values on descending branches. Column i corresponds to the i-th
        BMax level. Row j is the j-th B grid point (descending from BMax).
        Curves ordered from LARGEST BMax (col 0) to smallest (col K-1).
    Bplay : ndarray
        Play operator grid, descending from Bmax to -Bmax with step dB/2.

    Returns
    -------
    f_interps : list of callable
        Np interpolation functions. f_interps[k](B) returns H contribution.
    """
    row0, column0 = HdataLUT.shape
    row1 = column0        # number of BMax levels
    column1 = row0 - 1    # number of B intervals

    mu = np.zeros((row1, column1))

    # Step 1: First H-difference for each BMax level
    for ii in range(row1):
        mu[ii, 0] = HdataLUT[0, ii] - HdataLUT[1, ii]

    # Step 2: Double differences (nested finite differences)
    for ii in range(row1 - 1):
        n_jj = column1 - 2 * (ii + 1)  # MATLAB: column1 - 2*ii
        for jj in range(1, n_jj + 1):
            mu[ii, jj] = (HdataLUT[jj, ii] - HdataLUT[jj + 1, ii]
                          - (HdataLUT[jj - 1, ii + 1]
                             - HdataLUT[jj, ii + 1]))

    # Step 3: Adjust last non-zero element of each row
    mu_temp = mu.sum(axis=1)
    for ii in range(row1):
        col = column1 - 2 * ii - 1
        mu[ii, col] = (HdataLUT[col, ii] - HdataLUT[col + 1, ii]
                       - mu_temp[ii])

    # Step 4: Expand mu (row1 x column1) -> mu2 (2*row1 x column1)
    # Each mu row is duplicated into 2 rows with halved values.
    # Triangular fill: column jj starts at row ii0=jj.
    mu2 = np.zeros((row1 * 2, column1))
    ii0 = 0
    for jj in range(column1):
        for ii in range(row1 * 2 - ii0):
            mu2[ii + ii0, jj] = mu[ii // 2, jj] / 2.0
        ii0 += 1

    # Step 5: Reverse rows
    mu2 = mu2[::-1, :].copy()

    # Step 6: Cumulative sum -> descending branch
    Hdown = -np.cumsum(mu2, axis=0)
    Hup = -Hdown[::-1, :].copy()

    # Step 7: Assemble [Hup; zeros; Hdown]
    Hfunc = np.vstack([Hup, np.zeros((1, column1)), Hdown])

    # Step 8: Create interpolants (ascending x for interp1d)
    Bplay_asc = Bplay[::-1].copy()
    f_interps = []
    for jj in range(Hfunc.shape[1]):
        Hfunc_asc = Hfunc[::-1, jj].copy()
        f_interps.append(interp1d(Bplay_asc, Hfunc_asc, kind='linear',
                                  bounds_error=False,
                                  fill_value=(Hfunc_asc[0], Hfunc_asc[-1])))

    return f_interps


def play_hysteron(f_interps, Bx, By, eta, px0, py0):
    """Evaluate vector play hysterons for one time step.

    Python translation of PlayHysteron.m.

    Parameters
    ----------
    f_interps : list of callable
        Shape function interpolants from compute_shape_functions().
    Bx, By : float
        Current B field components [T].
    eta : ndarray, shape (Np,)
        Play thresholds for each operator [T].
    px0, py0 : ndarray, shape (Np,)
        Previous play output states. Modified in-place.

    Returns
    -------
    hx, hy : float
        H field components [A/m].
    """
    Np = len(eta)

    # Update play outputs
    B_p0 = np.sqrt((Bx - px0)**2 + (By - py0)**2)
    denom = np.maximum(eta, B_p0)
    # Avoid division by zero for eta=0 operators
    mask = denom > 0
    px0_new = px0.copy()
    py0_new = py0.copy()
    px0_new[mask] = Bx - eta[mask] * (Bx - px0[mask]) / denom[mask]
    py0_new[mask] = By - eta[mask] * (By - py0[mask]) / denom[mask]
    # For eta=0: p_k tracks B exactly
    px0_new[~mask] = Bx
    py0_new[~mask] = By

    # First operator (eta=0) always tracks B
    px0_new[0] = Bx
    py0_new[0] = By

    p0 = np.sqrt(px0_new**2 + py0_new**2)

    # Evaluate shape functions and compute H
    hx = 0.0
    hy = 0.0
    for jj in range(Np):
        if p0[jj] == 0:
            continue
        f0 = float(f_interps[jj](p0[jj]))
        hx += f0 * px0_new[jj] / p0[jj]
        hy += f0 * py0_new[jj] / p0[jj]

    # Update state in-place
    px0[:] = px0_new
    py0[:] = py0_new

    return hx, hy


def load_binput_data(mat_file):
    """Load B-input descending branch data from .mat file.

    Parameters
    ----------
    mat_file : str
        Path to B_input.mat file (from CASE_02.m).

    Returns
    -------
    dict with keys:
        'BH': list of (H_array, B_array) tuples for each BMax level
        'dB': float, B step size [T]
        'BMax': ndarray, BMax levels [T] (ascending)
        'HdataLUT': ndarray, shape (nB_max, K), H lookup table
        'Bplay': ndarray, play grid (descending from Bmax to -Bmax)
        'eta': ndarray, play thresholds [T]
    """
    import scipy.io as sio

    d = sio.loadmat(mat_file)
    BH_raw = d['BH']
    dB = float(d['dB'].item())
    BMax = d['BMax'].flatten()

    K = BH_raw.shape[1]  # number of BMax levels

    # Extract BH curves
    BH_curves = []
    for i in range(K):
        H = BH_raw[0, i]['H'].flatten()
        B = BH_raw[0, i]['B'].flatten()
        BH_curves.append((H, B))

    # Build HdataLUT: shape (max_nB, K)
    max_nB = len(BH_curves[0][0])  # largest curve has most points
    HdataLUT = np.zeros((max_nB, K))
    for i in range(K):
        H, _ = BH_curves[i]
        HdataLUT[:len(H), i] = H

    # Play grid and thresholds
    Bmax = BMax[-1]  # largest BMax value
    Bplay = np.arange(Bmax, -Bmax - dB / 4, -dB / 2)
    # Ensure symmetric and exact
    n_half = int(round(Bmax / (dB / 2)))
    Bplay = np.linspace(Bmax, -Bmax, 2 * n_half + 1)

    eta = np.arange(0, Bmax - dB / 4, dB / 2)  # [0, dB/2, dB, ..., Bmax-dB/2]

    return {
        'BH': BH_curves,
        'dB': dB,
        'BMax': BMax,
        'K': K,
        'HdataLUT': HdataLUT,
        'Bplay': Bplay,
        'eta': eta,
        'Bmax': Bmax,
    }


def run_play_model(f_interps, eta, Bx_traj, By_traj):
    """Run B-input play model over a trajectory.

    Parameters
    ----------
    f_interps : list of callable
        Shape function interpolants from compute_shape_functions().
    eta : ndarray, shape (Np,)
        Play thresholds [T].
    Bx_traj, By_traj : ndarray, shape (n_steps,)
        B field trajectory components [T].

    Returns
    -------
    Hx, Hy : ndarray, shape (n_steps,)
        H field trajectory [A/m].
    """
    Np = len(eta)
    n_steps = len(Bx_traj)

    px0 = np.zeros(Np)
    py0 = np.zeros(Np)

    Hx = np.zeros(n_steps)
    Hy = np.zeros(n_steps)

    for nt in range(n_steps):
        hx, hy = play_hysteron(f_interps, Bx_traj[nt], By_traj[nt],
                               eta, px0, py0)
        Hx[nt] = hx
        Hy[nt] = hy

    return Hx, Hy


def play_to_egger(f_interps, eta, Bmax, n_table=501, ndim=2):
    """Convert B-input play model to Egger energy model.

    The play model shape functions f_k are directly usable as energy model
    shape functions because:
    - Play: H = sum_k f_k(|p_k|) * p_k/|p_k|,  p_k = play(B, eta_k)
    - Energy: dU_k/d|J_k| = f_k(|J_k|),  J_k = p_k,  chi_k = eta_k

    Parameters
    ----------
    f_interps : list of callable
        Shape function interpolants from compute_shape_functions().
    eta : ndarray, shape (Np,)
        Play thresholds [T] (become pinning strengths chi_k).
    Bmax : float
        Maximum B [T].
    n_table : int
        Number of points per shape function table.
    ndim : int
        Spatial dimension.

    Returns
    -------
    EggerInverse instance.
    """
    Np = len(eta)
    chi = eta.copy()  # pinning = play threshold

    f_k_tables = []
    for k in range(Np):
        # Maximum |J_k| = Bmax - eta_k (play output range)
        r_max_k = Bmax - eta[k]
        if r_max_k <= 0:
            r_max_k = 1e-6  # degenerate operator

        # Sample shape function on uniform r grid
        r_grid = np.linspace(0, r_max_k, n_table)
        f_grid = np.array([float(f_interps[k](r)) for r in r_grid])

        f_k_tables.append((r_grid, f_grid))

    return EggerInverse(Np, chi, f_k_tables, ndim=ndim)
