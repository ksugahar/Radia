"""c4v_symmetry.py — D_4h symmetry decomposition for HDivDivFreeHex basis.

For the A1 cuboid (a=b≠c) the eddy-current problem has D_4h symmetry.  The
HDivDivFreeHex basis (5 sub-blocks A/B/C/D/E indexed by (block, i, j, k))
maps closed under each of the 16 group operations as a permutation + ±
phase.  This module:

  Phase 1: build M(g) ∈ ℝ^{n × n} for every g ∈ D_4h (sparse permutation
           with ±1 entries) on a given HDivDivFreeHexBasisMP.
  Phase 2: assemble the character-projected operators P_Γ = (d_Γ/|G|) Σ_g
           χ_Γ(g) M(g) for each irrep Γ; extract U_Γ via eigendecomposition.
  Phase 3: assemble the symmetry-adapted transformation C; block-
           diagonalize K, M, b.
  Phase 4: verification utilities (trace, dim sum, B-direction selection).

Basis parity reference
----------------------
  L[m] = L_{m+2}(ξ) (integrated Legendre)  parity in ξ:  (-1)^m
  L'[m] = dL_{m+2}/dξ                       parity in ξ:  (-1)^{m+1}
  P[n] = P_n(ξ) (standard Legendre)         parity in ξ:  (-1)^n
  (so P_{j+1}(ξ) has parity (-1)^{j+1})

Block components (vx, vy, vz) from extract_tau_A1_square_mp.HDivDivFreeHexBasisMP:
  A_{ijk}: vx = L[i](ξ) P[j+1](η) L'[k](ζ),    vy = 0,
           vz = -L'[i](ξ) P[j+1](η) L[k](ζ)
  B_{ijk}: vx = 0,
           vy = P[i+1](ξ) L[j](η) L'[k](ζ),
           vz = -P[i+1](ξ) L'[j](η) L[k](ζ)
  C_{0jk}: vx = 0,
           vy = L[j](η) L'[k](ζ),
           vz = -L'[j](η) L[k](ζ)
  D_{i0k}: vx = -L[i](ξ) L'[k](ζ),
           vy = 0,
           vz = L'[i](ξ) L[k](ζ)
  E_{ij0}: vx = L[i](ξ) L'[j](η),
           vy = -L'[i](ξ) L[j](η),
           vz = 0
"""
from __future__ import annotations

from enum import IntEnum
from typing import Callable, Dict, List, Tuple

import numpy as np


# Group element labels --------------------------------------------------------

D4H_ELEMENTS: Tuple[str, ...] = (
    # D_4 (no horizontal reflection):
    "E",        # identity
    "C4",       # rotation by 90°  about z
    "C2",       # rotation by 180° about z
    "C43",      # rotation by 270° about z (= C4 inverse)
    "sv_x",     # reflection ξ → -ξ  (vertical mirror through y-axis)
    "sv_y",     # reflection η → -η  (vertical mirror through x-axis)
    "sd_p",     # diagonal reflection ξ ↔ η  ("σ_d through x+y")
    "sd_m",     # anti-diagonal reflection ξ ↔ -η  ("σ_d through x-y")
    # × σ_h:
    "sh",
    "S4",       # C4 * σ_h
    "i",        # C2 * σ_h  (inversion through origin)
    "S43",      # C43 * σ_h
    "Cv_x",     # sv_x * σ_h  (C_2 axis along x in plane)
    "Cv_y",     # sv_y * σ_h  (C_2 axis along y in plane)
    "Cd_p",     # sd_p * σ_h
    "Cd_m",     # sd_m * σ_h
)
assert len(D4H_ELEMENTS) == 16


# D_4h irreducible representations (10 irreps, sum of dim² = 16) ---------------
# Reference: standard D_4h character table.  Conjugacy classes:
#   E, 2C4, C2, 2C2', 2C2'', i, 2S4, σ_h, 2σ_v, 2σ_d
#
# For our 16-element list above:
#   class "E":   {E}
#   class "2C4": {C4, C43}
#   class "C2":  {C2}
#   class "2C2'": {Cv_x, Cv_y}        (axes through edges in xy plane)
#   class "2C2''": {Cd_p, Cd_m}       (axes through corners in xy plane)
#   class "i":   {i}
#   class "2S4": {S4, S43}
#   class "σ_h": {sh}
#   class "2σ_v": {sv_x, sv_y}
#   class "2σ_d": {sd_p, sd_m}
D4H_CHARACTER_TABLE: Dict[str, Dict[str, int]] = {
    # irrep: { class: χ }
    "A1g": {"E": 1, "2C4": 1,  "C2": 1,  "2C2'": 1,  "2C2''": 1,
            "i":  1, "2S4":  1, "sh": 1, "2sv":  1, "2sd": 1},
    "A2g": {"E": 1, "2C4": 1,  "C2": 1,  "2C2'": -1, "2C2''": -1,
            "i":  1, "2S4":  1, "sh": 1, "2sv": -1, "2sd": -1},
    "B1g": {"E": 1, "2C4": -1, "C2": 1,  "2C2'": 1,  "2C2''": -1,
            "i":  1, "2S4": -1, "sh": 1, "2sv":  1, "2sd": -1},
    "B2g": {"E": 1, "2C4": -1, "C2": 1,  "2C2'": -1, "2C2''": 1,
            "i":  1, "2S4": -1, "sh": 1, "2sv": -1, "2sd": 1},
    "Eg":  {"E": 2, "2C4": 0,  "C2": -2, "2C2'": 0,  "2C2''": 0,
            "i":  2, "2S4":  0, "sh": -2,"2sv":  0, "2sd": 0},
    "A1u": {"E": 1, "2C4": 1,  "C2": 1,  "2C2'": 1,  "2C2''": 1,
            "i": -1, "2S4": -1, "sh": -1,"2sv": -1, "2sd": -1},
    "A2u": {"E": 1, "2C4": 1,  "C2": 1,  "2C2'": -1, "2C2''": -1,
            "i": -1, "2S4": -1, "sh": -1,"2sv":  1, "2sd": 1},
    "B1u": {"E": 1, "2C4": -1, "C2": 1,  "2C2'": 1,  "2C2''": -1,
            "i": -1, "2S4":  1, "sh": -1,"2sv": -1, "2sd": 1},
    "B2u": {"E": 1, "2C4": -1, "C2": 1,  "2C2'": -1, "2C2''": 1,
            "i": -1, "2S4":  1, "sh": -1,"2sv":  1, "2sd": -1},
    "Eu":  {"E": 2, "2C4": 0,  "C2": -2, "2C2'": 0,  "2C2''": 0,
            "i": -2, "2S4":  0, "sh":  2,"2sv":  0, "2sd": 0},
}
IRREP_DIM: Dict[str, int] = {
    "A1g": 1, "A2g": 1, "B1g": 1, "B2g": 1, "Eg": 2,
    "A1u": 1, "A2u": 1, "B1u": 1, "B2u": 1, "Eu": 2,
}

# Map from element label → its conjugacy class label
ELEMENT_CLASS: Dict[str, str] = {
    "E":     "E",
    "C4":    "2C4",   "C43":   "2C4",
    "C2":    "C2",
    "Cv_x":  "2C2'",  "Cv_y":  "2C2'",
    "Cd_p":  "2C2''", "Cd_m":  "2C2''",
    "i":     "i",
    "S4":    "2S4",   "S43":   "2S4",
    "sh":    "sh",
    "sv_x":  "2sv",   "sv_y":  "2sv",
    "sd_p":  "2sd",   "sd_m":  "2sd",
}


# Group action on the basis ---------------------------------------------------

class HexBasisBlock(IntEnum):
    A = 0
    B = 1
    C = 2
    D = 3
    E = 4


def _C4_action(key: Tuple[int, int, int, int]) -> Tuple[Tuple[int, int, int, int], int]:
    """C4 (90° rotation about z): (ξ,η,ζ) → (-η,ξ,ζ), (vx,vy,vz) → (-vy,vx,vz).

    Derivation: substituting (η,-ξ,ζ) into A_{ijk} and applying the vector
    transform yields
        C4·A_{i,j,k} = (-1)^{j+1} B_{j,i,k}
        C4·B_{i,j,k} = (-1)^{j+1} A_{j,i,k}
        C4·C_{0,j,k} = (-1)^j    D_{j,0,k}
        C4·D_{i,0,k} = -1        C_{0,i,k}
        C4·E_{i,j,0} = (-1)^j    E_{j,i,0}
    """
    b, i, j, k = key
    if b == HexBasisBlock.A:
        return (HexBasisBlock.B, j, i, k), (-1) ** (j + 1)
    if b == HexBasisBlock.B:
        return (HexBasisBlock.A, j, i, k), (-1) ** (j + 1)
    if b == HexBasisBlock.C:
        return (HexBasisBlock.D, j, 0, k), (-1) ** j
    if b == HexBasisBlock.D:
        return (HexBasisBlock.C, 0, i, k), -1
    if b == HexBasisBlock.E:
        return (HexBasisBlock.E, j, i, 0), (-1) ** j
    raise ValueError(f"unknown block {b}")


def _sv_x_action(key: Tuple[int, int, int, int]) -> Tuple[Tuple[int, int, int, int], int]:
    """σ_v(x=0 plane): ξ → -ξ, (vx,vy,vz) → (-vx,vy,vz).

    Computes parity (-1)^{(parity of vx in ξ) + 1} (the +1 from vector flip).
    For each block, vx and vz share the same parity (consistent eigenvector).
    """
    b, i, j, k = key
    if b == HexBasisBlock.A:
        # vx = L[i] P[j+1] L'[k], parity of vx in ξ = (-1)^i;
        # then σ_v flips vx sign → overall phase (-1)^{i+1}.
        return key, (-1) ** (i + 1)
    if b == HexBasisBlock.B:
        # vy = P[i+1] L[j] L'[k], parity in ξ = (-1)^{i+1};
        # σ_v leaves vy alone → phase (-1)^{i+1}.
        return key, (-1) ** (i + 1)
    if b == HexBasisBlock.C:
        # No ξ dependence, vx=0, vy and vz invariant in ξ.
        # σ_v flips vx (which is 0) → phase +1.
        return key, +1
    if b == HexBasisBlock.D:
        # vx = -L[i] L'[k], parity in ξ = (-1)^i; σ_v flips vx → (-1)^{i+1}.
        return key, (-1) ** (i + 1)
    if b == HexBasisBlock.E:
        # vx = L[i] L'[j], parity in ξ = (-1)^i; σ_v flips vx → (-1)^{i+1}.
        return key, (-1) ** (i + 1)
    raise ValueError(f"unknown block {b}")


def _sv_y_action(key: Tuple[int, int, int, int]) -> Tuple[Tuple[int, int, int, int], int]:
    """σ_v(y=0 plane): η → -η, (vx,vy,vz) → (vx,-vy,vz).

    Symmetric counterpart to _sv_x: swap roles of i↔j and Block A↔B,
    Block C↔D (sort of).
    """
    b, i, j, k = key
    if b == HexBasisBlock.A:
        # vx = L[i] P[j+1] L'[k], parity in η = (-1)^{j+1};
        # σ_v(y) leaves vx alone → phase (-1)^{j+1}.
        return key, (-1) ** (j + 1)
    if b == HexBasisBlock.B:
        # vy = P[i+1] L[j] L'[k], parity in η = (-1)^j;
        # σ_v(y) flips vy → (-1)^{j+1}.
        return key, (-1) ** (j + 1)
    if b == HexBasisBlock.C:
        # vy = L[j] L'[k], parity in η = (-1)^j; σ_v(y) flips vy → (-1)^{j+1}.
        return key, (-1) ** (j + 1)
    if b == HexBasisBlock.D:
        # No η dependence; vx and vz invariant in η; σ_v(y) leaves vx,vz alone.
        return key, +1
    if b == HexBasisBlock.E:
        # vy = -L'[i] L[j], parity in η = (-1)^j; σ_v(y) flips vy → (-1)^{j+1}.
        return key, (-1) ** (j + 1)
    raise ValueError(f"unknown block {b}")


def _sh_action(key: Tuple[int, int, int, int]) -> Tuple[Tuple[int, int, int, int], int]:
    """σ_h (z = c/2 plane in physical coords, ζ → -ζ),
       (vx,vy,vz) → (vx,vy,-vz).

    Block A/B/C/D have ζ dependence; Block E does not.  After substitution
    and applying the vector flip on vz, the common phase is:
        A/B/C/D  →  (-1)^{k+1}
        E        →  +1
    """
    b, i, j, k = key
    if b == HexBasisBlock.E:
        return key, +1
    # Otherwise — A, B, C, D depend on ζ via L[k] / L'[k]:
    return key, (-1) ** (k + 1)


def _compose(*actions: Callable) -> Callable:
    """Compose group action callables left-to-right: f1 then f2 then ..."""
    def composed(key):
        kk, phase = key, 1
        for act in actions:
            kk, p = act(kk)
            phase *= p
        return kk, phase
    return composed


# Build all 16 D_4h actions from a generating set (C4, sv_x, sh)
GROUP_ACTIONS: Dict[str, Callable] = {}
GROUP_ACTIONS["E"]     = lambda k: (k, 1)
GROUP_ACTIONS["C4"]    = _C4_action
GROUP_ACTIONS["C2"]    = _compose(_C4_action, _C4_action)
GROUP_ACTIONS["C43"]   = _compose(_C4_action, _C4_action, _C4_action)
GROUP_ACTIONS["sv_x"]  = _sv_x_action
GROUP_ACTIONS["sv_y"]  = _sv_y_action
# σ_d(x+y) plane: ξ ↔ η reflection.  Implemented via σ_v(x) ∘ C4 (or C43 ∘ σ_v(x)).
# Choose σ_d_p ≡ σ_v(x) ∘ C4 (the formula gives ξ ↔ -η or ξ ↔ η depending on
# the convention).  We use: σ_d_p · v(ξ,η,ζ) = (vy(η,ξ,ζ), vx(η,ξ,ζ), vz(η,ξ,ζ)).
GROUP_ACTIONS["sd_p"]  = _compose(_C4_action, _sv_y_action)
GROUP_ACTIONS["sd_m"]  = _compose(_C4_action, _sv_x_action)
GROUP_ACTIONS["sh"]    = _sh_action
GROUP_ACTIONS["S4"]    = _compose(_C4_action, _sh_action)
GROUP_ACTIONS["S43"]   = _compose(_C4_action, _C4_action, _C4_action, _sh_action)
GROUP_ACTIONS["i"]     = _compose(_C4_action, _C4_action, _sh_action)
GROUP_ACTIONS["Cv_x"]  = _compose(_sv_x_action, _sh_action)
GROUP_ACTIONS["Cv_y"]  = _compose(_sv_y_action, _sh_action)
GROUP_ACTIONS["Cd_p"]  = _compose(GROUP_ACTIONS["sd_p"], _sh_action)
GROUP_ACTIONS["Cd_m"]  = _compose(GROUP_ACTIONS["sd_m"], _sh_action)


# Index helpers ---------------------------------------------------------------

def build_dof_keys(p: int) -> List[Tuple[int, int, int, int]]:
    """Mirror of HDivDivFreeHexBasisMP._build_index_map."""
    keys: List[Tuple[int, int, int, int]] = []
    for i in range(p):
        for j in range(p):
            for k in range(p):
                keys.append((HexBasisBlock.A, i, j, k))
    for i in range(p):
        for j in range(p):
            for k in range(p):
                keys.append((HexBasisBlock.B, i, j, k))
    for j in range(p):
        for k in range(p):
            keys.append((HexBasisBlock.C, 0, j, k))
    for i in range(p):
        for k in range(p):
            keys.append((HexBasisBlock.D, i, 0, k))
    for i in range(p):
        for j in range(p):
            keys.append((HexBasisBlock.E, i, j, 0))
    return keys


# Phase 1 ---------------------------------------------------------------------

def group_action_matrix(g: str, p: int) -> np.ndarray:
    """Return M_g ∈ ℝ^{n×n} such that  g · φ_α = Σ_β M_g[β, α] φ_β.

    Sparse permutation + sign (each column has exactly one nonzero ±1).
    Stored as a dense numpy array for downstream linear-algebra simplicity
    at p ≤ 7 (n ≤ 833).
    """
    keys = build_dof_keys(p)
    n = len(keys)
    M = np.zeros((n, n), dtype=np.int8)
    key_to_index = {k: i for i, k in enumerate(keys)}
    action = GROUP_ACTIONS[g]
    for alpha, key in enumerate(keys):
        new_key, phase = action(key)
        if new_key not in key_to_index:
            raise RuntimeError(
                f"group action {g} mapped {key} → {new_key} not in basis (p={p})"
            )
        beta = key_to_index[new_key]
        M[beta, alpha] = phase
    return M


def build_all_group_matrices(p: int) -> Dict[str, np.ndarray]:
    """Return {g: M_g} for all 16 D_4h elements at given order p."""
    return {g: group_action_matrix(g, p) for g in D4H_ELEMENTS}


# Phase 2 ---------------------------------------------------------------------

def character(irrep: str, g: str) -> int:
    """Look up χ_Γ(g) using the conjugacy-class character table."""
    return D4H_CHARACTER_TABLE[irrep][ELEMENT_CLASS[g]]


def projector(irrep: str, M_group: Dict[str, np.ndarray]) -> np.ndarray:
    """P_Γ = (d_Γ / |G|) Σ_g χ_Γ(g) M(g).  |G| = 16."""
    d = IRREP_DIM[irrep]
    G = 16
    n = next(iter(M_group.values())).shape[0]
    P = np.zeros((n, n), dtype=np.float64)
    for g, M_g in M_group.items():
        chi = character(irrep, g)
        if chi == 0:
            continue
        P = P + (chi * M_g.astype(np.float64))
    P *= (d / G)
    return P


def extract_irrep_basis(P: np.ndarray, *, tol: float = 1e-10) -> np.ndarray:
    """Given a projector P with eigenvalues in {0, 1}, return a basis U_Γ
    (n × d_Γ block) for the +1 eigenspace via eigendecomposition.
    """
    # P is symmetric for unitary group representation on real basis.
    w, V = np.linalg.eigh(0.5 * (P + P.T))
    keep = w > 1.0 - tol
    return V[:, keep]


# Phase 3 ---------------------------------------------------------------------

def block_diagonalize_basis(p: int) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """Assemble the global symmetry-adapted change-of-basis C ∈ ℝ^{n×n}.

    Returns
    -------
    C : ndarray (n × n)
        Each column is a symmetry-adapted basis vector belonging to one irrep.
    irrep_blocks : dict {irrep_name: ndarray (n × n_Γ)}
        Per-irrep basis vectors so that  C = [U_A1g, U_A2g, U_B1g, ...].
    """
    M_group = build_all_group_matrices(p)
    irrep_blocks: Dict[str, np.ndarray] = {}
    total_dim = 0
    for irrep in IRREP_DIM.keys():
        P = projector(irrep, M_group)
        U = extract_irrep_basis(P)
        if U.shape[1] == 0:
            continue
        irrep_blocks[irrep] = U
        total_dim += U.shape[1]
    n = M_group["E"].shape[0]
    assert total_dim == n, f"total irrep dims {total_dim} != n_dofs {n}"
    C = np.hstack(list(irrep_blocks.values()))
    return C, irrep_blocks


# Phase 4 ---------------------------------------------------------------------

def verify_group_matrices(p: int, *, atol: float = 1e-10) -> Dict[str, bool]:
    """Sanity checks: M(g)·M(h) = M(gh) for a sample of products."""
    M_group = build_all_group_matrices(p)
    checks = {}
    # Identity:
    n = M_group["E"].shape[0]
    checks["E = I"] = np.allclose(M_group["E"], np.eye(n))
    # C4^2 = C2:
    checks["C4^2 = C2"] = np.allclose(
        M_group["C4"] @ M_group["C4"], M_group["C2"], atol=atol)
    # C4^4 = E:
    checks["C4^4 = E"] = np.allclose(
        M_group["C4"] @ M_group["C4"] @ M_group["C4"] @ M_group["C4"],
        np.eye(n), atol=atol)
    # sv_x · sv_x = E:
    checks["sv_x^2 = E"] = np.allclose(
        M_group["sv_x"] @ M_group["sv_x"], np.eye(n), atol=atol)
    # sh^2 = E:
    checks["sh^2 = E"] = np.allclose(
        M_group["sh"] @ M_group["sh"], np.eye(n), atol=atol)
    # sh commutes with C4:
    checks["[sh, C4] = 0"] = np.allclose(
        M_group["sh"] @ M_group["C4"],
        M_group["C4"] @ M_group["sh"], atol=atol)
    # i = C2 · sh:
    checks["i = C2·sh"] = np.allclose(
        M_group["i"], M_group["C2"] @ M_group["sh"], atol=atol)
    return checks


def verify_irrep_dimensions(p: int) -> Dict[str, int]:
    """Report dim of each irrep in the basis decomposition."""
    M_group = build_all_group_matrices(p)
    out = {}
    for irrep in IRREP_DIM.keys():
        P = projector(irrep, M_group)
        U = extract_irrep_basis(P)
        out[irrep] = U.shape[1]
    return out


# Phase 3 — block-diagonal K, M, b in symmetry-adapted basis -----------------

def transform_to_irrep_block(
    K: np.ndarray, M: np.ndarray, b: np.ndarray,
    U_irrep: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Project K, M, b onto the column space spanned by U_irrep.

    Returns (K_block, M_block, b_block) of shape (d_irrep, d_irrep) and
    (d_irrep,) respectively, where d_irrep = U_irrep.shape[1].

    Inputs may be numpy float arrays OR mpmath matrices (treated by .tolist()
    coercion).  For mpmath inputs, the caller should convert beforehand.
    """
    Kb = U_irrep.T @ K @ U_irrep
    Mb = U_irrep.T @ M @ U_irrep
    bb = U_irrep.T @ b
    return Kb, Mb, bb


def block_diagonalize_KMb(
    K: np.ndarray, M: np.ndarray, b: np.ndarray, p: int, *,
    target_irreps: List[str] = None,
) -> Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Block-diagonalize K, M, b across all (or selected) irreps.

    Returns a dict {irrep_name: (K_block, M_block, b_block)}.

    target_irreps=None projects onto all irreps.  Otherwise only the listed
    irreps (e.g. ["A2g"] for B∥z) are processed.
    """
    _, irrep_blocks = block_diagonalize_basis(p)
    if target_irreps is None:
        target_irreps = list(irrep_blocks.keys())
    out = {}
    for irrep in target_irreps:
        if irrep not in irrep_blocks:
            continue
        U = irrep_blocks[irrep]
        out[irrep] = transform_to_irrep_block(K, M, b, U)
    return out


def get_irrep_basis(p: int, irrep: str) -> np.ndarray:
    """Return U_Γ ∈ ℝ^{n×d_Γ} for the requested irrep at order p."""
    _, blocks = block_diagonalize_basis(p)
    if irrep not in blocks:
        raise ValueError(f"irrep {irrep} not present in basis at p={p}")
    return blocks[irrep]


# Phase 4 (B-direction selection rule) ----------------------------------------

def b_irrep_assignment(b_dir: Tuple[float, float, float]) -> List[str]:
    """Return the irrep(s) the source vector b(B) belongs to.

    For applied B along z (B_z ẑ): A_ext = ½B_z (-y, x, 0) is axial → A_2g.
    For B along x (B_x x̂): A_ext = ½B_x (0, -z, y) is planar polar → E_u.
    For B along y (B_y ŷ): A_ext = ½B_y ( z, 0, -x) is planar polar → E_u.
    """
    Bx, By, Bz = b_dir
    irreps = []
    if abs(Bz) > 0:
        irreps.append("A2g")
    if abs(Bx) > 0 or abs(By) > 0:
        irreps.append("Eu")
    return irreps


if __name__ == "__main__":
    import argparse
    import sys
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    ap = argparse.ArgumentParser()
    ap.add_argument("--p", type=int, default=3,
                    help="Polynomial order (basis size = 2p^3+3p^2)")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    p = args.p
    keys = build_dof_keys(p)
    n = len(keys)
    print(f"HDivDivFreeHex order p={p}, n_dofs = {n} (= 2p³+3p² = {2*p**3+3*p**2})")
    print()

    print("=== Phase 1 sanity (group multiplication closure) ===")
    checks = verify_group_matrices(p)
    for k, v in checks.items():
        mark = "✓" if v else "✗"
        print(f"  {mark}  {k}")
    if not all(checks.values()):
        raise SystemExit("Phase 1 group action check FAILED")

    print()
    print("=== Phase 2 irrep dimension decomposition ===")
    dims = verify_irrep_dimensions(p)
    total = 0
    print(f"  {'Irrep':>6}  {'dim_Γ':>5}  {'mult':>4}  {'block_size':>10}")
    for irrep, d in dims.items():
        mult = d // IRREP_DIM[irrep]
        print(f"  {irrep:>6}  {IRREP_DIM[irrep]:>5}  {mult:>4}  {d:>10}")
        total += d
    print(f"  {'TOTAL':>6}  {'':>5}  {'':>4}  {total:>10}  (expected {n})")
    assert total == n

    print()
    print("=== Phase 3 block-diagonalize a random symmetric matrix ===")
    # Test: build random symmetric K, M, b; verify that for B∥z the
    # response only lives in A2g.
    np.random.seed(0)
    K_test = np.random.randn(n, n);  K_test = 0.5 * (K_test + K_test.T)
    M_test = np.eye(n) + 0.1 * np.random.randn(n, n)
    M_test = 0.5 * (M_test + M_test.T)
    # b mimicking A_z source: random vector projected onto A2g.
    U_A2g = get_irrep_basis(p, "A2g")
    b_test_random = np.random.randn(n)
    b_test_A2g_only = U_A2g @ U_A2g.T @ b_test_random   # project to A2g
    blocks = block_diagonalize_KMb(K_test, M_test, b_test_A2g_only, p)
    print(f"  {'irrep':>6}  {'dim':>4}  {'|b_irrep|':>12}")
    for irrep, (Kb, Mb, bb) in blocks.items():
        print(f"  {irrep:>6}  {Kb.shape[0]:>4}  {np.linalg.norm(bb):>12.4e}")
    print()
    print("  → All non-A2g blocks should have |b| ≈ 0 (B∥z selection rule).")

    print()
    print("=== Phase 4 B-direction selection rule ===")
    for bd in [(0, 0, 1), (1, 0, 0), (0, 1, 0), (1, 1, 1)]:
        print(f"  B={bd}  →  excited irreps: {b_irrep_assignment(bd)}")
