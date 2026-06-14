"""Debug K[26, 24] discrepancy: compare CuPy K vs C++ K at single-entry level."""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path("S:/Radia/01_GitHub/packages/radia-vim/src")))

from hex_vim_cupy import evaluate_basis
from hex_vim_cupy_kassembly import assemble_K_cupy
import radia_vim


def main():
    p = 2
    a, b, c = 5e-3, 2e-3, 1e-3
    basis = radia_vim.HDivDivFreeHexBasis(p)
    n = basis.n_dofs

    n_th, n_ph, n_rh, n_c = 8, 12, 8, 4

    K_cpu = assemble_K_cupy(p, a, b, c, n_th, n_ph, n_rh, n_c,
                             use_gpu=False, verbose=False)
    K_cpp = radia_vim.assemble_K_bare(basis, a, b, c,
                                      n_omega_theta=n_th, n_omega_phi=n_ph,
                                      n_rho=n_rh, n_c=n_c, verbose=0)

    # Compute single C++ entry directly
    val_cpp = radia_vim.compute_K_entry_spherical_duffy(
        basis, 26, 24, a, b, c, n_th, n_ph, n_rh, n_c)
    print(f"C++ compute_K_entry(26, 24) direct = {val_cpp:.6e}")
    print(f"C++ assemble_K K[26, 24]           = {K_cpp[26, 24]:.6e}")
    print(f"CuPy K_cpu[26, 24]                  = {K_cpu[26, 24]:.6e}")
    print()
    val_cpp_24_26 = radia_vim.compute_K_entry_spherical_duffy(
        basis, 24, 26, a, b, c, n_th, n_ph, n_rh, n_c)
    print(f"C++ compute_K_entry(24, 26) direct = {val_cpp_24_26:.6e}")
    print(f"C++ K[24, 26]                       = {K_cpp[24, 26]:.6e}")
    print(f"CuPy K_cpu[24, 26]                  = {K_cpu[24, 26]:.6e}")
    print()
    print("Note: assemble_K_bare uses upper triangle then mirrors, so")
    print("K[26, 24] = compute_K_entry(24, 26) (i<j convention).")
    print()

    # Sample a single point and evaluate basis
    pts_test = [(0.3, 0.4, 0.5), (0.5, 0.5, 0.5), (0.6, 0.7, 0.8)]
    print("--- Block E basis values at sample points ---")
    for x, y, z in pts_test:
        v24 = basis.evaluate(24, x, y, z)
        v26 = basis.evaluate(26, x, y, z)
        # CuPy basis
        x_arr = np.asarray([x])
        y_arr = np.asarray([y])
        z_arr = np.asarray([z])
        Phi_py = evaluate_basis(x_arr, y_arr, z_arr, p)  # (1, n, 3)
        v24_py = Phi_py[0, 24, :]
        v26_py = Phi_py[0, 26, :]
        print(f"  pt=({x}, {y}, {z}):")
        print(f"    [24] cpp=({v24[0]:.4e}, {v24[1]:.4e}, {v24[2]:.4e})")
        print(f"    [24] py =({v24_py[0]:.4e}, {v24_py[1]:.4e}, {v24_py[2]:.4e})")
        print(f"    [26] cpp=({v26[0]:.4e}, {v26[1]:.4e}, {v26[2]:.4e})")
        print(f"    [26] py =({v26_py[0]:.4e}, {v26_py[1]:.4e}, {v26_py[2]:.4e})")


if __name__ == "__main__":
    main()
