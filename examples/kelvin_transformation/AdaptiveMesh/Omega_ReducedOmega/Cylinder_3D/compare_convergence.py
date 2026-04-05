"""
3D円筒問題の収束履歴比較グラフ
3つの手法 x 3つのorder (p=2, 3, 4) を1x3のサブプロットで比較
Note: 円筒問題には解析解が存在しないため、誤差推定値の収束を比較
"""
import os
import numpy as np
import scipy.io as sio
import matplotlib
import matplotlib.pyplot as plt

# Times New Roman font settings
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['xtick.minor.size'] = 3
plt.rcParams['ytick.minor.size'] = 3

script_dir = os.path.dirname(os.path.abspath(__file__))

# Data structure: {order: {method: (ndof, error)}}
data = {}

methods = {
    'Refine_with_zz_estimator': ('Adaptive (ZZ Estimator)', 'bo-', 2, 6),
    'maxh': ('Maxh-based Remesh', 'rs-', 1.5, 5),
    'Refine_all_elements': ('Uniform (All Elements)', 'g^--', 1.5, 8),
}

orders = [2, 3, 4]

# Load data
for order in orders:
    data[order] = {}
    for method_dir, (label, marker, lw, ms) in methods.items():
        if method_dir == 'maxh':
            mat_file = os.path.join(script_dir, f"order={order}", method_dir,
                                    "Cylinder_3D_maxh_with_Kelvin.mat")
        else:
            mat_file = os.path.join(script_dir, f"order={order}", method_dir,
                                    "Cylinder_3D_adaptive_with_Kelvin.mat")

        if os.path.exists(mat_file):
            mat_data = sio.loadmat(mat_file)
            ndof = mat_data['ndof'].flatten()
            error = mat_data['error'].flatten()
            # Filter out NaN values
            valid = ~np.isnan(error)
            data[order][method_dir] = (ndof[valid], error[valid])
            print(f"Loaded: order={order}, {method_dir}: {len(ndof[valid])} points")
        else:
            print(f"Not found: {mat_file}")

# Create 1x3 subplot
fig, axes = plt.subplots(1, 3, figsize=(15, 5), dpi=150)

for idx, order in enumerate(orders):
    ax = axes[idx]

    # Plot each method
    for method_dir, (label, marker, lw, ms) in methods.items():
        if method_dir in data[order]:
            ndof, error = data[order][method_dir]
            if len(ndof) > 0:
                ax.loglog(ndof, error, marker, linewidth=lw, markersize=ms, label=label)

    # Reference line: O(N^{-p/3}) for 3D (instead of N^{-p/2} for 2D)
    ndof_line = np.array([1e2, 1e5])

    # Draw theoretical slope for current order only
    ref_styles = {2: ('k--', '$O(N^{-2/3})$'),
                  3: ('k--', '$O(N^{-1})$'),
                  4: ('k--', '$O(N^{-4/3})$')}
    if 'maxh' in data[order] and len(data[order]['maxh'][0]) > 0:
        # Use final point of metric-based for this order as reference
        N_ref = data[order]['maxh'][0][-1]  # final DOF
        err_ref = data[order]['maxh'][1][-1]  # final error
        err_line = err_ref * (N_ref / ndof_line) ** (order / 3)
        style, label = ref_styles[order]
        ax.loglog(ndof_line, err_line, style, linewidth=1, alpha=0.5, label=label)

    ax.set_xlabel('DOFs', fontsize=12)
    ax.set_ylabel('Error Estimator', fontsize=12)
    ax.set_title(f'$p = {order}$', fontsize=14, fontweight='bold')
    ax.set_xlim(1e2, 1e5)
    ax.set_ylim(1e-4, 1e-2)  # Fixed scale for all orders
    ax.grid(True, alpha=0.3)
    ax.legend(loc='lower left', fontsize=9)
    ax.tick_params(which='both', direction='in', top=True, right=True)

plt.suptitle('3D Cylinder (1/8 Model): Convergence Comparison (Omega-formulation + Kelvin)',
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()

# Save
output_file = os.path.join(script_dir, "convergence_comparison.png")
plt.savefig(output_file, dpi=150, bbox_inches='tight')
print(f"\nSaved: {output_file}")

# Open the file
os.startfile(output_file)

#plt.show()
