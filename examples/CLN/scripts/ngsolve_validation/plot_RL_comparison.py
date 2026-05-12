"""Compare R_n, L_n between NGSolve Kameari (tree-cotree+GS) and Mathematica Cauer-I.

Note: these two methods produce DIFFERENT ladder forms (Krylov-Kameari vs Cauer-I),
so coefficient-by-coefficient agreement is NOT expected. The Z(s) frequency response
should match in the limits, however, when accumulated through their respective
network topologies.
"""
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

OUT_DIR = Path(__file__).parent.parent

# NGSolve Kameari (tree-cotree + Gram-Schmidt, h=0.20mm, order 2)
ngsolve_R = np.array([
    1.7241,    # R_0 = 1/(sigma V) directly
    0.6026,
    3.2119,
    9.8078,
    28.148,
    133.70,
    758.03,
    2166.4,
    8714.7,
    4.763e4,
    2.161e5,
    7.459e5,
])
ngsolve_L_nH = np.array([
    4666.5,
    33320,
    23808,
    26269,
    38275,
    26999,
    24623,
    31642,
    28172,
    23202,
    20466,
    20678,
])

# Mathematica Cauer-I (125 modes = {1,3,5,7,9}^3, exact rational + symbolic Pi)
math_R = np.array([
    1.951168559832186,
    7.447179062839928,
    11.25055650340680,
    17.96982314267895,
    35.87057476477225,
    36.56093930906552,
    157.6774050444724,
    223.1254976399821,
    474.8206047116640,
    996.8382693594679,
    3131.460955952874,
    5139.408668127379,
])
math_L_nH = np.array([
    8142.175028417007,
    9379.606521589173,
    18836.21054209553,
    34469.37169432091,
    33201.28986845171,
    60548.45900026239,
    112820.0239315821,
    155831.9783195020,
    365899.4460713181,
    1.061614832939471e6,
    1.990968173418074e6,
    4.827025600338959e6,
])

n_ngs = np.arange(len(ngsolve_R))
n_math = np.arange(len(math_R))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

ax1.semilogy(n_ngs, ngsolve_R, "o-", label="NGSolve Kameari (tree-cotree+GS)", markersize=8)
ax1.semilogy(n_math, math_R, "s--", label="Mathematica Cauer-I (125 modes)", markersize=8)
ax1.set_xlabel("stage $n$")
ax1.set_ylabel(r"$R_n\ [\Omega]$")
ax1.set_title("Resistance ladder: NGSolve Kameari vs Mathematica Cauer-I")
ax1.grid(True, which="both", alpha=0.3)
ax1.legend(fontsize=9)

ax2.semilogy(n_ngs, ngsolve_L_nH, "o-", label="NGSolve Kameari (tree-cotree+GS)", markersize=8)
ax2.semilogy(n_math, math_L_nH, "s--", label="Mathematica Cauer-I (125 modes)", markersize=8)
ax2.set_xlabel("stage $n$")
ax2.set_ylabel(r"$L_n$ [nH]")
ax2.set_title("Inductance ladder: NGSolve Kameari vs Mathematica Cauer-I")
ax2.grid(True, which="both", alpha=0.3)
ax2.legend(fontsize=9)

fig.suptitle("5×2×1 mm Cu cuboid Case A — circuit constants comparison\n"
             "(Krylov-Kameari and Cauer-I are different ladder forms; matching only at Stage 0 R_DC limit)",
             fontsize=10)
fig.tight_layout()

out = OUT_DIR / "RL_comparison_ngsolve_vs_mathematica.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"Wrote {out}")
