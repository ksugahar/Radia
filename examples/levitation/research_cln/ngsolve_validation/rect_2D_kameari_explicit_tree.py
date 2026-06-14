"""Historical field-recurrence experiment retired.

Only the independent geometry, mesh and linear-algebra helpers remain:
build_rect, build_explicit_tree_mask.
There is no reduction experiment entry point or validation result in this file.

Historical method attribution: Kameari.
"""
import json
from collections import deque
from math import pi
from pathlib import Path

from netgen.geom2d import SplineGeometry
from ngsolve import (
    Mesh, H1, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, grad, dx, BitArray, BND,
    Integrate, TaskManager, ngsglobals,
)

mu0 = 4 * pi * 1e-7
sigma = 5.8e7
ax, ay = 5e-3, 2e-3
H_MESH = 0.05e-3   # match rect_2D_kameari.py for bit-identity check
ORDER = 3
N_STAGES = 12


def build_rect():
    geo = SplineGeometry()
    # No bcs assigned — we will mask explicitly
    geo.AddRectangle((0, 0), (ax, ay), bcs=["b", "r", "t", "l"])
    return geo


def build_explicit_tree_mask(fes, mesh):
    """Build H1 'tree' by masking all boundary vertices.

    In 2D, the boundary of the conductor is a closed cycle. Masking all
    boundary vertices = adding the boundary cycle to the spanning tree =
    full Dirichlet. This is the natural 2D analog of 3D BFS-tree masking
    on HCurl.
    """
    free = BitArray(fes.FreeDofs())
    n_masked = 0
    # Mask all DoFs (vertex + edge interior + ...) of boundary elements
    for el in mesh.Elements(BND):
        for d in fes.GetDofNrs(el):
            if d >= 0 and free[d]:
                free[d] = False
                n_masked += 1
    return free, n_masked







# CLN experiment entry point retired; independent helpers above are retained.
