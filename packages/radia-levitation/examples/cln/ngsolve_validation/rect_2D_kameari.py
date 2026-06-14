"""Historical field-recurrence experiment retired.

Only the independent geometry, mesh and linear-algebra helpers remain:
build_rect.
There is no reduction experiment entry point or validation result in this file.

Historical method attribution: Kameari, Tanimoto.
"""
import json
from math import pi
from pathlib import Path

from netgen.geom2d import SplineGeometry
from ngsolve import (
    Mesh, H1, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, grad, dx,
    Integrate, TaskManager, ngsglobals,
)

mu0 = 4 * pi * 1e-7
sigma = 5.8e7
ax, ay = 5e-3, 2e-3
H_MESH = 0.05e-3   # 50 µm — finer mesh for higher precision
ORDER = 3          # cubic elements
N_STAGES = 12


def build_rect():
    geo = SplineGeometry()
    geo.AddRectangle((0, 0), (ax, ay), bcs=["bot", "right", "top", "left"])
    return geo







# Historical method attributions retained: Tanimoto.

# CLN experiment entry point retired; independent helpers above are retained.
