# Beam Tracking Examples

Particle trajectory computation in magnetic fields for accelerator physics applications.

## Overview

The `beam_tracking` module provides NumPy-based interface for charged particle trajectory computation. It wraps the core Radia functions (`FldPtcTrj`, `FldFocPot`, etc.) with numpy array output for easy data analysis.

## Examples

### demo_particle_trajectory.py

Basic particle trajectory computation through a dipole magnet.

- Creates a simple H-type dipole magnet
- Computes electron trajectory using Runge-Kutta integration
- Demonstrates batch trajectory computation for phase space scans
- Visualizes trajectory with matplotlib

```bash
python demo_particle_trajectory.py
```

### demo_undulator_trajectory.py

Electron trajectory in a planar undulator.

- Creates a simplified planar undulator (periodic magnet array)
- Computes electron oscillation through the undulator
- Analyzes trajectory to extract K parameter
- Compares with theoretical predictions

```bash
python demo_undulator_trajectory.py
```

## Usage

```python
import radia as rad
from radia import beam_tracking as bt
import numpy as np

# Create magnetic field source (coordinates in meters)
magnet = rad.ObjRecMag([0, 0, 0], [0.1, 0.01, 0.02], [0, 0, 954930])

# Compute trajectory (returns numpy array)
trajectory = bt.particle_trajectory(
    magnet,
    energy_gev=6.0,
    initial_conditions=[0, 0, 0, 0],  # [x0, x'0, z0, z'0]
    s_range=[-500, 500],
    n_points=1001
)

# trajectory.shape = (1001, 5)
# Columns: [s, x, x', z, z']
s = trajectory[:, 0]   # Longitudinal position [mm]
x = trajectory[:, 1]   # Horizontal position [mm]
xp = trajectory[:, 2]  # Horizontal angle [rad]
z = trajectory[:, 3]   # Vertical position [mm]
zp = trajectory[:, 4]  # Vertical angle [rad]
```

## Available Functions

| Function | Description |
|----------|-------------|
| `particle_trajectory()` | Compute single particle trajectory |
| `batch_trajectory()` | Compute multiple trajectories (phase space scan) |
| `focusing_potential()` | Compute focusing potential along a line |
| `shim_signature()` | Compute field variation from element displacement |
| `focusing_kick_periodic()` | Compute 2nd order kick matrices for periodic fields |

## Physical Model

The trajectory computation uses the relativistic Lorentz force equation:

```
dp/dt = q(v x B)
```

Integration is performed using Runge-Kutta methods (RK4 or adaptive RK5).

**Scope**: Magnetic field only (no electric field). Energy is conserved during tracking.

## Future Extensions (2026 Roadmap)

- Hamiltonian formulation with vector potential A
- Symplectic integrators (Forest-Ruth, Yoshida)
- GPU parallel batch tracking
