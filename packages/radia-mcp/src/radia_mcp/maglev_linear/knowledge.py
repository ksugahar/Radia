"""Maglev + linear drive knowledge."""


# Authoritative topic enum for the dispatcher tool (wired into
# `<short>_topics()` via common.register_topics_tool).
TOPICS: dict[str, str] = {
    "overview": "(see knowledge file for details)",
    "maglev": "(see knowledge file for details)",
    "linear": "(see knowledge file for details)",
    "all": "return \"\n\n\".join([OVERVIEW, MAGLEV, LINEAR])",
}

OVERVIEW = r"""
# Magnetic levitation + linear drive landscape

## Levitation principles

| Method | Force source | Examples | Stability |
|--------|--------------|----------|-----------|
| **EMS** (Electromagnetic Suspension) | Attractive (controlled current) | Transrapid, JR maglev (low-speed) | Active (unstable, needs control) |
| **EDS** (Electrodynamic Suspension) | Repulsive (induced eddy in track) | SCMaglev (Chuo Shinkansen) | Self-stable above threshold velocity |
| **Superconducting EDS** | Repulsive (SC coil on vehicle) | SCMaglev | Self-stable |
| **Permanent magnet repulsion** | Repulsive (Halbach array) | Halbach maglev, Inductrack | Self-stable above threshold |
| **Diamagnetic** | Repulsive (bismuth, etc.) | physics demos | Stable but small force |
| **Hybrid (PM + EM)** | Bias by PM, control by EM | Modern industrial | Active near zero current |

## Linear drive types

| Type | Principle | Example |
|------|-----------|---------|
| **LIM** (Linear Induction Motor) | Stator AC, secondary plate eddy current | airport trams, JR Linimo |
| **LSM** (Linear Synchronous Motor) | Stator AC, rotor PMs (or SC coils) | SCMaglev, MagLev launch |
| **DC linear** | DC excitation | small actuators |
| **Linear stepper** | Pulse-fed | precision positioning |

## Lab focus

- ★ Bearingless motor + WPT (cross-link to `radia_mcp.motor` + `radia_mcp.wpt`)
- Maglev research (磁気浮上 papers)
- Linear induction motor (lab Sugahara group has some papers)
"""


MAGLEV = r"""
# Magnetic levitation (磁気浮上)

## EMS — Electromagnetic Suspension

DC electromagnet attracts ferromagnetic track from above (vehicle
beneath rail). Control system maintains ~10 mm gap.

```
F_lift = B² A / (2 μ_0)    (Maxwell stress on iron plate)
```
where A = pole face area, B = pole face field.

## EDS — Electrodynamic Suspension

Induced eddy currents in track repel moving permanent magnet (or
superconducting coil) on vehicle.

```
F_lift ~ v² / (v² + v_crit²)    (saturates at high v)
```
Self-stable but needs auxiliary wheels at low speed (< v_crit ≈ 100 km/h).

## SCMaglev (Chuo Shinkansen, ★ Japan)

Superconducting EDS with Nb-Ti coils on vehicle, ground coils as track:
- Operating: 500 km/h
- Levitation gap: 10 cm (much larger than EMS)
- Active power steering via 8-shape coils

## Halbach array maglev (Inductrack)

PM Halbach array on vehicle + passive track coils.
- No on-board power for levitation
- Demonstrated by LLNL (US)

## Lab papers (07_磁気浮上/)

W:/.../99_アプリケーション/07_磁気浮上/:
- 永久磁石磁気浮上 (PM levitation)
- 渦電流で磁気浮上 (eddy current levitation)
- Characterization of Eddy-Current Position Sensing for Magnetic Levitation
- Winding Optimization of a Magnetic Levitation Induction Motor (Multipolar Suspension)
- 関西大学 papers (lab collaboration)

## Lab specialty: bearingless motor + WPT

```
Stator coils provide BOTH:
  - Magnetic levitation force (radial bearing replacement)
  - Wireless power transfer to rotor
Rotor has:
  - Receiver coil + rectifier
  - Motor windings (drive torque)
  - No physical bearings
```

References (multiple lab papers):
- 整流コイルを用いたベアリングレスモータ
- 磁界共鳴給電によるロータ磁化方式ベアリングレスモータ
- 非接触給電を用いたベアリングレスモータの回転機構
- 非接触給電を用いたベアリングレスモータの開発
- ワイヤレスホイールモータ (Fujimoto, JSAE)

## Cross-references

- `radia_mcp.wpt.applications.robot_bearingless` — bearingless motor + WPT detail
- `radia_mcp.motor` — motor design knowledge
- `radia_mcp.team_benchmark.force_motion.problem_28` — Electrodynamic Levitation benchmark
"""


LINEAR = r"""
# Linear drives (リニアドライブ)

## Linear Induction Motor (LIM)

Slotted primary stator with traveling magnetic wave; secondary is a
conductive plate (aluminum + iron back).

```
F_thrust = (V²/R) * s / (1 + s²)    (slip-thrust relation, simplified)
```
where V = induced voltage, R = secondary resistance, s = slip.

Issues:
- **End effect** — finite primary length distorts field at entry/exit
- **Air gap** — much larger than rotary motor (5-30 mm vs <1 mm)
- **Efficiency** — lower (50-70%) than rotary IM due to end effects

Reference:
- Basic Consideration of End Effect Compensator of Linear Induction Motor for Transit
  (in W:/.../06_非破壊検査/ misc folder, related to LIM end effects)

## Linear Synchronous Motor (LSM)

Stator excited with 3-phase AC; mover has permanent magnets or SC coils.
Synchronous = speed proportional to AC frequency.

```
v = 2 * τ_pole * f    (synchronous speed, τ_pole = pole pitch)
```

Used in:
- SCMaglev (rotor = SC coils)
- High-precision positioning (rotor = NdFeB array)
- Magnetic launch systems

## Lab papers (09_リニアドライブ/)

W:/.../99_アプリケーション/09_リニアドライブ/ (32 files, 555 MB):
- Various linear drive papers, transit applications, end-effect studies

## Cross-references

- `radia_mcp.motor` — analogous rotary motor analysis
- `radia_mcp.accelerator` — linear accelerator (different physics)
- `radia_mcp.team_benchmark.special.motors_30b_34` — IM analytic references
"""


def get_knowledge(topic: str = "overview") -> str:
    """Dispatch maglev / linear drive topics.

    Topics:
        overview     - Levitation + linear drive landscape (DEFAULT)
        maglev       - Magnetic levitation (EMS, EDS, SCMaglev, Halbach, bearingless)
        linear       - Linear induction motor (LIM), linear synchronous motor (LSM)
        all          - Everything
    """
    topic = topic.lower().strip()
    if topic in ("overview", "intro"):
        return OVERVIEW
    if topic in ("maglev", "levitation", "磁気浮上"):
        return MAGLEV
    if topic in ("linear", "lim", "lsm", "リニア"):
        return LINEAR
    if topic == "all":
        return "\n\n".join([OVERVIEW, MAGLEV, LINEAR])
    return f"Unknown topic '{topic}'. Available: overview, maglev, linear, all."
