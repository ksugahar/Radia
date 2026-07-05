"""WPT applications: dynamic EV charging, robot, bearingless motor, lab-specific."""

DYNAMIC_EV = r"""
# Dynamic EV charging (走行中給電)

Charge EV while it's moving (no stop required). Key challenges:
- Continuous coupling as vehicle moves over primary array
- Rapid k variation → impedance control bandwidth needed
- Power factor maintenance under varying load
- Stray field outside vehicle (safety)

## Primary configurations

| Architecture | When |
|--------------|------|
| **Long-track single primary** | Old OLEV (Korea), shorter coil pad |
| **Segmented (multiple short coils)** | Modern EV dynamic, switched as vehicle passes |
| **Continuous DD coil array** | Highest power, complex switching |

## Control challenges

| Challenge | Solution |
|-----------|----------|
| k(t) varies fast | LCC compensation (k-insensitive constant current) |
| Receiver must hand-off between primary segments | Position-aware switching, overlap windows |
| Power factor drops during transition | Feed-forward Q control |
| Reactive power surges | DC-link voltage buffer |

## Lab papers (04_dynamic_charging in 03_IEICE/)

7 papers covering:
- Dynamic primary with feed-forward control
- Multi-leg inverter with real-time reflected load detection
- Segmental supply edge position detection
- Sensor-less position detection (cross-cited with FOD)

## References

public-safe curated corpus (7 papers)
public-safe curated corpus (10 papers, dual-use)
"""


ROBOT_BEARINGLESS = r"""
# Robot WPT + bearingless motor (★ lab specialty)

## Multi-joint arm robot WPT

Powering wrist tools through multiple rotating joints without slip
rings or trailing cables.

Key lab papers (99_その他/論文_ロボット/):
- 多関節アームロボットへの電力・信号同時伝送システム
- ディスクリピーター for arm robot (WPT2016-23, 26, 39)
- 製造業ロボットにおける省配線・ワイヤレス給電の動向
- シリアル機構の解析

Approach: chain of LCC-LCC stages, one per joint axis. Each stage
allows rotation around its axis while transferring power.

## Bearingless motor + WPT

A bearingless motor has the rotor supported by **magnetic levitation
forces** instead of physical bearings. Combine with WPT:
- Stator: stationary coil for both levitation AND power
- Rotor: receiver coil + rectifier + load (sensors, motor windings)

Lab paper series:
- 整流コイルを用いたベアリングレスモータ
- 磁界共鳴給電によるロータ磁化方式ベアリングレスモータ
- 非接触給電を用いたベアリングレスモータの回転機構
- 非接触給電を用いたベアリングレスモータの開発
- ワイヤレスホイールモータ (Fujimoto JSAE)

Cross-references:
- `radia_mcp.motor` — motor design knowledge
- `radia_mcp.electromagnet` — Hantila + Play hysteresis for stator-rotor coupling

## References

public-safe curated corpus (6 papers)
public-safe curated corpus (5 papers)
"""


WPT_LAB_LINEAGE = r"""
# Sugahara lab WPT paper lineage

## Sugahara lab WPT-relevant authors

| Author | Topic | Years |
|--------|-------|-------|
| 田中 (Tanaka) | RNA + WPT | 2014-2017 (マグ研 papers) |
| 羽根 (Hane) | PEEC + WPT | 2017-2019 |
| 久田 (Hisada) | 修士論文_本審査 久田 |  |
| 吉田 (Yoshida) | 博士論文_本審査 吉田 |  |

(All in public-safe curated corpus)

## Lab core methods for WPT analysis

1. **PEEC for primary/secondary L, M, k** → `radia_mcp.peec`
2. **Carstensen AC copper loss** for Litz wire R → `radia_mcp.peec.carstensen`
3. **RNA + PEEC coupling** for transformer-like core models →
   `radia_mcp.bem.mmm_msc.rna_mmm`
4. **MOR / equivalent circuit** for time-domain control →
   `radia_mcp.mor.systematic` (CLN, Lanczos, PRIMA)
5. **FOD via search coil arrays** (research direction)

## Implementation pattern (lab production)

```python
import radia as rad
from radia.peec_matrices import PyPEECBuilder

# Build primary + secondary as PEEC filaments
builder = PyPEECBuilder()
# ... add coils ...
topology = builder.build_topology()

# Compute L, M, k vs frequency
from radia.peec_topology import PEECCircuitSolver
solver = PEECCircuitSolver(topology)
Y = solver.compute_admittance_matrix(freq=85e3)
# Extract L, M, k from Y

# AC copper loss
from radia.carstensen_loss import compute_carstensen_loss
loss = compute_carstensen_loss(coil_geometry, freq, current)

# Couple to ferrite via Radia MMM
ferrite = rad.ObjHexahedron(verts, [0, 0, 0])  # initial M=0
rad.MatApl(ferrite, rad.MatSatIsoTab(BH_DATA))
# Combined PEEC + MMM solve...
```

## References

public-safe curated corpus (90 papers, lab WPT IEICE corpus)
public-safe curated corpus (10+ files, lab RNA-WPT lineage)
public-safe curated corpus (44 papers, Litz wire physics)
"""


def get_applications_knowledge(topic: str = "dynamic_ev") -> str:
    """Dispatch WPT applications topics.

    Topics:
        dynamic_ev               - EV running, segmented primary (DEFAULT)
        robot_bearingless        - ★ lab specialty (multi-joint robot, bearingless motor)
        pcb_lab_lineage          - Sugahara lab paper lineage + production methods
        all                      - Everything
    """
    topic = topic.lower().strip()
    if topic in ("dynamic_ev", "dynamic", "ev", "moving"):
        return DYNAMIC_EV
    if topic in ("robot_bearingless", "robot", "bearingless"):
        return ROBOT_BEARINGLESS
    if topic in ("pcb_lab_lineage", "lab", "lineage", "sugahara"):
        return WPT_LAB_LINEAGE
    if topic == "all":
        return "\n\n".join([DYNAMIC_EV, ROBOT_BEARINGLESS, WPT_LAB_LINEAGE])
    return (f"Unknown topic '{topic}'. Available: dynamic_ev, "
            "robot_bearingless, pcb_lab_lineage, all.")
