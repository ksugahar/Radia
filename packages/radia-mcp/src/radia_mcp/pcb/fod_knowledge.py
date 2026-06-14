"""Foreign Object Detection (FOD) — Sugahara lab core research area."""

OVERVIEW = r"""
# Foreign Object Detection (FOD) — ★ lab core research

When a metal object (coin, tool, screw) falls between primary and
secondary coils during WPT charging, it heats via induction and may
catch fire. **FOD is mandated by IEC 61980-3** before/during charging.

Lab has 53 papers in W:/.../05_ワイヤレス給電/10_異物検出/ across
6+ detection method families.

## Detection method families (lab classification)

| # | Family | Principle | Detection target |
|---|--------|-----------|------------------|
| 1 | **位置検知方式** (position) | Detect change in receiver position | Vehicle alignment |
| 2 | **センサレス方式** (sensorless) | Monitor primary V, I, freq, Q | Metal object in gap |
| 3 | **センサ等の方式** (sensor) | Search coil array | Metal & living object |
| 4 | **透磁率発熱** | Detect change in eff. permeability | Magnetic objects |
| 5 | **特許 (Qualcomm, TDK)** | Patented methods | Various |
| 6 | **IEC 61980-3 標準** | Standardized test items | Compliance |

Each family has trade-offs in sensitivity, false-alarm rate, cost.

## Cross-reference

- W:/.../05_ワイヤレス給電/10_異物検出/01_位置検知方式/ (10 papers)
- W:/.../05_ワイヤレス給電/10_異物検出/02_センサレス方式/ (15 papers)
- W:/.../05_ワイヤレス給電/10_異物検出/03_センサ等の方式/ (10 papers)
- W:/.../05_ワイヤレス給電/10_異物検出/04_透磁率発熱/ (2 papers)
- W:/.../05_ワイヤレス給電/10_異物検出/05_特許/ (4 patents)
- W:/.../05_ワイヤレス給電/10_異物検出/06_標準_IEC/ (3 IEC docs)
- W:/.../05_ワイヤレス給電/10_異物検出/07_レビュー/ (3 review papers)
"""


SENSORLESS = r"""
# Sensorless FOD methods (monitor primary V/I/freq)

Most cost-effective family — no extra hardware. Monitor primary-side
parameters (V, I, freq, phase, Q, efficiency) and detect changes
that indicate metal in the gap.

## Sub-methods

| Method | Signal | Sensitivity | False alarm |
|--------|--------|-------------|-------------|
| **Output端 位相** (output phase) | Phase of V_out vs I_in | High | Medium |
| **効率** (efficiency) | Input/output power ratio | Medium | Low (with load known) |
| **周波数特性** (frequency sweep) | Multiple freqs, detect Z shift | High | Low |
| **Q-factor monitoring** | Loaded Q before charging | High | Medium |
| **電圧電流** (V/I direct) | Detect anomalous I waveform | Low | High |

## Frequency sweep method (most rigorous)

Before main charging, sweep frequency 70-100 kHz and measure
input impedance Z_in(f). Build a baseline map. During operation,
compare to baseline:
- Coupled Z_in shifts due to metal → trigger FOD

References (lab):
- 04_ダイヘン_異物検出_2 周波数異物検知 (WPT2015-37)
- センサレス方式_周波数特性_Detection of Metal Obstacles WPT EV
- センサレス方式_周波数特性_Foreign Metal Detection by Coil Impedance for EV

## Two-frequency method

Drive at two frequencies simultaneously (e.g. main 85 kHz + probe 100 kHz).
- Main: power transfer
- Probe: sensitive to metal Z change

Pros: continuous monitoring, no interruption.
Cons: probe signal needs careful design to avoid efficiency loss.

## Efficiency-based detection

If P_in increases without P_out changing → loss from metal heating.
Simple but late (object already heating).

## Cross-reference

- `pcb_efficiency_safety('coupling_k')` — k changes due to metal
"""


SENSOR_BASED = r"""
# Sensor-based FOD methods (search coil array)

Add **search coil array** between primary and secondary. Each coil
senses local magnetic field; metal disturbs the pattern.

## Coil array geometries

| Pattern | When | Pro / Con |
|---------|------|-----------|
| **Non-overlapping** ★ | EV (lab) | Simple decode, less crosstalk |
| **Overlapping** | Phone Qi | Higher resolution, more wires |
| **Symmetric pair** | Detection only | Cancels main field signal |
| **Hexagonal** | Large area | Coverage uniformity |

## Lab papers (10_異物検出/03_センサ等の方式/)

- センサ等の方式_サーチコイル_Dual-Purpose Non-overlapping Coil Sets ★
  - Same array detects BOTH metal AND vehicle position
  - Critical insight for EV (one sensor system instead of two)
- センサ等の方式_サーチコイル_Self-Inductance-Based Metal Object Detection
  - Use mistuned LC resonance to null out main field, sense metal
- センサ等の方式_カメラ・レーダー_Video-based Foreign Object Debris
  - Optical/radar alternative (more expensive, more reliable for living objects)
- センサ等の方式_静電容量_Living Object Detection Based on Comb Pattern Capacitive Sensor
  - Detect humans/animals via capacitance change (high-permittivity body)

## Living Object Detection (LOD)

Subset of sensor-based FOD focused on humans/animals (safety-critical):
- Capacitive (skin permittivity)
- Thermal (body heat)
- Camera (visual)
- Doppler radar (movement)

## Decision tree

```
What to detect?
├── Metal object only (FOD)
│   ├── Cost-sensitive (Qi) → sensorless efficiency
│   ├── EV high-power → search coil array (non-overlapping)
│   └── Industrial AGV → sensor combined with positioning
│
└── Living object (LOD) ★ safety
    ├── Capacitive sensor (lab + production)
    ├── Camera (more cost) — also catches metal
    └── Doppler radar — pet detection
```

## Reference

W:/.../05_ワイヤレス給電/10_異物検出/03_センサ等の方式/ (10 papers)
"""


PERMEABILITY_PATENTS = r"""
# Permeability-change FOD + patent methods

## Permeability method (透磁率による発熱)

For magnetic metal objects (iron, nickel), permeability change of
the gap region affects:
- Primary inductance L_1 (small, ~0.1%)
- Coupling k (small but measurable)
- Heating rate (proportional to μ' σ for ferromagnetic metals)

Trade-off: very sensitive to MAGNETIC objects (iron screw) but NOT
to non-magnetic conductors (aluminum can, copper coin). Combine with
other methods for full coverage.

W:/.../05_ワイヤレス給電/10_異物検出/04_透磁率発熱/ (2 papers)
- Comparative Study of Metal Obstacle Variations in Disturbing WPT

## Patent methods (5_特許/)

- **Qualcomm 異物検出特許** — multi-coil differential sensing
- **TDK 異物検出** — saturable reactor-based pre-charge check

## IEC 61980-3 test items (06_標準_IEC/)

| Test item | What it checks | Pass criteria |
|-----------|----------------|---------------|
| 項目1 | Steel washer detection | < N% efficiency loss |
| 項目2 | Coin detection | < N° C temp rise |
| 項目3 | Aluminum foil detection | Auto-shutdown < 2 s |

Lab is targeting IEC 61980-3 compliance for production WPT systems.

## Review papers (07_レビュー/)

- A review of foreign object detection (FOD) for inductive power
- ワイヤレス電力伝送における異物検出と周囲技術
- 非接触給電における金属遺物の発熱量の最大値に関する研究

## Lab Sugahara FOD research direction

- **Multi-method fusion**: combine sensorless freq sweep + search coil
- **Living object detection**: capacitive + thermal
- **Real-time online**: continuous monitoring during charging
- **Daihen 共研, Panasonic 共研** active collaborations (W:/.../08_ダイヘン共研/)

## Cross-reference

- `pcb_efficiency_safety('safety_human_exposure')` — LOD aspect
- `radia_mcp.peec` — search coil L computation
- `radia_mcp.ih` — analogous metal heating physics
"""


def get_fod_knowledge(topic: str = "overview") -> str:
    """Dispatch FOD topics.

    Topics:
        overview                 - 6 FOD method families (DEFAULT)
        sensorless               - V/I/freq monitoring methods
        sensor_based             - Search coil array + camera + radar
        permeability_patents     - Permeability change + patents + IEC 61980-3
        all                      - Everything
    """
    topic = topic.lower().strip()
    if topic in ("overview", "all_families"):
        return OVERVIEW
    if topic in ("sensorless", "sensorless_methods"):
        return SENSORLESS
    if topic in ("sensor_based", "sensor", "search_coil"):
        return SENSOR_BASED
    if topic in ("permeability_patents", "permeability", "patents", "iec"):
        return PERMEABILITY_PATENTS
    if topic == "all":
        return "\n\n".join([OVERVIEW, SENSORLESS, SENSOR_BASED,
                            PERMEABILITY_PATENTS])
    return (f"Unknown topic '{topic}'. Available: overview, sensorless, "
            "sensor_based, permeability_patents, all.")
