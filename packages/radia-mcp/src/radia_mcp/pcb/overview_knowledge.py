"""WPT overview: landscape, decision tree, lab focus."""

LANDSCAPE = r"""
# Wireless Power Transfer (WPT) landscape

## WPT regimes (by frequency / distance / power)

| Regime | Frequency | Distance | Efficiency | Power | Typical use |
|--------|-----------|----------|------------|-------|-------------|
| **Inductive (Qi)** | 100-200 kHz | < 5 mm | 70-90% | < 15 W | Phone charging |
| **Magnetic resonant** ★ | 85 kHz (SAE J2954) | 5-30 cm | 85-95% | 1-22 kW | EV charging (stationary) |
| **Dynamic resonant** | 85 kHz | < 30 cm + moving | 75-90% | 10-100 kW | EV running, AGV |
| **Capacitive** | 1-10 MHz | < 5 cm | 70-90% | < 100 W | Tablet, small EV |
| **Microwave (rectenna)** | 2.45 / 5.8 / 28 GHz | meters - km | 30-60% | < 1 kW | Drone, satellite |
| **Laser / optical** | THz | km | 20-50% | research | Space-based |

★ = main lab focus (kHz magnetic resonant for EV / IH / robot WPT)

## Lab focus areas (production research)

1. **Coil design** for kHz resonant — Litz wire, planar / cylindrical
2. **Foreign Object Detection (FOD)** — critical for safety, multiple methods
3. **Dynamic charging** — moving load (EV, robot)
4. **Bearingless motor + WPT** — rotor-side power without slip rings
5. **Standards compliance** — IEC 61980, SAE J2954

## Cross-references

- `radia_mcp.peec` — coil inductance via PEEC filaments (kQ product)
- `radia_mcp.ih` — analogous SIBC / eddy current (workpiece = receiver coil)
- `radia_mcp.bem` — open-boundary BEM for coil field evaluation
- `radia_mcp.matrix_solvers` — coupled coil + receiver system solvers
- `radia_mcp.mor` — model order reduction for transient WPT simulation
"""


DECISION_TREE = r"""
# Decision tree: which WPT topic for which design question

```
1. What's the design question?
   ├── Coil geometry / sizing
   │   → pcb_coil_compensation('coil_design')
   │   → pcb_efficiency_safety('q_factor_kQ')
   │
   ├── Compensation topology (SS / LCC / LCL)
   │   → pcb_coil_compensation('compensation_topology')
   │
   ├── Resonant frequency / matching
   │   → pcb_coil_compensation('resonance_matching')
   │
   ├── Efficiency analysis
   │   → pcb_efficiency_safety('q_factor_kQ')
   │   → pcb_efficiency_safety('coupling_k')
   │
   ├── Foreign Object Detection (FOD) ★ lab core
   │   → pcb_fod('overview') — 5+ detection method families
   │
   ├── Dynamic charging (moving EV / robot)
   │   → pcb_applications('dynamic_charging_ev')
   │
   ├── Safety / human exposure
   │   → pcb_efficiency_safety('safety_human_exposure')
   │   → pcb_efficiency_safety('standards_IEC_SAE')
   │
   ├── Alternative: capacitive WPT
   │   → pcb_alternatives('capacitive')
   │
   ├── Alternative: microwave rectenna
   │   → pcb_alternatives('microwave_rectenna')
   │
   └── Bearingless motor + WPT (lab specialty)
       → pcb_applications('bearingless_motor')

2. Implementation
   ├── Coil inductance / mutual L → radia_mcp.peec MCP
   ├── Field evaluation → radia_mcp.bem / radia_mcp.radia_ngsolve
   ├── Eddy current in metal load → radia_mcp.ih (SIBC)
   └── Loss in Litz wire → radia_mcp.peec.carstensen_ac_copper_loss
```

## Lab-specific decision points

| Situation | Recommendation |
|-----------|----------------|
| Designing 85 kHz EV stationary charging | SS compensation + ferrite plate + Litz wire |
| 100 kHz robot WPT with rotation | LCC compensation + thin coil + repeater |
| 6.78 MHz IoT WPT | Higher Q metamaterial coil + smaller capacitor |
| Safety verification | IEC 61980-1 induction current limit |
| FOD for metal AGV environment | Symmetric coil array + Δ I/V detection |
"""


HISTORY = r"""
# WPT history / lab genealogy

```
1864  Maxwell — theoretical foundation
1899  Tesla — early experiments (Wardenclyffe Tower)
1960s Microwave WPT theory (Brown, Glaser)
1980s Inductive charging for medical implants
1990s Electric toothbrush, AGV inductive charging (Daifuku)
2007  MIT WiTricity demo (resonant magnetic coupling, 60 W @ 2m)
2009  Qi standard (Wireless Power Consortium)
2011  Smartphone Qi adoption begins
2014  IEC TS 61980-1 (general EV WPT)
2016  Toyota Prius PHV adopts WPT charging
2018  SAE J2954 (EV stationary WPT) released
2020  Dynamic charging trials (Stellantis, OLEV Korea)
2024  IEC 61980-3 (dynamic) FDIS
```

## Lab Sugahara WPT lineage (selected)

```
~2010  Daifuku conveyor non-contact power (early industrial)
2014  田中 (lab member) WPT papers begin
2015-2018  田中 / 羽根 (Tanaka, Hane) — マグ研 papers on RNA + WPT
2016  ディスクリピーター for arm robot (WPT2016-23, 26, 39)
2018+  Bearingless motor + WPT integration (multiple lab papers)
2020+  FOD research (Daihen 共研 + Pana共研)
```

## Cross-references

- `radia_mcp.peec` — Tanaka-Hane PEEC + RNA coupling lineage
- `radia_mcp.ih` — Daifuku conveyor lineage (induction heating sibling)
- `radia_mcp.motor` — bearingless motor MCP (lab specialty)
"""


def get_overview_knowledge(topic: str = "decision_tree") -> str:
    """Dispatch WPT overview topics.

    Topics:
        landscape       - WPT regimes (inductive/resonant/capacitive/microwave)
        decision_tree   - Which topic for which design question (DEFAULT)
        history         - WPT history + Sugahara lab lineage
        all             - Everything
    """
    topic = topic.lower().strip()
    if topic in ("landscape", "overview", "regimes"):
        return LANDSCAPE
    if topic in ("decision_tree", "decision", "choose", "tree"):
        return DECISION_TREE
    if topic in ("history", "genealogy", "lineage"):
        return HISTORY
    if topic == "all":
        return "\n\n".join([LANDSCAPE, DECISION_TREE, HISTORY])
    return (f"Unknown topic '{topic}'. Available: landscape, decision_tree, "
            "history, all.")
