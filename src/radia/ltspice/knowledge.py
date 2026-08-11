"""Small circuit-knowledge helpers for first-pass SPICE design.

These helpers intentionally contain general engineering rules of thumb, not
private corpus excerpts or lab-internal references.  They are meant to produce
reasonable simulation seeds that users must still validate in SPICE.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any


@dataclass(frozen=True)
class BuckSeed:
    """First-pass asynchronous buck-converter sizing result."""

    vin_v: float
    vout_v: float
    iout_a: float
    fsw_hz: float
    duty: float
    period_s: float
    ton_s: float
    load_ohm: float
    inductance_h: float
    ripple_current_a: float
    capacitance_f: float
    ripple_voltage_target_v: float
    switch_ron_ohm: float
    diode_drop_v: float
    inductor_dcr_ohm: float
    capacitor_esr_ohm: float

    def to_dict(self) -> dict[str, float]:
        return asdict(self)

    def to_netlist(self) -> str:
        """Return an LTspice-ready open-loop transient netlist."""
        return "\n".join([
            "* Asynchronous PWM buck seed",
            f"VIN vin 0 DC {self.vin_v:.6g}",
            f"VGATE gate 0 PULSE(0 10 0 20n 20n {_eng(self.ton_s, 's')} {_eng(self.period_s, 's')})",
            "SMAIN vin sw gate 0 SW_MAIN",
            f"L1 sw lx {_eng(self.inductance_h, 'H')}",
            f"RL1 lx out {_eng(self.inductor_dcr_ohm, 'ohm')}",
            f"RCOUT out cout {_eng(self.capacitor_esr_ohm, 'ohm')}",
            f"COUT cout 0 {_eng(self.capacitance_f, 'F')}",
            f"RLOAD out 0 {_eng(self.load_ohm, 'ohm')}",
            "D1 0 sw DFREE",
            f".model SW_MAIN SW(Ron={_eng(self.switch_ron_ohm, 'ohm')} Roff=10Meg Vt=4 Vh=0.5)",
            ".model DFREE D(Is=20n Rs=35m N=1.2 Cjo=120p Tt=30n)",
            f".tran 0 {_eng(max(500.0 * self.period_s, 5e-3), 's')} 0 {_eng(self.period_s / 80.0, 's')} startup",
            ".save V(vin) V(sw) V(out) V(cout) I(L1) I(VIN)",
            ".end",
        ])


def buck_seed(
    vin_v: float,
    vout_v: float,
    iout_a: float,
    fsw_hz: float = 100_000.0,
    ripple_fraction: float = 0.25,
    ripple_voltage_fraction: float = 0.01,
    switch_ron_ohm: float = 0.05,
    diode_drop_v: float = 0.55,
) -> BuckSeed:
    """Size a first-pass open-loop asynchronous buck-converter seed.

    The duty estimate uses average inductor voltage with conduction-loss and
    diode-drop terms.  The inductor is sized from a target peak-to-peak ripple
    fraction, and the capacitor is sized from the triangular ripple estimate.
    """
    vin_v = float(vin_v)
    vout_v = float(vout_v)
    iout_a = max(float(iout_a), 1e-6)
    fsw_hz = max(float(fsw_hz), 1.0)
    ripple_fraction = min(max(float(ripple_fraction), 0.05), 0.8)
    ripple_voltage_fraction = min(max(float(ripple_voltage_fraction), 0.001), 0.2)

    load_ohm = max(vout_v / iout_a, 1e-3)
    period_s = 1.0 / fsw_hz
    inductor_dcr_ohm = max(0.02, min(0.25, 0.08 * (1.0 / max(iout_a, 0.2))))
    duty = (
        (vout_v + iout_a * inductor_dcr_ohm + diode_drop_v)
        / max(vin_v - iout_a * switch_ron_ohm + diode_drop_v, 1e-6)
    )
    duty = min(max(duty, 0.03), 0.9)
    ripple_current_a = max(iout_a * ripple_fraction, 0.05)
    inductance_h = max((vin_v - vout_v) * duty / (ripple_current_a * fsw_hz), 1e-9)
    ripple_voltage_target_v = max(vout_v * ripple_voltage_fraction, 0.01)
    capacitance_f = max(ripple_current_a / (8.0 * fsw_hz * ripple_voltage_target_v), 1e-9)
    capacitor_esr_ohm = max(0.01, min(0.15, ripple_voltage_target_v / max(ripple_current_a, 1e-6)))

    return BuckSeed(
        vin_v=vin_v,
        vout_v=vout_v,
        iout_a=iout_a,
        fsw_hz=fsw_hz,
        duty=duty,
        period_s=period_s,
        ton_s=duty * period_s,
        load_ohm=load_ohm,
        inductance_h=inductance_h,
        ripple_current_a=ripple_current_a,
        capacitance_f=capacitance_f,
        ripple_voltage_target_v=ripple_voltage_target_v,
        switch_ron_ohm=switch_ron_ohm,
        diode_drop_v=diode_drop_v,
        inductor_dcr_ohm=inductor_dcr_ohm,
        capacitor_esr_ohm=capacitor_esr_ohm,
    )


def circuit_knowledge(topic: str = "") -> dict[str, Any]:
    """Return compact public circuit-design rules by topic."""
    topic_l = topic.lower()
    rules = {
        "buck": [
            "Start from volt-second balance, then include diode drop and conduction losses before simulation.",
            "Choose inductor ripple deliberately; 20-40% of load current is a useful simulation seed range.",
            "Judge the settled window separately from startup; average over the whole transient can mislead.",
            "Do not sign off an open-loop PWM seed. Add feedback, current limit, thermal checks, and layout parasitics.",
        ],
        "switching": [
            "Save switch-node voltage, inductor current, input current, and output ripple in every transient run.",
            "Sweep duty ratio and parasitics; small losses move the operating point in low-voltage outputs.",
            "Check startup overshoot separately from steady-state ripple.",
        ],
        "conversion": [
            "Round-trip checks need topology comparison, not only component count.",
            "Controlled switches have four electrical pins plus a model reference.",
            "Subcircuits and device models live in different namespaces.",
        ],
        "opamp": [
            "A first ideal op-amp seed is useful for topology, but replace it with a finite-bandwidth model before judging stability.",
            "Always include source impedance, load impedance, and supply rails in the simulation seed.",
        ],
    }
    if "buck" in topic_l or "降圧" in topic_l:
        selected = rules["buck"] + rules["switching"]
    elif "switch" in topic_l or "pwm" in topic_l or "power" in topic_l:
        selected = rules["switching"]
    elif "convert" in topic_l or "asc" in topic_l or "netlist" in topic_l:
        selected = rules["conversion"]
    elif "opamp" in topic_l or "op-amp" in topic_l:
        selected = rules["opamp"]
    else:
        selected = rules["buck"][:2] + rules["conversion"][:2]
    return {"topic": topic, "rules": selected}


def _eng(value: float, unit: str) -> str:
    suffixes = [
        (1e9, "G"),
        (1e6, "Meg"),
        (1e3, "k"),
        (1.0, ""),
        (1e-3, "m"),
        (1e-6, "u"),
        (1e-9, "n"),
        (1e-12, "p"),
    ]
    abs_value = abs(value)
    for scale, suffix in suffixes:
        if abs_value >= scale or scale == 1e-12:
            number = value / scale
            if unit in {"ohm", "F", "H", "Hz", "s"}:
                return f"{number:.4g}{suffix}"
            return f"{number:.4g}{suffix}{unit}"
    return f"{value:.4g}"
