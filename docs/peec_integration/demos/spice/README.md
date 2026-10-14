# spice — SPICE / Verilog-A Export

PRIMA model order reduction with SPICE netlist and Verilog-A export.

## Files

| File | Description |
|------|-------------|
| `demo_veriloga_export.py` | Verilog-A export demo |
| `demo_dowell_spice.py` | Dowell skin-effect SPICE model |
| `prima_with_dowell_correction.py` | PRIMA(DC) + Dowell correction verification |
| `dowell_to_prima.py` | Trial of an empirical RL ladder circuit for the Dowell formula (element values from DC values and an ad hoc scaling; neither a PRIMA projection nor a derived continued fraction; the file name is historical) |

Numerical PRIMA/Dowell evidence is maintained under
`validation_test/peec_integration/verification/`.

## SPICE Models

| File | Description |
|------|-------------|
| `dowell_skin.sp` | Dowell skin-effect SPICE netlist |
| `skin_effect.sp` | Skin effect model |
| `wire_full.sp` | Full wire model |

## Verilog-A Models

| File | Description |
|------|-------------|
| `demo_cole_cole_cap.va` | Cole-Cole capacitor model |
| `demo_debye_cap.va` | Debye capacitor model |
| `demo_dowell_skin.va` | Dowell skin-effect model |
| `demo_multi_debye.va` | Multi-Debye relaxation model |
| `demo_peec_segment.va` | PEEC segment model |

## Usage

```bash
python demo_veriloga_export.py
```
