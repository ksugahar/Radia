# spice — SPICE / Verilog-A Export

PRIMA model order reduction with SPICE netlist and Verilog-A export.

## Files

| File | Description |
|------|-------------|
| `demo_prima_spice_export.py` | Export PRIMA reduced model as SPICE netlist or Verilog-A |
| `demo_peec_prima_reduction.py` | PEEC + PRIMA integration demo |
| `demo_veriloga_export.py` | Verilog-A export demo |
| `demo_dowell_spice.py` | Dowell skin-effect SPICE model |
| `prima_with_dowell_correction.py` | PRIMA(DC) + Dowell correction verification |
| `dowell_to_prima.py` | Dowell formula to PRIMA ladder conversion |

## SPICE Models

| File | Description |
|------|-------------|
| `dowell_skin.sp` | Dowell skin-effect SPICE netlist |
| `skin_effect.sp` | Skin effect model |
| `wire_full.sp` | Full wire model |
| `wire_prima.sp` | PRIMA reduced wire model |
| `wire_prima_skin.sp` | PRIMA wire with skin effect |

## Verilog-A Models

| File | Description |
|------|-------------|
| `demo_cole_cole_cap.va` | Cole-Cole capacitor model |
| `demo_debye_cap.va` | Debye capacitor model |
| `demo_dowell_skin.va` | Dowell skin-effect model |
| `demo_multi_debye.va` | Multi-Debye relaxation model |
| `demo_peec_segment.va` | PEEC segment model |
| `wire_prima.va` | PRIMA wire model |

## Usage

```bash
# SPICE netlist export
python demo_prima_spice_export.py

# Verilog-A export
python demo_prima_spice_export.py --verilog-a

# Set Lanczos order
python demo_prima_spice_export.py --lanczos 10
```
