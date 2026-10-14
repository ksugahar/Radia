"""Convert a mm-declared STEP file to a m-declared STEP file with
identical physical geometry (CARTESIAN_POINT values rescaled by 1/1000
so the M-declared output reads at the same physical size).

Why we need this: the radia PEEC-BEM parser (`filaments_from_step` in
calc_peec_bem.py) reads STEP coordinates as raw m-values per CLAUDE.md
"Unit System Policy".  A STEP declaring MILLIMETRE with (68.021, ...)
is misread as 68 m, scrambling helix detection (n_filaments collapses
to a degenerate fallback).  Run this tool once per CAD-source file to
get a m-declared sister STEP for the panel.

Implementation note: `Interface_Static.SetCVal_s("write.step.unit","M")`
ONLY changes the LENGTH_UNIT declaration in the output header — it does
NOT scale the geometry (verified empirically with OCP 7.7).  We therefore
apply an explicit ×0.001 transform via `BRepBuilderAPI_Transform` before
write.  The output is physically identical to the input.

Usage:
    python tools/step_mm_to_m.py input.stp output.step
"""
import sys
import os

if len(sys.argv) != 3:
    print("Usage: python step_mm_to_m.py input.step output.step")
    sys.exit(2)

src, dst = sys.argv[1], sys.argv[2]
if not os.path.isfile(src):
    print(f"FAIL: input not found: {src}")
    sys.exit(1)

from OCP.STEPControl import (STEPControl_Reader, STEPControl_Writer,
                              STEPControl_AsIs)
from OCP.IFSelect import IFSelect_RetDone
from OCP.Interface import Interface_Static
from OCP.gp import gp_Trsf
from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform

reader = STEPControl_Reader()
status = reader.ReadFile(src)
if status != IFSelect_RetDone:
    print(f"FAIL: STEPControl_Reader.ReadFile returned {status}")
    sys.exit(1)
reader.TransferRoots()
shape = reader.OneShape()

# OCP / OpenCascade keeps geometry internally in mm.  When we write a
# STEP with `write.step.unit = M`, the writer is supposed to scale by
# 1/1000 on output, but in practice the unit declaration is changed
# without a scale (verified empirically).  Apply explicit ×0.001 to
# the geometry so the output is physically correct in m.
trsf = gp_Trsf()
trsf.SetScaleFactor(0.001)
shape = BRepBuilderAPI_Transform(shape, trsf, True).Shape()

# Probe input unit declaration (informational only).
header_unit = "?"
try:
    with open(src, "r", errors="ignore") as f:
        head = f.read(8192)
    if "METRE" in head and "MILLIMETRE" not in head:
        header_unit = "METRE"
    elif "MILLIMETRE" in head:
        header_unit = "MILLIMETRE"
except Exception:
    pass
print(f"Input  STEP: {src}")
print(f"  declared LENGTH_UNIT: {header_unit}")

# Force writer to declare LENGTH_UNIT = M.  Per OCAF write convention.
Interface_Static.SetCVal_s("write.step.unit", "M")

writer = STEPControl_Writer()
status = writer.Transfer(shape, STEPControl_AsIs)
if status != IFSelect_RetDone:
    print(f"FAIL: STEPControl_Writer.Transfer returned {status}")
    sys.exit(1)
status = writer.Write(dst)
if status != IFSelect_RetDone:
    print(f"FAIL: STEPControl_Writer.Write returned {status}")
    sys.exit(1)

# Verify output declared unit.
out_unit = "?"
try:
    with open(dst, "r", errors="ignore") as f:
        head = f.read(8192)
    if "METRE" in head and "MILLIMETRE" not in head:
        out_unit = "METRE"
    elif "MILLIMETRE" in head:
        out_unit = "MILLIMETRE"
except Exception:
    pass

print(f"Output STEP: {dst}  ({os.path.getsize(dst)/(1024*1024):.2f} MB)")
print(f"  declared LENGTH_UNIT: {out_unit}")

# Cross-check: pick first CARTESIAN_POINT in input vs output.
def first_nonzero_pt(path):
    """Return first non-origin CARTESIAN_POINT in the file."""
    with open(path, "r", errors="ignore") as f:
        for line in f:
            if "CARTESIAN_POINT" not in line or "(" not in line:
                continue
            idx0 = line.find("(", line.find("CARTESIAN_POINT"))
            blob = line[idx0:idx0+120].strip()
            # skip the (0,0,0) origin convention point
            if "(0.,0.,0.)" in blob.replace(" ", ""):
                continue
            if "(0.0,0.0,0.0)" in blob.replace(" ", ""):
                continue
            return blob
    return "?"

# Re-detect declared unit by scanning the whole file (header may be
# anywhere in the first ~200 KB for medium-size STEP files).
def declared_unit(path):
    with open(path, "r", errors="ignore") as f:
        head = f.read(2_000_000)
    has_mm = "MILLIMETRE" in head or "MILLIMETER" in head
    has_m  = ("SI_UNIT($,.METRE.)" in head
              or ("LENGTH_UNIT" in head and "METRE" in head and not has_mm))
    if has_m and not has_mm:
        return "METRE"
    if has_mm:
        return "MILLIMETRE"
    return "?"

# Patch script-level header_unit + out_unit using the broader scanner.
header_unit = declared_unit(src)
out_unit = declared_unit(dst)
print(f"  input  declared LENGTH_UNIT: {header_unit}")
print(f"  output declared LENGTH_UNIT: {out_unit}")

print(f"  input  first non-zero pt: {first_nonzero_pt(src)}")
print(f"  output first non-zero pt: {first_nonzero_pt(dst)}")
