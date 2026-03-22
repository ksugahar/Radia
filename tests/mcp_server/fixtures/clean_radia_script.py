"""Clean Radia script -- should produce zero lint findings.

Used by --selftest to verify no false positives.
"""
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))
import radia as rad


def main():
    rad.UtiDelAll()

    # Create magnet (meters, A/m)
    pm = rad.ObjRecMag([0, 0, 0], [0.01, 0.01, 0.01], [0, 0, 954930])

    # Background field with callable
    bkg = rad.ObjBckg(lambda p: [0, 0, 0.1])

    # Compute field
    B = rad.Fld(pm, "b", [0.02, 0, 0])
    print(f"B = {B}")

    rad.UtiDelAll()


if __name__ == "__main__":
    main()
