"""
Generate coil STEP file from parameters (subprocess, no Cubit needed).

Called from Cubit GUI via subprocess:
    python generate_coil.py --type racetrack --straight 0.1 --radius 0.02 --output coil.step

Or with a coil script:
    python generate_coil.py --script my_coil.py --output coil.step

IMPORTANT: NGSolve/OCC must be available (Python 3.12, not Cubit's 3.10).
Outputs JSON to stdout (last line) with coil info.
"""

import argparse
import json
import math
import os
import sys


def build_racetrack(straight, radius, width=None, height=None, current=1.0,
                    center_z=0.0, center_y=0.0):
    """Build a racetrack coil (2 straights + 2 semicircular arcs).

    Args:
        straight: Straight section length [m]
        radius: Arc radius [m]
        width: Cross-section width [m] (optional, for solid display)
        height: Cross-section height [m] (optional)
        current: Current [A]
        center_z: Z offset of coil center [m]
        center_y: Y offset of coil center [m]

    Returns:
        CoilBuilder instance
    """
    _this_dir = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.dirname(_this_dir))

    from radia.coil_builder import CoilBuilder

    cb = CoilBuilder(current=current)
    cb.set_start([0, center_y - straight / 2, center_z])

    # Cross-section: default to radius/5 square if not specified
    if not width:
        width = radius / 5
    if not height:
        height = radius / 5
    cb.set_cross_section(width=width, height=height)

    cb.add_straight(length=straight, tilt=0)
    cb.add_arc(radius=radius, arc_angle=180, tilt=0)
    cb.add_straight(length=straight, tilt=0)
    cb.add_arc(radius=radius, arc_angle=180, tilt=0)

    return cb


def build_from_script(script_path):
    """Build coil from a Python script that defines build_coil().

    The script must define:
        def build_coil():
            cb = CoilBuilder(current=...)
            cb.set_start(...)
            ...
            return cb

    Args:
        script_path: Path to .py script

    Returns:
        CoilBuilder instance
    """
    import importlib.util
    spec = importlib.util.spec_from_file_location("coil_script", script_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    if not hasattr(mod, "build_coil"):
        raise ValueError(f"{script_path} must define build_coil() function")

    return mod.build_coil()


def main():
    parser = argparse.ArgumentParser(
        description="Generate coil STEP file for Cubit visualization")
    parser.add_argument("--type", default="racetrack",
                        choices=["racetrack"],
                        help="Coil type (default: racetrack)")
    parser.add_argument("--script", help="Custom coil script (.py with build_coil())")
    parser.add_argument("--straight", type=float, default=0.1,
                        help="Straight section length [m]")
    parser.add_argument("--radius", type=float, default=0.02,
                        help="Arc radius [m]")
    parser.add_argument("--width", type=float, default=0,
                        help="Cross-section width [m] (0=wire)")
    parser.add_argument("--height", type=float, default=0,
                        help="Cross-section height [m] (0=wire)")
    parser.add_argument("--current", type=float, default=1.0,
                        help="Current [A]")
    parser.add_argument("--center-y", type=float, default=0.0,
                        help="Coil center Y offset [m]")
    parser.add_argument("--center-z", type=float, default=0.0,
                        help="Coil center Z offset [m]")
    parser.add_argument("--output", "-o", default="coil.step",
                        help="Output STEP file path")
    args = parser.parse_args()

    if args.script:
        cb = build_from_script(args.script)
    else:
        cb = build_racetrack(
            straight=args.straight, radius=args.radius,
            width=args.width if args.width > 0 else None,
            height=args.height if args.height > 0 else None,
            current=args.current,
            center_y=args.center_y, center_z=args.center_z)

    cb.write_step(args.output)

    # Wire segments for Biot-Savart
    segs, current = cb.to_wire_segments()

    result = {
        "output": args.output,
        "n_segments": len(cb.segments),
        "n_wire_segments": len(segs),
        "current": current,
    }

    print(json.dumps(result))


if __name__ == "__main__":
    main()
