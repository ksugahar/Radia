"""
demo_box.py — minimal end-to-end demo for the
build123d -> Netgen -> Gmsh pipeline.

Produces (in ./runs/):
    demo_box.brep            build123d export
    demo_box.msh             Netgen mesh as Gmsh v2.2 .msh
    demo_box_post.msh        same mesh + a dummy scalar view attached
    demo_box.json            machine-readable run record

Run:
    python demo_box.py
"""

from build123d import Box, Cylinder
from _pipeline import run_pipeline, save_record


def main():
    # Plate 100 x 100 x 20 with a through-hole (r = 15).  All units are
    # the build123d default (mm).  Netgen maxh is in the same units.
    plate = Box(100, 100, 20) - Cylinder(radius=15, height=20)
    plate.label = "plate_with_hole"

    record = run_pipeline(plate, out_dir="./runs", label="demo_box", maxh=8.0)

    print(f"[demo_box] status = {record['status']}")
    for stage in record["stages"]:
        print(f"  {stage['stage']:4s} ok={stage.get('ok')}")
    save_record(record, "./runs")


if __name__ == "__main__":
    main()
