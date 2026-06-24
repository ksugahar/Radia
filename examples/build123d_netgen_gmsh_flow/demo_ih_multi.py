"""
demo_ih_multi.py — multi-region demo: workpiece + coil + air.

Demonstrates the run_pipeline_multi helper:
    - build 3 labeled build123d Parts
    - hand them as (part, name) tuples in strict order
    - pipeline emits a .step (Compound) + .msh (named physical groups)
      + a per-element 'region_id' view ready for colored rendering

Expected post-stage output in runs/ih_multi*_post.msh:
    physical groups:
        (dim=3, tag=100001) name='workpiece'
        (dim=3, tag=100002) name='coil'
        (dim=3, tag=100003) name='air'
    view 'ih_multi_region_id' = 1.0 / 2.0 / 3.0 per element

Run:
    python demo_ih_multi.py
"""

from build123d import Cylinder, Pos
from _pipeline import run_pipeline_multi, save_record


def main():
    # Workpiece: solid steel cylinder (r=15, h=20) centered at origin
    workpiece = Cylinder(radius=15.0, height=20.0)
    workpiece.label = "workpiece"

    # Coil: ring (r_in=22, r_out=28, h=10) sitting above the workpiece
    coil = Pos(0, 0, 25.0) * (
        Cylinder(radius=28.0, height=10.0)
        - Cylinder(radius=22.0, height=10.0)
    )
    coil.label = "coil"

    # Air domain: large cylinder minus the two solid regions
    air = Cylinder(radius=60.0, height=80.0) - workpiece - coil
    air.label = "air"

    regions = [
        (workpiece, "workpiece"),
        (coil,      "coil"),
        (air,       "air"),
    ]

    rec = run_pipeline_multi(regions, out_dir="./runs",
                             label="ih_multi", maxh=6.0)

    print(f"[ih_multi] status = {rec['status']}")
    for s in rec["stages"]:
        if s["stage"] == "cad":
            print("  cad:  ", [r["name"] for r in s["regions"]])
        elif s["stage"] == "mesh":
            print(f"  mesh: materials={s['materials']}  "
                  f"nv={s['nv']}  ne={s['ne']}")
        elif s["stage"] == "post":
            print(f"  post: n_regions={s['n_regions']}")
            for r in s["regions"]:
                print(f"    [{r['index']}] {r['name']:12s} "
                      f"phys_tag={r['physical_tag']}  "
                      f"elements={r['n_elem']}")
    save_record(rec, "./runs")


if __name__ == "__main__":
    main()
