"""
Translation of 02_sweep_revolve.jou to build123d (algebra mode).

Verdict: TRANSLATE.
Rationale: revolve + fillet are directly supported by build123d and covered
by the cubit_rosetta mapping (revolve surface X about z angle 360 ->
revolve(sketch, axis=Axis.Z, revolution_arc=360); modify curve ... blend radius R
-> fillet(edges, radius=R)).

Line-by-line correspondence (source .jou -> build123d):
  .jou L4   reset                                               -> (no-op)
  .jou L7   create surface rectangle width 5 height 10 zplane   -> Rectangle(5, 10) placed on Plane.XZ
  .jou L8   move surface 1 x 20 y 0 z 0                         -> Pos(20, 0, 0) * sketch
  .jou L11  revolve surface 1 about z angle 360                 -> revolve(sketch, axis=Axis.Z, revolution_arc=360)
  .jou L14  modify curve in volume 1 blend radius 1             -> fillet(ring.edges(), radius=1)
  .jou L16  rename volume 1 "ring_bobbin"                       -> ring_bobbin.label = "ring_bobbin"

Note on the sketch plane: Cubit's "zplane" means the plane normal to Z
(i.e. the XY plane) for surface creation; the subsequent move+revolve
about Z produces an axisymmetric ring of radius 20, width 5 (radial),
height 10 (along Z). For revolve about Z in build123d, the sketch must
live in a plane that contains the Z axis, i.e. Plane.XZ. We therefore
place the Rectangle on Plane.XZ with width=5 (X) and height=10 (Z), then
translate it by +20 in X so its inner edge is at r=17.5, outer at r=22.5.
(The Cubit script works the same way: "zplane" creates a horizontal
rectangle, but revolving it about Z while it sits offset from the axis
produces the ring; in build123d we express the same final geometry more
directly by sketching on Plane.XZ.)
"""

from build123d import (
    Axis,
    Plane,
    Pos,
    Rectangle,
    fillet,
    revolve,
)


def build() -> "Part":
    # Rectangular profile in XZ plane (.jou L7 + L8):
    #   width = 5 (radial, along X), height = 10 (axial, along Z),
    #   centered at x=20.
    profile = Pos(20, 0, 0) * (Plane.XZ * Rectangle(5, 10))

    # Full 360 deg revolve about Z (.jou L11).
    ring = revolve(profile, axis=Axis.Z, revolution_arc=360)

    # Fillet every edge with radius 1 (.jou L14). For a full 360 ring
    # the edge set is the top and bottom outer/inner circles (4 edges).
    ring_bobbin = fillet(ring.edges(), radius=1)

    # Label (.jou L16).
    ring_bobbin.label = "ring_bobbin"
    return ring_bobbin


if __name__ == "__main__":
    part = build()
    bb = part.bounding_box()
    print(f"label:       {part.label}")
    print(f"is_valid:    {part.is_valid}")
    print(f"volume:      {part.volume:.4f}")
    print(f"bbox.min:    ({bb.min.X:.3f}, {bb.min.Y:.3f}, {bb.min.Z:.3f})")
    print(f"bbox.max:    ({bb.max.X:.3f}, {bb.max.Y:.3f}, {bb.max.Z:.3f})")
    print(f"bbox.size:   ({bb.size.X:.3f}, {bb.size.Y:.3f}, {bb.size.Z:.3f})")
