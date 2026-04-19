"""
Translation of 01_pure_cad.jou to build123d (algebra mode).

Verdict: TRANSLATE.
Rationale: purely CAD-subset .jou (primitives + boolean + move + rename),
all verbs are covered by the cubit_rosetta / build123d_crossref mapping tables.

Line-by-line correspondence (source .jou -> build123d):
  .jou L3  reset                                        -> (no-op; fresh Python module)
  .jou L6  create brick x 100 y 50 z 20                 -> Box(100, 50, 20)
  .jou L9  create cylinder radius 10 height 20          -> Cylinder(radius=10, height=20)
  .jou L10 move volume 2 x 30 y 0 z 0                   -> Pos(30, 0, 0) * hole
  .jou L11 subtract volume 2 from volume 1              -> plate - hole_r
  .jou L13 create cylinder radius 10 height 20          -> Cylinder(radius=10, height=20)
  .jou L14 move volume 3 x -30 y 0 z 0                  -> Pos(-30, 0, 0) * hole
  .jou L15 subtract volume 3 from volume 1              -> plate - hole_l
  .jou L18 rename volume 1 "bracket"                    -> bracket.label = "bracket"
"""

from build123d import Box, Cylinder, Pos


def build() -> "Part":
    # Base plate (.jou L6)
    plate = Box(100, 50, 20)

    # Right mounting hole (.jou L9-L11)
    hole_r = Pos(30, 0, 0) * Cylinder(radius=10, height=20)
    plate = plate - hole_r

    # Left mounting hole (.jou L13-L15)
    hole_l = Pos(-30, 0, 0) * Cylinder(radius=10, height=20)
    bracket = plate - hole_l

    # Label (.jou L18)
    bracket.label = "bracket"
    return bracket


if __name__ == "__main__":
    part = build()
    bb = part.bounding_box()
    print(f"label:       {part.label}")
    print(f"is_valid:    {part.is_valid}")
    print(f"volume:      {part.volume:.4f}")
    print(f"bbox.min:    ({bb.min.X:.3f}, {bb.min.Y:.3f}, {bb.min.Z:.3f})")
    print(f"bbox.max:    ({bb.max.X:.3f}, {bb.max.Y:.3f}, {bb.max.Z:.3f})")
    print(f"bbox.size:   ({bb.size.X:.3f}, {bb.size.Y:.3f}, {bb.size.Z:.3f})")
