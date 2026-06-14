"""coil_spec.py

JSON-serializable specification for coil designs used in the LLM-
driven optimization loop.

A CoilSpec describes a coil as plain data (profile + path + current)
independent of build123d / CoilBuilder. build_coil_from_spec() turns
it back into an executable CoilBuilder. This keeps the optimizer's
"proposer" side (which may be an LLM emitting JSON) decoupled from
the geometry-authoring side.

MVP supports:
  - profile_type='circle' with r_wire [m]  -> CircleProfile
  - profile_type='rect'   with wire_w, wire_h [m] -> RectProfile
  - path = single arc (R_major [m], arc_deg in (0, 360]) starting at
    (R_major, 0, 0) and going CCW in the XY plane

Extending to multi-segment paths, loft profiles, non-planar / helical
paths is a natural next step -- just add fields and corresponding
branches in build_coil_from_spec.
"""

import dataclasses
import json
from typing import Optional


@dataclasses.dataclass
class CoilSpec:
    profile_type: str = 'circle'           # 'circle' | 'rect'
    r_wire: Optional[float] = None         # [m] for 'circle'
    wire_w: Optional[float] = None         # [m] for 'rect'
    wire_h: Optional[float] = None         # [m] for 'rect'
    R_major: float = 30e-3                 # [m] arc radius
    arc_deg: float = 355.0                 # [deg] arc span, (0, 360]
    current: float = 1.0                   # [A]

    def to_dict(self) -> dict:
        d = dataclasses.asdict(self)
        return {k: v for k, v in d.items() if v is not None}

    def to_json(self, indent: Optional[int] = None) -> str:
        return json.dumps(self.to_dict(), indent=indent)

    @classmethod
    def from_dict(cls, d: dict) -> 'CoilSpec':
        known = {f.name for f in dataclasses.fields(cls)}
        return cls(**{k: v for k, v in d.items() if k in known})

    @classmethod
    def from_json(cls, s: str) -> 'CoilSpec':
        return cls.from_dict(json.loads(s))

    def validate(self):
        if self.profile_type == 'circle':
            if self.r_wire is None or self.r_wire <= 0:
                raise ValueError("circle profile needs r_wire > 0")
        elif self.profile_type == 'rect':
            if (self.wire_w is None or self.wire_w <= 0 or
                    self.wire_h is None or self.wire_h <= 0):
                raise ValueError("rect profile needs wire_w > 0 and wire_h > 0")
        else:
            raise ValueError(f"unknown profile_type: {self.profile_type!r}")
        if self.R_major <= 0:
            raise ValueError("R_major must be > 0")
        if not (0 < self.arc_deg <= 360):
            raise ValueError("arc_deg must be in (0, 360]")
        # Geometric sanity: inner radius > 0
        if self.profile_type == 'circle':
            if self.r_wire >= self.R_major:
                raise ValueError(
                    f"r_wire ({self.r_wire*1e3:.2f} mm) >= R_major "
                    f"({self.R_major*1e3:.2f} mm): arc inner radius <= 0")
        else:
            if self.wire_w / 2.0 >= self.R_major:
                raise ValueError(
                    f"wire_w/2 >= R_major: arc inner radius <= 0")


def build_coil_from_spec(spec: CoilSpec):
    """Build a CoilBuilder from a CoilSpec. Raises ValueError on invalid."""
    from radia.coil_builder import CoilBuilder
    from radia.coil_profile import CircleProfile, RectProfile

    spec.validate()

    if spec.profile_type == 'circle':
        profile = CircleProfile(r=spec.r_wire)
    else:
        profile = RectProfile(w=spec.wire_w, h=spec.wire_h)

    return (CoilBuilder(current=spec.current)
            .set_start([spec.R_major, 0.0, 0.0])
            .set_profile(profile)
            .add_arc(radius=spec.R_major, arc_angle=spec.arc_deg))
