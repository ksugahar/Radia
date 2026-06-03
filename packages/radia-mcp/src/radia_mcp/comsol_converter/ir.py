"""Neutral intermediate representation (IR) for comsol-converter.

COMSOL <-> IR <-> NGSolve. The IR is solver-agnostic: each side maps only
to/from this hub, avoiding brittle N x M direct pairings. v1 scope =
low-frequency magnetics (magnetostatic A-formulation, B = curl A), which
covers the lab-core TEAM problems (6 linear / 13 nonlinear / 20 force).

Design notes
------------
* The IR carries *physics intent*, not solver syntax: a domain has a
  material (linear mu_r / nonlinear B-H / permanent magnet), a source is a
  coil (N, I) or a current density, a boundary condition is named by its
  physical meaning (magnetic insulation, PMC, symmetry, infinity).
* ``mesh.symmetry_fraction`` is first-class on purpose: the single biggest
  COMSOL trip-wire (interop comsol_lab_tips 'coil_current_scaling') is that a
  1/8 model needs the post-processed force/flux scaled by 8. Translators must
  carry this explicitly so the NGSolve and COMSOL sides agree.
* ``provenance`` records where the IR came from (round-trip auditing).
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any, Optional

IR_VERSION = "0.1"

# Formulations the IR currently models (extensible; v1 = magnetostatics).
FORMULATIONS = ("magnetostatic_A",)

# Boundary-condition vocabulary (physical meaning, solver-agnostic).
BC_KINDS = (
    "magnetic_insulation",        # n x A = 0  (COMSOL Magnetic Insulation; NGSolve Dirichlet on A_t)
    "perfect_magnetic_conductor", # n x H = 0  (COMSOL PMC; NGSolve natural BC)
    "symmetry_normalB0",          # B.n = 0 symmetry plane
    "symmetry_tangentialH0",      # H x n = 0 symmetry plane
    "infinity",                   # far-field decay (COMSOL infinite element / NGSolve Kelvin)
    "continuity",                 # interior interface
)
MATERIAL_KINDS = ("air", "linear", "nonlinear_BH", "permanent_magnet")
SOURCE_KINDS = ("coil", "current_density")
OUTPUT_KINDS = ("B_at_point", "flux", "force_maxwell", "energy", "torque")


@dataclass
class Material:
    name: str
    domain: str                                  # geometry region this applies to
    kind: str = "air"                            # see MATERIAL_KINDS
    mu_r: float = 1.0                            # linear relative permeability
    bh_curve: Optional[str] = None               # name/ref of a B-H table (nonlinear iron)
    conductivity: float = 0.0                    # S/m (0 = non-conducting magnetostatic)
    br: float = 0.0                              # PM remanence [T]
    magnetization_dir: Optional[list[float]] = None  # PM easy axis (unit vector)


@dataclass
class Source:
    name: str
    kind: str                                    # "coil" | "current_density"
    domain: str
    n_turns: float = 1.0
    current: float = 0.0                         # per-turn current [A] (total NI = n_turns*current)
    j_density: Optional[list[float]] = None      # [Jx,Jy,Jz] A/m^2 for current_density
    direction: Optional[list[float]] = None      # coil winding direction (unit vector)


@dataclass
class BoundaryCondition:
    boundary: str
    kind: str                                    # see BC_KINDS


@dataclass
class Output:
    name: str
    kind: str                                    # see OUTPUT_KINDS
    target: Optional[str] = None                 # body / boundary / point identifier
    point: Optional[list[float]] = None          # [x,y,z] for B_at_point


@dataclass
class Study:
    kind: str = "stationary"                     # "stationary" | "frequency"
    frequency: float = 0.0
    nonlinear: bool = False


@dataclass
class Mesh:
    element_order: int = 1
    max_size: Optional[float] = None
    symmetry_fraction: float = 1.0               # 0.125 = 1/8 model -> scale outputs by 8 (lab tip)


@dataclass
class ModelIR:
    name: str
    formulation: str = "magnetostatic_A"
    length_unit: str = "m"
    geometry: dict[str, Any] = field(default_factory=dict)   # OCC/Netgen spec or source-native ref
    materials: list[Material] = field(default_factory=list)
    sources: list[Source] = field(default_factory=list)
    boundary_conditions: list[BoundaryCondition] = field(default_factory=list)
    mesh: Mesh = field(default_factory=Mesh)
    study: Study = field(default_factory=Study)
    outputs: list[Output] = field(default_factory=list)
    provenance: dict[str, Any] = field(default_factory=dict)  # {source_tool, source_file, notes}

    def to_dict(self) -> dict:
        return asdict(self)

    @staticmethod
    def from_dict(d: dict) -> "ModelIR":
        return ModelIR(
            name=d["name"],
            formulation=d.get("formulation", "magnetostatic_A"),
            length_unit=d.get("length_unit", "m"),
            geometry=d.get("geometry", {}),
            materials=[Material(**m) for m in d.get("materials", [])],
            sources=[Source(**s) for s in d.get("sources", [])],
            boundary_conditions=[BoundaryCondition(**b) for b in d.get("boundary_conditions", [])],
            mesh=Mesh(**d["mesh"]) if d.get("mesh") else Mesh(),
            study=Study(**d["study"]) if d.get("study") else Study(),
            outputs=[Output(**o) for o in d.get("outputs", [])],
            provenance=d.get("provenance", {}),
        )


def team20_example() -> ModelIR:
    """A representative IR for TEAM problem 20 (3-D static force).

    Coil between a steel pole and yoke, DC excitation, nonlinear steel;
    quantity of interest = lifting force on the pole (Maxwell stress).
    Geometry is left as a named ref (the geometry layer is built per-side);
    this example pins down the physics/material/source/BC/output intent.
    """
    return ModelIR(
        name="TEAM20_static_force",
        formulation="magnetostatic_A",
        length_unit="mm",
        geometry={"ref": "team20", "note": "center pole + yoke + coil + air box"},
        materials=[
            Material(name="air", domain="air", kind="air", mu_r=1.0),
            Material(name="pole", domain="pole", kind="nonlinear_BH",
                     bh_curve="team20_steel"),
            Material(name="yoke", domain="yoke", kind="nonlinear_BH",
                     bh_curve="team20_steel"),
        ],
        sources=[Source(name="coil", kind="coil", domain="coil",
                        n_turns=1000.0, current=3.0)],  # 3000 AT
        boundary_conditions=[
            BoundaryCondition(boundary="outer", kind="magnetic_insulation"),
            BoundaryCondition(boundary="sym_x", kind="symmetry_normalB0"),
            BoundaryCondition(boundary="sym_y", kind="symmetry_normalB0"),
        ],
        mesh=Mesh(element_order=2, symmetry_fraction=0.25),  # 1/4 model
        study=Study(kind="stationary", nonlinear=True),
        outputs=[
            Output(name="lifting_force", kind="force_maxwell", target="pole"),
            Output(name="B_center", kind="B_at_point", point=[0.0, 0.0, 0.0]),
        ],
        provenance={"source_tool": "example", "team_problem": 20},
    )
