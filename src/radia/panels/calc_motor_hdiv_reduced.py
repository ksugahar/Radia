"""Stage-2 CLI for the planar HDiv-VIM reduced reluctance-motor path."""
from __future__ import annotations

import argparse
import math
import os
import sys

import numpy as np

_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)

from calc_common import calc_main, progress, setup_paths  # noqa: E402


def solve_motor_hdiv_reduced(
    vol_file: str,
    mu_r: float,
    h_amplitude: float,
    field_angle_deg: float,
    rotor_angle_start_deg: float,
    rotor_angle_stop_deg: float,
    rotor_angle_steps: int,
    maxwell_radius: float,
    stack_length: float,
    energy_delta_deg: float,
    circle_points: int,
    center_x: float,
    center_y: float,
    eta: float,
    order: int,
    msh_output: str = "",
):
    setup_paths()
    if not vol_file or not os.path.isfile(vol_file):
        raise FileNotFoundError(f"rotor .vol not found: {vol_file!r}")
    if rotor_angle_steps < 1:
        raise ValueError("rotor-angle-steps must be >= 1")

    import ngsolve as ng
    from radia.motor_hdiv import HDivReducedMotor

    mesh = ng.Mesh(vol_file)
    field_angle = math.radians(field_angle_deg)
    field_global = h_amplitude * np.array(
        [math.cos(field_angle), math.sin(field_angle)])
    rotor_angles = np.radians(np.linspace(
        rotor_angle_start_deg, rotor_angle_stop_deg, rotor_angle_steps))

    progress("MOTOR-HDIV", f"mesh={vol_file}, elements={mesh.ne}")
    with ng.TaskManager():
        motor = HDivReducedMotor(
            mesh, mu_r, stack_length=stack_length,
            center=(center_x, center_y), eta=eta, order=order)
        progress(
            "MOTOR-HDIV",
            f"Gram built once: ndof={motor.body.ndof}, charges={motor.body.n_charge}")
        result = motor.sweep(
            rotor_angles,
            field_global,
            maxwell_radius=maxwell_radius,
            circle_points=circle_points,
            energy_delta_angle=math.radians(energy_delta_deg),
        )
        gmsh_file = ""
        if msh_output:
            from radia.gmsh_post_export import GmshPostExport

            state = motor.solve_angle(float(rotor_angles[-1]), field_global)
            magnetization = ng.GridFunction(motor.body.fes)
            magnetization.vec.FV().NumPy()[:] = state.coefficients
            magnetization_3d = ng.CoefficientFunction(
                (magnetization[0], magnetization[1], 0.0)
            )
            gmsh_file = os.path.abspath(msh_output)
            os.makedirs(os.path.dirname(gmsh_file), exist_ok=True)
            post = GmshPostExport(mesh)
            post.add_vector_field(
                "M_local_A_per_m", magnetization_3d, cell_data=True
            )
            post.write(gmsh_file)
    result["vol_file"] = os.path.abspath(vol_file)
    result["field_angle_deg"] = float(field_angle_deg)
    result["gmsh_file"] = gmsh_file
    if gmsh_file:
        result["gmsh_field_rotor_angle_deg"] = float(rotor_angle_stop_deg)
    return result


def build_argparser():
    parser = argparse.ArgumentParser(
        description="Planar BDM1/BDM2 HDiv-VIM reduced reluctance-motor sweep")
    parser.add_argument("--vol", required=True, help="rotor-only 2D Netgen .vol mesh")
    parser.add_argument("--mu-r", type=float, default=1000.0,
                        help="linear rotor relative permeability")
    parser.add_argument("--H-amplitude", type=float, default=80000.0,
                        help="global stator-equivalent H amplitude [A/m]")
    parser.add_argument("--field-angle-deg", type=float, default=0.0,
                        help="global applied-H angle [deg]")
    parser.add_argument("--rotor-angle-start-deg", type=float, default=-45.0)
    parser.add_argument("--rotor-angle-stop-deg", type=float, default=45.0)
    parser.add_argument("--rotor-angle-steps", type=int, default=7)
    parser.add_argument("--maxwell-radius", type=float, default=0.05,
                        help="air-gap integration-circle radius [m]")
    parser.add_argument("--stack-length", type=float, default=0.05, help="axial length [m]")
    parser.add_argument("--energy-delta-deg", type=float, default=0.25,
                        help="coenergy central-difference half-step [deg]")
    parser.add_argument("--circle-points", type=int, default=1440)
    parser.add_argument("--center-x", type=float, default=0.0)
    parser.add_argument("--center-y", type=float, default=0.0)
    parser.add_argument("--eta", type=float, default=2.0,
                        help="charge-Gram admissibility parameter")
    parser.add_argument("--order", type=int, choices=(1, 2), default=1,
                        help="HDiv order: BDM1 or BDM2")
    parser.add_argument(
        "--msh-output", default="",
        help="GMSH .msh v4.1 output for the final rotor-angle magnetization",
    )
    return parser


def main():
    parser = build_argparser()

    def run(args):
        return solve_motor_hdiv_reduced(
            args.vol, args.mu_r, args.H_amplitude, args.field_angle_deg,
            args.rotor_angle_start_deg, args.rotor_angle_stop_deg,
            args.rotor_angle_steps, args.maxwell_radius, args.stack_length,
            args.energy_delta_deg, args.circle_points,
            args.center_x, args.center_y, args.eta, args.order,
            args.msh_output)

    calc_main(run, parser)


if __name__ == "__main__":
    main()
