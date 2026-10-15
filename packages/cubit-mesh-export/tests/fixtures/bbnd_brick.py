"""Licensed regeneration input: an independently authored SI unit cube.

Run with the headless Cubit console and the candidate CUBIT_PLUGIN_DIR.
BBND_OUTPUT chooses the output path; BBND_CASE optionally selects a negative case.
"""
import os
import cubit

for command in (
    "reset", "create brick x 1 y 1 z 1", "volume all scheme tetmesh",
    "volume all size 0.5", "mesh volume all", "block 1 add tet all",
    'block 1 name "domain"', "sideset 1 add surface all",
    'sideset 1 name "boundary"',
):
    cubit.cmd(command)

# Select the vertical edge by its geometric centre, independent of Cubit IDs.
curves = [cid for cid in cubit.parse_cubit_list("curve", "all")
          if all(abs(a-b) < 1e-10 for a, b in zip(
              cubit.get_center_point("curve", cid), (-0.5, -0.5, 0.0)))]
assert len(curves) == 1, curves
curve = curves[0]
cubit.cmd(f"sideset 100 add curve {curve}")
cubit.cmd('sideset 100 name "test_edge"')
case = os.environ.get("BBND_CASE", "valid")
if case == "conflict":
    cubit.cmd(f"sideset 101 add curve {curve}")
    cubit.cmd('sideset 101 name "other_edge"')
elif case == "unrepresented":
    cubit.cmd("create vertex 2 0 0")
    first = max(cubit.parse_cubit_list("vertex", "all"))
    cubit.cmd("create vertex 3 0 0")
    second = max(cubit.parse_cubit_list("vertex", "all"))
    cubit.cmd(f"create curve vertex {first} {second}")
    extra = max(cubit.parse_cubit_list("curve", "all"))
    cubit.cmd(f"sideset 101 add curve {extra}")
    cubit.cmd('sideset 101 name "unrepresented_edge"')
elif case != "valid":
    raise ValueError(case)
output = os.environ["BBND_OUTPUT"].replace("\\", "/")
order = int(os.environ.get("BBND_ORDER", "2"))
cubit.cmd(f'export netgen "{output}" order {order} overwrite')
