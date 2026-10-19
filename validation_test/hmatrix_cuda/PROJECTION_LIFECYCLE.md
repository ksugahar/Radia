# Headless projection ownership

The standalone projection bridge creates NGSolve GridFunctions with the
explicit `novisual` flag. It uses NGSolve's FE and interpolation operations
unchanged and does not register a Netgen visualization object.

In the pinned NGSolve v6.2.2607 implementation, GridFunction's `visual` flag
defaults to true. Its concrete constructor calls `Visualize`, which creates
`VisualizeCoefficientFunction` and registers it with `Ng_SetSolutionData`.
That visualization object strongly owns MeshAccess. Its GridFunction alias
does not keep the original GridFunction owner alive. This accounts for the
observed surviving mesh with expired FE, GridFunction and target owners.
See the upstream [GridFunction implementation](https://github.com/NGSolve/ngsolve/blob/v6.2.2607/comp/gridfunction.cpp)
and [visualization owner declarations](https://github.com/NGSolve/ngsolve/blob/v6.2.2607/comp/gridfunction.hpp).

The diagnostic build compares mesh-only, the actual Targets CF, FE
construction/update, GridFunction construction/update and interpolation.
`checkProjectionOwnerPhases` retains default-visual construction as a
negative control and requires no surviving mesh in the corresponding
`novisual` cases. These controls deliberately register visualization data;
they are not the baseline for a subsequent memory gate.

Run the production create/project/destroy memory gate in a fresh MATLAB
process with the shared initialized-device protocol, unchanged limits and
independent CPU coefficient checks. Weak-owner results alone do not close
that gate. No global Netgen registry is cleared: it may contain another
consumer's visualization data. Long-lived context reuse remains useful but
does not waive correct destruction or justify an unbounded creation loop.
