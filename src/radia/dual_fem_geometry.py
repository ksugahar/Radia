"""Mapping checks for the tetrahedral Kelvin benchmark meshes."""

from itertools import product


def audit_curved_kelvin_mesh(mesh, offset, radius, *, deformation_order=None):
    """Check sampled tetrahedron orientation and the entire polynomial seam.

    Compare each mapped determinant with its own affine orientation. Match
    triangle vertices through the periodic identifications and compare their
    maps at a degree-k unisolvent barycentric lattice, including the edges.
    """
    import ngsolve as ng
    import numpy as np
    from radia.kelvin_identify_ngsolve import _identified_point_pairs

    order = deformation_order if deformation_order is not None else mesh.GetCurveOrder()
    tetra_lattice = [
        (i / order, j / order, k / order)
        for i, j, k in product(range(order + 1), repeat=3)
        if i + j + k <= order
    ]
    reference_rule = ng.IntegrationRule(tetra_lattice, [1.] * len(tetra_lattice))
    sample_order = 2 * order + 2
    volume_rule = ng.IntegrationRule(ng.ET.TET, sample_order)
    jacobian = ng.specialcf.JacobianMatrix(3)
    failed = {}
    minimum = float("inf")
    for element in mesh.Elements(ng.VOL):
        if element.type != ng.ET.TET:
            raise ValueError("Kelvin benchmark requires tetrahedra")
        vertices = np.array([tuple(mesh[v].point) for v in element.vertices])
        orientation = np.sign(np.linalg.det(vertices[:3] - vertices[3]))
        trafo = mesh.GetTrafo(element)
        samples = np.concatenate([
            np.asarray(jacobian(trafo(rule))).reshape(-1, 3, 3)
            for rule in (reference_rule, volume_rule)
        ])
        with np.errstate(divide="ignore", invalid="ignore"):
            scaled = orientation * np.linalg.det(samples) / np.prod(
                np.linalg.norm(samples, axis=1), axis=1)
        if not np.all(np.isfinite(scaled)) or np.min(scaled) <= 0:
            failed[element.mat] = failed.get(element.mat, 0) + 1
        if np.all(np.isfinite(scaled)):
            minimum = min(minimum, float(np.min(scaled)))

    inner, outer = ([
        face for face in mesh.Elements(ng.BND) if face.mat == label
    ] for label in ("kelvin_int", "kelvin_ext"))
    if not inner or not outer or any(face.type != ng.ET.TRIG for face in inner + outer):
        raise ValueError("Kelvin seam requires two triangulated boundaries")
    inner_vertices = {v.nr for face in inner for v in face.vertices}
    outer_vertices = {v.nr for face in outer for v in face.vertices}
    identification = {}
    for first, second in _identified_point_pairs(mesh.ngmesh.GetIdentifications()):
        for source, target in ((first - 1, second - 1), (second - 1, first - 1)):
            if source in inner_vertices and target in outer_vertices:
                if source in identification and identification[source] != target:
                    raise ValueError("Ambiguous periodic Kelvin vertex")
                identification[source] = target
    if identification.keys() != inner_vertices or set(identification.values()) != outer_vertices:
        raise ValueError("Incomplete periodic Kelvin vertex pairing")
    outer_faces = {frozenset(v.nr for v in face.vertices): face for face in outer}
    if len(inner) != len(outer_faces) or len(outer) != len(outer_faces):
        raise ValueError("Different Kelvin seam face counts")

    bary = np.array([
        (i, j, order - i - j)
        for i in range(order + 1) for j in range(order + 1 - i)
    ], dtype=float) / order
    weights = [1.] * len(bary)
    rule = ng.IntegrationRule(bary[:, :2].tolist(), weights)
    xyz = ng.CoefficientFunction((ng.x, ng.y, ng.z))
    max_seam_error = 0.0
    for inner_face in inner:
        paired = [identification[v.nr] for v in inner_face.vertices]
        outer_face = outer_faces.pop(frozenset(paired), None)
        if outer_face is None:
            raise ValueError("Periodic Kelvin vertices do not pair whole faces")
        permutation = [paired.index(v.nr) for v in outer_face.vertices[:2]]
        other_rule = ng.IntegrationRule(bary[:, permutation].tolist(), weights)
        left = np.asarray(xyz(mesh.GetTrafo(inner_face)(rule)))
        right = np.asarray(xyz(mesh.GetTrafo(outer_face)(other_rule)))
        max_seam_error = max(max_seam_error, float(np.max(np.linalg.norm(
            right - left - offset, axis=1))))
    return {
        "failed_tetrahedra_by_material": failed,
        "min_signed_scaled_jacobian": minimum if np.isfinite(minimum) else None,
        "jacobian_integration_order": sample_order,
        "kelvin_vertex_pairs": len(identification),
        "kelvin_face_pairs": len(inner),
        "kelvin_curved_match_max_m": max_seam_error,
        "kelvin_curved_match_relative": max_seam_error / radius,
    }

def curve_analytic_ellipsoid(mesh, axes, offset, radius, order):
    """Curve an affine mesh by one shared H1 displacement.

    Interpolate radial CAD projections on the three interfaces, then extend
    the boundary displacement harmonically to the adjacent tetrahedra.
    """
    import ngsolve as ng

    if mesh.GetCurveOrder() != 1:
        raise ValueError("Analytic curving requires an affine input mesh")
    xyz = ng.CoefficientFunction((ng.x, ng.y, ng.z))
    exterior = xyz - ng.CoefficientFunction(offset)
    iron_radius = ng.sqrt(sum((xyz[i] / axis)**2 for i, axis in enumerate(axes)))
    displacement = mesh.BoundaryCF({
        "iron_air": (1 / iron_radius - 1) * xyz,
        "kelvin_int": (radius / ng.Norm(xyz) - 1) * xyz,
        "kelvin_ext": (radius / ng.Norm(exterior) - 1) * exterior,
    }, default=ng.CoefficientFunction((0., 0., 0.)))
    boundary = "iron_air|kelvin_int|kelvin_ext"
    space = ng.VectorH1(mesh, order=order, dirichlet=boundary)
    deformation = ng.GridFunction(space)
    deformation.Set(displacement, definedon=mesh.Boundaries(boundary))
    trial, test = space.TnT()
    energy = ng.BilinearForm(space, symmetric=True)
    energy += ng.InnerProduct(ng.grad(trial), ng.grad(test)) * ng.dx
    energy.Assemble()
    deformation.vec.data -= energy.mat.Inverse(space.FreeDofs(), inverse="sparsecholesky") * (
        energy.mat * deformation.vec)
    mesh.SetDeformation(deformation)


def ellipsoid_boundary_residuals(mesh, axes, offset, radius, order):
    """RMS dimensionless level-set residual against each analytic surface."""
    import ngsolve as ng
    from math import sqrt

    xyz = (ng.x, ng.y, ng.z)
    levelsets = {
        "iron_air": sum((x / a)**2 for x, a in zip(xyz, axes)) - 1,
        "kelvin_int": sum((x / radius)**2 for x in xyz) - 1,
        "kelvin_ext": sum(((x - o) / radius)**2 for x, o in zip(xyz, offset)) - 1,
    }
    return {
        label: sqrt(ng.Integrate(levelset**2, mesh, ng.BND,
                                 definedon=mesh.Boundaries(label), order=order)
                    / ng.Integrate(1, mesh, ng.BND,
                                   definedon=mesh.Boundaries(label), order=order))
        for label, levelset in levelsets.items()
    }
