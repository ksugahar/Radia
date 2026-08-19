"""Section-optics inverse design: an optics specification becomes a field difference.

Specifying the exit state of a magnet is a badly posed design problem.  It
is a handful of numbers against a map with hundreds of free coefficients,
so the solution set is a large manifold and any optimizer wanders along
it.  Specifying what one SECTION of the magnet must do optically is well
posed, because the map is built by ordered composition of per-element
maps and therefore

    M = M_after . M_S . M_before

holds exactly for any contiguous run ``S`` of elements.  ``M_S`` depends
only on the field inside ``S``, so a specification written on it is local
in the field and rich in content -- a whole transfer block rather than a
few exit numbers.  That identity is what this module rests on, and
:func:`section_composition_defect` measures it rather than assuming it.

The design step solves

    minimize   1/2 |W^(1/2) d|^2     subject to   J d = requested

where ``J = d(spec)/d(parameters)`` and ``W`` is a PHYSICAL field metric,
not a coefficient norm.  The metric is where the designer's intent enters
quantitatively, and it is not a formality: with a plain coefficient
2-norm the recovered difference cancels violently between components --
mathematically minimal, physically absurd, and unreachable by iron.

Two representations are supported, and they are not interchangeable.

``multipole_*``
    The per-segment transverse multipole profile through decapole, with
    the diagonal aperture-edge metric of :func:`aperture_field_metric`.
    This is the standard accelerator description and its Jacobian is
    analytic, because the map builder already propagates forward-mode
    tangents with respect to every multipole it consumes.  It fails on
    strongly curved, strongly graded magnets: at field index 15 the
    builder's own analytic-versus-native consistency gate refuses the
    profile outright.

``chain_*``
    The CanonicalHCurl chain's reduced coefficients, with the dense
    physical field norm ``A^T Omega A`` built from the chain's own
    evaluation operator.  This is the representation the three-route
    tracking certificate is taken in, and it is the one to use when the
    multipole route refuses.  Its Jacobian is by central differences --
    the chain's segment-polynomial path carries no forward-mode tangents
    -- which costs two map builds per coefficient and is still seconds.

Both feed the same solver, :func:`constrained_minimum_norm_step`.

Three things this module reports rather than hides.

*Inert knobs are dropped.*  A specification built from linear-map entries
is controlled by the dipole and quadrupole alone; the sextupole and above
move ``T/U/V`` and leave ``R`` untouched, so their Jacobian columns carry
no information.  A minimum-norm solution still assigns them something,
and under an aperture metric their inverse weights are of order ``a^-8``,
so whatever noise the Jacobian carries there is divided by a very small
number before it lands in the coefficients that downstream code consumes.
Whether there is noise to divide depends on how the Jacobian was
obtained: the analytic multipole route returns exact zeros, while a
DIFFERENCED Jacobian -- which is what the chain route must use -- returns
round-off.  Holding those columns at zero makes the recovered difference
independent of that noise, and the filter measures each column's
authority in the whitened metric the inversion actually minimizes.  The
failure that motivated it was a differenced Jacobian producing high-order
coefficients large enough for the map builder's own analytic-versus
-native consistency gate to refuse the profile.

*The planes may be locked together.*  :func:`focusing_coupling` measures
how nearly the horizontal and vertical focal powers are one knob.  A
straight gradient slab sits at ``|cos| = 0.9999``, so holding one plane
while changing the other costs about seventy times the field; a curved
exit fringe measures worse, not better, because its curvature is frozen
into the reference orbit and offers no independent freedom.  Asking such
a section to decouple the planes is not a hard inverse problem but very
nearly an impossible one, and the multiplier says so numerically instead
of returning an enormous field difference with no explanation.

*The linearization is checked, not trusted.*
:func:`verify_multipole_linearization` applies the recovered difference to
the actual profile and rebuilds the map from scratch.  A ratio far from
one means the step left the linear regime.

Finally, a warning that cost a full design run: a section is a GEOMETRIC
region and must be pinned to a fixed arc-length interval.  Defining it by
a field-dependent threshold -- "from where the body field last exceeded
90 % of its peak" -- makes the boundary move when the design changes the
field, and the before/after comparison then measures a change of
definition rather than of physics.  :func:`snap_breaks_to_section` builds
element breaks that represent a fixed interval exactly.
"""
from __future__ import annotations

import numpy as np

from radia.accelerator_lie_topopt import (
    FOURTH_ORDER_MULTIPOLE_COMPONENTS,
    _fourth_order_lie_map_from_vector_potential_polynomials,
    fourth_order_lie_map_from_multipoles,
)

__all__ = [
    "COMPONENT_COUNT",
    "NORMAL_BLOCKS",
    "NORMAL_DEGREES",
    "SPEC_NAMES",
    "aperture_field_metric",
    "chain_bend_row",
    "chain_field_metric",
    "chain_field_operator",
    "chain_section_spec",
    "chain_section_spec_jacobian",
    "chain_section_transfer",
    "constrained_minimum_norm_step",
    "focusing_coupling",
    "minimum_norm_field_difference",
    "multipole_section_map",
    "multipole_section_spec",
    "multipole_section_spec_jacobian",
    "multipole_section_spec_jacobian_by_differences",
    "multipole_section_transfer",
    "scatter_multipole_difference",
    "section_composition_defect",
    "snap_breaks_to_section",
    "verify_multipole_linearization",
]

COMPONENT_COUNT = len(FOURTH_ORDER_MULTIPOLE_COMPONENTS)
#: Normal dipole/quadrupole/sextupole/octupole/decapole blocks.  The skew
#: blocks vanish identically for a midplane-symmetric magnet, so they are
#: not design variables and are held at zero.
NORMAL_BLOCKS = (0, 1, 3, 5, 7)
NORMAL_DEGREES = (0, 1, 2, 3, 4)
#: The specification vector produced by the ``*_section_spec`` functions.
SPEC_NAMES = ("horizontal_focal_power", "vertical_focal_power",
              "integrated_bend")


# --------------------------------------------------------------------------
# representation-independent core
# --------------------------------------------------------------------------
def constrained_minimum_norm_step(jacobian, requested, metric, *,
                                  controllability_tolerance=1.0e-8):
    """Smallest ``metric``-norm step satisfying ``jacobian @ step = requested``.

    ``metric`` is any symmetric positive-definite matrix, so one solver
    serves both representations: a diagonal aperture-field weight for the
    multipole profile, a dense ``A^T Omega A`` for the chain coefficients.

    The step itself is computed in whitened coordinates ``u = L^T d``
    where ``metric = L L^T``, in which the objective is a plain 2-norm.

    Rows are equilibrated first.  The specification entries differ by
    orders of magnitude in units (1/m, 1/m, T m), and without
    equilibration the conditioning of the reduced Gram matrix is set by
    that unit choice rather than by the physics.

    The inert-knob filter measures each column on the EQUILIBRATED
    JACOBIAN, not in the whitened metric.  Whether a column carries
    information is a property of the Jacobian alone; the metric only sets
    what the column costs.  Mixing the two rejects the wrong columns: an
    aperture metric can give a decapole column an inverse weight of
    ``a^-8``, which lifts pure round-off well above a relative threshold
    that was meant to catch exactly that round-off.  Comparing raw column
    norms asks the right question -- is this entry signal or noise -- and
    is independent of how the design chose to price the column.

    Returns ``(step, kept)``; ``kept`` is the boolean mask of columns that
    were treated as design variables.
    """
    jacobian = np.asarray(jacobian, dtype=float)
    scale = np.linalg.norm(jacobian, axis=1)
    if not np.all(np.isfinite(scale)) or np.any(scale <= 0.0):
        raise ValueError(
            "a specification functional has an identically zero Jacobian "
            "row: it cannot be controlled by this section")
    scaled = jacobian / scale[:, None]
    target = np.asarray(requested, dtype=float).reshape(-1) / scale
    lower = np.linalg.cholesky(np.asarray(metric, dtype=float))
    authority = np.linalg.norm(scaled, axis=0)
    kept = authority >= controllability_tolerance * float(np.max(authority))
    if not np.any(kept):
        raise ValueError("no design variable in this section controls the "
                         "specification")
    whitened = np.linalg.solve(lower, scaled.T).T
    active = whitened[:, kept]
    coordinates = np.zeros(jacobian.shape[1])
    coordinates[kept] = active.T @ np.linalg.solve(active @ active.T, target)
    return np.linalg.solve(lower.T, coordinates), kept


def minimum_norm_field_difference(jacobian, requested, weight, *,
                                  controllability_tolerance=1.0e-8):
    """:func:`constrained_minimum_norm_step` with a diagonal metric."""
    return constrained_minimum_norm_step(
        jacobian, requested, np.diag(np.asarray(weight, dtype=float)),
        controllability_tolerance=controllability_tolerance)


def focusing_coupling(jacobian, metric, kept=None):
    """How locked together the horizontal and vertical focal powers are.

    Returns ``(cosine, cost_multiplier)``.  The cosine is taken in the
    metric the inversion actually minimizes, so it answers a design
    question rather than a linear-algebra one: at ``|cos| = 1`` the two
    focal powers cannot be moved independently at all, and the price of
    holding one while changing the other is ``1 / sqrt(1 - cos^2)`` times
    the field the unconstrained change would need.

    ``metric`` may be the diagonal weight vector or the full matrix.
    """
    matrix = np.asarray(jacobian, dtype=float)
    metric = np.asarray(metric, dtype=float)
    if kept is not None:
        matrix = matrix[:, kept]
        metric = metric[kept] if metric.ndim == 1 else metric[np.ix_(kept,
                                                                    kept)]
    if metric.ndim == 1:
        whitened = matrix / np.sqrt(metric)[None, :]
    else:
        whitened = np.linalg.solve(np.linalg.cholesky(metric), matrix.T).T
    horizontal, vertical = whitened[0], whitened[1]
    denominator = np.linalg.norm(horizontal) * np.linalg.norm(vertical)
    cosine = float(horizontal @ vertical / denominator) if denominator > 0.0 \
        else 0.0
    residual = max(1.0 - cosine**2, 0.0)
    multiplier = float("inf") if residual <= 0.0 else 1.0 / np.sqrt(residual)
    return cosine, multiplier


def snap_breaks_to_section(breaks, section_start):
    """Move the element break nearest ``section_start`` exactly onto it.

    A section is a geometric region of the magnet, so before/after
    comparisons must use the same arc-length interval.  Element breaks are
    normally placed by a field-driven grading rule, which shifts them when
    the design changes the field; snapping restores an exact interval
    boundary while leaving the grading otherwise intact.

    Returns ``(breaks, index, moved_m)``.
    """
    breaks = np.array(breaks, dtype=float)
    index = int(np.argmin(np.abs(breaks - float(section_start))))
    if index in (0, breaks.size - 1):
        raise ValueError("the section start falls on a chain end")
    moved = float(abs(breaks[index] - section_start))
    breaks[index] = float(section_start)
    if not np.all(np.diff(breaks) > 0.0):
        raise ValueError("snapping the section start broke break ordering")
    return breaks, index, moved


def section_composition_defect(before, section, whole):
    """``max |R_S . R_before - R_whole|``, the composition identity residual.

    The whole design model rests on ``M = M_after . M_S . M_before`` being
    exact rather than approximate.  Measuring it costs three matrix
    products and turns a structural assumption into a reported number.
    """
    composed = np.asarray(section, dtype=float) @ np.asarray(before,
                                                             dtype=float)
    return float(np.max(np.abs(composed - np.asarray(whole, dtype=float))))


# --------------------------------------------------------------------------
# multipole representation
# --------------------------------------------------------------------------
def multipole_section_map(response, lengths, curvature, rigidity, begin, end):
    """Fourth-order map of the contiguous segment run ``[begin, end)``."""
    blocks = np.asarray(response, dtype=float).reshape(COMPONENT_COUNT, -1)
    return fourth_order_lie_map_from_multipoles(
        blocks[:, begin:end].reshape(-1), np.asarray(lengths)[begin:end],
        rigidity,
        reference_curvature_per_m=np.asarray(curvature)[begin:end])


def multipole_section_transfer(response, lengths, curvature, rigidity,
                               begin, end):
    """Linear transfer matrix of the segment run ``[begin, end)``."""
    return multipole_section_map(response, lengths, curvature, rigidity,
                                 begin, end).factorization.R


def multipole_section_spec(response, lengths, curvature, rigidity,
                           begin, end):
    """The section's ``[R_S[1,0], R_S[3,2], integral B_y ds]``.

    The integrated bend is in tesla metre and is what must be held fixed
    when only the focusing is being changed: the construction is expressed
    about a fixed reference orbit whose curvature is a frozen geometric
    input, and holding the bend is what keeps that orbit valid to first
    order.  It does not forbid a local dipole bump that integrates to
    zero; only a closed loop that re-tracks the orbit can catch that.
    """
    matrix = multipole_section_transfer(response, lengths, curvature,
                                        rigidity, begin, end)
    blocks = np.asarray(response, dtype=float).reshape(COMPONENT_COUNT, -1)
    lengths = np.asarray(lengths, dtype=float)
    return np.array([matrix[1, 0], matrix[3, 2],
                     float(np.sum(blocks[0][begin:end] * lengths[begin:end]))])


def _section_columns(begin, end, blocks_used):
    return [(int(block), begin + local)
            for block in blocks_used for local in range(end - begin)]


def multipole_section_spec_jacobian(response, lengths, curvature, rigidity,
                                    begin, end, *,
                                    blocks_used=NORMAL_BLOCKS):
    """Analytic ``d(spec)/d(multipole)`` for the normal blocks in the section.

    The map builder propagates forward-mode tangents with respect to every
    multipole it consumes, so one map build yields the whole Jacobian.
    The bend row is exact by inspection: only the dipole block integrates
    to bend, with weight ``ds``.

    Returns ``(jacobian, index)`` where ``index`` lists the
    ``(component_block, segment)`` pair each column belongs to.
    """
    result = multipole_section_map(response, lengths, curvature, rigidity,
                                   begin, end)
    tangent = result.factorization.R_jacobian
    count = end - begin
    expected = COMPONENT_COUNT * count
    if tangent.shape != (expected, 6, 6):
        raise RuntimeError(
            f"map Jacobian carries {tangent.shape[0]} parameters, expected "
            f"{expected} = {COMPONENT_COUNT} components x {count} segments")
    lengths = np.asarray(lengths, dtype=float)
    index = _section_columns(begin, end, blocks_used)
    columns = [
        [tangent[block * count + segment - begin, 1, 0],
         tangent[block * count + segment - begin, 3, 2],
         float(lengths[segment]) if block == 0 else 0.0]
        for block, segment in index]
    return np.column_stack(columns), np.asarray(index)


def multipole_section_spec_jacobian_by_differences(
        response, lengths, curvature, rigidity, begin, end, *,
        relative_step=1.0e-5, blocks_used=NORMAL_BLOCKS):
    """Central differences, kept to check the analytic Jacobian."""
    base = np.asarray(response, dtype=float).reshape(COMPONENT_COUNT, -1)
    index = _section_columns(begin, end, blocks_used)
    columns = []
    for block, segment in index:
        scale = float(np.max(np.abs(base[block][begin:end])))
        step = relative_step * (scale if scale > 0.0 else 1.0)
        perturbed = base.copy()
        perturbed[block, segment] += step
        plus = multipole_section_spec(perturbed.reshape(-1), lengths,
                                      curvature, rigidity, begin, end)
        perturbed[block, segment] -= 2.0 * step
        minus = multipole_section_spec(perturbed.reshape(-1), lengths,
                                       curvature, rigidity, begin, end)
        columns.append((plus - minus) / (2.0 * step))
    return np.column_stack(columns), np.asarray(index)


def aperture_field_metric(index, lengths, half_aperture):
    """Diagonal weight making the objective an aperture-edge field norm.

    Multipole ``b_n`` contributes ``b_n x^n`` to ``B_y`` at transverse
    offset ``x``, so weighting by ``(a^n)^2 ds`` measures the squared field
    difference at the aperture edge integrated along the section.  That is
    a physical quantity, unlike a bare coefficient 2-norm whose components
    do not even share units.
    """
    degree = dict(zip(NORMAL_BLOCKS, NORMAL_DEGREES))
    lengths = np.asarray(lengths, dtype=float)
    return np.asarray([
        (float(half_aperture) ** degree[int(block)]) ** 2
        * float(lengths[segment]) for block, segment in index])


def scatter_multipole_difference(difference, index, shape):
    """Place a recovered column vector back into a full response array."""
    full = np.zeros(shape, dtype=float)
    for value, (block, segment) in zip(difference, index):
        full[int(block), int(segment)] = value
    return full


def verify_multipole_linearization(response, lengths, curvature, rigidity,
                                   begin, end, difference, requested):
    """Apply the recovered difference and re-measure the specification.

    The Jacobian is exact only to first order.  This applies the actual
    field difference to the actual profile, rebuilds the section map from
    scratch, and reports how much of the request was really delivered.  A
    ratio far from one means the step left the linear regime, not that the
    inversion was solved badly.

    Returns ``(before, after, delivered, ratio)``.
    """
    base = np.asarray(response, dtype=float).reshape(COMPONENT_COUNT, -1)
    before = multipole_section_spec(base.reshape(-1), lengths, curvature,
                                    rigidity, begin, end)
    after = multipole_section_spec((base + difference).reshape(-1), lengths,
                                   curvature, rigidity, begin, end)
    delivered = after - before
    asked = np.asarray(requested, dtype=float)
    active = np.abs(asked) > 0.0
    ratio = np.ones_like(asked)
    ratio[active] = delivered[active] / asked[active]
    return before, after, delivered, ratio


# --------------------------------------------------------------------------
# CanonicalHCurl chain representation
# --------------------------------------------------------------------------
def _element_offsets(chain):
    return np.concatenate(([0], np.cumsum(
        [element.dimension for element in chain.elements])))


def chain_field_operator(chain, probe_s, nodes_x, nodes_y):
    """Rows mapping chain coefficients to ``(g bx, g by, bs)`` on a grid.

    Built from the chain's own element row builder, so the operator is the
    exact linearization of the field the design may change.  That is what
    makes ``A^T Omega A`` a physical field norm rather than a coefficient
    norm.  Returns ``(operator, grid_s, offsets)``.
    """
    grid_s, grid_x, grid_y = (array.reshape(-1) for array in np.meshgrid(
        np.asarray(probe_s, dtype=float), np.asarray(nodes_x, dtype=float),
        np.asarray(nodes_y, dtype=float), indexing="ij"))
    index, zeta = chain._locate(grid_s)
    offsets = _element_offsets(chain)
    count = grid_s.size
    rows = np.zeros((3 * count, offsets[-1]))
    for element_index, element in enumerate(chain.elements):
        mask = index == element_index
        if not np.any(mask):
            continue
        gbx, gby, bs_columns = element.b_row_columns(
            grid_x[mask] / element.half_width_m,
            grid_y[mask] / element.half_height_m, zeta[mask])
        where = np.flatnonzero(mask)
        columns = slice(offsets[element_index], offsets[element_index + 1])
        rows[where, columns] = gbx
        rows[count + where, columns] = gby
        rows[2 * count + where, columns] = bs_columns
    return rows @ chain._reduced, grid_s, offsets


def chain_field_metric(operator, grid_s, section_bounds, *,
                       outside_weight=1.0e2, ridge=1.0e-10):
    """Physical field norm that localizes the change to the section.

    The section cannot be isolated by freezing coefficients outside it:
    the chain reduction imposes interface continuity, so the elements are
    tied together and a hard freeze would be infeasible rather than local.
    Locality is imposed by making a field change outside the section
    ``outside_weight`` times dearer than one inside.
    """
    begin, end = (float(value) for value in section_bounds)
    inside = (np.asarray(grid_s) >= begin) & (np.asarray(grid_s) < end)
    omega = np.concatenate([np.where(inside, 1.0, float(outside_weight))] * 3)
    operator = np.asarray(operator, dtype=float)
    metric = operator.T @ (omega[:, None] * operator)
    metric += float(ridge) * float(np.trace(metric)) / metric.shape[0] \
        * np.eye(metric.shape[0])
    return metric, inside


def chain_bend_row(chain, monitor_s):
    """Exact ``d(integral B_y ds)/d(coefficients)`` on the design orbit."""
    monitor_s = np.asarray(monitor_s, dtype=float)
    quadrature = np.gradient(monitor_s)
    index, zeta = chain._locate(monitor_s)
    offsets = _element_offsets(chain)
    row = np.zeros(offsets[-1])
    for element_index, element in enumerate(chain.elements):
        mask = index == element_index
        if not np.any(mask):
            continue
        zeros = np.zeros(int(np.sum(mask)))
        _gbx, gby, _bs = element.b_row_columns(zeros, zeros, zeta[mask])
        row[offsets[element_index]:offsets[element_index + 1]] += \
            quadrature[mask] @ gby
    return row @ chain._reduced


def chain_section_transfer(chain, coefficients, members, rigidity, *,
                           reference_orbit_tolerance=2.0e-2):
    """Linear transfer matrix of the chain elements selected by ``members``.

    ``members`` is any slice or index array over the chain's elements.
    ``coefficients`` are the REDUCED coefficients; the chain's stored fit
    is temporarily replaced and always restored.
    """
    from dataclasses import replace

    saved = chain._fit
    if saved is None:
        raise RuntimeError("fit the chain before evaluating a section map")
    chain._fit = replace(
        saved, coefficients=chain._reduced @ np.asarray(coefficients,
                                                        dtype=float))
    try:
        ay, a_s, lengths, curvatures = chain.lie_element_spoly_arrays(degree=5)
    finally:
        chain._fit = saved
    result = _fourth_order_lie_map_from_vector_potential_polynomials(
        ay[members], a_s[members], lengths[members], float(rigidity),
        reference_curvature_per_m=curvatures[members],
        longitudinal_component="covariant",
        reference_orbit_tolerance=float(reference_orbit_tolerance),
        parameter_jacobians=False)
    return result.transfer.factorization.R


def chain_section_spec(chain, coefficients, members, rigidity, bend_row, **kw):
    """The section's ``[R_S[1,0], R_S[3,2], integral B_y ds]`` from the chain."""
    matrix = chain_section_transfer(chain, coefficients, members, rigidity,
                                    **kw)
    return np.array([float(matrix[1, 0]), float(matrix[3, 2]),
                     float(np.asarray(bend_row) @ np.asarray(coefficients))])


def chain_section_spec_jacobian(chain, coefficients, members, rigidity,
                                bend_row, *, relative_step=1.0e-6, **kw):
    """``d(spec)/d(coefficients)`` by central differences.

    The chain's segment-polynomial path carries no forward-mode tangents,
    so unlike the multipole route this is a difference quotient.  It costs
    two map builds per coefficient, which is seconds for a production
    chain, and the bend row is supplied exactly rather than differenced.
    """
    coefficients = np.asarray(coefficients, dtype=float)
    bend_row = np.asarray(bend_row, dtype=float)
    step = float(relative_step) * float(np.max(np.abs(coefficients)))
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("cannot choose a finite difference step from these "
                         "coefficients")
    jacobian = np.zeros((3, coefficients.size))
    jacobian[2] = bend_row
    for column in range(coefficients.size):
        plus, minus = coefficients.copy(), coefficients.copy()
        plus[column] += step
        minus[column] -= step
        upper = chain_section_transfer(chain, plus, members, rigidity, **kw)
        lower = chain_section_transfer(chain, minus, members, rigidity, **kw)
        jacobian[0, column] = (upper[1, 0] - lower[1, 0]) / (2.0 * step)
        jacobian[1, column] = (upper[3, 2] - lower[3, 2]) / (2.0 * step)
    return jacobian
