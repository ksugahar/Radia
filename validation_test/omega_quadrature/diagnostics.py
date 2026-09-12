"""Independent energy contractions and embedding diagnostics for mixed Omega."""
import ngsolve as ng
import numpy as np

MU0 = 4e-7 * np.pi

def norm_report(squared, reference):
    if not np.isfinite(squared) or not np.isfinite(reference):
        raise ValueError("nonfinite embedding norm")
    if squared < -1e-24 or reference < -1e-24:
        raise ValueError("negative squared embedding norm")
    absolute = float(np.sqrt(max(0.0, squared)))
    return {"absolute": absolute,
            "relative": absolute / np.sqrt(reference) if reference > 0 else None}


def constraint_violation(low, high, mesh):
    fes = high["fes"]
    embedded = ng.GridFunction(fes)
    embedded.vec[:] = 0
    reports = {}
    for index, name, materials in (
            (0, "phi_reduced", "air|kelvin"), (1, "phi_total", "iron")):
        source = low["solution"].components[index]
        target = embedded.components[index]
        target.Set(source)
        delta = target - source
        gradient_delta = ng.grad(target) - ng.grad(source)
        region = mesh.Materials(materials)
        norms = {}
        for label, error, ref, selector, vb in (
                ("value", delta * delta, source * source, region, ng.VOL),
                ("gradient", ng.InnerProduct(gradient_delta, gradient_delta),
                 ng.InnerProduct(ng.grad(source), ng.grad(source)), region, ng.VOL),
                ("interface_trace", delta * delta, source * source,
                 mesh.Boundaries("iron_air_interface"), ng.BND)):
            norms[label] = norm_report(
                float(ng.Integrate(error, mesh, vb, definedon=selector, order=16)),
                float(ng.Integrate(ref, mesh, vb, definedon=selector, order=16)))
        constrained = ~np.array(list(target.space.FreeDofs()), dtype=bool)
        values = target.vec.FV().NumPy()
        norms["constrained_coefficient_max_abs"] = (
            float(np.max(np.abs(values[constrained]))) if constrained.any() else None)
        norms["constrained_dof_count"] = int(constrained.sum())
        norms["constraint_scope"] = "all nonfree DOFs, including inactive and zero gauge"
        reports[name] = norms
    r = embedded.vec.CreateVector()
    system = high["system"]
    r.data = system["linear_form"].vec - system["bilinear_form"].mat * embedded.vec
    span = fes.Range(2)
    mask = np.zeros(fes.ndof, dtype=bool)
    mask[span.start:span.stop] = True
    mask &= np.array(list(fes.FreeDofs()), dtype=bool)
    rhs = system["linear_form"].vec.FV().NumPy()
    absolute = float(np.linalg.norm(r.FV().NumPy()[mask]))
    denominator = float(np.linalg.norm(rhs[mask]))
    global_rhs = float(np.linalg.norm(rhs))
    return {
        "embedding_relative_error": {
            name: values["gradient"]["relative"] for name, values in reports.items()},
        "embedding_checks": reports,
        "constraint_violation_l2": absolute,
        "constraint_rhs_l2": denominator,
        "constraint_violation_relative": absolute / denominator if denominator else None,
        "constraint_global_rhs_relative": absolute / global_rhs if global_rhs else None,
        "interpretation": "diagnostic only; inspect value, trace and gauge before nesting"}


def audit_energy(result, mesh, h_source, source_ext, harmonic, bonus, q):
    """Audit the already re-solved system, never label this a new solve."""
    fes = result["fes"]
    trial, test = fes.TnT()
    solution = result["solution"]
    pr, pt = solution.components[:2]
    mu = result["mu_cf"]
    rules = {et: ng.IntegrationRule(et, q) for et in (ng.ET.TET, ng.ET.HEX, ng.ET.PRISM)}
    rows = []
    for rule_name, options in (("assembly_default", {"bonus_intorder": bonus}),
                               ("fixed_rule", {"intrules": rules})):
        terms = []
        for name, region, expr, value, bilinear in (
                ("quadratic_coupled", "air|kelvin",
                 mu * ng.grad(trial[0]) * ng.grad(test[0]),
                 mu * ng.InnerProduct(ng.grad(pr), ng.grad(pr)), True),
                ("quadratic_iron", "iron", mu * ng.grad(trial[1]) * ng.grad(test[1]),
                 mu * ng.InnerProduct(ng.grad(pt), ng.grad(pt)), True),
                ("load_air", "air", mu * h_source * ng.grad(test[0]),
                 mu * h_source * ng.grad(pr), False),
                ("load_kelvin", "kelvin", -mu * source_ext * ng.grad(test[0]),
                 -mu * source_ext * ng.grad(pr), False),
                ("load_iron", "iron", mu * harmonic * ng.grad(test[1]),
                 mu * harmonic * ng.grad(pt), False)):
            measure = ng.dx(definedon=mesh.Materials(region), **options)
            form = ng.BilinearForm(fes, symmetric=True) if bilinear else ng.LinearForm(fes)
            form += expr * measure
            form.Assemble()
            if bilinear:
                product = solution.vec.CreateVector()
                product.data = form.mat * solution.vec
                contracted = 0.5 * float(np.dot(solution.vec.FV().NumPy(), product.FV().NumPy()))
            else:
                contracted = float(np.dot(form.vec.FV().NumPy(), solution.vec.FV().NumPy()))
            row = {"term": name, "contracted": contracted}
            if rule_name == "fixed_rule":
                integrated = float(ng.Integrate(value * measure, mesh)) * (0.5 if bilinear else 1)
                row.update(integrated=integrated, difference=integrated - contracted)
            terms.append(row)
            del form
        half = sum(t["contracted"] for t in terms[:2])
        linear = sum(t["contracted"] for t in terms[2:])
        rows.append({"rule": rule_name, "terms": terms, "half_xAx": half,
                     "b_dot_x": linear, "energy": half - linear})
        print("  term audit", rule_name, "J=", half-linear, flush=True)
    result["assembled_energy"]["term_audit"] = rows
    result["assembled_energy"]["fixed_rule_order"] = q
    result["assembled_energy"]["default_reconstruction_difference"] = (
        rows[0]["energy"] - result["assembled_energy"]["energy"])
    h = result["H_cf"]
    fixed_W = 0.0
    for region, density in (
            ("air", MU0 * ng.InnerProduct(h-h_source, h+h_source)),
            ("iron", mu * ng.InnerProduct(h, h) - MU0 * ng.InnerProduct(h_source, h_source)),
            ("kelvin", mu * ng.InnerProduct(h-source_ext, h+source_ext))):
        fixed_W += 0.5 * float(ng.Integrate(
            density * ng.dx(definedon=mesh.Materials(region), intrules=rules), mesh))
    offset = 0.5 * float(ng.Integrate(
        (mu * ng.InnerProduct(harmonic, harmonic)
         - MU0 * ng.InnerProduct(h_source, h_source))
        * ng.dx(definedon=mesh.Materials("iron"), intrules=rules), mesh))
    result["assembled_energy"]["fixed_rule_identity"] = {
        "W": fixed_W, "J": rows[1]["energy"], "source_offset": offset,
        "W_minus_J_minus_offset": fixed_W - rows[1]["energy"] - offset}
    return result
