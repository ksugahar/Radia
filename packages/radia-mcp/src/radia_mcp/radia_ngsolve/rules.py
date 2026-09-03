"""
Lint rules for NGSolve Python scripts.

Each rule function receives (filepath, lines) and returns a list of findings.
A finding is a dict with keys: line, severity, rule, message.
"""

import re
import os
from typing import List, Dict


# ── NGSolve FEM rules ─────────────────────────────────────────

def check_hcurl_missing_nograds(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: Magnetostatic HCurl spaces should use nograds=True."""
    findings = []
    has_magnetostatic = any(
        kw in line.lower()
        for line in lines
        for kw in ['magnetostatic', 'curl(u)*curl(v)', 'curl-curl',
                    'vector potential', 'magnet']
    )
    if not has_magnetostatic:
        return findings

    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith('#'):
            continue
        if re.search(r'HCurl\s*\(', stripped):
            has_nograds = 'nograds=True' in stripped or 'nograds = True' in stripped
            is_complex = 'complex=True' in stripped or 'complex = True' in stripped
            if not has_nograds and not is_complex:
                findings.append({
                    'line': i,
                    'severity': 'HIGH',
                    'rule': 'hcurl-missing-nograds',
                    'message': (
                        'HCurl space for magnetostatics should use nograds=True '
                        'to remove gradient null space. Without it, the curl-curl '
                        'system is singular. Add nograds=True parameter.'
                    ),
                })
    return findings


def check_ngsolve_precond_after_assemble(filepath: str, lines: List[str]) -> List[Dict]:
    """MODERATE: BDDC preconditioner must be registered before assembly."""
    findings = []
    has_bddc = any('Preconditioner' in line and 'bddc' in line for line in lines)
    if not has_bddc:
        return findings

    assemble_lines = []
    precond_lines = []
    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith('#'):
            continue
        if '.Assemble()' in stripped:
            assemble_lines.append(i)
        if re.search(r'Preconditioner\s*\(.*bddc', stripped):
            precond_lines.append(i)

    for p_line in precond_lines:
        for a_line in assemble_lines:
            if a_line < p_line:
                findings.append({
                    'line': p_line,
                    'severity': 'MODERATE',
                    'rule': 'ngsolve-precond-after-assemble',
                    'message': (
                        'BDDC Preconditioner registered AFTER .Assemble() '
                        '(line {}). It must be registered BEFORE assembly '
                        'to access element matrices. Move Preconditioner() '
                        'before .Assemble().'.format(a_line)
                    ),
                })
                break
    return findings


def check_ngsolve_missing_trace_bem(filepath: str, lines: List[str]) -> List[Dict]:
    """CRITICAL: BEM on HDivSurface requires .Trace() on trial/test functions."""
    findings = []
    has_bem_hdiv = any(
        'HDivSurface' in line and ('LaplaceSL' in ''.join(lines) or
                                    'HelmholtzSL' in ''.join(lines))
        for line in lines
    )
    if not has_bem_hdiv:
        return findings

    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith('#'):
            continue
        if re.search(r'(?:LaplaceSL|HelmholtzSL)\s*\(', stripped):
            if '.Trace()' not in stripped:
                context = ' '.join(lines[max(0,i-2):min(len(lines),i+2)])
                if '.Trace()' not in context:
                    findings.append({
                        'line': i,
                        'severity': 'CRITICAL',
                        'rule': 'ngsolve-missing-trace-bem',
                        'message': (
                            'BEM operator (LaplaceSL/HelmholtzSL) on HDivSurface '
                            'requires .Trace() on trial/test functions. Without '
                            '.Trace(), boundary-edge DOFs get corrupted diagonal '
                            'entries, causing wildly wrong results. '
                            'Use: j_trial.Trace() * ds(...)'
                        ),
                    })
    return findings


def check_ngsolve_overwrite_xyz(filepath: str, lines: List[str]) -> List[Dict]:
    """MODERATE: Overwriting x/y/z coordinate variables in loops."""
    findings = []
    has_ngsolve = any(
        kw in line
        for line in lines
        for kw in ['from ngsolve', 'import ngsolve', 'from ngsolve import']
    )
    if not has_ngsolve:
        return findings

    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith('#'):
            continue
        match = re.match(r'for\s+(x|y|z)\s+in\s+', stripped)
        if match:
            var = match.group(1)
            findings.append({
                'line': i,
                'severity': 'MODERATE',
                'rule': 'ngsolve-overwrite-xyz',
                'message': (
                    f'Loop variable "{var}" overwrites NGSolve coordinate '
                    f'CoefficientFunction. After this loop, {var} will be a '
                    f'scalar, not a coordinate. Use a different variable name '
                    f'(e.g., "{var}i" or "{var}_val").'
                ),
            })
    return findings


def check_ngsolve_vec_assign(filepath: str, lines: List[str]) -> List[Dict]:
    """MODERATE: NGSolve vector assignment must use .data attribute."""
    findings = []
    has_ngsolve = any(
        kw in line
        for line in lines
        for kw in ['from ngsolve', 'import ngsolve', 'GridFunction']
    )
    if not has_ngsolve:
        return findings

    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith('#'):
            continue
        if re.search(r'\.vec\s*=\s*(?!.*\.data)', stripped):
            if '.vec.data' not in stripped and '.vec[' not in stripped:
                if 'CreateVector' not in stripped:
                    findings.append({
                        'line': i,
                        'severity': 'MODERATE',
                        'rule': 'ngsolve-vec-assign',
                        'message': (
                            'Direct assignment to .vec creates a symbolic '
                            'expression, not an evaluated result. '
                            'Use .vec.data = ... to evaluate and store, or '
                            '.vec[:] = ... for slice assignment.'
                        ),
                    })
    return findings


def check_ngsolve_dim2_occ(filepath: str, lines: List[str]) -> List[Dict]:
    """MODERATE: 2D OCC geometry requires dim=2 parameter."""
    findings = []
    has_2d_geom = any(
        kw in line
        for line in lines
        for kw in ['Rectangle', 'Circle', 'Face()', 'WorkPlane']
    )
    if not has_2d_geom:
        return findings

    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith('#'):
            continue
        if re.search(r'OCCGeometry\s*\(', stripped) and 'dim=' not in stripped:
            context = ' '.join(lines[max(0,i-5):i])
            if any(kw in context for kw in ['Rectangle', '.Face()', 'WorkPlane',
                                              'Circle']):
                findings.append({
                    'line': i,
                    'severity': 'MODERATE',
                    'rule': 'ngsolve-dim2-occ',
                    'message': (
                        'OCCGeometry with 2D shapes (Rectangle, Face) requires '
                        'dim=2 parameter. Without it, a 3D surface mesh is '
                        'generated instead of a 2D mesh. '
                        'Use: OCCGeometry(shape, dim=2)'
                    ),
                })
    return findings


def check_ngsolve_cg_on_saddle_point(filepath: str, lines: List[str]) -> List[Dict]:
    """MODERATE: CG used on A-Omega saddle-point system (should use GMRes/MinRes)."""
    findings = []
    has_saddle = False
    for line in lines:
        stripped = line.strip()
        if any(kw in stripped for kw in [
            'curl(N).Trace() * normal',
            'curl(A).Trace() * normal',
            'A_ReducedOmega',
            'fesA * fesOmega',
        ]):
            has_saddle = True
            break

    if not has_saddle:
        return findings

    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith('#'):
            continue
        if re.search(r'solvers\.CG\s*\(', stripped) or \
           re.search(r'CGSolver\s*\(', stripped):
            findings.append({
                'line': i,
                'severity': 'MODERATE',
                'rule': 'ngsolve-cg-on-saddle-point',
                'message': (
                    'CG solver used on a saddle-point system (A-Omega mixed '
                    'formulation). The system is INDEFINITE, so CG may diverge. '
                    'Use MinRes, GMRes, or a direct solver instead.'
                ),
            })
    return findings


def check_ngsolve_vectorh1_for_em(filepath: str, lines: List[str]) -> List[Dict]:
    """MODERATE: VectorH1 is wrong for electromagnetic fields (use HCurl/HDiv)."""
    findings = []
    has_em_context = any(
        kw in line.lower()
        for line in lines
        for kw in ['magnetostatic', 'maxwell', 'curl(u)', 'curl(v)',
                    'vector potential', 'electric field', 'magnetic',
                    'hcurl', 'hdiv', 'b-field', 'e-field']
    )
    if not has_em_context:
        return findings

    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith('#'):
            continue
        if re.search(r'VectorH1\s*\(', stripped):
            findings.append({
                'line': i,
                'severity': 'MODERATE',
                'rule': 'ngsolve-vectorh1-for-em',
                'message': (
                    'VectorH1 used in electromagnetic context. VectorH1 enforces '
                    'full C^0 continuity on ALL components, which is wrong for EM '
                    'fields. Use HCurl (tangential continuity for E, A) or '
                    'HDiv (normal continuity for B, J).'
                ),
            })
    return findings


def check_ngsolve_pinvit_no_projection(filepath: str, lines: List[str]) -> List[Dict]:
    """MODERATE: PINVIT/LOBPCG eigenvalue solver without gradient projection."""
    findings = []
    has_pinvit = any(
        kw in line
        for line in lines
        for kw in ['solvers.PINVIT', 'PINVIT(', 'solvers.LOBPCG', 'LOBPCG(']
    )
    if not has_pinvit:
        return findings

    has_projection = any(
        kw in line
        for line in lines
        for kw in ['CreateGradient', 'gradmat', 'grad_projection', 'proj @']
    )
    has_hcurl = any('HCurl' in line for line in lines)

    if has_hcurl and not has_projection:
        for i, line in enumerate(lines, 1):
            stripped = line.strip()
            if re.search(r'(?:PINVIT|LOBPCG)\s*\(', stripped):
                findings.append({
                    'line': i,
                    'severity': 'MODERATE',
                    'rule': 'ngsolve-pinvit-no-projection',
                    'message': (
                        'PINVIT/LOBPCG eigenvalue solver on HCurl space without '
                        'gradient projection. The curl-curl null space (gradient '
                        'fields) produces spurious zero eigenvalues. Use '
                        'fes.CreateGradient() to build projection matrix.'
                    ),
                })
    return findings


def check_eddy_current_missing_complex(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: Eddy current HCurl/H1 spaces must use complex=True."""
    findings = []
    has_eddy_context = any(
        kw in line.lower()
        for line in lines
        for kw in ['eddy', 'induction_heating', 'joule', 'freq', '2j *',
                    'a + grad(phi)', 'a+grad(phi)']
    )
    if not has_eddy_context:
        return findings

    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith('#'):
            continue
        if re.search(r'(?:HCurl|H1)\s*\(', stripped):
            has_complex = 'complex=True' in stripped or 'complex = True' in stripped
            is_thermal = any(
                kw in stripped.lower()
                for kw in ['temperature', 'thermal', 'heat']
            )
            if not has_complex and not is_thermal:
                context = ' '.join(lines[max(0,i-3):i+2])
                is_thermal_ctx = any(
                    kw in context.lower()
                    for kw in ['temperature', 'thermal', 'heat equation',
                                'rho_c', 'kappa', 'convection']
                )
                if not is_thermal_ctx:
                    findings.append({
                        'line': i,
                        'severity': 'HIGH',
                        'rule': 'eddy-current-missing-complex',
                        'message': (
                            'FE space in eddy current context without complex=True. '
                            'Frequency-domain eddy current analysis requires complex-'
                            'valued spaces. Add complex=True parameter.'
                        ),
                    })
    return findings


def check_joule_heat_missing_conj(filepath: str, lines: List[str]) -> List[Dict]:
    """MODERATE: Joule heat Q must use Conj(E), not E*E."""
    findings = []
    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith('#'):
            continue
        if re.search(r'InnerProduct\s*\(\s*E\s*,\s*E\s*\)', stripped):
            if 'Conj' not in stripped:
                findings.append({
                    'line': i,
                    'severity': 'MODERATE',
                    'rule': 'joule-heat-missing-conj',
                    'message': (
                        'Joule heat uses InnerProduct(E, E) instead of '
                        'InnerProduct(E, Conj(E)). For complex fields, '
                        'E*E gives a complex number, not real power. '
                        'Use: 0.5 * sigma * InnerProduct(E, Conj(E)).real'
                    ),
                })
    return findings


def check_ngsolve_kelvin_missing_bonus_intorder(filepath: str, lines: List[str]) -> List[Dict]:
    """MODERATE: Kelvin transform bilinear form without bonus_intorder."""
    findings = []
    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith('#'):
            continue
        if re.search(r'dx\s*\(\s*["\']Kelvin["\']', stripped) and \
           'bonus_intorder' not in stripped:
            findings.append({
                'line': i,
                'severity': 'MODERATE',
                'rule': 'ngsolve-kelvin-missing-bonus-intorder',
                'message': (
                    'Integration over Kelvin domain without bonus_intorder. '
                    'The spatially-varying Kelvin Jacobian requires higher '
                    'quadrature for accuracy. '
                    'Use: dx("Kelvin", bonus_intorder=4)'
                ),
            })
    return findings


# ── Shared PEEC/BEM rules (also in radia/rules.py) ───────────

# SHARED: also in radia/rules.py - keep synchronized
def check_bessel_jv_for_sibc(filepath: str, lines: List[str]) -> List[Dict]:
    """CRITICAL: Circular wire SIBC must use iv (modified Bessel), not jv."""
    findings = []
    has_sibc_context = any(
        kw in line.lower()
        for line in lines
        for kw in ['sibc', 'bessel', 'skin_depth', 'impedance_circular',
                    'bessel_impedance']
    )
    if not has_sibc_context:
        return findings

    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith('#'):
            continue
        if re.search(r'from\s+scipy\.special\s+import\s+.*\bjv\b', stripped):
            if not re.search(r'\biv\b', stripped):
                findings.append({
                    'line': i,
                    'severity': 'CRITICAL',
                    'rule': 'bessel-jv-not-iv',
                    'message': (
                        'Circular wire SIBC requires modified Bessel functions '
                        'iv (I0, I1), NOT regular jv (J0, J1). '
                        'jv gives wrong sign on internal inductance. '
                        'Use: from scipy.special import iv'
                    ),
                })
    return findings


# SHARED: also in radia/rules.py - keep synchronized
def check_peec_low_nseg(filepath: str, lines: List[str]) -> List[Dict]:
    """MODERATE: PEEC circular coil should use n_seg >= 32 for coupling accuracy."""
    findings = []
    has_peec_context = any(
        kw in line
        for line in lines
        for kw in ['FastHenryParser', 'PEECBuilder', 'peec', 'n_seg']
    )
    if not has_peec_context:
        return findings

    pattern = re.compile(r'\bn_seg\s*=\s*(\d+)')
    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith('#'):
            continue
        match = pattern.search(stripped)
        if match:
            n_seg = int(match.group(1))
            if n_seg < 32:
                findings.append({
                    'line': i,
                    'severity': 'MODERATE',
                    'rule': 'peec-low-nseg',
                    'message': (
                        f'n_seg={n_seg} may be too coarse for circular coil '
                        f'coupling accuracy. n_seg=16 gives ~26% Delta_L error; '
                        f'use n_seg>=64 for <10% error with magnetic cores.'
                    ),
                })
    return findings


# SHARED: also in radia/rules.py - keep synchronized
def check_efie_v_minus_sign(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: EFIE V term must have positive sign (Lenz's law)."""
    findings = []
    has_efie_context = any(
        kw in line.lower()
        for line in lines
        for kw in ['shieldbemsibc', 'efie', 'v_ll', 'maxwell_slp',
                    'maxwellsinglelayer']
    )
    if not has_efie_context:
        return findings

    patterns = [
        re.compile(r'[-]\s*(?:1j|1\.?j)\s*\*\s*(?:omega|w)\s*\*\s*(?:mu_0|MU_0|mu0)\s*\*\s*(?:V_LL|V_loop|self\._V_LL)'),
        re.compile(r'[-]\s*(?:mu_0|MU_0|mu0)\s*\*\s*(?:omega|w)\s*\*\s*(?:V_LL|V_loop|self\._V_LL)'),
    ]
    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith('#'):
            continue
        for pattern in patterns:
            if pattern.search(stripped):
                findings.append({
                    'line': i,
                    'severity': 'HIGH',
                    'rule': 'efie-v-minus-sign',
                    'message': (
                        'EFIE V_LL term has MINUS sign. The correct EFIE is: '
                        '(Zs*M_LL + jw*mu_0*V_LL)*I = -jw*b. '
                        'A minus sign on V violates Lenz\'s law '
                        '(A_scat would reinforce A_inc instead of opposing it).'
                    ),
                })
                break
    return findings


# SHARED: also in radia/rules.py - keep synchronized
def check_classical_efie_breakdown(filepath: str, lines: List[str]) -> List[Dict]:
    """MODERATE: Classical EFIE 1/kappa^2 causes low-frequency breakdown."""
    findings = []
    has_helmholtz = any(
        kw in line
        for line in lines
        for kw in ['HelmholtzSL', 'HelmholtzSingleLayer']
    )
    if not has_helmholtz:
        return findings

    patterns = [
        re.compile(r'1\s*/\s*(?:kappa|kap)\s*\*\*\s*2'),
        re.compile(r'1\.0\s*/\s*\(?(?:kappa|kap)\s*\*\*\s*2'),
        re.compile(r'1\.0\s*/\s*\(?(?:kappa|kap)\s*\*\s*(?:kappa|kap)'),
    ]
    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith('#'):
            continue
        for pattern in patterns:
            if pattern.search(stripped):
                findings.append({
                    'line': i,
                    'severity': 'MODERATE',
                    'rule': 'classical-efie-breakdown',
                    'message': (
                        'Classical EFIE with 1/kappa^2 scaling causes O(kappa^{-2}) '
                        'condition number blow-up at low frequency. '
                        'Use stabilized formulation: kappa^2 * V_kappa '
                        '(multiply V by kappa^2 instead of dividing). '
                        'Ref: Weggler stabilized EFIE.'
                    ),
                })
                break
    return findings


# SHARED: also in radia/rules.py - keep synchronized
def check_peec_p_over_jw(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: PEEC Loop-Star P/(jw) causes low-frequency breakdown."""
    findings = []
    has_peec_ls = any(
        kw in line
        for line in lines
        for kw in ['loop_star', 'Loop-Star', 'LoopStar', 'Z_SS', 'P_matrix',
                    'NGBEMPEECSolver', 'ngsbem_peec']
    )
    if not has_peec_ls:
        return findings

    patterns = [
        re.compile(r'(?:self\.)?P\s*/\s*\(?\s*(?:1j|1\.?j)\s*\*\s*(?:omega|w)'),
        re.compile(r'(?:self\.)?P\s*/\s*\(?\s*(?:omega|w)\s*\*\s*(?:1j|1\.?j)'),
        re.compile(r'Z_SS\s*=\s*.*(?:self\.)?P\s*/'),
    ]
    in_docstring = False
    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        if '"""' in stripped:
            count = stripped.count('"""')
            if count == 1:
                in_docstring = not in_docstring
            continue
        if in_docstring or stripped.startswith('#'):
            continue
        for pattern in patterns:
            if pattern.search(stripped):
                findings.append({
                    'line': i,
                    'severity': 'HIGH',
                    'rule': 'peec-p-over-jw',
                    'message': (
                        'PEEC Loop-Star P/(jw) causes low-frequency breakdown '
                        '(40-340% error). Use reformulated Schur complement: '
                        'precompute P^{-1}@M_LS, multiply by jw at runtime. '
                        'Or use stabilized mode with Weggler\'s k^2*V_0 block. '
                        'See NGBEMPEECSolver mode="full" or mode="stabilized".'
                    ),
                })
                break
    return findings


def check_axisymmetric_h1_over_r(filepath, lines):
    """Detect naive 1/r weighted H1 form for axisymmetric problems.

    The naive formulation `nu/r * grad(u)*grad(v) * dx` suffers from
    1/r singularity at the axis. The mixed formulation (H1 * VectorL2)
    with Bop(phi) = CF((-1/x*grad(phi)[1], 1/x*grad(phi)[0])) is
    recommended (Schoeberl).
    """
    findings = []
    patterns = [
        # nu/r * grad * grad  or  nu / x * grad * grad
        re.compile(r'(?:nu|reluct)\s*/\s*(?:r_coord|r|x)\s*\*\s*grad'),
        # 1/r * grad * grad  or  1/x * grad * grad
        re.compile(r'1\s*/\s*(?:r_coord|r|x)\s*\*\s*grad\s*\(\s*\w+\s*\)\s*\*\s*grad'),
    ]
    # Only trigger if the file looks axisymmetric
    has_axi_hint = False
    for line in lines:
        if re.search(r'(?:axi|r_coord|半径|radial|1/x\s*\*\s*grad)', line, re.IGNORECASE):
            has_axi_hint = True
            break
    if not has_axi_hint:
        return findings

    in_docstring = False
    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        if '"""' in stripped:
            count = stripped.count('"""')
            if count == 1:
                in_docstring = not in_docstring
            continue
        if in_docstring or stripped.startswith('#'):
            continue
        for pattern in patterns:
            if pattern.search(stripped):
                findings.append({
                    'line': i,
                    'severity': 'MODERATE',
                    'rule': 'axisymmetric-naive-1-over-r',
                    'message': (
                        'Naive 1/r weighted H1 formulation for axisymmetric problems '
                        'suffers from singularity at r=0. Consider using the mixed '
                        'formulation (H1 * VectorL2) with Bop(phi) = CF((-1/x*grad(phi)[1], '
                        '1/x*grad(phi)[0])). See ngsolve_usage topic "axisymmetric".'
                    ),
                })
                break
    return findings


# SHARED: also in radia/rules.py - keep synchronized
def check_ngsbem_volume_mesh(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: NGSBEM BEM requires surface-only mesh, not volume mesh."""
    findings = []
    has_bem = any(
        kw in line
        for line in lines
        for kw in ['LaplaceSL', 'HelmholtzSL', 'HDivSurface', 'bem_inductance']
    )
    if not has_bem:
        return findings

    has_box = False
    has_glue = False
    box_line = 0
    for i, line in enumerate(lines, 1):
        if 'Box(' in line:
            has_box = True
            box_line = i
        if 'Glue(' in line:
            has_glue = True

    if has_box and not has_glue:
        findings.append({
            'line': box_line,
            'severity': 'HIGH',
            'rule': 'ngsbem-volume-mesh',
            'message': (
                'BEM with Box() creates a volume mesh. HDivSurface on volume '
                'mesh includes interior edges -> singular SL matrix (cond ~1e17). '
                'Use: OCCGeometry(Glue(box.faces)) for surface-only mesh.'
            ),
        })
    return findings


def check_scattered_field_sibc_missing_nxH(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: Scattered-field SIBC needs BOTH RHS terms: sibc(A_inc) + n x H_inc."""
    findings = []
    # Only check files with scattered-field SIBC context
    has_scat_sibc = any(
        kw in line.lower()
        for line in lines
        for kw in ['scattered', 'a_scat', 'a_inc']
    )
    has_sibc = any('sibc' in line.lower() for line in lines)
    if not (has_scat_sibc and has_sibc):
        return findings

    # Look for SIBC RHS with A_inc but without n x H_inc
    has_sibc_A_inc_rhs = False
    has_nxH_term = False
    sibc_rhs_line = 0
    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith('#'):
            continue
        # Detect SIBC RHS: something like sibc_coeff * A_inc * v * ds
        if re.search(r'sibc.*A_inc.*ds|A_inc.*sibc.*ds', stripped, re.IGNORECASE):
            has_sibc_A_inc_rhs = True
            sibc_rhs_line = i
        # Detect n x H_inc term: nxH or n_cross_H or curl(A_inc)
        if re.search(r'nxH|n_cross_H|curl.*A_inc|H_inc.*ds', stripped, re.IGNORECASE):
            has_nxH_term = True

    if has_sibc_A_inc_rhs and not has_nxH_term:
        findings.append({
            'line': sibc_rhs_line,
            'severity': 'HIGH',
            'rule': 'scattered-sibc-missing-nxH',
            'message': (
                'Scattered-field SIBC RHS has -(jw/Z_s)*<A_inc, v> but is missing '
                'the second term -<n x H_inc, v>. Both terms are required. '
                'Missing this term causes a factor-of-3 error. '
                'Total-field formulation (J_source in volume) does not have this issue.'
            ),
        })
    return findings



# ALL_RULES is defined at the end of the file (after all function definitions)


def check_ngsbem_missing_curvaturesafety(filepath: str, lines: List[str]) -> List[Dict]:
    """MODERATE: OCC mesh for BEM should set curvaturesafety to prevent degenerate elements."""
    findings = []
    has_bem = any(
        kw in line
        for line in lines
        for kw in ['LaplaceSL', 'LaplaceDL', 'HelmholtzSL', 'HDivSurface']
    )
    if not has_bem:
        return findings

    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith('#'):
            continue
        if re.search(r'\.GenerateMesh\s*\(', stripped):
            if 'curvaturesafety' not in stripped:
                findings.append({
                    'line': i,
                    'severity': 'MODERATE',
                    'rule': 'ngsbem-missing-curvaturesafety',
                    'message': (
                        'GenerateMesh() for BEM without curvaturesafety parameter. '
                        'OCC meshing can produce degenerate elements on curved surfaces, '
                        'causing zero eigenvalues in BEM operators. '
                        'Add curvaturesafety=1: geo.GenerateMesh(maxh=h, curvaturesafety=1)'
                    ),
                })
    return findings


def check_ngsbem_taskmanager_reproducibility(filepath: str, lines: List[str]) -> List[Dict]:
    """LOW: TaskManager with BEM causes non-deterministic results."""
    findings = []
    has_bem = any(
        kw in line
        for line in lines
        for kw in ['LaplaceSL', 'LaplaceDL', 'HelmholtzSL']
    )
    if not has_bem:
        return findings

    has_taskmanager = False
    has_setnumthreads_1 = False
    tm_line = 0
    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith('#'):
            continue
        if re.search(r'TaskManager\s*\(', stripped):
            has_taskmanager = True
            tm_line = i
        if re.search(r'SetNumThreads\s*\(\s*1\s*\)', stripped):
            has_setnumthreads_1 = True

    if has_taskmanager and not has_setnumthreads_1:
        findings.append({
            'line': tm_line,
            'severity': 'LOW',
            'rule': 'ngsbem-taskmanager-nondeterministic',
            'message': (
                'TaskManager with BEM operators causes non-deterministic results '
                '(~1.4% fluctuation at 8 threads) due to parallel summation order. '
                'For reproducible results, remove TaskManager or add SetNumThreads(1). '
                'See ngsbem_inductance topic "known_limitations".'
            ),
        })
    return findings


def check_efie_sibc_use_scalar_bie(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: EFIE-SIBC is wrong for finite Z_s. Use Scalar BIE + SIBC instead."""
    findings = []
    has_efie = any(
        re.search(r'LaplaceSL.*HDivSurface|HDivSurface.*LaplaceSL', line)
        for line in lines
    )
    has_sibc_zs = any(
        re.search(r'Z_s|surface.impedance|sibc', line, re.IGNORECASE)
        for line in lines
    )
    if not (has_efie and has_sibc_zs):
        return findings

    # Check if already using ScalarBIESIBCSolver
    has_scalar_bie = any('ScalarBIESIBCSolver' in line for line in lines)
    if has_scalar_bie:
        return findings

    for i, line in enumerate(lines, 1):
        if re.search(r'Z_s.*LaplaceSL|LaplaceSL.*Z_s', line):
            findings.append({
                'line': i,
                'severity': 'HIGH',
                'rule': 'efie-sibc-use-scalar-bie',
                'message': (
                    'EFIE with HDivSurface + Z_s is wrong for finite Z_s '
                    '(SL eigenvalue R/3 causes factor-of-3 error). '
                    'Use ScalarBIESIBCSolver (H1 scalar potential BIE) instead: '
                    'from radia.bem_sibc_solver import ScalarBIESIBCSolver. '
                    'It uses existing LaplaceSL/DL with H1 and surface '
                    'Laplacian for SIBC coupling. Validated <0.1% on sphere.'
                ),
            })
            break
    return findings


# ============================================================
# Radia C++ API rules (from mcp-server-radia)
# ============================================================


def check_objbckg_no_lambda(filepath: str, lines: List[str]) -> List[Dict]:
    """CRITICAL: rad.ObjBckg() must receive a callable, not a list."""
    findings = []
    pattern = re.compile(r'rad\.ObjBckg\s*\(\s*\[')
    for i, line in enumerate(lines, 1):
        if pattern.search(line):
            findings.append({
                'line': i, 'severity': 'CRITICAL',
                'rule': 'objbckg-needs-callable',
                'message': 'rad.ObjBckg() requires a callable, not a list. '
                           'Use: rad.ObjBckg(lambda p: [Bx, By, Bz])',
            })
    return findings


def check_missing_utidelall(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: Scripts should call rad.UtiDelAll() for cleanup."""
    findings = []
    has_radia = any('import radia' in line or 'import rad' in line for line in lines)
    if not has_radia:
        return findings
    basename = os.path.basename(filepath)
    if basename.startswith('__') or basename in ('radia.py',):
        return findings
    has_main = any("__name__" in line and "__main__" in line for line in lines)
    if not has_main:
        return findings
    if not any('UtiDelAll' in line for line in lines):
        findings.append({
            'line': len(lines), 'severity': 'HIGH',
            'rule': 'missing-utidelall',
            'message': 'Script imports radia but does not call rad.UtiDelAll().',
        })
    return findings


def check_hardcoded_absolute_paths(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: No hardcoded absolute paths in sys.path."""
    findings = []
    patterns = [
        re.compile(r'sys\.path\.insert\s*\(\s*\d+\s*,\s*r?["\'][A-Za-z]:\\'),
        re.compile(r'(?:repo_root|work_dir|base_dir)\s*=\s*r?["\'][A-Za-z]:\\'),
    ]
    for i, line in enumerate(lines, 1):
        if line.strip().startswith('#'):
            continue
        for p in patterns:
            if p.search(line.strip()):
                findings.append({
                    'line': i, 'severity': 'HIGH',
                    'rule': 'hardcoded-absolute-path',
                    'message': 'Use os.path.dirname(__file__) for relative paths.',
                })
                break
    return findings


def check_removed_fldunits(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: rad.FldUnits() has been removed."""
    findings = []
    pattern = re.compile(r'(?:rad\.)?FldUnits\s*\(')
    for i, line in enumerate(lines, 1):
        if not line.strip().startswith('#') and pattern.search(line):
            findings.append({
                'line': i, 'severity': 'HIGH', 'rule': 'removed-fldunits',
                'message': 'rad.FldUnits() removed. Radia always uses meters.',
            })
    return findings


def check_removed_fldbatch(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: FldBatch/FldA/FldPhi removed."""
    findings = []
    pattern = re.compile(r'(?:rad\.)?(?:FldBatch|FldA|FldPhi)\s*\(')
    for i, line in enumerate(lines, 1):
        if not line.strip().startswith('#') and pattern.search(line):
            findings.append({
                'line': i, 'severity': 'HIGH', 'rule': 'removed-fldbatch',
                'message': 'Use rad.Fld(obj, field_type, points) with (N,3) for batch.',
            })
    return findings


def check_removed_solver_apis(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: Reject retired solver APIs and stale SolverConfig keywords."""
    findings = []
    kernel_apis = [
        'SetHACApKParams', 'SetHMatrixEpsilon', 'SetBiCGSTABTol',
        'GetBiCGSTABTol', 'SetHMatrixFieldEval',
    ]
    nonlinear_apis = [
        'SetRelaxParam', 'GetRelaxParam', 'SetNewtonMethod', 'GetNewtonMethod',
        'SetNewtonDamping', 'GetNewtonDampingStats',
    ]
    kernel_pattern = re.compile(
        r'(?:rad\.)?(?:' + '|'.join(kernel_apis) + r')\s*\('
    )
    nonlinear_pattern = re.compile(
        r'(?:rad\.)?(?:' + '|'.join(nonlinear_apis) + r')\s*\('
    )
    retired_config_pattern = re.compile(
        r'(?:rad\.)?SolverConfig\s*\([^)]*\b'
        r'(?:hacapk_eps|hacapk_leaf|hacapk_eta|hmatrix_eps|bicgstab_tol)\s*='
    )
    for i, line in enumerate(lines, 1):
        if line.strip().startswith('#'):
            continue
        if kernel_pattern.search(line) or retired_config_pattern.search(line):
            findings.append({
                'line': i, 'severity': 'HIGH', 'rule': 'removed-solver-api',
                'message': (
                    'Legacy compact-magnetostatic/Krylov controls are retired. '
                    'Configure HDiv, PEEC, or BEM compression through that solver API.'
                ),
            })
        elif nonlinear_pattern.search(line):
            findings.append({
                'line': i, 'severity': 'HIGH', 'rule': 'removed-solver-api',
                'message': 'Use rad.SolverConfig(**kwargs) and rad.GetSolverConfig().',
            })
    return findings


def check_docstring_hardcoded_mm(filepath: str, lines: List[str]) -> List[Dict]:
    """MODERATE: Docstrings should not hardcode 'in mm'."""
    findings = []
    if 'src/radia' not in filepath.replace('\\', '/') and 'src\\radia' not in filepath:
        return findings
    in_docstring = False
    for i, line in enumerate(lines, 1):
        if '"""' in line.strip():
            if line.strip().count('"""') == 1:
                in_docstring = not in_docstring
        if in_docstring and re.search(r'\bpoint.*\bin mm\b', line.strip()):
            findings.append({
                'line': i, 'severity': 'MODERATE', 'rule': 'docstring-hardcoded-mm',
                'message': 'Use "in constructor length units" for unit-agnostic API.',
            })
    return findings


def check_build_release_path(filepath: str, lines: List[str]) -> List[Dict]:
    """LOW: build/Release path import may be outdated."""
    findings = []
    pattern = re.compile(r"sys\.path\.insert.*['\"].*build/Release['\"]")
    for i, line in enumerate(lines, 1):
        if pattern.search(line):
            findings.append({
                'line': i, 'severity': 'LOW', 'rule': 'build-release-path',
                'message': 'Prefer src/radia which has the latest binaries.',
            })
    return findings


def check_scattered_eddy_missing_a0(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: Scattered-field eddy/Joule loss must include the background A0 in E.

    In a scattered-field harmonic eddy solve the FE unknown gfA is the SCATTERED
    vector potential; the total electric field is E = -jw (A0 + gfA + grad(Phi)).
    Computing the Joule loss from the scattered gfA alone (omitting the background
    A0) overestimates the loss by ~10x (observed in an independent induction-
    heating cross-validation). Coil-source (total-field) problems carry no
    background A0 and correctly pass A0=None -- those files do not define an A0
    and are therefore not flagged.
    """
    findings = []
    text = '\n'.join(lines)
    # The helper's own definition legitimately references gfA without A0 in the
    # "A0 is None" branch; never flag the definition file itself.
    if 'def joule_loss_density' in text:
        return findings
    low = text.lower()
    if not any(k in low for k in ['joule', 'eddy', 'induction_heating', 'loss_density']):
        return findings
    # Only a scattered/background-field problem needs A0 added; coil/total-field
    # sources correctly use A0=None and must NOT be flagged. Gate on a background
    # marker so those files are excluded.
    has_background = (
        any(k in low for k in ['background', 'applied field', 'incident', 'scattered'])
        or any(re.search(r'\bB0\b', l) for l in lines)
    )
    if not has_background:
        return findings

    in_docstring = False
    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        # Skip triple-quoted blocks -- knowledge/doc text shows these patterns as prose.
        if stripped.count('"""') == 1 or stripped.count("'''") == 1:
            in_docstring = not in_docstring
            continue
        if in_docstring or stripped.startswith('#'):
            continue
        code = stripped.split('#', 1)[0]   # ignore inline comments
        # Unambiguous signature: the joule_loss_density() helper called WITHOUT A0=.
        # (Hand-rolled E = -jw(gfA + grad(Phi)) is deliberately NOT flagged: it is
        # correct for total-field/coil-source formulations where gfA is already the
        # total potential, which is statically indistinguishable from scattered.)
        if re.search(r'joule_loss_density\s*\(', code):
            call_ctx = ' '.join(l.split('#', 1)[0] for l in lines[i - 1:i + 2])
            if not re.search(r'\bA0\s*=', call_ctx):
                findings.append({
                    'line': i, 'severity': 'HIGH',
                    'rule': 'scattered-eddy-missing-a0',
                    'message': (
                        'joule_loss_density() called without A0= but this file sets a '
                        'background/applied field. In a scattered-field eddy solve the '
                        'FE unknown gfA is the SCATTERED potential, so the total '
                        'E = -jw (A0 + gfA + grad(Phi)); omitting A0 overestimates the '
                        'Joule loss by ~10x. Pass A0=<background> (A0=None ONLY for '
                        'coil/total-field sources, where gfA is already total).'
                    ),
                })
    return findings


# NOTE: the "direct solver on 3D HCurl hits the ~84k order-2 ceiling" footgun is
# NOT a static lint -- "uses a direct solver" is correct the vast majority of the
# time (every 2D-H1 and small 3D solve), so a static check floods correct library
# code with false positives. It is instead enforced as a RUNTIME guard in solve.py
# (`_direct_inverse`), which turns the cryptic "bad array new length" overflow into
# an actionable "switch to CG+BDDC" message. Right guard, right layer.


def check_gridfunction_set_definedon_in_loop(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: GridFunction.Set(..., definedon=...) inside a loop.

    NGSolve's GridFunction.Set() ZEROS the entire vector before projecting onto the
    definedon region. Looping over boundaries/materials and calling Set() once per
    region therefore leaves ONLY the LAST region's data -- every prior region is wiped
    (e.g. a capacitor whose conductor potentials are set in a loop collapses to V=0,
    energy 0). Set ALL regions in ONE call via mesh.BoundaryCF / mesh.MaterialCF.
    """
    findings = []
    if not any(('ngsolve' in l) or ('GridFunction' in l) for l in lines):
        return findings
    set_pat = re.compile(r'\b(\w+)\.Set\s*\(.*\bdefinedon\s*=')
    loop_pat = re.compile(r'^(\s*)(?:for|while)\b')
    loop_stack = []          # (indent, lineno) of currently-open loops
    in_docstring = False
    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.count('"""') == 1 or stripped.count("'''") == 1:
            in_docstring = not in_docstring
            continue
        if in_docstring or stripped.startswith('#'):
            continue
        indent = len(line) - len(line.lstrip())
        while loop_stack and indent <= loop_stack[-1][0]:
            loop_stack.pop()
        code = stripped.split('#', 1)[0]
        sm = set_pat.search(code)
        if sm and loop_stack:
            findings.append({
                'line': i, 'severity': 'HIGH',
                'rule': 'ngsolve-set-definedon-in-loop',
                'message': (
                    f"'{sm.group(1)}.Set(..., definedon=...)' is inside the loop at "
                    f"line {loop_stack[-1][1]}. GridFunction.Set() ZEROS the whole "
                    "vector first, so each iteration wipes the previous region -- only "
                    "the LAST survives (e.g. conductor potentials collapse to 0). Set "
                    "ALL regions in ONE call: gf.Set(mesh.BoundaryCF({...}, default=0), "
                    "definedon=mesh.Boundaries('a|b'))."
                ),
            })
        lm = loop_pat.match(line)
        if lm:
            loop_stack.append((len(lm.group(1)), i))
    return findings


def check_curve_after_vol_import(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: mesh.Curve() after importing a .vol FLATTENS the baked-in curving to facets.

    A .vol exported with a curve order (Cubit ``export netgen "<f>" order N``) stores
    the high-order curved-node data, and ``ngsolve.Mesh("<f>.vol")`` loads it
    faithfully. But an imported .vol carries NO CAD geometry, so ``mesh.Curve(k)`` --
    which rebuilds the curved nodes by projecting them onto the ATTACHED geometry --
    has nothing to project onto and recomputes them FLAT, collapsing the surface to
    the inscribed polytope. (For HEX this is the only curving path at all: NGSolve
    never builds hex from geometry, so the curving must come from the export.)

    Verified empirically on NGSolve 6.2.2604 (a Cubit order-2 curved hex ball):
    ``Curve(1) == Curve(2) == Curve(3)`` all degrade the DtN eigenvalue floor
    1.6e-5 -> 5.3e-3 (~320x) and the Gamma surface deviation 5.7e-5 -> 2.3e-2 (~400x,
    surface area collapsing to 96.7% of 4*pi = the flat-facet area) -- IDENTICAL for
    every k, because the requested order is irrelevant once the geometry is gone.
    Worse, ``GetCurveOrder()`` then reports the NEW order while the geometry is
    actually flat, so the loss is SILENT. (Contrast: on an OCCGeometry mesh, which
    HAS faces to project onto, Curve(k) correctly IMPROVES the surface -- so this
    footgun is specific to geometry-less imported meshes, hence the .vol gate.)

    Do NOT .Curve() an imported .vol: set the FE order via ``H1(mesh, order=k)``; to
    change the GEOMETRY order, re-export the .vol at a higher order.
    """
    findings = []
    msg_sep = (
        "'{var}.Curve()' is called on a mesh imported from a .vol. An imported .vol "
        "carries no CAD geometry, so mesh.Curve() rebuilds the curved nodes FLAT and "
        "collapses the surface to facets -- it does NOT add curving (and NGSolve has "
        "no hex-from-geometry path, so hex curving must come from the export). "
        "Verified on NGSolve 6.2.2604: Curve(k) degrades the DtN floor ~320x "
        "(1.6e-5 -> 5.3e-3) for every k, while GetCurveOrder() misreports success -- "
        "a SILENT accuracy trap. Remove the .Curve(); set the FE order with "
        "H1(mesh, order=k); bake the geometry order into the .vol on export "
        "(Cubit: export netgen ... order N)."
    )
    msg_chain = (
        '.Curve() chained onto Mesh(".vol"). An imported .vol has no CAD geometry, so '
        "mesh.Curve() rebuilds the curved nodes FLAT (surface -> facets); it does not "
        "add curving and silently degrades accuracy (~320x on a curved hex ball, "
        "NGSolve 6.2.2604), with GetCurveOrder() misreporting success. Drop the "
        ".Curve(); set order via H1(mesh, order=k); bake the geometry order into the "
        ".vol on export."
    )

    pathvar_pat = re.compile(r'(\w+)\s*=\s*[rfb]*["\'][^"\']*\.vol\b', re.IGNORECASE)
    chained_pat = re.compile(r'Mesh\s*\(\s*[^)]*\.vol\b[^)]*\)\s*\.\s*Curve\s*\(', re.IGNORECASE)
    vol_literal_pat = re.compile(
        r'(?:(\w+)\s*=\s*)?(?:\w+\.)?Mesh\s*\(\s*[^)]*\.vol\b', re.IGNORECASE)
    mesh_of_var_pat = re.compile(r'(?:(\w+)\s*=\s*)?(?:\w+\.)?Mesh\s*\(\s*(\w+)\s*\)')
    load_pat = re.compile(r'(\w+)\s*\.\s*Load\s*\(\s*[^)]*\.vol\b', re.IGNORECASE)
    curve_pat = re.compile(r'\b(\w+)\s*\.\s*Curve\s*\(')

    path_vars = set()
    vol_mesh_vars = set()

    def _iter_code():
        in_doc = False
        for idx, ln in enumerate(lines, 1):
            s = ln.strip()
            if s.count('"""') == 1 or s.count("'''") == 1:
                in_doc = not in_doc
                continue
            if in_doc or s.startswith('#'):
                continue
            yield idx, s.split('#', 1)[0]

    # pass 1: identify .vol-imported mesh variables (+ flag chained Mesh(".vol").Curve())
    for i, code in _iter_code():
        pv = pathvar_pat.search(code)
        if pv:
            path_vars.add(pv.group(1))
        if chained_pat.search(code):
            findings.append({'line': i, 'severity': 'HIGH',
                             'rule': 'ngsolve-curve-after-vol-import', 'message': msg_chain})
            continue
        ml = vol_literal_pat.search(code)
        if ml:
            if ml.group(1):
                vol_mesh_vars.add(ml.group(1))
            continue
        mv = mesh_of_var_pat.search(code)
        if mv and mv.group(2) in path_vars and mv.group(1):
            vol_mesh_vars.add(mv.group(1))
        lo = load_pat.search(code)
        if lo:
            vol_mesh_vars.add(lo.group(1))

    if not vol_mesh_vars:
        return findings

    # pass 2: flag <var>.Curve() for any var that was loaded from a .vol
    for i, code in _iter_code():
        cm = curve_pat.search(code)
        if cm and cm.group(1) in vol_mesh_vars:
            findings.append({'line': i, 'severity': 'HIGH',
                             'rule': 'ngsolve-curve-after-vol-import',
                             'message': msg_sep.format(var=cm.group(1))})
    return findings


def check_coil_scalar_potential_as_lift(filepath: str, lines: List[str]) -> List[Dict]:
    """CRITICAL: a coil's magnetic SCALAR potential must not drive a reduced-Omega FEM.

    A current loop's magnetic scalar potential Omega_s (= -rad.Fld(coil, 'phi'),
    the solid-angle potential) is MULTIVALUED -- it jumps by the ampere-turns NI
    across the loop's spanning surface.  For a coil that links the iron magnetic
    circuit (any real electromagnet), using Omega_s as a Dirichlet LIFT on the
    iron-air interface (``gf.Set(Omega_s, BND, ...)``) is ill-posed: the iron
    magnetises but the air/gap field comes out ~100% wrong, MESH-INDEPENDENTLY.
    (The bug hides because the only validations were UNIFORM sources, whose
    Omega_s = H0*z is single-valued.)

    The production accelerator-magnet pipeline (``calc_accel_magnet.py``) does it
    correctly: the source enters as the single-valued VECTOR field
    ``H_s = rad.RadiaField(coil, 'h')`` in a VOLUME integral,

        a += mu * grad(u) * grad(v) * dx
        f += mu * H_s * grad(v) * dx          # H = H_s - grad(Omega)

    Likewise ``ScalarPotentialSolver.solve_total_reduced_potential`` (the
    Omega_s-lift two-scalar method) is broken for coil sources -- prefer
    ``solve_single_potential`` / ``solve_nonlinear_newton`` (vector H_s).
    """
    findings = []
    # The owner of the (deprecated/broken) method legitimately defines and
    # dispatches to it -- do not self-flag that library file.
    if os.path.basename(filepath) == 'scalar_potential_solver.py':
        return findings
    text = '\n'.join(lines)
    has_phi = bool(re.search(r"(?:RadiaField|rad\.Fld)\s*\([^)]*['\"]phi['\"]", text))
    calls_total_reduced = 'solve_total_reduced_potential' in text
    if not (has_phi or calls_total_reduced):
        return findings

    # Track variables bound to a coil scalar potential: direct + one-level aliases.
    phi_assign = re.compile(
        r"^(\w+)\s*=\s*[^#\n]*(?:RadiaField|rad\.Fld)\s*\([^)]*['\"]phi['\"]")
    phi_vars = set()
    in_doc = False
    for line in lines:
        s = line.strip()
        if s.count('"""') == 1 or s.count("'''") == 1:
            in_doc = not in_doc
            continue
        if in_doc or s.startswith('#'):
            continue
        m = phi_assign.search(s.split('#', 1)[0])
        if m:
            phi_vars.add(m.group(1))
    # one-level aliases:  Om = -1.0 * phi   /   Omega_s = -phi
    if phi_vars:
        alias = re.compile(
            r"^(\w+)\s*=\s*[-\d.\s()*+]*\b(" + '|'.join(map(re.escape, phi_vars)) + r")\b")
        in_doc = False
        for line in lines:
            s = line.strip()
            if s.count('"""') == 1 or s.count("'''") == 1:
                in_doc = not in_doc
                continue
            if in_doc or s.startswith('#'):
                continue
            m = alias.search(s.split('#', 1)[0])
            if m:
                phi_vars.add(m.group(1))

    # GridFunctions filled FROM a coil scalar potential (`gf.Set(<...phi...>)`):
    # projecting phi into a GridFunction is fine, but the GridFunction then
    # carries the multivalued potential, so a LATER `.Set(gf, BND, ...)` lift is
    # the bug.  Track those GridFunctions too (the projection line has no BND, so
    # it is not itself flagged).
    set_from_phi = re.compile(r"^(\w+)\.Set\s*\(([^)]*)\)")
    in_doc = False
    for line in lines:
        s = line.strip()
        if s.count('"""') == 1 or s.count("'''") == 1:
            in_doc = not in_doc
            continue
        if in_doc or s.startswith('#'):
            continue
        code = s.split('#', 1)[0]
        m = set_from_phi.search(code)
        if m and 'BND' not in code and 'Boundaries(' not in code:
            arg = m.group(2)
            if re.search(r"(?:RadiaField|rad\.Fld)\s*\([^)]*['\"]phi['\"]", arg) \
                    or any(re.search(r'\b' + re.escape(pv) + r'\b', arg) for pv in phi_vars):
                phi_vars.add(m.group(1))

    in_doc = False
    for i, line in enumerate(lines, 1):
        s = line.strip()
        if s.count('"""') == 1 or s.count("'''") == 1:
            in_doc = not in_doc
            continue
        if in_doc or s.startswith('#'):
            continue
        code = s.split('#', 1)[0]
        # (a) a coil scalar potential used as a boundary Dirichlet lift
        is_bnd_set = ('.Set(' in code and ('BND' in code or 'Boundaries(' in code))
        if is_bnd_set:
            uses_phi = bool(re.search(r"(?:RadiaField|rad\.Fld)\s*\([^)]*['\"]phi['\"]", code)) \
                or any(re.search(r'\b' + re.escape(pv) + r'\b', code) for pv in phi_vars)
            if uses_phi:
                findings.append({
                    'line': i,
                    'severity': 'CRITICAL',
                    'rule': 'coil-scalar-potential-as-lift',
                    'message': (
                        "A coil magnetic SCALAR potential (rad.Fld/RadiaField 'phi') is "
                        "used as a Dirichlet boundary LIFT (.Set(..., BND, ...)). The "
                        "coil scalar potential is MULTIVALUED (solid angle, jumps by NI) "
                        "for any current loop linking iron, so this lift is ill-posed -> "
                        "~100% wrong air/gap field (mesh-independent). Drive the reduced-"
                        "Omega FEM with the VECTOR source H_s = rad.RadiaField(coil,'h') in "
                        "a VOLUME integral (H = H_s - grad(Omega)), as calc_accel_magnet.py "
                        "does: f += mu * H_s * grad(v) * dx."
                    ),
                })
        # (b) the broken Omega_s-lift two-scalar method, called (not defined)
        if 'solve_total_reduced_potential' in code and not code.lstrip().startswith('def '):
            findings.append({
                'line': i,
                'severity': 'HIGH',
                'rule': 'coil-scalar-potential-as-lift',
                'message': (
                    "solve_total_reduced_potential() (Omega_s-lift two-scalar) is broken "
                    "for coil sources: it couples the coil via its MULTIVALUED scalar "
                    "potential Omega_s on the iron boundary -> ~100% wrong gap field "
                    "(validated only on uniform sources). Use solve_single_potential / "
                    "solve_nonlinear_newton, which use the vector H_s = RadiaField(coil,'h')."
                ),
            })
    return findings


# ── Radia API-drift / example-bitrot rules (2026-06, examples->docs sweep) ──
# Mirrors the bug_patterns catalog entries radia-magnetization-tesla-not-apm,
# rad-solve-return-4-tuple, radia-constructor-arity-drift.

def check_solve_result_bad_index(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: rad.Solve returns a 4-tuple [residual, _, _, iterations]."""
    findings = []
    solve_vars = set()
    for line in lines:
        if line.strip().startswith('#'):
            continue
        m = re.search(r'(\w+)\s*=\s*(?:rad|rd)?\.?Solve\s*\(', line)
        if m:
            solve_vars.add(m.group(1))
    if not solve_vars:
        return findings
    # only [4] is always wrong (out of range); [3] is the legitimate iteration count
    idx_pat = re.compile(r'(\w+)\s*\[\s*4\s*\]')
    for i, line in enumerate(lines, 1):
        if line.strip().startswith('#'):
            continue
        for m in idx_pat.finditer(line):
            if m.group(1) in solve_vars:
                findings.append({
                    'line': i, 'severity': 'HIGH', 'rule': 'solve-result-bad-index',
                    'message': (
                        'rad.Solve returns a 4-tuple [residual, _, _, iterations] '
                        '(indices 0-3). Index [4] is out of range. Use [0] for the '
                        'convergence residual and [3] for the iteration count.'
                    ),
                })
    return findings


def check_objarccur_missing_axis(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: ObjArcCur needs 8 args incl. man_auto + axis ('man','z')."""
    findings = []
    for i, line in enumerate(lines, 1):
        s = line.strip()
        if s.startswith('#') or 'ObjArcCur(' not in s:
            continue
        # only judge COMPLETE single-line calls (the call closes on this line);
        # multi-line calls may carry 'man','z' on a later line -> skip them
        tail = s[s.index('ObjArcCur(') + len('ObjArcCur('):]
        if ')' not in tail:
            continue
        # a correct call carries a quoted axis ('x'/'y'/'z')
        if not re.search(r"""['"][xyz]['"]""", s):
            findings.append({
                'line': i, 'severity': 'HIGH', 'rule': 'objarccur-missing-axis',
                'message': (
                    "ObjArcCur takes 8 args: (center, radii, phi, h, nseg, "
                    "man_auto, axis, j). Old 6-arg calls omit man_auto + axis -> "
                    "TypeError. Add e.g. 'man','z' before the current density."
                ),
            })
    return findings


def check_matsatisofrm_multiple_lists(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: MatSatIsoFrm takes ONE nested list, not 3 positional lists."""
    findings = []
    pat = re.compile(r'MatSatIsoFrm\(\s*\[[^\[\]]*\]\s*,\s*\[')
    for i, line in enumerate(lines, 1):
        if not line.strip().startswith('#') and pat.search(line):
            findings.append({
                'line': i, 'severity': 'HIGH', 'rule': 'matsatisofrm-multiple-lists',
                'message': (
                    'MatSatIsoFrm takes a SINGLE nested list '
                    '[[ksi1,ms1],[ksi2,ms2],...], not 3 positional lists -> '
                    'TypeError. Wrap them: MatSatIsoFrm([[..],[..],[..]]).'
                ),
            })
    return findings


def check_magnetization_looks_tesla(filepath: str, lines: List[str]) -> List[Dict]:
    """LOW: magnetization vector looks like Br (Tesla); Radia uses A/m."""
    findings = []
    ctor = re.compile(r'Obj(?:CylMag|Hexahedron|Wedge|RecMag)\s*\(')
    vec = re.compile(r'\[\s*(-?\d+\.?\d*)\s*,\s*(-?\d+\.?\d*)\s*,\s*(-?\d+\.?\d*)\s*\]')
    for i, line in enumerate(lines, 1):
        s = line.strip()
        if s.startswith('#') or not ctor.search(s):
            continue
        for m in vec.finditer(s):
            comps = [abs(float(m.group(k))) for k in (1, 2, 3)]
            nonzero = [c for c in comps if c > 0]
            # one nonzero component in the typical remanence band 0.8-1.6 T
            if len(nonzero) == 1 and 0.8 <= nonzero[0] <= 1.6:
                findings.append({
                    'line': i, 'severity': 'LOW', 'rule': 'magnetization-looks-tesla',
                    'message': (
                        'Magnetization {} looks like Br in Tesla. Radia uses M in '
                        'A/m: M = Br/mu_0 (Br=1.0 T -> 7.96e5 A/m). Verify (a '
                        'Tesla value gives ~zero field).'.format(m.group(0))
                    ),
                })
    return findings


def check_cp932_stdout_rewrap(filepath: str, lines: List[str]) -> List[Dict]:
    """LOW: codecs.getwriter(sys.stdout.buffer) breaks Jupyter / is unneeded."""
    findings = []
    for i, line in enumerate(lines, 1):
        s = line.strip()
        if s.startswith('#'):
            continue
        if 'codecs.getwriter' in s and re.search(r'sys\.std(?:out|err)\.buffer', s):
            findings.append({
                'line': i, 'severity': 'LOW', 'rule': 'cp932-stdout-rewrap',
                'message': (
                    'Re-wrapping sys.stdout via codecs.getwriter(...buffer) breaks '
                    'in Jupyter (OutStream has no .buffer) and is an unneeded cp932 '
                    'workaround. Remove it; keep output ASCII instead.'
                ),
            })
    return findings


# All rules in execution order
ALL_RULES = [
    # Radia API-drift / example-bitrot
    check_solve_result_bad_index,
    check_objarccur_missing_axis,
    check_matsatisofrm_multiple_lists,
    check_magnetization_looks_tesla,
    check_cp932_stdout_rewrap,
    # NGSolve FEM rules
    check_hcurl_missing_nograds,
    check_ngsolve_precond_after_assemble,
    check_ngsolve_missing_trace_bem,
    check_ngsolve_overwrite_xyz,
    check_ngsolve_vec_assign,
    check_ngsolve_dim2_occ,
    check_ngsolve_cg_on_saddle_point,
    check_ngsolve_vectorh1_for_em,
    check_ngsolve_pinvit_no_projection,
    check_eddy_current_missing_complex,
    check_joule_heat_missing_conj,
    check_scattered_eddy_missing_a0,
    check_gridfunction_set_definedon_in_loop,
    check_ngsolve_kelvin_missing_bonus_intorder,
    check_curve_after_vol_import,
    check_axisymmetric_h1_over_r,
    check_coil_scalar_potential_as_lift,
    # Shared PEEC/BEM rules
    check_bessel_jv_for_sibc,
    check_peec_low_nseg,
    check_efie_v_minus_sign,
    check_classical_efie_breakdown,
    check_peec_p_over_jw,
    # BEM mesh rules
    check_ngsbem_volume_mesh,
    # FEM-SIBC rules
    check_scattered_field_sibc_missing_nxH,
    # Scalar BIE recommendation
    check_efie_sibc_use_scalar_bie,
    # BEM mesh/assembly pitfalls
    check_ngsbem_missing_curvaturesafety,
    check_ngsbem_taskmanager_reproducibility,
    # Radia C++ API rules
    check_objbckg_no_lambda,
    check_missing_utidelall,
    check_hardcoded_absolute_paths,
    check_removed_fldunits,
    check_removed_fldbatch,
    check_removed_solver_apis,
    check_docstring_hardcoded_mm,
    check_build_release_path,
]
