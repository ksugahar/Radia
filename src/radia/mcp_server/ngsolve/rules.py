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


# All rules in execution order
ALL_RULES = [
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
    check_ngsolve_kelvin_missing_bonus_intorder,
    check_axisymmetric_h1_over_r,
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
]
