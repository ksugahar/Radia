"""
Lint rules for Radia Python scripts.

Each rule function receives (filepath, lines) and returns a list of findings.
A finding is a dict with keys: line, severity, rule, message.
"""

import re
import os
from typing import List, Dict


def check_objbckg_no_lambda(filepath: str, lines: List[str]) -> List[Dict]:
    """CRITICAL: rad.ObjBckg() must receive a callable, not a list."""
    findings = []
    pattern = re.compile(r'rad\.ObjBckg\s*\(\s*\[')
    for i, line in enumerate(lines, 1):
        if pattern.search(line):
            findings.append({
                'line': i,
                'severity': 'CRITICAL',
                'rule': 'objbckg-needs-callable',
                'message': (
                    'rad.ObjBckg() requires a callable (lambda or function), '
                    'not a list. Use: rad.ObjBckg(lambda p: [Bx, By, Bz])'
                ),
            })
    return findings


def check_missing_utidelall(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: Scripts should call rad.UtiDelAll() for cleanup."""
    findings = []
    # Only check files that import radia
    has_radia = any('import radia' in line or 'import rad' in line for line in lines)
    if not has_radia:
        return findings

    # Skip __init__.py and module files (they shouldn't call UtiDelAll)
    basename = os.path.basename(filepath)
    if basename.startswith('__') or basename in ('radia.py',):
        return findings

    # Only check example scripts (files with if __name__ == '__main__')
    has_main = any("__name__" in line and "__main__" in line for line in lines)
    if not has_main:
        return findings

    has_utidelall = any('UtiDelAll' in line for line in lines)
    if not has_utidelall:
        findings.append({
            'line': len(lines),
            'severity': 'HIGH',
            'rule': 'missing-utidelall',
            'message': (
                'Script imports radia but does not call rad.UtiDelAll(). '
                'Add cleanup before script exit.'
            ),
        })
    return findings


def check_hardcoded_absolute_paths(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: No hardcoded absolute paths in sys.path or path assignments."""
    findings = []
    # Patterns for hardcoded Windows absolute paths used as imports/paths
    patterns = [
        re.compile(r'sys\.path\.insert\s*\(\s*\d+\s*,\s*r?["\'][A-Za-z]:\\'),
        re.compile(r'(?:repo_root|work_dir|base_dir)\s*=\s*r?["\'][A-Za-z]:\\'),
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
                    'rule': 'hardcoded-absolute-path',
                    'message': (
                        'Hardcoded absolute path detected. '
                        'Use os.path.dirname(__file__) for relative paths.'
                    ),
                })
                break
    return findings


def check_removed_fldunits(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: rad.FldUnits() has been removed. Radia always uses meters."""
    findings = []
    pattern = re.compile(r'(?:rad\.)?FldUnits\s*\(')
    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith('#'):
            continue
        if pattern.search(stripped):
            findings.append({
                'line': i,
                'severity': 'HIGH',
                'rule': 'removed-fldunits',
                'message': (
                    'rad.FldUnits() has been removed. Radia always uses meters. '
                    'Delete this call.'
                ),
            })
    return findings


def check_removed_fldbatch(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: FldBatch/FldA/FldPhi removed. Use unified Fld() with batch points."""
    findings = []
    pattern = re.compile(r'(?:rad\.)?(?:FldBatch|FldA|FldPhi)\s*\(')
    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith('#'):
            continue
        if pattern.search(stripped):
            findings.append({
                'line': i,
                'severity': 'HIGH',
                'rule': 'removed-fldbatch',
                'message': (
                    'FldBatch/FldA/FldPhi are removed. '
                    'Use rad.Fld(obj, field_type, points) with shape (N,3) for batch.'
                ),
            })
    return findings


def check_removed_solver_apis(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: Old solver parameter APIs removed. Use SolverConfig()."""
    findings = []
    old_apis = [
        'SetHACApKParams', 'SetHMatrixEpsilon',
        'SetBiCGSTABTol', 'GetBiCGSTABTol',
        'SetRelaxParam', 'GetRelaxParam',
        'SetNewtonMethod', 'GetNewtonMethod',
        'SetNewtonDamping', 'GetNewtonDampingStats',
        'SetHMatrixFieldEval',
    ]
    pattern = re.compile(
        r'(?:rad\.)?(?:' + '|'.join(old_apis) + r')\s*\('
    )
    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith('#'):
            continue
        if pattern.search(stripped):
            findings.append({
                'line': i,
                'severity': 'HIGH',
                'rule': 'removed-solver-api',
                'message': (
                    'Old solver parameter API removed. '
                    'Use rad.SolverConfig(**kwargs) and rad.GetSolverConfig() instead.'
                ),
            })
    return findings


def check_docstring_hardcoded_mm(filepath: str, lines: List[str]) -> List[Dict]:
    """MODERATE: Public API docstrings should not hardcode 'in mm'."""
    findings = []
    # Only check src/radia modules
    if 'src/radia' not in filepath.replace('\\', '/') and 'src\\radia' not in filepath:
        return findings

    # Look for "in mm" in docstrings (inside triple-quoted strings)
    in_docstring = False
    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        if '"""' in stripped:
            count = stripped.count('"""')
            if count == 1:
                in_docstring = not in_docstring
            # If count == 2, it's an inline docstring - stays outside
        if in_docstring and re.search(r'\bpoint.*\bin mm\b', stripped):
            findings.append({
                'line': i,
                'severity': 'MODERATE',
                'rule': 'docstring-hardcoded-mm',
                'message': (
                    'Docstring hardcodes "in mm". '
                    'Use "in constructor length units" for unit-agnostic API.'
                ),
            })
    return findings


def check_build_release_path(filepath: str, lines: List[str]) -> List[Dict]:
    """LOW: build/Release path import may be outdated."""
    findings = []
    pattern = re.compile(r"sys\.path\.insert.*['\"].*build/Release['\"]")
    for i, line in enumerate(lines, 1):
        if pattern.search(line):
            findings.append({
                'line': i,
                'severity': 'LOW',
                'rule': 'build-release-path',
                'message': (
                    'Import from build/Release may be outdated. '
                    'Prefer src/radia which has the latest binaries.'
                ),
            })
    return findings


# SHARED: also in ngsolve/rules.py - keep synchronized
def check_bessel_jv_for_sibc(filepath: str, lines: List[str]) -> List[Dict]:
    """CRITICAL: Circular wire SIBC must use iv (modified Bessel), not jv."""
    findings = []
    # Check if file uses SIBC/impedance AND imports jv (not iv)
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
        # Check for import of jv from scipy.special in SIBC context
        if re.search(r'from\s+scipy\.special\s+import\s+.*\bjv\b', stripped):
            # Check it's not also importing iv
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


# SHARED: also in ngsolve/rules.py - keep synchronized
def check_peec_low_nseg(filepath: str, lines: List[str]) -> List[Dict]:
    """MODERATE: PEEC circular coil should use n_seg >= 32 for coupling accuracy."""
    findings = []
    # Only check files that look like PEEC scripts
    has_peec_context = any(
        kw in line
        for line in lines
        for kw in ['FastHenryParser', 'PEECBuilder', 'peec', 'n_seg']
    )
    if not has_peec_context:
        return findings

    # Look for n_seg=<small number> patterns
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


# SHARED: also in ngsolve/rules.py - keep synchronized
def check_efie_v_minus_sign(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: EFIE V term must have positive sign (Lenz's law)."""
    findings = []
    # Only check files with EFIE/BEM/shield context
    has_efie_context = any(
        kw in line.lower()
        for line in lines
        for kw in ['shieldbemsibc', 'efie', 'v_ll', 'maxwell_slp',
                    'maxwellsinglelayer']
    )
    if not has_efie_context:
        return findings

    # Look for minus sign before V_LL or mu_0 * V in EFIE system assembly
    # Patterns that indicate wrong sign:
    #   - Zs * M_LL - jw * mu_0 * V_LL  (minus before jw*mu_0*V)
    #   - M_LL - ...V_LL
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


# SHARED: also in ngsolve/rules.py - keep synchronized
def check_classical_efie_breakdown(filepath: str, lines: List[str]) -> List[Dict]:
    """MODERATE: Classical EFIE 1/kappa^2 causes low-frequency breakdown."""
    findings = []
    # Only check files using ngsolve.bem/HelmholtzSL
    has_helmholtz = any(
        kw in line
        for line in lines
        for kw in ['HelmholtzSL', 'HelmholtzSingleLayer']
    )
    if not has_helmholtz:
        return findings

    # Look for 1/kappa**2 or (1/kappa^2) patterns
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


# SHARED: also in ngsolve/rules.py - keep synchronized
def check_peec_p_over_jw(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: PEEC Loop-Star P/(jw) causes low-frequency breakdown."""
    findings = []
    # Only check files with PEEC/Loop-Star context
    has_peec_ls = any(
        kw in line
        for line in lines
        for kw in ['loop_star', 'Loop-Star', 'LoopStar', 'Z_SS', 'P_matrix',
                    'NGBEMPEECSolver', 'ngsbem_peec']
    )
    if not has_peec_ls:
        return findings

    # Detect P/(jw) or P/(1j*omega) patterns in actual code (not comments/docstrings)
    patterns = [
        re.compile(r'(?:self\.)?P\s*/\s*\(?\s*(?:1j|1\.?j)\s*\*\s*(?:omega|w)'),
        re.compile(r'(?:self\.)?P\s*/\s*\(?\s*(?:omega|w)\s*\*\s*(?:1j|1\.?j)'),
        re.compile(r'Z_SS\s*=\s*.*(?:self\.)?P\s*/'),
    ]
    in_docstring = False
    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        # Track docstring state
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




def check_msc_wrong_sign(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: MSC system matrix must be (1/chi + N), not (-1/chi - N)."""
    findings = []
    # Only check files with MSC/MMMBuilder context
    has_msc = any(
        kw in line
        for line in lines
        for kw in ['MMMSolver', 'MMMBuilder', 'inv_chi', 'K_msc', 'msc_matvec']
    )
    if not has_msc:
        return findings

    patterns = [
        re.compile(r'-\s*np\.diag\s*\(\s*inv_chi\s*\)\s*-\s*N'),
        re.compile(r'-\s*self\._inv_chi\s*\*\s*x\s*-\s*self\._N'),
        re.compile(r'A\s*=\s*-N'),
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
                    'rule': 'msc-wrong-sign',
                    'message': (
                        'MSC system matrix has wrong sign. '
                        'Correct: K = diag(1/chi) + N (positive). '
                        'The old convention (-1/chi - N) produces inverted magnetization.'
                    ),
                })
                break
    return findings


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

    # Detect OCCGeometry(Box(...)) without Glue
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
                'mesh includes interior edges → singular SL matrix (cond ~1e17). '
                'Use: OCCGeometry(Glue(box.faces)) for surface-only mesh.'
            ),
        })
    return findings


def check_msc_eval_face_center(filepath: str, lines: List[str]) -> List[Dict]:
    """MODERATE: MSC eval point should be midpoint, not face center."""
    findings = []
    has_msc = any(
        kw in line
        for line in lines
        for kw in ['eval_point', 'eval_pts', 'MMMBuilder', 'Yano']
    )
    if not has_msc:
        return findings

    # Detect face_center without midpoint averaging with element center
    pattern = re.compile(r'eval_point.*=.*face_center\b(?!.*elem)')
    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith('#'):
            continue
        if pattern.search(stripped):
            findings.append({
                'line': i,
                'severity': 'MODERATE',
                'rule': 'msc-eval-face-center',
                'message': (
                    'MSC evaluation point should be Yano-Sugahara midpoint: '
                    '(face_center + element_center) / 2, not face_center alone. '
                    'Using face_center makes the self-interaction term singular.'
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


# All rules in execution order
ALL_RULES = [
    # Radia rules
    check_objbckg_no_lambda,
    check_missing_utidelall,
    check_hardcoded_absolute_paths,
    check_removed_fldunits,
    check_removed_fldbatch,
    check_removed_solver_apis,
    check_docstring_hardcoded_mm,
    check_build_release_path,
    check_bessel_jv_for_sibc,
    check_peec_low_nseg,
    check_efie_v_minus_sign,
    check_classical_efie_breakdown,
    check_peec_p_over_jw,
    # MSC / multi-level rules
    check_msc_wrong_sign,
    check_ngsbem_volume_mesh,
    check_msc_eval_face_center,
    # FEM-SIBC rules
    check_scattered_field_sibc_missing_nxH,
]
