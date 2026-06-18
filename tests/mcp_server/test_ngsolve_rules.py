"""Tests for NGSolve lint rules.

Each test provides intentionally bad code and verifies the rule fires,
then provides correct code and verifies no findings.
"""

import pytest
from radia_mcp.radia_ngsolve.rules import (
    ALL_RULES,
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
    check_bessel_jv_for_sibc,
    check_peec_low_nseg,
    check_efie_v_minus_sign,
    check_classical_efie_breakdown,
    check_peec_p_over_jw,
)


def _run(rule_func, code: str, filepath: str = "test_script.py"):
    lines = code.splitlines()
    return rule_func(filepath, lines)


# ============================================================
# NGSolve FEM rules
# ============================================================

class TestHCurlMissingNograds:
    def test_detects_hcurl_without_nograds(self):
        code = (
            '# magnetostatic solver\n'
            'fes = HCurl(mesh, order=2)\n'
        )
        findings = _run(check_hcurl_missing_nograds, code)
        assert len(findings) == 1

    def test_allows_nograds_true(self):
        code = (
            '# magnetostatic solver\n'
            'fes = HCurl(mesh, order=2, nograds=True)\n'
        )
        assert _run(check_hcurl_missing_nograds, code) == []

    def test_allows_complex(self):
        code = (
            '# magnetostatic solver\n'
            'fes = HCurl(mesh, order=2, complex=True)\n'
        )
        assert _run(check_hcurl_missing_nograds, code) == []


class TestNgsolvePrecondAfterAssemble:
    def test_detects_late_precond(self):
        code = (
            'a.Assemble()\n'
            'pre = Preconditioner(a, "bddc")\n'
        )
        findings = _run(check_ngsolve_precond_after_assemble, code)
        assert len(findings) == 1

    def test_allows_precond_before_assemble(self):
        code = (
            'pre = Preconditioner(a, "bddc")\n'
            'a.Assemble()\n'
        )
        assert _run(check_ngsolve_precond_after_assemble, code) == []


class TestNgsolveMissingTraceBem:
    def test_detects_missing_trace(self):
        code = (
            'fes = HDivSurface(mesh, order=1)\n'
            'u, v = fes.TnT()\n'
            'a = BilinearForm(LaplaceSL(3) * u * v * ds)\n'
        )
        findings = _run(check_ngsolve_missing_trace_bem, code)
        assert len(findings) == 1
        assert findings[0]['severity'] == 'CRITICAL'

    def test_allows_trace(self):
        code = (
            'fes = HDivSurface(mesh, order=1)\n'
            'u, v = fes.TnT()\n'
            'a = BilinearForm(LaplaceSL(3) * u.Trace() * v.Trace() * ds)\n'
        )
        assert _run(check_ngsolve_missing_trace_bem, code) == []


class TestNgsolveOverwriteXyz:
    def test_detects_for_x_in(self):
        code = (
            'from ngsolve import *\n'
            'for x in range(10):\n'
            '    pass\n'
        )
        findings = _run(check_ngsolve_overwrite_xyz, code)
        assert len(findings) == 1
        assert 'overwrites' in findings[0]['message']

    def test_allows_for_xi_in(self):
        code = (
            'from ngsolve import *\n'
            'for xi in range(10):\n'
            '    pass\n'
        )
        assert _run(check_ngsolve_overwrite_xyz, code) == []

    def test_skips_non_ngsolve(self):
        code = 'for x in range(10):\n    pass\n'
        assert _run(check_ngsolve_overwrite_xyz, code) == []


class TestNgsolveVecAssign:
    def test_detects_bad_assign(self):
        code = (
            'from ngsolve import *\n'
            'gfu.vec = f.vec\n'
        )
        findings = _run(check_ngsolve_vec_assign, code)
        assert len(findings) == 1

    def test_allows_data_assign(self):
        code = (
            'from ngsolve import *\n'
            'gfu.vec.data = f.vec\n'
        )
        assert _run(check_ngsolve_vec_assign, code) == []


class TestNgsolveDim2Occ:
    def test_detects_missing_dim2(self):
        code = (
            'shape = Rectangle(1, 1).Face()\n'
            'geo = OCCGeometry(shape)\n'
        )
        findings = _run(check_ngsolve_dim2_occ, code)
        assert len(findings) == 1

    def test_allows_dim2(self):
        code = (
            'shape = Rectangle(1, 1).Face()\n'
            'geo = OCCGeometry(shape, dim=2)\n'
        )
        assert _run(check_ngsolve_dim2_occ, code) == []


class TestNgsolveCgOnSaddlePoint:
    def test_detects_cg_on_saddle(self):
        code = (
            'fes = fesA * fesOmega\n'
            'gfu = solvers.CG(a.mat, rhs.vec)\n'
        )
        findings = _run(check_ngsolve_cg_on_saddle_point, code)
        assert len(findings) == 1

    def test_allows_gmres(self):
        code = (
            'fes = fesA * fesOmega\n'
            'gfu = solvers.GMRes(a.mat, rhs.vec)\n'
        )
        assert _run(check_ngsolve_cg_on_saddle_point, code) == []


class TestNgsolveVectorH1ForEm:
    def test_detects_vectorh1_em(self):
        code = (
            '# magnetostatic simulation\n'
            'fes = VectorH1(mesh, order=2)\n'
        )
        findings = _run(check_ngsolve_vectorh1_for_em, code)
        assert len(findings) == 1

    def test_allows_non_em(self):
        code = (
            '# structural mechanics\n'
            'fes = VectorH1(mesh, order=2)\n'
        )
        assert _run(check_ngsolve_vectorh1_for_em, code) == []


class TestNgsolvePinvitNoProjection:
    def test_detects_no_projection(self):
        code = (
            'fes = HCurl(mesh, order=2)\n'
            'evals, evecs = solvers.PINVIT(a.mat, m.mat, pre, num=10)\n'
        )
        findings = _run(check_ngsolve_pinvit_no_projection, code)
        assert len(findings) == 1

    def test_allows_with_projection(self):
        code = (
            'fes = HCurl(mesh, order=2)\n'
            'gradmat = fes.CreateGradient()\n'
            'evals, evecs = solvers.PINVIT(a.mat, m.mat, pre, num=10)\n'
        )
        assert _run(check_ngsolve_pinvit_no_projection, code) == []


class TestEddyCurrentMissingComplex:
    def test_detects_missing_complex(self):
        code = (
            '# eddy current simulation\n'
            'fes = HCurl(mesh, order=2)\n'
        )
        findings = _run(check_eddy_current_missing_complex, code)
        assert len(findings) == 1
        assert findings[0]['severity'] == 'HIGH'

    def test_allows_complex_true(self):
        code = (
            '# eddy current simulation\n'
            'fes = HCurl(mesh, order=2, complex=True)\n'
        )
        assert _run(check_eddy_current_missing_complex, code) == []

    def test_skips_thermal(self):
        code = (
            '# eddy current simulation\n'
            'fes_T = H1(mesh, order=2)  # temperature\n'
        )
        assert _run(check_eddy_current_missing_complex, code) == []


class TestJouleHeatMissingConj:
    def test_detects_missing_conj(self):
        code = 'Q = sigma * InnerProduct(E, E)'
        findings = _run(check_joule_heat_missing_conj, code)
        assert len(findings) == 1

    def test_allows_conj(self):
        code = 'Q = 0.5 * sigma * InnerProduct(E, Conj(E)).real'
        assert _run(check_joule_heat_missing_conj, code) == []


class TestKelvinMissingBonusIntorder:
    def test_detects_missing_bonus(self):
        code = 'a += mu * curl(u) * curl(v) * dx("Kelvin")'
        findings = _run(check_ngsolve_kelvin_missing_bonus_intorder, code)
        assert len(findings) == 1

    def test_allows_bonus(self):
        code = 'a += mu * curl(u) * curl(v) * dx("Kelvin", bonus_intorder=4)'
        assert _run(check_ngsolve_kelvin_missing_bonus_intorder, code) == []


# ============================================================
# Shared PEEC/BEM rules (also tested in test_radia_rules.py)
# ============================================================

class TestBesselJvForSibc:
    def test_detects_jv_in_sibc_context(self):
        code = (
            '# SIBC calculation\n'
            'from scipy.special import jv\n'
            'Z = jv(0, ka) / jv(1, ka)\n'
        )
        findings = _run(check_bessel_jv_for_sibc, code)
        assert len(findings) == 1
        assert findings[0]['severity'] == 'CRITICAL'


class TestPeecLowNseg:
    def test_detects_low_nseg(self):
        code = (
            '# PEEC coil model\n'
            'n_seg = 16\n'
        )
        findings = _run(check_peec_low_nseg, code)
        assert len(findings) == 1


class TestEfieVMinusSign:
    def test_detects_minus_sign(self):
        code = (
            '# ShieldBEMSIBC solver\n'
            'Z_LL = Zs * M_LL - 1j * omega * mu_0 * V_LL\n'
        )
        findings = _run(check_efie_v_minus_sign, code)
        assert len(findings) == 1


class TestClassicalEfieBreakdown:
    def test_detects_one_over_kappa_squared(self):
        code = (
            'from ngbem import HelmholtzSL\n'
            'P = 1 / kappa ** 2 * V_SL\n'
        )
        findings = _run(check_classical_efie_breakdown, code)
        assert len(findings) == 1


class TestPeecPOverJw:
    def test_detects_p_over_jw(self):
        code = (
            '# Loop-Star PEEC\n'
            'Z_SS = self.P / (1j * omega)\n'
        )
        findings = _run(check_peec_p_over_jw, code)
        assert len(findings) == 1


# ============================================================
# Meta-tests
# ============================================================

class TestAllRulesList:
    def test_all_rules_count(self):
        # Tripwire: bump when adding/removing a rule in
        # radia_ngsolve.rules.ALL_RULES (NOT covered by ci_preflight's
        # pre-push gate -- only this build-test pytest runs it).
        assert len(ALL_RULES) == 36

    def test_all_rules_callable(self):
        for rule in ALL_RULES:
            assert callable(rule)

    def test_clean_code_no_findings(self):
        """Clean code should produce zero findings across all rules."""
        clean = (
            'import numpy as np\n'
            'x = np.linspace(0, 1, 100)\n'
        )
        lines = clean.splitlines()
        for rule in ALL_RULES:
            findings = rule("clean.py", lines)
            assert findings == [], f"{rule.__name__} fired on clean code"
