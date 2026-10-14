from __future__ import annotations

import sys
import types

from radia_mcp.radia_ngsolve.solve import laminated_mu_eff


def test_laminated_mu_eff_delegates_to_radia_domain_api(monkeypatch):
    calls = []

    def domain_api(*args):
        calls.append(args)
        return 3.0 - 4.0j

    radia_package = types.ModuleType("radia")
    radia_package.__path__ = []
    lamination_module = types.ModuleType("radia.lamination")
    lamination_module.laminated_mu_eff = domain_api
    monkeypatch.setitem(sys.modules, "radia", radia_package)
    monkeypatch.setitem(sys.modules, "radia.lamination", lamination_module)

    result = laminated_mu_eff(1200.0, 2.1e6, 314.0, 0.35e-3, 0.96)

    assert result == 3.0 - 4.0j
    assert calls == [(1200.0, 2.1e6, 314.0, 0.35e-3, 0.96)]
