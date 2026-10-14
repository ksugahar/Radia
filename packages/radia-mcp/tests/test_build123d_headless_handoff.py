"""The Radia-owned optional Cubit consumer must remain headless."""

import inspect

from radia_mcp.build123d import server


def test_build123d_cubit_handoff_names_the_batch_session():
    source = inspect.getsource(server)
    assert 'CubitSession.get(mode="batch")' in source
    assert 'CubitSession.get(mode="gui")' not in source
