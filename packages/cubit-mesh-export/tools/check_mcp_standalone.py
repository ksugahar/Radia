"""Installed-wheel stdio acceptance in a venv containing no Radia products."""
import asyncio
import importlib.util
import json
from pathlib import Path
import sys
import tempfile

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client


async def check():
    assert importlib.util.find_spec('radia') is None, 'Radia must be absent'
    assert importlib.util.find_spec('radia_mcp') is None, 'radia-mcp must be absent'
    assert importlib.util.find_spec('cae_mcp_core') is None, 'Retired foundation must be absent'
    import cubit_mesh_export
    assert Path(cubit_mesh_export.__file__).resolve().is_relative_to(Path(sys.prefix).resolve())
    with tempfile.TemporaryDirectory(prefix='cubit-mcp-wheel-', dir='C:/temp') as work:
        params = StdioServerParameters(command=sys.executable,
            args=['-I', '-m', 'cubit_mesh_export.mcp.server'], cwd=work)
        async with stdio_client(params) as (read, write):
            async with ClientSession(read, write) as client:
                initialized = await client.initialize()
                assert 'Never launch or attach to the Cubit GUI' in initialized.instructions
                names = {t.name for t in (await client.list_tools()).tools}
                assert {'cubit_status', 'cubit_exec', 'cubit_import_journal',
                        'cubit_validation_catalog'}.issubset(names)
                assert not {'cubit_show', 'open_in_cubit'}.intersection(names)
                result = await client.call_tool('cubit_status', {})
                assert not result.isError
                status = json.loads(result.content[0].text)
                provenance = status['runtime_provenance']
                assert provenance['distribution']['name'] == 'cubit-mesh-export'
                assert provenance['distribution']['version'] == cubit_mesh_export.__version__
                assert provenance['module'] == 'cubit_mesh_export.mcp.server'
                assert status['runtime_contract']['complete']
                print(json.dumps({'passed': True, 'tools': len(names),
                                  'provenance': provenance}, ensure_ascii=False))


if __name__ == '__main__':
    asyncio.run(check())
