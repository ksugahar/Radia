import asyncio
import hashlib
import json
import os
from pathlib import Path
import sys
import time
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

async def main():
    source = Path(__file__).resolve().parents[2] / 'packages/radia-mcp/src'
    env = os.environ.copy()
    env['PYTHONPATH'] = str(source)
    params = StdioServerParameters(command=sys.executable, args=['-m', 'radia_mcp.mathematica.server'], env=env)
    out = {'interpreter': sys.executable, 'source': str(source), 'tools_sha256': hashlib.sha256((source/'radia_mcp/mathematica/tools.py').read_bytes()).hexdigest(), 'checks': []}
    async with stdio_client(params) as (read, write):
        async with ClientSession(read, write) as session:
            await session.initialize()
            catalog = await session.list_tools()
            out['tools'] = [t.name for t in catalog.tools]
            for name, args, expected in [
                ('mathematica_evaluate', {'code': '{$Version, $ProcessorCount, 1+1}', 'timeout': 30}, '2'),
                ('mathematica_simplify', {'expression': 'Sin[x]^2+Cos[x]^2', 'timeout': 30}, '1'),
                ('mathematica_check_identity', {'lhs': 'D[Sin[x],x]', 'rhs': 'Cos[x]', 'timeout': 30}, 'True'),
            ]:
                start = time.perf_counter()
                result = await session.call_tool(name, args)
                content = [c.text for c in result.content if c.type == 'text']
                assert not result.isError, (name, content)
                payload = json.loads(content[0])
                if name == 'mathematica_evaluate':
                    assert payload['exit_code'] == 0 and not payload['timed_out'] and '14.2' in payload['result'], payload
                elif name == 'mathematica_check_identity':
                    assert payload['ok'] and payload['identical'] is True, payload
                else:
                    assert payload['ok'] and payload['value'].strip() == expected, payload
                out['checks'].append({'tool': name, 'seconds': time.perf_counter()-start, 'result': payload})
    out['passed'] = True
    target = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).with_name('result.json')
    target.write_text(json.dumps(out, indent=2), encoding='utf-8')
    print(json.dumps({'passed': True, 'tool_count': len(out['tools']), 'checks': len(out['checks']), 'result': str(target)}))

asyncio.run(main())

