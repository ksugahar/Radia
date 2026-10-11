"""Exercise a real owned headless Cubit session through the stdio MCP protocol."""
import argparse, asyncio, json, os, sys
from pathlib import Path
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

async def main():
    root = Path(__file__).resolve().parents[2]
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env['PYTHONPATH'] = str(root/'packages/cubit-mesh-export/src')
    if args.plugin_dir:
        env['CUBIT_PLUGIN_DIR'] = str(args.plugin_dir.resolve())
    env['CUBIT_MCP_TEMP'] = str(work)
    params = StdioServerParameters(command=sys.executable,args=['-m','cubit_mesh_export.mcp.server'],env=env)
    results = {}
    async with stdio_client(params) as (read,write):
        async with ClientSession(read,write) as session:
            await session.initialize()
            catalog = await session.list_tools()
            results['tool_count'] = len(catalog.tools)
            async def call(name,args):
                result = await session.call_tool(name,args)
                if result.isError:
                    raise RuntimeError(f'{name} failed: {result}')
                payload = json.loads(next(c.text for c in result.content if c.type=='text'))
                results[name] = payload
                return payload
            try:
                sample = root/'packages/cubit-mesh-export/src/cubit_mesh_export/cubit_gui/solver_ready_sample.jou'
                mesh = work/'mcp.vol'
                payload = await call('cubit_exec',{'commands':[f'play "{sample.as_posix()}"',f'export netgen "{mesh.as_posix()}" order 2 overwrite'],'timeout_s':120})
                if not payload['all_ok'] or payload['gui_started'] is not False:
                    raise RuntimeError(f'Headless export failed: {payload}')
                payload = await call('cubit_check_vol',{'vol_path':str(mesh),'strict_labels':True,'quality':True,'threshold_pct':0.05,'report_json':str(work/'check-vol.json')})
                if payload.get('passed') is not True:
                    raise RuntimeError(f'Mesh acceptance failed: {payload}')
                status = await call('cubit_session_status',{})
                if status['gui_started'] is not False or status['owned'] is not True:
                    raise RuntimeError(f'Unexpected session ownership or mode: {status}')
            finally:
                shutdown = await call('cubit_session_shutdown',{})
                if shutdown.get('stopped') != 'owned-child':
                    raise RuntimeError(f'Owned session shutdown failed: {shutdown}')
    results['passed'] = True
    (work/'result.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    print(json.dumps({'passed':True,'tool_count':results['tool_count'],'result':str(work/'result.json')}))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work', type=Path, required=True, help='Owned output directory outside the checkout')
    parser.add_argument('--plugin-dir', type=Path, help='Directory containing the candidate plugin DLL/CCM')
    args = parser.parse_args()
    asyncio.run(main())

