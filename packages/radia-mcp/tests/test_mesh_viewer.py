"""Guards: literal .vol paths, checked Netgen runtime, asynchronous GUI launch.
Assumes: mocked native APIs here; real current-user GUI evidence is separate.
"""
import asyncio
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace as NS
import pytest
import radia_mcp
from radia_mcp.mesh_viewer import launcher, mcp as adapter


@pytest.fixture(autouse=True)
def netgen_metadata_pin(monkeypatch):
    monkeypatch.setattr(launcher.importlib.metadata, 'requires',
                        lambda name: ['netgen-mesher==6.2.2608'])


@pytest.mark.parametrize("extension",[".msh",".py",".geo"])
def test_other_extensions_are_rejected_before_spawn(tmp_path,monkeypatch,extension):
    path=tmp_path/("mesh"+extension);path.write_text("data",encoding="utf-8")
    monkeypatch.setattr(adapter.subprocess,"Popen",lambda *a,**k:pytest.fail("unexpected process"))
    with pytest.raises(ValueError,match=".vol"):
        asyncio.run(adapter.open_mesh_viewer(str(path)))


def test_missing_file_is_rejected(tmp_path):
    with pytest.raises(FileNotFoundError):launcher.validate_mesh_path(str(tmp_path/"missing.vol"))


def test_wrong_netgen_version_is_not_a_fallback(tmp_path,monkeypatch):
    monkeypatch.setattr(launcher.importlib.metadata,"version",lambda _:"6.2.2607")
    with pytest.raises(RuntimeError,match="6.2.2608 required"):
        launcher.runtime_identity(tmp_path/"mesh.vol")


def test_identity_binds_the_actual_native_module(tmp_path,monkeypatch):
    module=tmp_path/"__init__.py";module.write_bytes(b"module")
    native=tmp_path/"libngpy.pyd";native.write_bytes(b"native payload")
    fake=NS(__version__="6.2.2608",__file__=str(module),libngpy=NS(__file__=str(native)))
    monkeypatch.setitem(sys.modules,"netgen",fake)
    monkeypatch.setattr(launcher.importlib.metadata,"version",lambda _:"6.2.2608")
    identity=launcher.runtime_identity(tmp_path/"mesh.vol")
    assert identity["runtime"]["native_module"]["sha256"]==hashlib.sha256(native.read_bytes()).hexdigest()
    fake.__version__="6.2.2607"
    with pytest.raises(RuntimeError,match="module/distribution"):
        launcher.runtime_identity(tmp_path/"mesh.vol")


@pytest.mark.parametrize("close_style",["window","native-menu","child-only"])
def test_filename_is_one_literal_tcl_argument_and_window_closes_once(tmp_path,monkeypatch,close_style):
    path=tmp_path/"日本語 空白 {} [] $ &;.vol";path.write_text("mesh3d",encoding="utf-8")
    calls=[];protocol={};destroyed=[];bindings={}
    def destroy():
        destroyed.append(True)
        if "<Destroy>" in bindings:bindings["<Destroy>"](NS(widget=win))
    def mainloop():
        if close_style=="window":protocol["WM_DELETE_WINDOW"]()
        elif close_style=="native-menu":destroy()
        elif "<Destroy>" in bindings:bindings["<Destroy>"](NS(widget=NS()))
    win=NS(tk=NS(call=lambda *args:calls.append(args)),
           bind=lambda name,fn,add=None:bindings.update({name:fn}),
           protocol=lambda name,fn:protocol.update({name:fn}),title=lambda title:None,
           destroy=destroy,mainloop=mainloop)
    monkeypatch.setitem(sys.modules,"gmsh",None)  # .vol must not import the optional Gmsh runtime.
    monkeypatch.setitem(sys.modules,"netgen",NS(__path__=[]))
    monkeypatch.setitem(sys.modules,"netgen.gui",NS(win=win))
    launcher.open_mesh(path)
    assert calls[0]==("Ng_LoadMesh",path.as_posix())
    assert len(destroyed)==1
    assert path.read_text(encoding="utf-8")=="mesh3d"


def test_check_mode_does_not_enter_gui(tmp_path,monkeypatch,capsys):
    path=tmp_path/"mesh.vol";path.write_text("mesh3d",encoding="utf-8")
    monkeypatch.setattr(launcher,"require_interactive_desktop",lambda:{"session_id":1})
    monkeypatch.setattr(launcher,"runtime_identity",lambda p:{"backend":"netgen"})
    monkeypatch.setattr(launcher,"open_mesh",lambda p:pytest.fail("unexpected GUI"))
    assert launcher.main([str(path),"--check"])==0
    assert json.loads(capsys.readouterr().out)=={"backend":"netgen","session_id":1}


def test_non_windows_session_fails_before_gui(monkeypatch):
    monkeypatch.setattr(launcher.sys,"platform","linux")
    with pytest.raises(RuntimeError,match="interactive Windows"):
        launcher.require_interactive_desktop()


@pytest.mark.parametrize("returncode",[0,1])
def test_async_preflight_and_detached_launch(tmp_path,monkeypatch,returncode):
    path=tmp_path/"日本語 空白 [$];.vol";path.write_text("mesh3d",encoding="utf-8")
    calls=[];heartbeat=[]
    async def communicate():
        await asyncio.sleep(.01)
        return json.dumps({"backend":"netgen"}).encode(),b"bad runtime"
    async def probe(*args,**kwargs):
        calls.append(("probe",args,kwargs));return NS(returncode=returncode,communicate=communicate)
    def spawn(*args,**kwargs):calls.append(("GUI",args,kwargs));return NS(pid=123)
    for name,value in [("CREATE_NO_WINDOW",0x08000000),("DETACHED_PROCESS",8),("CREATE_NEW_PROCESS_GROUP",0x200)]:
        monkeypatch.setattr(adapter.subprocess,name,value,raising=False)
    monkeypatch.setattr(adapter,"require_interactive_desktop",lambda:{"session_id":1})
    monkeypatch.setattr(adapter.asyncio,"create_subprocess_exec",probe)
    monkeypatch.setattr(adapter.subprocess,"Popen",spawn)
    async def scenario():
        task=asyncio.create_task(adapter.open_mesh_viewer(str(path)))
        await asyncio.sleep(.001);heartbeat.append(True)
        return await task
    if returncode:
        with pytest.raises(RuntimeError,match="bad runtime"):asyncio.run(scenario())
        assert [row[0] for row in calls]==["probe"]
    else:
        result=asyncio.run(scenario())
        assert result["pid"]==123 and result["status"]=="launched"
        assert result["argv"][-1]==str(path.resolve())
        assert calls[1][2]["creationflags"]==0x208 and "shell" not in calls[1][2]
        assert adapter.SOURCE_ROOT == Path(radia_mcp.__file__).resolve().parent.parent
        assert calls[1][2]["env"]["PYTHONPATH"] == str(adapter.SOURCE_ROOT)
    assert heartbeat==[True]


def test_registered_mcp_adapter_is_async():
    from mcp.server.fastmcp import FastMCP
    server=FastMCP("netgen-viewer-contract");adapter.register_mesh_viewer_tool(server)
    tools=asyncio.run(server.list_tools())
    assert [tool.name for tool in tools]==["open_mesh_viewer"]
    assert tools[0].annotations.readOnlyHint is False


def test_netgen_pin_follows_installed_radia_metadata(tmp_path, monkeypatch):
    monkeypatch.setattr(launcher.importlib.metadata, 'requires',
                        lambda name: ['netgen-mesher==6.2.2700'])
    monkeypatch.setattr(launcher.importlib.metadata, 'version', lambda name: '6.2.2608')
    with pytest.raises(RuntimeError, match='6.2.2700 required by radia metadata'):
        launcher.runtime_identity(tmp_path / 'mesh.vol')


@pytest.mark.parametrize('requirements', [[], ['netgen-mesher>=6.2'],
                                        ['netgen-mesher==6.2.1', 'netgen-mesher==6.2.2']])
def test_ambiguous_or_unpinned_netgen_metadata_fails(tmp_path, monkeypatch, requirements):
    monkeypatch.setattr(launcher.importlib.metadata, 'requires', lambda name: requirements)
    with pytest.raises(RuntimeError, match='one exact netgen-mesher pin'):
        launcher.runtime_identity(tmp_path / 'mesh.vol')
