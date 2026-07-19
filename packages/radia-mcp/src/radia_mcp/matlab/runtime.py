import hashlib,json,os,re,shutil
from pathlib import Path
HERE=Path(__file__).parent; EXT=HERE/"extensions/radia-matlab-tools.json"; MATLAB=HERE/"matlab"
PROFILES={"existing":("existing","desktop"),"auto_nodesktop":("auto","nodesktop"),"new_nodesktop":("new","nodesktop")}
def matlab_extension_path():
    if not EXT.is_file(): raise FileNotFoundError(EXT)
    return EXT
def matlab_extension_contract():
    raw=matlab_extension_path().read_bytes(); data=json.loads(raw); tools=data.get("tools",[]); sig=data.get("signatures",{}); files=list((MATLAB/"+radia_mcp_matlab").glob("*.m")); names={p.stem for p in files}; errors=[]
    if len(tools)!=43 or len(sig)!=43: errors.append("expected 43 tools/signatures")
    if len(files)!=86: errors.append("expected 86 MATLAB functions")
    if any(not t["name"].startswith("matlab_") for t in tools): errors.append("invalid tool prefix")
    for p in files:
        text=p.read_text(encoding="utf-8")
        if "acoustic" in text.lower() or "fembem" in text.lower(): errors.append(f"{p.name}: namespace leak")
        errors += [f"{p.name}: missing {d}" for d in re.findall(r"radia_mcp_matlab\.([A-Za-z]\w*)",text) if d not in names]
    return {"schema":"radia-mcp.matlab-extension/v1","ok":not errors,"status":"ok" if not errors else "error","runtime_owner":"MathWorks MATLAB MCP Core Server","extension_file":str(EXT),"matlab_root":str(MATLAB),"sha256":hashlib.sha256(raw).hexdigest(),"tool_count":len(tools),"tool_names":[t["name"] for t in tools],"signature_count":len(sig),"matlab_function_count":len(files),"errors":errors}
def matlab_official_server_config(profile="existing",*,include_generic_extension=False):
    if profile not in PROFILES: raise ValueError(profile)
    session,display=PROFILES[profile]; args=[f"--matlab-session-mode={session}"]; files=[]; setup=""
    if include_generic_extension:
        c=matlab_extension_contract()
        if not c["ok"]: raise RuntimeError(c["errors"])
        files=[c["extension_file"]]; args += [f"--extension-file={c['extension_file']}"]; setup=f"addpath('{c['matlab_root']}');"
    if display=="nodesktop": args.append("--matlab-display-mode=nodesktop")
    cmd=os.getenv("RADIA_MATLAB_MCP_SERVER") or next((shutil.which(x) for x in ["matlab-mcp-core-server","matlab-mcp-core-server-win64.exe"] if shutil.which(x)),"matlab-mcp-core-server")
    return {"schema":"radia-mcp.matlab-server-config/v1","status":"ok","runtime_owner":"MathWorks MATLAB MCP Core Server","integration_owner":"radia-mcp.matlab","command_id":"matlab-mcp-core-server","command":cmd,"profile":profile,"args":args,"extension_files":files,"matlab_setup_code":setup}
def matlab_radia_acoustic_interface_contract(): return {"runtime_owner":"MathWorks MATLAB MCP Core Server","generic_operations_owner":"radia_mcp.matlab","generic_matlab_package":"radia_mcp_matlab","education_solver_owner":"radia_mcp.acoustic_fembem"}
