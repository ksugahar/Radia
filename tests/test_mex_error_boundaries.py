from pathlib import Path
import importlib.util
import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('mex_boundary', ROOT/'packages/radia-simulink/tools/check_mex_error_boundaries.py')
CHECK = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHECK)

GOOD = '''void mexFunction(int nl,mxArray** out,int nr,const mxArray** in){
char errorMessage[2048]={0};
if(!radia_mex::checked(dispatch,nl,out,nr,in,errorMessage,sizeof(errorMessage))){
mexErrMsgIdAndTxt("test:failure","%s",errorMessage);}}'''

def test_all_owned_mex_error_exits_have_pod_gateway():
    assert CHECK.audit(ROOT/'packages/radia-simulink/native') == []

@pytest.mark.parametrize('local', ['std::string value;', 'std::unique_ptr<int> value;', 'Owner<int> value;', 'auto callback=[]{};'])
def test_gateway_rejects_nontrivial_or_unrecognized_locals(local):
    assert CHECK.inspect_text(GOOD.replace('char errorMessage', local+'char errorMessage'))

def test_catch_scope_and_helper_error_exits_are_rejected():
    assert CHECK.inspect_text('void mexFunction(){try{run();}catch(const std::exception& e){mexErrMsgIdAndTxt("test:failure","%s",e.what());}}')
    assert CHECK.inspect_text('void helper(){mexErrMsgTxt("failure");}')

def test_literals_and_comments_do_not_invent_error_exits():
    assert CHECK.inspect_text('// mexErrMsgTxt("comment")\nconst char* x="mexErrMsgTxt({)";') == []
    assert CHECK.inspect_text(GOOD) == []

def test_new_source_file_cannot_escape_scan(tmp_path):
    (tmp_path/'new.cu').write_text('void helper(){Owner x;mexErrMsgTxt("failure");}')
    assert CHECK.audit(tmp_path)
