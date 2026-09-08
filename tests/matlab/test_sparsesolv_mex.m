function tests = test_sparsesolv_mex
tests = functiontests(localfunctions);
end

function setupOnce(t)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'matlab'));
radia.setup(Force=true);
folder = tempname('C:\temp');
mkdir(folder);
t.TestData.folder = folder;
helper = fullfile(root,'tests','matlab','sparsesolv_python_reference.py');
[status, output] = radia.internal.runPythonProcess('python "' + string(helper) + '" "' + string(folder) + '"');
assert(status == 0, 'Python oracle failed: %s', output);
t.TestData.ref = jsondecode(fileread(fullfile(folder,'reference.json')));
end

function teardownOnce(t)
if isfolder(t.TestData.folder)
    rmdir(t.TestData.folder,'s');
end
end

function testPythonParityAndLifetime(t)
r = t.TestData.ref;
metrics = struct([]);
for i = 1:numel(r.cases)
    c = r.cases(i);
    fprintf('Checking %s\n',c.name);
    mesh = radia.ngsolve.Mesh.create(r.mesh);
    space = radia.ngsolve.FESpace.create(mesh,"hcurl",1,NoGrads=true);
    form = radia.ngsolve.BilinearForm.create(space,"mass");
    aux = form.matrix();
    complexCase = startsWith(string(c.name),"complex");
    if complexCase
        cs = radia.ngsolve.FESpace.create(mesh,"hcurl",1,NoGrads=true,Complex=true);
        weight = radia.ngsolve.CoefficientFunction.constant(1+0.5i);
        cf = radia.ngsolve.BilinearForm.createFromCoefficient(cs,"mass",weight);
        a = cf.matrix();
    else
        a = aux;
    end
    if endsWith(string(c.name),"ams")
        pre = radia.sparsesolv.AMS(aux,space,Complex=complexCase);
    else
        pre = radia.sparsesolv.IC(a);
    end
    rhs = a.vector();
    verifyEqual(t,pre.IsComplex,complexCase);
    rhs.setValues(unpack(c.rhs,complexCase));
    applied = pre.matvec(rhs);
    verifyEqual(t,applied.values(),unpack(c.applied,complexCase),RelTol=1e-10,AbsTol=1e-10);
    inv = radia.sparsesolv.COCR(a,pre,Tolerance=1e-11);
    dense = sparse(a);
    delete(pre); delete(a); delete(aux); delete(form); delete(space); delete(mesh);
    if complexCase
        delete(cf); delete(cs); delete(weight);
    end
    solution = inv.matvec(rhs);
    verifyEqual(t,solution.values(),unpack(c.solution,complexCase),RelTol=1e-8,AbsTol=1e-8);
    verifyLessThan(t,norm(dense*solution.values()-rhs.values())/norm(rhs.values()),1e-9);
    metrics(i).name = c.name;
    metrics(i).ndof = rhs.Size;
    metrics(i).preconditioner_relative_error = norm(applied.values()-unpack(c.applied,complexCase))/norm(applied.values());
    metrics(i).solution_relative_error = norm(solution.values()-unpack(c.solution,complexCase))/norm(solution.values());
    metrics(i).relative_residual = norm(dense*solution.values()-rhs.values())/norm(rhs.values());
    stale = inv.nativeHandle();
    delete(inv);
    verifyError(t,@() radia.internal.callMex('ngsolve.matrix.info',stale),'radia:mex:Exception');
    delete(solution); delete(applied); delete(rhs);
end
path = getenv('RADIA_SPARSESOLV_METRICS');
if ~isempty(path)
    fid = fopen(path,'w');
    assert(fid >= 0,'Cannot write validation metrics');
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid,'%s',jsonencode(struct('cases',metrics,'oracle',rmfield(r,'cases'))));
end
end

function testRejectWrongSpace(t)
mesh = radia.ngsolve.Mesh.create(t.TestData.ref.mesh);
cleanup = onCleanup(@() delete(mesh));
space = radia.ngsolve.FESpace.create(mesh,"h1",1);
spaceCleanup = onCleanup(@() delete(space));
form = radia.ngsolve.BilinearForm.create(space,"mass");
formCleanup = onCleanup(@() delete(form));
a = form.matrix();
aCleanup = onCleanup(@() delete(a));
verifyError(t,@() radia.sparsesolv.AMS(a,space),'radia:mex:Exception');
end

function values = unpack(value,isComplex)
values = value.real;
if isComplex
    values = complex(values,value.imag);
end
end

function testExplicitPythonFallback(t)
value = radia.python.sparsesolv("has_compact_ams");
verifyEqual(t,value.backend,"python-fallback");
verifyTrue(t,logical(value.value));
end

function testElectromagnetBatchEntries(t)
r = radia.python.electromagnetValidation("static_electromagnet_three_engine_contract");
verifyEqual(t,r.backend,"python-fallback");
verifyNotEmpty(t,r.value);
r = radia.python.esrfExamples("get_esrf_example_spec",{int32(1)});
verifyEqual(t,double(r.value.number),1);
r = radia.python.staticElectromagnet("StaticElectromagnetMixedDomain", ...
    {{'air'}, {'iron','kelvin'}, {'iron'}});
verifyEqual(t,string(r.value.ground_boundary),"GND");
end
