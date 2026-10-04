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
pythonExecutable = string(getenv('RADIA_PYTHON_EXECUTABLE'));
if strlength(pythonExecutable) == 0
    pythonExecutable = "python";
end
command = '"' + pythonExecutable + '" -u -X faulthandler "' + ...
    string(helper) + '" "' + string(folder) + '" 2>&1';
[status, output] = radia.internal.runPythonProcess(command);
assert(status == 0, 'Python oracle failed (exit %d): %s', status, output);
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
    if strcmp(c.name,'ams')
        % A real AMS matvec must overwrite every value in a caller-owned
        % matrix vector. This exercises the parallel-vector output type used
        % by the MEX gateway and catches allocator NaNs leaking into A*x.
        poisoned = a.vector();
        poisoned.setValues(nan(rhs.Size,1));
        pre.matvecInto(rhs,poisoned);
        verifyEqual(t,poisoned.values(),unpack(c.applied,false), ...
            RelTol=1e-10,AbsTol=1e-10);
        delete(poisoned);
    end
    inv = radia.sparsesolv.COCR(a,pre,Tolerance=1e-11);
    dense = sparse(a);
    delete(pre); delete(a); delete(aux); delete(form); delete(space); delete(mesh);
    if complexCase
        delete(cf); delete(cs); delete(weight);
    end
    repeatCount = 1;
    if strcmp(c.name,'ams')
        % The real AMS path owns TaskManager-parallel work vectors. Exercise
        % repeated direct matvec calls so an omitted gateway TaskManager region
        % cannot intermittently return the zero initial guess.
        repeatCount = 3;
    end
    for repeatIndex = 1:repeatCount
        solution = inv.matvec(rhs);
        verifyEqual(t,solution.values(),unpack(c.solution,complexCase),RelTol=1e-8,AbsTol=1e-8);
        verifyLessThan(t,norm(dense*solution.values()-rhs.values())/norm(rhs.values()),1e-9);
        if repeatIndex < repeatCount
            delete(solution);
        end
    end
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

function [space, form, a, cleanup] = iccgSystem(t, complexCase)
% H1 order-2 Dirichlet Laplace system matching sparsesolv_python_reference.py.
mesh = radia.ngsolve.Mesh.create(t.TestData.ref.mesh);
space = radia.ngsolve.FESpace.create(mesh,"h1",2,Dirichlet=".*",Complex=complexCase);
if complexCase
    weight = radia.ngsolve.CoefficientFunction.constant(1+0.5i);
    form = radia.ngsolve.BilinearForm.createFromCoefficient(space,"stiffness",weight);
else
    weight = [];
    form = radia.ngsolve.BilinearForm.create(space,"stiffness");
end
a = form.matrix();
cleanup = onCleanup(@() cellfun(@delete, {a, form, space, weight, mesh}));
end

function options = iccgOptions(name)
switch name
    case "default", options = {};
    case "fixed_unscaled", options = {'Shift',1.2,'AutoShift',false,'DiagonalScaling',false};
    case "iteration_limit", options = {'MaxIterations',4};
    case "strict_stagnation", options = {'DivergenceThreshold',1.0,'DivergenceCount',0};
    case "abmc", options = {'UseABMC',true};
    otherwise, error('unknown ICCG case %s', name);
end
end

function testHermitianICCG(t)
mesh = radia.ngsolve.Mesh.create(t.TestData.ref.mesh);
space = radia.ngsolve.FESpace.create(mesh,"h1",2,Dirichlet=".*",Complex=true);
form = radia.ngsolve.BilinearForm.create(space,"stiffness");
a = form.matrix();
rhs = a.vector();
cleanup = onCleanup(@() cellfun(@delete,{rhs,a,form,space,mesh})); %#ok<NASGU>
free = logical(radia.internal.callMex('ngsolve.fespace.free_dofs',space.nativeHandle()));
b = exp(0.3i*(1:rhs.Size)'); b(~free) = 0;
rhs.setValues(b);
A = sparse(a);
reference = A(free,free)\b(free);
for abmc = [false,true]
    [x,info] = radia.sparsesolv.ICCG(a,rhs,Conjugate=true,UseABMC=abmc,Tolerance=1e-11);
    verifyTrue(t,info.converged);
    values = x.values();
    verifyLessThan(t,norm(values(free)-reference)/norm(reference),1e-9);
    verifyLessThan(t,norm(A(free,free)*values(free)-b(free))/norm(b(free)),1e-10);
    delete(x);
end
[~,~,bad,badCleanup] = iccgSystem(t,true); %#ok<ASGLU>
zero = bad.vector(); zero.setZero();
zeroCleanup = onCleanup(@() delete(zero)); %#ok<NASGU>
verifyError(t,@() radia.sparsesolv.ICCG(bad,zero,Conjugate=true),'radia:mex:Exception');
end

function testICCGPythonParity(t)
% The MEX entry runs the same C++ ICCG as the Python module; results agree to
% thread-order rounding, with the same option meaning and result fields.
cases = t.TestData.ref.iccg;
for i = 1:numel(cases)
    c = cases(i);
    fprintf('Checking ICCG %s (complex=%d)\n', c.name, c.complex);
    [space, ~, a, cleanup] = iccgSystem(t, logical(c.complex)); %#ok<ASGLU>
    rhs = a.vector();
    rhs.setValues(unpack(c.rhs, logical(c.complex)));
    opts = iccgOptions(string(c.name));
    [x, info] = radia.sparsesolv.ICCG(a, rhs, opts{:}, SaveResidualHistory=true);
    verifyEqual(t, info.converged, logical(c.converged));
    verifyLessThanOrEqual(t, abs(info.iterations - c.iterations), 1);
    verifyEqual(t, info.actual_shift, c.actual_shift, AbsTol=1e-12);
    verifyEqual(t, numel(info.residual_history), info.iterations + 1);
    verifyEqual(t, info.residual_history(1), c.residual_history(1), RelTol=1e-12);
    k = min(numel(info.residual_history), numel(c.residual_history));
    verifyEqual(t, info.residual_history(1:min(k,5)), c.residual_history(1:min(k,5)), RelTol=1e-8);
    ref = unpack(c.solution, logical(c.complex));
    verifyLessThan(t, norm(x.values() - ref) / max(norm(ref), realmin), 1e-6);
    % true_residual is the original free-DOF residual of the returned x
    free = logical(radia.internal.callMex('ngsolve.fespace.free_dofs', space.nativeHandle()));
    A = sparse(a);
    r = rhs.values() - A * x.values();
    verifyEqual(t, info.true_residual, norm(r(free)) / norm(rhs.values()), RelTol=1e-6, AbsTol=1e-14);
    verifyEqual(t, info.final_residual, info.residual_history(info.best_iteration + 1), RelTol=0);
    delete(x); delete(rhs); clear cleanup
end
end

function testICCGContractAndLifecycle(t)
[~, ~, a, cleanup] = iccgSystem(t, false); %#ok<ASGLU>
rhs = a.vector();
rhs.setValues(cos(0.3 * (1:rhs.Size)'));
before = radia.apiInfo();
[x, info] = radia.sparsesolv.ICCG(a, rhs);
verifyTrue(t, info.converged);
verifyEqual(t, info.best_iteration, info.iterations);
verifyEmpty(t, info.residual_history);  % default SaveResidualHistory=false, as in Python
% an initial guess that already solves the system returns before the IC factor
[x2, info2] = radia.sparsesolv.ICCG(a, rhs, InitialGuess=x, Tolerance=1e-6);
verifyTrue(t, info2.converged);
verifyEqual(t, info2.iterations, 0);
verifyEqual(t, info2.actual_shift, 0);
verifyEqual(t, x2.values(), x.values());
% an exactly zero right-hand side returns x = 0
zero = a.vector(); zero.setZero();
[x0, info0] = radia.sparsesolv.ICCG(a, zero, InitialGuess=x);
verifyTrue(t, info0.converged);
verifyEqual(t, info0.iterations, 0);
verifyEqual(t, x0.values(), zeros(rhs.Size, 1));
% the iteration limit returns the best iterate, and history has iterations + 1 entries
[x4, info4] = radia.sparsesolv.ICCG(a, rhs, MaxIterations=2, SaveResidualHistory=true);
verifyFalse(t, info4.converged);
verifyEqual(t, info4.iterations, 2);
verifyEqual(t, info4.final_residual, min(info4.residual_history));
% option and input validation fail loudly
verifyError(t, @() radia.sparsesolv.ICCG(a, rhs, Shift=0.5), 'MATLAB:validators:mustBeGreaterThanOrEqual');
% Direct native calls request both outputs so they reach option validation,
% not the arity check; the complete struct is the positive control.
good = struct('tolerance',1e-8,'max_iterations',0,'shift',1,'auto_shift',true, ...
    'diagonal_scaling',true,'save_best_result',true,'save_residual_history',false, ...
    'divergence_check',true,'divergence_threshold',10,'divergence_count',10, ...
    'use_abmc',false,'abmc_block_size',4,'abmc_num_colors',4,'abmc_reorder_spmv',false, ...
    'abmc_use_rcm',false,'conjugate',false);
[hc, infoc] = rawIccg(a.nativeHandle(), rhs.nativeHandle(), uint64(0), good);
verifyTrue(t, infoc.converged);
radia.internal.callMex('ngsolve.vector.destroy', hc);
verifyError(t, @() rawIccg(a.nativeHandle(), rhs.nativeHandle(), uint64(0), ...
    rmfield(good, 'divergence_count')), 'radia:mex:Exception');          % missing option
withExtra = good; withExtra.unknown_option = true;
verifyError(t, @() rawIccg(a.nativeHandle(), rhs.nativeHandle(), uint64(0), withExtra), ...
    'radia:mex:Exception');                                              % unknown option
negative = good; negative.max_iterations = -1;
verifyError(t, @() rawIccg(a.nativeHandle(), rhs.nativeHandle(), uint64(0), negative), ...
    'radia:mex:Exception');                                              % invalid value
lowShift = good; lowShift.shift = 0.5;
verifyError(t, @() rawIccg(a.nativeHandle(), rhs.nativeHandle(), uint64(0), lowShift), ...
    'radia:mex:Exception');                                              % native shift check
verifyError(t, @() radia.internal.callMex('sparsesolv.iccg', a.nativeHandle(), ...
    rhs.nativeHandle(), uint64(0), good), 'radia:mex:Exception');        % one output
[otherSpace, ~, other, otherCleanup] = iccgSystem(t, true); %#ok<ASGLU>
otherRhs = other.vector();
verifyError(t, @() rawIccg(a.nativeHandle(), otherRhs.nativeHandle(), uint64(0), good), ...
    'radia:mex:Exception');                                              % complex rhs, real matrix
short = radia.ngsolve.Vector.fromNativeHandle( ...
    radia.internal.callMex('ngsolve.matrix.vector', a.nativeHandle()));
staleRhs = short.nativeHandle();
delete(short);
verifyError(t, @() rawIccg(a.nativeHandle(), staleRhs, uint64(0), good), ...
    'radia:mex:Exception');                                              % stale rhs handle
verifyError(t, @() rawIccg(a.nativeHandle(), rhs.nativeHandle(), staleRhs, good), ...
    'radia:mex:Exception');                                              % stale initial guess
delete(otherRhs); clear otherCleanup
nanRhs = a.vector(); nanRhs.setValues(nan(rhs.Size, 1));
verifyError(t, @() radia.sparsesolv.ICCG(a, nanRhs), 'radia:mex:Exception');
% the solution vector keeps its matrix alive; every returned handle is released
values = x.values();
delete(a);
verifyEqual(t, x.values(), values);
stale = x.nativeHandle();
delete(nanRhs); delete(zero); delete(x0); delete(x4); delete(x2); delete(x);
verifyError(t, @() radia.internal.callMex('ngsolve.vector.info', stale), 'radia:mex:Exception');
after = radia.apiInfo();
verifyEqual(t, after.handle_count, before.handle_count - 1);  % only the deleted matrix
delete(rhs);
clear cleanup
end

function [h, info] = rawIccg(varargin)
[h, info] = radia.internal.callMex('sparsesolv.iccg', varargin{:});
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

function testRejectHigherOrderHCurl(t)
% Order 2 with NoGrads=true keeps one gradient column per vertex, so only an
% explicit order or structure check rejects it. Python raises RuntimeError for
% the same space (test_ams_rejects_higher_order_space); the routes must agree.
mesh = radia.ngsolve.Mesh.create(t.TestData.ref.mesh);
cleanup = onCleanup(@() delete(mesh));
space = radia.ngsolve.FESpace.create(mesh,"hcurl",2,NoGrads=true);
spaceCleanup = onCleanup(@() delete(space));
form = radia.ngsolve.BilinearForm.create(space,"mass");
formCleanup = onCleanup(@() delete(form));
a = form.matrix();
aCleanup = onCleanup(@() delete(a));
verifyError(t,@() radia.sparsesolv.AMS(a,space),'radia:mex:Exception');
verifyError(t,@() radia.sparsesolv.AMS(a,space,Complex=true),'radia:mex:Exception');
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
