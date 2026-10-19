function result=diagnoseNearRootNonlinear(mesh,cpu,preconditioner)
% Probe globalization near an independently computed CPU solution, not a speed test.
[reference,~]=radia.gpu.readHDivH5(cpu);
if nargin<3,preconditioner='mass-polynomial-graph';end
assert(ismember(preconditioner,{'mass-polynomial','mass-polynomial-graph'}));
opts={'NonlinearTolerance',1e-12,'LinearTolerance',1e-12,'LinearMaxIterations',4000,'MaxIterations',8,'MaterialBudgetBytes',2^30,'CacheBytes',ceil(1.6*2^30),'MaterialChunkElements',128,'InexactNewton',true,'Preconditioner',preconditioner};
obj=radia.gpu.PreparedHDivNonlinear(mesh,2,[0 0 1e5],reference.table,opts{:});owner=onCleanup(@()delete(obj));
m=reference.coefficients(:);direction=sin((1:numel(m))');direction=direction/norm(direction);cases={};
for relative=[0 1e-10 -1e-10 1e-9 -1e-9]
    row=struct('initial_relative_perturbation',relative,'completed',false);
    try
        s=obj.solve(1,m+relative*norm(m)*direction);row.completed=true;
        row.nonlinear_relative_residual=s.nonlinear_relative_residual;
        row.cpu_error=norm(s.coefficients-m)/norm(m);row.history=s.history;
        assert(row.cpu_error<=1e-8 && s.nonlinear_relative_residual<=1e-12);
    catch ex
        row.error=ex.message;row.identifier=ex.identifier;
        marker='diagnostics=';where=strfind(ex.message,marker);
        if ~isempty(where),row.diagnostics=jsondecode(ex.message(where(1)+numel(marker):end));end
        cases{end+1}=row;break
    end
    cases{end+1}=row;
end
liveBefore=radia_cuda_hmatrix('live_count');clear owner obj;wait(gpuDevice);
assert(radia_cuda_hmatrix('live_count')==0);
accepted=numel(cases)==5 && all(cellfun(@(row)row.completed,cases));
result=struct('accepted',accepted,'preconditioner',preconditioner,'diagnostic_max_iterations',8,'scope','Bounded near-root regression; independent CPU initial state with declared tiny perturbations; no performance acceptance','cases',{cases},'live_before_release',liveBefore,'live_count',radia_cuda_hmatrix('live_count'));
end
