function report = topology_information_gate(seed)
arguments
    seed (1,1) double {mustBeInteger,mustBeNonnegative} = 227
end
rng(seed,"twister"); a=1.2;b=-.3;c=.8; tropicalLeft=a+min(b,c); tropicalRight=min(a+b,a+c);
points=[-.1 0;.1 0;3 0;3.2 0]; epsValues=[.15 .3 2 3.5]; components=zeros(size(epsValues));
for k=1:numel(epsValues), components(k)=componentCount(points,epsValues(k)); end
p=[.3;-.2]; r2=sum(p.^2); sphere=[2*p/(1+r2);(r2-1)/(1+r2)]; pBack=sphere(1:2)/(1-sphere(3)); chartError=norm(pBack-p);
prob=.37; fisher=1/(prob*(1-prob)); badProb=1.2; badRejected=badProb<=0||badProb>=1;
checks=struct("tropical_distributive",abs(tropicalLeft-tropicalRight)<1e-12, ...
    "components_monotone",all(diff(components)<=0),"manifold_chart_roundtrip",chartError<1e-12, ...
    "fisher_metric_positive",fisher>0,"invalid_probability_rejected",badRejected);
ids=["DH8_VXXoKz4","5x4r_deONQM","CiPHQVkaXNI","UhQ80ajAQHY"];
topics=["tropical_geometry","topological_data_analysis","manifold_intrinsic","information_geometry"];
lessons=repmat(struct("video_id","","topic","","public_url",""),4,1);
for k=1:4, lessons(k)=struct("video_id",ids(k),"topic",topics(k),"public_url","https://www.youtube.com/watch?v="+ids(k)); end
report=struct("schema","radia-mcp.matlab.topology-information.v1","lessons",lessons, ...
    "seed",seed,"units",struct("coordinate","1","metric","1"),"component_counts",components, ...
    "chart_roundtrip_error",chartError,"fisher_metric",fisher,"geometry_is_ground_truth",false, ...
    "promotion_requires_forward_solver",true,"checks",checks,"ok",all(structfun(@logical,checks)));
end
function n=componentCount(points,epsilon)
adj=squareformLocal(points)<=epsilon; seen=false(size(points,1),1);n=0;
for i=1:numel(seen),if seen(i),continue,end,n=n+1;queue=i;seen(i)=true;while ~isempty(queue),j=queue(1);queue(1)=[];nei=find(adj(j,:)&~seen');seen(nei)=true;queue=[queue nei];end,end
end
function d=squareformLocal(p), n=size(p,1);d=zeros(n);for i=1:n,for j=1:n,d(i,j)=norm(p(i,:)-p(j,:));end,end,end
