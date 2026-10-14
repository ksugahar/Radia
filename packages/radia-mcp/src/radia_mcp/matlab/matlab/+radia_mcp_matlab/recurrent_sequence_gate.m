function report = recurrent_sequence_gate(seed)
arguments
    seed (1,1) double {mustBeInteger,mustBeNonnegative} = 163
end
rng(seed,"twister"); t=30; d=4; x=randn(t,d); wf=.35*randn(d); wb=.35*randn(d); hf=zeros(t,d); hb=zeros(t,d);
for i=1:t, prev=zeros(1,d); if i>1,prev=hf(i-1,:);end, hf(i,:)=tanh(x(i,:)+prev*wf); end
for i=t:-1:1, prev=zeros(1,d); if i<t,prev=hb(i+1,:);end, hb(i,:)=tanh(x(i,:)+prev*wb); end
bidir=[hf hb]; forget=1./(1+exp(-x)); inputGate=1-forget; cell=zeros(1,d);
for i=1:t, cell=forget(i,:).*cell+inputGate(i,:).*tanh(x(i,:)); end
update=1./(1+exp(-x(end,:))); old=hf(end-1,:); candidate=tanh(x(end,:)); gru=update.*old+(1-update).*candidate;
exploding=1.2.^(1:t);
checks=struct("rnn_finite",all(isfinite(hf),"all"),"bilstm_shape",isequal(size(bidir),[t 2*d]), ...
    "lstm_gates_bounded",all(forget>0&forget<1,"all"),"cell_finite",all(isfinite(cell)), ...
    "gru_convex_update",all(gru>=min(old,candidate)-eps&gru<=max(old,candidate)+eps), ...
    "exploding_negative_control",exploding(end)>100);
ids=["O1PCh_aaprE","oxygME2UBFc","K8ktkhAEuLM","IcCIu5Gx6uA","NJdrYvYgaPM"];
topics=["bilstm","lstm","gru","rnn_usage","rnn_meaning"];
lessons=repmat(struct("video_id","","topic","","public_url",""),5,1);
for k=1:5, lessons(k)=struct("video_id",ids(k),"topic",topics(k),"public_url","https://www.youtube.com/watch?v="+ids(k)); end
report=struct("schema","radia-mcp.matlab.recurrent-sequence.v1","lessons",lessons, ...
    "seed",seed,"units",struct("sequence_feature","1"),"sequence_length",t, ...
    "bidirectional_shape",size(bidir),"exploding_control_final",exploding(end), ...
    "sequence_output_is_ground_truth",false,"promotion_requires_forward_solver",true, ...
    "checks",checks,"ok",all(structfun(@logical,checks)));
end
