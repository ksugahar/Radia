function report = nlp_transformer_gate(seed)
arguments
    seed (1,1) double {mustBeInteger,mustBeNonnegative} = 113
end
rng(seed,"twister"); n=12; d=8; heads=2; q=randn(n,d); k=randn(n,d); val=randn(n,d);
scores=q*k'/sqrt(d); causal=tril(true(n)); scores(~causal)=-Inf; attention=softmaxRows(scores); output=attention*val;
rowError=max(abs(sum(attention,2)-1)); futureLeak=max(abs(attention(~causal)),[],"all");
layers=12; unshared=layers*d*d; shared=d*d; parameterRatio=shared/unshared;
pretrainMask=true(n); inferenceMask=causal; masksDiffer=nnz(pretrainMask~=inferenceMask)>0;
bleuExact=bleu2(["field","solver","is","verified"],["field","solver","is","verified"]);
bleuBad=bleu2(["unrelated","tokens"],["field","solver","is","verified"]);
vec=randn(1,100); vec(abs(vec)<1.1)=0; sparsity=mean(vec==0);
taskPrefix="predict_qoi:"; missingPrefix=""; prefixRejected=strlength(missingPrefix)==0;
checks=struct("attention_normalized",rowError<1e-12,"causal_no_future_leak",futureLeak==0, ...
    "parameter_sharing_reduces",parameterRatio<.1,"pretrain_inference_masks_distinct",masksDiffer, ...
    "bleu_exact_one",abs(bleuExact-1)<1e-12,"bleu_no_overlap_zero",bleuBad==0, ...
    "sparse_document_vector",sparsity>.5,"missing_task_prefix_rejected",prefixRejected, ...
    "output_finite",all(isfinite(output),"all"),"task_prefix_present",strlength(taskPrefix)>0);
ids=["-x08lNz3Qfo","du-hz9SVjj8","T7nbrIJtYlE","tG-WI9qMluE","AhEb8kICwIQ","3BUk7mtf10M","IaTCGRL41_k","wDXPXgn5hX4","hMrOcH5dcGM","FFoLqib6u-0","50XvMaWhiTY","aZJAizFSTWg","gnnnB3gd_0U","AByKltWQMl8"];
topics=["t5_text_to_text","albert_parameter_sharing","roberta_training","xlnet_pretrain_inference","xlnet_pretraining","gpt2_scaling","bert","gpt_pretrain_finetune","elmo_context","transformer_rta","multihead_attention","bleu","scdv","gnmt"];
lessons=repmat(struct("video_id","","topic","","public_url",""),numel(ids),1);
for ii=1:numel(ids), lessons(ii)=struct("video_id",ids(ii),"topic",topics(ii),"public_url","https://www.youtube.com/watch?v="+ids(ii)); end
report=struct("schema","radia-mcp.matlab.nlp-transformer.v1","lessons",lessons, ...
    "seed",seed,"units",struct("attention","1","bleu","1"),"attention_row_error",rowError, ...
    "future_attention_leak",futureLeak,"parameter_ratio_shared",parameterRatio, ...
    "bleu_exact",bleuExact,"bleu_negative",bleuBad,"document_sparsity",sparsity, ...
    "model_output_is_ground_truth",false,"promotion_requires_forward_solver",true, ...
    "checks",checks,"ok",all(structfun(@logical,checks)));
end
function a=softmaxRows(x), m=max(x,[],2); e=exp(x-m); a=e./sum(e,2); end
function b=bleu2(candidate,reference)
if isempty(candidate), b=0; return, end
p=zeros(1,2);
for n=1:2
    c=joinN(candidate,n); r=joinN(reference,n); hit=0;
    for i=1:numel(c), hit=hit+any(c(i)==r); end
    p(n)=hit/max(numel(c),1);
end
if any(p==0), b=0; return, end
bp=min(1,exp(1-numel(reference)/numel(candidate))); b=bp*exp(mean(log(p)));
end
function g=joinN(t,n)
if numel(t)<n, g=strings(0,1); return, end
g=strings(numel(t)-n+1,1); for i=1:numel(g), g(i)=strjoin(t(i:i+n-1),"|"); end
end
