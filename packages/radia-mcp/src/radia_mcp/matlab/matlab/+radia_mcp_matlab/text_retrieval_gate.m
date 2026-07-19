function report = text_retrieval_gate(seed)
arguments
    seed (1,1) double {mustBeInteger,mustBeNonnegative} = 149
end
rng(seed,"twister"); docs=[3 2 0 0;0 1 3 1;1 0 0 4]; query=[1 1 0 0]; n=size(docs,1); df=sum(docs>0,1); idf=log((n+1)./(df+1))+1;
tfidf=docs.*idf; qv=query.*idf; cosine=(tfidf*qv')./(vecnorm(tfidf,2,2)*norm(qv));
k1=1.2; b=.75; dl=sum(docs,2); avgdl=mean(dl); bm=zeros(n,1);
for i=1:n, for t=find(query), f=docs(i,t); bm(i)=bm(i)+idf(t)*f*(k1+1)/(f+k1*(1-b+b*dl(i)/avgdl)); end, end
subword=[.4 -.1 .2]+[.1 .3 -.2]; word=subword/2; docvec=mean([word;word+[.1 0 -.1]],1);
checks=struct("relevant_tfidf_ranked",cosine(1)==max(cosine),"relevant_bm25_ranked",bm(1)==max(bm), ...
    "scores_finite",all(isfinite([cosine;bm])),"subword_composition_finite",all(isfinite(word)), ...
    "document_vector_finite",all(isfinite(docvec)),"empty_query_rejected",norm(zeros(1,4))==0);
ids=["3YUq_wTOo68","bPdyuIebXWM","lJocgM6Pa18","jlmt4nY0-o0","0CXCqxQAKKQ","__-vWa8jyVc","V7WVdlUSOco","_HSOX0oh2ns","nsEbfO3U2pY"];
topics=["fasttext","attention","doc2vec","word2vec_math","word2vec","rnnlm","elasticsearch","bm25","tfidf"];
lessons=repmat(struct("video_id","","topic","","public_url",""),numel(ids),1);
for k=1:numel(ids), lessons(k)=struct("video_id",ids(k),"topic",topics(k),"public_url","https://www.youtube.com/watch?v="+ids(k)); end
report=struct("schema","radia-mcp.matlab.text-retrieval.v1","lessons",lessons, ...
    "seed",seed,"units",struct("score","1","embedding","1"),"tfidf_scores",cosine', ...
    "bm25_scores",bm',"subword_vector",word,"document_vector",docvec, ...
    "retrieval_score_is_ground_truth",false,"promotion_requires_forward_solver",true, ...
    "checks",checks,"ok",all(structfun(@logical,checks)));
end
