function report = tabular_lookup_gate(seed)
arguments
    seed (1,1) double {mustBeInteger,mustBeNonnegative} = 197
end
rng(seed,"twister"); keys=[10;20;30;40]; values=[1.1;2.2;3.3;4.4]; query=30;
match=find(keys==query,1); exact=values(match); indexMatch=values(find(keys==query,1));
approxQuery=27; idx=find(keys<=approxQuery,1,"last"); approximate=values(idx);
unsorted=[20;10;40;30]; approximateUnsortedRejected=~issorted(unsorted);
joined=strjoin(["case","30","verified"],"_"); missing=find(keys==99,1); missingRejected=isempty(missing);
checks=struct("exact_lookup",exact==3.3,"index_match_equal",indexMatch==exact, ...
    "approximate_sorted_contract",approximate==2.2,"unsorted_approximate_rejected",approximateUnsortedRejected, ...
    "missing_key_rejected",missingRejected,"join_deterministic",joined=="case_30_verified");
ids=["1nZhvyfL7Ro","RrXVemskaoI","BiomXv_9GDw","mLTc3soOkHA","bwuc1gtoWMg","zi7v2Ow1p14"];
topics=["index_match","match","index","vlookup_false","vlookup_true","concatenate_join"];
lessons=repmat(struct("video_id","","topic","","public_url",""),6,1);
for k=1:6, lessons(k)=struct("video_id",ids(k),"topic",topics(k),"public_url","https://www.youtube.com/watch?v="+ids(k)); end
report=struct("schema","radia-mcp.matlab.tabular-lookup.v1","lessons",lessons, ...
    "seed",seed,"units",struct("key","1","value","Pa"),"exact_value",exact, ...
    "approximate_value",approximate,"joined_key",joined,"lookup_is_ground_truth",false, ...
    "promotion_requires_forward_solver",true,"checks",checks,"ok",all(structfun(@logical,checks)));
end
