function report = tabular_transform_gate(seed)
arguments
    seed (1,1) double {mustBeInteger,mustBeNonnegative} = 199
end
rng(seed,"twister"); text="case_030_pressure"; parts=split(text,"_"); left=extractBefore(text,5); right=extractBetween(text,10,17); mid=extractBetween(text,6,8); pos=strfind(text,"pressure");
values=[1 2 3 4]; weights=[.5 1 1.5 2]; sumValue=sum(values); sumProduct=sum(values.*weights); mask=values>=3; sumIf=sum(values(mask)); cumulative=cumsum(values);
missingFind=isempty(strfind(text,"temperature")); badSplit=split("","_");
checks=struct("split_fields",numel(parts)==3,"left_right",left=="case"&&right=="pressure", ...
    "mid_extract",mid=="030","len_find",strlength(text)==17&&pos==10,"sum",sumValue==10, ...
    "sumproduct",sumProduct==15,"sumifs",sumIf==7,"cumulative",isequal(cumulative,[1 3 6 10]), ...
    "missing_find_rejected",missingFind,"empty_split_detected",all(strlength(badSplit)==0));
ids=["3d-8nsF4eNA","x6G7nh7ENEk","0arrBAxw-zg","6_JKWGsg5YM","CLzv3R8bdag","oVHBSRrLrxw","TJcbEdz8xII","R9Qv08AKNro","_4jAL3-yqt4"];
topics=["split","mid","len_find","left_right","sum_summary","sumproduct","sumifs","cumulative_sum","sum"];
lessons=repmat(struct("video_id","","topic","","public_url",""),9,1);
for k=1:9, lessons(k)=struct("video_id",ids(k),"topic",topics(k),"public_url","https://www.youtube.com/watch?v="+ids(k)); end
report=struct("schema","radia-mcp.matlab.tabular-transform.v1","lessons",lessons, ...
    "seed",seed,"units",struct("value","Pa","key","string"),"parts",parts',"sum",sumValue, ...
    "sumproduct",sumProduct,"conditional_sum",sumIf,"tabular_result_is_ground_truth",false, ...
    "promotion_requires_forward_solver",true,"checks",checks,"ok",all(structfun(@logical,checks)));
end
