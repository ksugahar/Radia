function result=optunaFloatError(actual,expected)
% Real binary64 diagnostics; signed zeros compare equal, NaNs are separate.
actual=actual(:); expected=expected(:);
finite=isfinite(actual)&isfinite(expected);
a=actual(finite); b=expected(finite);
delta=abs(a-b);
relative=delta./abs(b);
relative(b==0)=NaN; % Undefined reference-relative error, never silently zero.
ka=ordered(a); kb=ordered(b);
ulps=max(ka,kb)-min(ka,kb);
result=struct('elements',numel(actual),'finite_elements',numel(a), ...
    'unequal_finite',nnz(a~=b),'nonfinite_mismatches', ...
    nnz(~finite & ~(actual==expected | (isnan(actual)&isnan(expected)))), ...
    'zero_reference_mismatches',nnz(b==0&a~=b), ...
    'max_absolute',0,'max_relative',NaN,'max_ulp','0', ...
    'absolute_pair',[],'relative_pair',[],'ulp_pair',[]);
if isempty(a), return; end
[result.max_absolute,i]=max(delta); result.absolute_pair=[a(i),b(i)];
if any(~isnan(relative))
    [result.max_relative,i]=max(relative,[],'omitnan');
    result.relative_pair=[a(i),b(i)];
end
[u,i]=max(ulps); result.max_ulp=char(string(u));
result.ulp_pair=[a(i),b(i)];
end

function key=ordered(values)
bits=reshape(typecast(values,'uint64'),[],1);
negative=bitget(bits,64)~=0;
key=bitset(bits,64);
key(negative)=bitcmp(bits(negative))+uint64(1);
end
