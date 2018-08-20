function s=mergestruct(s1,s2)
% -----------------------------------------------------
% -- Fast mmWave Ray Tracing Simulator (v0.2)
% -- 2018 (c) junquan.deng@aalto.fi
% -----------------------------------------------------
if(~isstruct(s1) || ~isstruct(s2))
    error('input parameters contain non-struct');
end
if(length(s1)>1 || length(s2)>1)
    error('can not merge struct arrays');
end
fn=fieldnames(s2);
s=s1;
for i=1:length(fn)              
    s=setfield(s,fn{i},getfield(s2,fn{i}));
end

