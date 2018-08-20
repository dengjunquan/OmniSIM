function val=jsonopt(key,default,varargin)
% -----------------------------------------------------
% -- Fast mmWave Ray Tracing Simulator (v0.2)
% -- 2018 (c) junquan.deng@aalto.fi
% -----------------------------------------------------
val=default;
if(nargin<=2) return; end
opt=varargin{1};
if(isstruct(opt))
    if(isfield(opt,key))
       val=getfield(opt,key);
    elseif(isfield(opt,lower(key)))
       val=getfield(opt,lower(key));
    end
end
