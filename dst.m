function [ res ] = dst( ps, pf )
res = distance(ps, pf);
if ps(1)>pf(1)
    res = -1*res;
end;

end

