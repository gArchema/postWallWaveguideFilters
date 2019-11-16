function [res]=angle2P(pf, ps)
    res = atan2((ps(2)-pf(2)), (ps(1)-pf(1)));
%     if res > pi
%         res = res - 2*pi;
%     end;
end