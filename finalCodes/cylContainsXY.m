function [ res ] = cylContainsXY( cylCoord, x, y )
    res = -1;
    for pos = 1:size(cylCoord, 1)
 %       display(x); display(y); display(cylCoord(pos, 1)); display(cylCoord(pos, 2));
        d = distance(x, y, cylCoord(pos, 1), cylCoord(pos, 2));
        if d < cylCoord(pos, 3)
          %  display(1)
            res = pos;
            break;
        end;
    end;
        

end

