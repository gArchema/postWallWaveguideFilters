for p=1:L*N+D+LV*NVW+LH*NHW
    xcur = cylCoord(p, 1);
    ycur = cylCoord(p, 2);
    rcur = cylCoord(p, 3);
    if xcur > xstart && xcur < xfinish && ycur > ystart && ycur < yfinal
        Expos = xresolution*abs(xcur-xstart)/rangex;
        Eypos = yrezolution*abs(ycur-ystart)/rangey;
        mx_delta = floor(rcur*xresolution/rangex);
        my_delta = floor(rcur*yrezolution/rangey);
        for mx_pos = -mx_delta:mx_delta
            for my_pos = -my_delta:my_delta
                change_x_pos = Expos+mx_pos;
                change_y_pos = Eypos+my_pos;
                if change_x_pos>0 && change_x_pos<rangex && change_y_pos>0 && change_y_pos<rangey
                    E(change_x_pos, change_y_pos) = 0;
                end;
            end;
        end;
    end;
end;

figure;
imagesc(x, y, abs(real(E)));
colorbar;
hold on;
axis([xstart xstart+rangex ystart yfinal]);
viscircles([cylCoord(:,1), cylCoord(:,2)], cylCoord(:, 3), 'LineStyle', ':', 'LineWidth', 1/70);
colormap 'hot'
toc