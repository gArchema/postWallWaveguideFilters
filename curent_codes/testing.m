%calculating field
tic
xresolution = 600;
yrezolution = xresolution/2;
ystart = min(cylCoord(:, 2))+1.5*dx;
yfinal = max(cylCoord(:, 2))-1.5*dx;
rangey = abs(ystart)+abs(yfinal);
rangex = rangey*2;
xstart = min(cylCoord(:, 1))+4*dx;
xfinish = xstart+rangex;
x = linspace(xstart, xstart+rangex, xresolution);
y = linspace(ystart, yfinal, yrezolution);

[X, Y] = meshgrid(x, y);
E = exp(1i*kE*distance(xs, ys, X, Y).*cos(atan2(Y-ys, X-xs)-fii));
for p=1:L*N+D+LV*NVW+LH*NHW
    for m = -M:M
        E = E + A((p-1)*(2*M+1) + m+M+1)*besselh(m, 2, kE*distance(cylCoord(p, 1), cylCoord(p, 2), X, Y)).*exp(1i*m*atan2(Y - cylCoord(p, 2), X - cylCoord(p, 1)));
    end;
end;


%% internal field calculation <3


for p=1:L*N+D+LV*NVW+LH*NHW
    xcur = cylCoord(p, 1);
    ycur = cylCoord(p, 2);
    rcur = cylCoord(p, 3);
    if xcur+rcur >= xstart && xcur-rcur <= xfinish && ycur+rcur >= ystart && ycur-rcur <= yfinal
        Expos = xresolution*abs(xcur-xstart)/rangex;
        Eypos = yrezolution*abs(ycur-ystart)/rangey;
        mx_delta = round(rcur*xresolution/rangex);
        my_delta = round(rcur*yrezolution/rangey);
        for mx_pos = -mx_delta:mx_delta
            for my_pos = -my_delta:my_delta
                change_x_pos = ceil(Expos+mx_pos);
                change_y_pos = ceil(Eypos+my_pos);
                cur_dist = rcur*sqrt(mx_pos^2+my_pos^2)/mx_delta;
                if cur_dist <= rcur
                    if change_x_pos>0 && change_x_pos<xresolution && change_y_pos>0 && change_y_pos<yrezolution
                        E(change_y_pos, change_x_pos) = 0;
                        if rcur == aD
                            des = distance([xcur, ycur], [x(change_x_pos), y(change_y_pos)]);
                            fes = angle2P([xcur, ycur], [x(change_x_pos), y(change_y_pos)]);
                            for m = -M:M
                                b = T((p-1)*(2*M+1)+m+M+1, 1)*A((p-1)*(2*M+1)+m+M+1, 1);
                                E(change_y_pos, change_x_pos) = E(change_y_pos, change_x_pos) + b*besselj(m, cylCoord(p, 4)*des)*exp(1i*m*fes);
                            end;
                        end;
                    end;
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

    
    