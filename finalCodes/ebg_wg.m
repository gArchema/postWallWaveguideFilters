%init
tic
format long;
dielE = 1;              %permittivity of environment
dielD = 11.56;              %permittivity of Dielectric cylinders in structure
dielDWall = 100000;             %permittivity of Dielectric wall
dx = 1e-6;              %period along x
dy = dx;                %period along y
aD = 0.18*dx;               %radius of dielectric cylinder
aP = 0.65e-3;               %radius of pec cylinder
dyhv = 2*aD;             %period of vertical and horizontal wall's cylinders
wg = 2*dy;               %width of waveguide
fii = pi;                %incident angle to opposite direction of x axis
M=3;
L=2;                    %amount of layers
N = 20*2;               %amount of cylinders
LV=1;                   %amount of layers for vertical wall
NVW = 12*2;             %amount of cylinders for vertical wall
LH = 1;                 %amount of layers for horizontal wall
NHW = floor(dx*N/2/(2*aD)+3)*2;      %amount of cylinders for horizontal wall
cylCoord = ones(L*N+LV*NVW+LH*NHW, 5);         %imformative matrix for cylinders: 1st column x coordinates, 2 - y coordinates, 3 - radius, 4 - wavenumber, 5 - permittivity
D = input('enter amount of defects: ');        %amount of defect cylinders
normFreq = 0.425;                              %normalized frequency
freq = physconst('LightSpeed')*normFreq/dx;    %frequency
kE = 2*pi*freq*sqrt(dielE)/physconst('LightSpeed');         %wavenumber of environment
kD = 2*pi*freq*sqrt(dielD)/physconst('LightSpeed');         %wavenumber of Dielectric cylinders in structure
kDW = 2*pi*freq*sqrt(dielDWall)/physconst('LightSpeed');    %wavenumber of Dielectric wall
xs = 0; ys = 0;                                         %source coordinates

%creating coordinate matrix.
for ll = 1:LV
    distFSource = 2*dx;
    cvy= wg/2 + (LV-ll)*aD;
    cvxi = xs + 2*dx + 4*(ll-1)*aD;
    for p = (ll-1)*NVW+1:2:(ll-1)*NVW+NVW-1
        cylCoord(p, :) = [cvxi, cvy, aD, kDW, dielDWall];
        cylCoord(p+1, :) = [cvxi, -cvy, aD, kDW, dielDWall];
        cvy = cvy+dyhv;
    end;
end;

for l = 1:L
    cx = cvxi + dx;
    cy = (l-1)*dy + wg/2;
    for p = LV*NVW + (l-1)*N + 1:2:LV*NVW + (l-1)*N + N-1
        cylCoord(p, :) = [cx, cy, aD, kD, dielD];
        cylCoord(p+1, :) = [cx, -cy, aD, kD, dielD];
        cx = cx+dx;
    end;
end;

for lh = 1:LH
    cxh = cvxi + 2*aD;
    cyh = wg/2 + (NVW/2-1)*2*aD-(lh-1)*aD*2;
    for p = LV*NVW + L*N + (lh-1)*NHW + 1:2:LV*NVW + L*N + (lh-1)*NHW + NHW - 1
        cylCoord(p, :) = [cxh, cyh, aD, kDW, dielDWall];
        cylCoord(p+1, :) = [cxh, -cyh, aD, kDW, dielDWall];
        cxh = cxh+2*aD;
    end;
end;
    
for p = L*N+LV*NVW+LH*NHW+1:L*N+LV*NVW+LH*NHW+D
    r = input('enter r: ');
    xd = input('enter x: ');
    yd = input('enter y: ');
    cylCoord(p, :) = [xd*1e-3, yd*1e-3, r*1e-3];
end;

% %find A matrix
C = zeros((L*N+D+LV*NVW+LH*NHW)*(2*M+1), (L*N+D+LV*NVW+LH*NHW)*(2*M+1));
F = zeros((L*N+D+LV*NVW+LH*NHW)*(2*M+1), 1);
T = zeros((L*N+D+LV*NVW+LH*NHW)*(2*M+1), 1);
for p=1:L*N+D+LV*NVW+LH*NHW
    df = distance([xs, ys], [cylCoord(p, 1), cylCoord(p,2)]);
    ff = angle2P([xs, ys], [cylCoord(p, 1), cylCoord(p,2)]);
    for m = -M:M
        fc = (1i)^(m)*exp(-1i*m*fii)*exp(1i*kE*df*cos(ff-fii));
        F((p-1)*(2*M+1) + m+M+1, 1) = fc;
    end;
    for q=1:L*N+D+LV*NVW+LH*NHW
        dc = distance([cylCoord(p, 1), cylCoord(p, 2)], [cylCoord(q, 1), cylCoord(q, 2)]);
        fc = angle2P([cylCoord(p, 1), cylCoord(p, 2)], [cylCoord(q, 1), cylCoord(q, 2)]);
        for m = -M:M
            for n = -M:M
                if p~=q
                    c = -besselh(m-n, 2, kE*dc)*exp(-1*(1i)*(m-n)*fc);
                elseif m==n
                     if cylCoord(p, 3) == aP
                          c = -besselh(m, 2, kE*cylCoord(p, 3))/besselj(m, kE*cylCoord(p, 3));
                          t = 0;
                     else
                         syms rk real; 
                         h1 = besselh(m, 2, kE*cylCoord(p, 3));
                         ht = diff(besselj(m, rk)-1i*bessely(m, rk));
                         hd1 = subs(ht, rk, kE*cylCoord(p, 3));
                         j1 = besselj(m, kE*cylCoord(p, 3));
                         j2 = besselj(m, cylCoord(p, 4)*cylCoord(p, 3));
                         jt = diff(besselj(m, rk));
                         jd1 = subs(jt, rk, kE*cylCoord(p, 3));
                         jd2 = subs(jt, rk, cylCoord(p, 4)*cylCoord(p, 3));
                         c = -(cylCoord(p, 5)*kE*h1*jd2-dielE*cylCoord(p, 4)*hd1*j2)/(cylCoord(p, 5)*kE*j1*jd2-dielE*cylCoord(p, 4)*jd1*j2);
                         t = dielE*cylCoord(p, 4)*(h1*jd1-hd1*j1)/(dielE*cylCoord(p, 4)*jd1*j2-cylCoord(p, 5)*kE*j1*jd2);
                     end;
                     T((p-1)*(2*M+1)+m+M+1, 1) = t;
                else
                    c=0;
                end;
                C((p-1)*(2*M+1) + m+M+1, (q-1)*(2*M+1) + n+M+1) = c;
            end;
        end;
    end;
end;
A = C\F;
toc
tic

%% calculating field
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

% external field calculation
for p=1:L*N+D+LV*NVW+LH*NHW
    for m = -M:M
        E = E + A((p-1)*(2*M+1) + m+M+1)*besselh(m, 2, kE*distance(cylCoord(p, 1), cylCoord(p, 2), X, Y)).*exp(1i*m*atan2(Y - cylCoord(p, 2), X - cylCoord(p, 1)));
    end;
end;

% internal field calculation
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
    
    

    
    
    