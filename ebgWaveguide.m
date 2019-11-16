%init
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
L=1;                    %amount of layers
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
xs = 0; ystart = 0;                                         %source coordinates

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
    df = distance([xs,ystart], [cylCoord(p, 1), cylCoord(p,2)]);
    ff = angle2P([xs, ystart], [cylCoord(p, 1), cylCoord(p,2)]);
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
        
%calculating field
xresolution = 50;
yrezolution = 30;
ystart = min(cylCoord(:, 2))+1.5*dx;
yfinal = max(cylCoord(:, 2))-1.5*dx;
rangey = abs(ystart)+abs(yfinal);
rangex = rangey*2;
xstart = min(cylCoord(:, 1))+4*dx;
x = linspace(xstart, xstart+rangex, xresolution);
y = linspace(ystart, yfinal, yrezolution);
E = zeros(size(y, 2), size(x, 2));
for yp = 1:size(y, 2)/2
    for xp = 1:size(x, 2)
        if cylContainsXY(cylCoord, x(xp), y(yp))==-1% && (y(yp)<wg/2 && y(yp)>-wg/2)
            dei = distance([xs, ystart], [x(xp), y(yp)]);
            fei = angle2P([xs, ystart], [x(xp), y(yp)]);
            ei = exp(1i*kE*dei*cos(fei-fii));
            es = 0;
            eesc = 0;
            for p=1:L*N+D+LV*NVW+LH*NHW
                des = distance([cylCoord(p, 1), cylCoord(p, 2)], [x(xp), y(yp)]);
                fes = angle2P([cylCoord(p, 1), cylCoord(p, 2)], [x(xp), y(yp)]);
                for m = -M:M
                    es = es + A((p-1)*(2*M+1) + m+M+1)*besselh(m, 2, kE*des)*exp(1i*m*fes);
                end;
            end;
            E(yp, xp) = es+ei;
            E(size(y,2)-yp, xp) = es+ei;
        else
            curCyl = cylContainsXY(cylCoord, x(xp), y(yp));
            if(cylCoord(curCyl, 3) == aD)
                des = distance([cylCoord(curCyl, 1), cylCoord(curCyl, 2)], [x(xp), y(yp)]);
                fes = angle2P([cylCoord(curCyl, 1), cylCoord(curCyl, 2)], [x(xp), y(yp)]);
                eint = 0;
                for m = -M:M
                  %  binJn = F((curCyl-1)*(2*M+1)+m+M+1)*besselj(m, kE*des)+A((curCyl-1)*(2*M+1)+m+M+1)*besselh(m, 2, kE*des);
                    b = T((curCyl-1)*(2*M+1)+m+M+1, 1)*A((curCyl-1)*(2*M+1)+m+M+1, 1);
                    eint = eint + b*besselj(m, cylCoord(curCyl, 4)*des)*exp(1i*m*fes);
                end;
                E(yp, xp) = eint;
                E(size(y,2)-yp, xp) = eint;
            else
                E(yp, xp) = 0;
                E(size(y,2)-yp, xp) = 0;
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








    
    
    