%init
format long;
dielE = 1;
dielDK =100000;
aP = 0.65e-3;
dx = 1e-6;
aD = 0.18*dx;
dy = dx;
fii = pi;
M=3;
LL=1;
NK = 5*2;
LH = 1;
NH = 10*2;
aditation = 0;
for ll = 1:LL
    aditation = aditation+1+aditation;
end;
cylCoord = ones(LL*NK+LH*NH, 5);
f = 10; %input('enter frequency: ');
freq = physconst('LightSpeed')*0.33/dx;
kE = 2*pi*freq*sqrt(dielE)/physconst('LightSpeed');
kDK = 2*pi*freq*sqrt(dielDK)/physconst('LightSpeed');
xs = 0; ystart = 0;
%creating coordinate matrix.

for ll = 1:LL
    cvy=aD+ (ll-1)*aD;
    cvxi = 2*dx + 4*(ll-1)*aD;
    for p = 1 + (ll-1)*NK:2:NK-1 + (ll-1)*NK
        cylCoord(p, :) = [cvxi, cvy, aD, kDK, dielDK];
        cylCoord(p+1, :) = [cvxi, -cvy, aD, kDK, dielDK];
    %     cylCoord(p+2, :) = [cvxf, cvy, aP];
    %     cylCoord(p+3, :) = [cvxf, -cvy, aP];
        cvy = cvy+2*aD;
    end;
end;

for lh = 1:LH
    cxh = cvxi + 2*aD;
    cyh = (NK/2-1)*2*aD-(lh-1)*aD*2;
    for p = LL*NK + (lh-1)*NH + 1:2:LL*NK + (lh-1)*NH + NH - 1
        cylCoord(p, :) = [cxh, cyh, aD, kDK, dielDK];
        cylCoord(p+1, :) = [cxh, -cyh, aD, kDK, dielDK];
        cxh = cxh+2*aD;
    end;
end;

% cylCoord(LL*NK+LH*NH+1, :) = [cvxi, aD, aD, kDK, dielDK];
% cylCoord(LL*NK+LH*NH+2, :) = [cvxi, -aD, aD, kDK, dielDK];
% cylCoord(LL*NK+LH*NH+3, :) = [cvxi-4*aD, 0, aD, kDK, dielDK];

%plot(cylCoord(:, 1), cylCoord(:, 2), 'o');


C = zeros((LL*NK+LH*NH)*(2*M+1), (LL*NK+LH*NH)*(2*M+1));
F = zeros((LL*NK+LH*NH)*(2*M+1), 1);
T = zeros((LL*NK+LH*NH)*(2*M+1), 1);
for p=1:LL*NK+LH*NH
    df = distance([xs,ystart], [cylCoord(p, 1), cylCoord(p,2)]);
    ff = angle2P([xs, ystart], [cylCoord(p, 1), cylCoord(p,2)]);
    for m = -M:M
        fc = (1i)^(m)*exp(-1i*m*fii)*exp(1i*kE*df*cos(ff-fii));
        F((p-1)*(2*M+1) + m+M+1, 1) = fc;
    end;
    for q=1:LL*NK+LH*NH
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
        

        
        
        
%incident field and seccondary field;
ystart = min(cylCoord(:, 2))-aD;
yfinal = max(cylCoord(:, 2))+aD;
rangey = abs(ystart)+abs(yfinal);
rangex = rangey*2;
xstart = min(cylCoord(:, 1))-2*aD;
x = linspace(xstart, xstart+rangex, 70);
y = linspace(ystart, yfinal, 50);
E = zeros(size(y, 2), size(x, 2));
Ei = zeros(size(y, 2), size(x, 2));
Es = zeros(size(y, 2), size(x, 2));
Eint = zeros(size(y, 2), size(x, 2));
figure;
for yp = 1:size(y, 2)/2
    for xp = 1:size(x, 2)
        if cylContainsXY(cylCoord, x(xp), y(yp))==-1% && (y(yp)<wg/2 && y(yp)>-wg/2)
            dei = distance([xs, ystart], [x(xp), y(yp)]);
            fei = angle2P([xs, ystart], [x(xp), y(yp)]);
            ei = exp(1i*kE*dei*cos(fei-fii));
            es = 0;
            eesc = 0;
            for p=1:LL*NK+LH*NH
                des = distance([cylCoord(p, 1), cylCoord(p, 2)], [x(xp), y(yp)]);
                fes = angle2P([cylCoord(p, 1), cylCoord(p, 2)], [x(xp), y(yp)]);
                for m = -M:M
                    es = es + A((p-1)*(2*M+1) + m+M+1)*besselh(m, 2, kE*des)*exp(1i*m*fes);
                end;
            end;
            E(yp, xp) = es+ei;
            E(size(y,2)-yp+1, xp) = es+ei;
            Ei(yp, xp) = ei;
            Es(yp, xp) = es;
        else
            curCyl = cylContainsXY(cylCoord, x(xp), y(yp));
            if(cylCoord(curCyl, 3) == aD)
                des = distance([cylCoord(curCyl, 1), cylCoord(curCyl, 2)], [x(xp), y(yp)]);
                fes = angle2P([cylCoord(curCyl, 1), cylCoord(curCyl, 2)], [x(xp), y(yp)]);
                eint = 0;
                for m = -M:M
                    binJn = F((curCyl-1)*(2*M+1)+m+M+1)*besselj(m, kE*des)+A((curCyl-1)*(2*M+1)+m+M+1)*besselh(m, 2, kE*des);
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
imagesc(x, y, abs(real(E)));
colorbar;
hold on;
axis([xstart xstart+rangex ystart yfinal]);
viscircles([cylCoord(:,1), cylCoord(:,2)], cylCoord(:, 3), 'LineStyle', ':', 'LineWidth', 1/70)
% figure;
% imagesc(x, y, (real(Es)));
% colorbar;
% figure;
% imagesc(x, y, (real(Ei)));
% colorbar;
% figure;
% imagesc(x, y, (real(Eesc)));
% colorbar;