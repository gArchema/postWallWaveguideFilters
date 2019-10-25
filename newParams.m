%init
format long;
dielE = 1;
dielD = 11.56;
dielDK = 39;
dx = 1e-6;
dy = dx;
aD = 0.18*dx;
aP = 0.65e-3;
dyh = 2*aD;
wg = 2*dy;
fii = pi;
M=3;
L=2;
N = 20*2;
LL=3;
NK = 20*2;
LH = 1;
NH = floor(dx*N/2/(2*aD)+3)*2;
cylCoord = ones(L*N+LL*NK+LH*NH, 5);
D = input('enter amount of defects: ');
f = 10; %input('enter frequency: ');
freq = physconst('LightSpeed')*0.44/dx;
kE = 2*pi*freq*sqrt(dielE)/physconst('LightSpeed');
kD = 2*pi*freq*sqrt(dielD)/physconst('LightSpeed');
kDK = 2*pi*freq*sqrt(dielDK)/physconst('LightSpeed');
xs = 0; ystart = 0;
%creating coordinate matrix.

for ll = 1:LL
    distFSource = 2*dx;
    cvy= wg/2 + (LL-ll)*aD;
    cvxi = xs + 2*dx + 4*(ll-1)*aD;
 %   cvxf = distFSource + LL*dx + N/2*dx;
    for p = (ll-1)*NK+1:2:(ll-1)*NK+NK-1
        cylCoord(p, :) = [cvxi, cvy, aD, kDK, dielDK];
        cylCoord(p+1, :) = [cvxi, -cvy, aD, kDK, dielDK];
%         cylCoord(p+2, :) = [cvxf, cvy, aD, kDK];
%         cylCoord(p+3, :) = [cvxf, -cvy, aD, kDK];
        cvy = cvy+dyh;
    end;
end;

for l = 1:L
    cx = cvxi + dx;
    cy = (l-1)*dy + wg/2;
    for p = LL*NK + (l-1)*N + 1:2:LL*NK + (l-1)*N + N-1
        cylCoord(p, :) = [cx, cy, aD, kD, dielD];
        cylCoord(p+1, :) = [cx, -cy, aD, kD, dielD];
%         cylCoord(p+2, :) = [cx, cy+(NK/2-1)*dyh, aD, kDK];
%         cylCoord(p+3, :) = [cx, -(cy+(NK/2-1)*dyh), aD, kDK];
        cx = cx+dx;
    end;
end;

for lh = 1:LH
    cxh = cvxi + 2*aD;
    cyh = wg/2 + (NK/2-1)*2*aD-(lh-1)*aD*2;
    for p = LL*NK + L*N + (lh-1)*NH + 1:2:LL*NK + L*N + (lh-1)*NH + NH - 1
        cylCoord(p, :) = [cxh, cyh, aD, kDK, dielDK];
        cylCoord(p+1, :) = [cxh, -cyh, aD, kDK, dielDK];
        cxh = cxh+2*aD;
    end;
end;



%plot(cylCoord(:, 1), cylCoord(:, 2), 'o');
    
for p = L*N+LL*NK+LH*NH+1:L*N+LL*NK+LH*NH+D
    r = input('enter r: ');
    xd = input('enter x: ');
    yd = input('enter y: ');
    cylCoord(p, :) = [xd*1e-3, yd*1e-3, r*1e-3];
end;

% %find A matrix
C = zeros((L*N+D+LL*NK+LH*NH)*(2*M+1), (L*N+D+LL*NK+LH*NH)*(2*M+1));
F = zeros((L*N+D+LL*NK+LH*NH)*(2*M+1), 1);
T = zeros((L*N+D+LL*NK+LH*NH)*(2*M+1), 1);
for p=1:L*N+D+LL*NK+LH*NH
    df = distance([xs,ystart], [cylCoord(p, 1), cylCoord(p,2)]);
    ff = angle2P([xs, ystart], [cylCoord(p, 1), cylCoord(p,2)]);
    for m = -M:M
        fc = (1i)^(m)*exp(-1i*m*fii)*exp(1i*kE*df*cos(ff-fii));
        F((p-1)*(2*M+1) + m+M+1, 1) = fc;
    end;
    for q=1:L*N+D+LL*NK+LH*NH
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
ystart = min(cylCoord(:, 2))+1.5*dx;
yfinal = max(cylCoord(:, 2))-1.5*dx;
rangey = abs(ystart)+abs(yfinal);
rangex = rangey*2;
xstart = min(cylCoord(:, 1))+4*dx;
x = linspace(xstart, xstart+rangex, 600);
y = linspace(ystart, yfinal, 300);
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
            for p=1:L*N+D+LL*NK+LH*NH
                des = distance([cylCoord(p, 1), cylCoord(p, 2)], [x(xp), y(yp)]);
                fes = angle2P([cylCoord(p, 1), cylCoord(p, 2)], [x(xp), y(yp)]);
                for m = -M:M
                    es = es + A((p-1)*(2*M+1) + m+M+1)*besselh(m, 2, kE*des)*exp(1i*m*fes);
                end;
            end;
            E(yp, xp) = es+ei;
            E(size(y,2)-yp, xp) = es+ei;
            Ei(yp, xp) = ei;
            Es(yp, xp) = es;
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







    
    
    