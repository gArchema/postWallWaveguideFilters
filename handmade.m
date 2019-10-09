%init
dielE = 9.8;
dielD = 1;
aP = 0.65e-3;
aD = 1.5e-3;
dx = 3.3e-3;
dxh = 2.45e-3;
dy = 4.09e-3;
wg = 8.91e-3;
dyv=wg/4;
fii = pi;
M=5;
L=3;
N = 22*2;
Nh = 0*floor(N/2*dx/dxh*2);
NK = 30*2;
cylCoord = ones(L*N+NK+Nh, 3);
D = input('enter amount of defects: ');
f = 10; %input('enter frequency: ');
freq = f*1e9;
kE = 2*pi*freq*sqrt(dielE)/physconst('LightSpeed');
kD = 2*pi*freq*sqrt(dielD)/physconst('LightSpeed');
xs = 0*dx; ystart = 0;
%creating coordinate matrix.
for l = 1:L
    cx = xs + 2*aP+aD;
    cy = (l-1)*dy + wg/2;
    for p = (l-1)*N + 1:2:(l-1)*N + N-1
        cylCoord(p, :) = [cx, cy, aD];
        cylCoord(p+1, :) = [cx, -cy, aD];
%         cylCoord(p+2, :) = [cx, cy+(NK/4-1)*dyh, aP];
%         cylCoord(p+3, :) = [cx, -(cy+(NK/4-1)*dyh), aP];
        cx = cx+dx;
    end;
end;

cvy=wg/2;
cvxi = xs;
cvxf = cvxi+(N/2+1)*dx;
for p = L*N+1:2:L*N+NK-1
    cylCoord(p, :) = [cvxi, cvy, aP];
    cylCoord(p+1, :) = [cvxi, -cvy, aP];
%     cylCoord(p+2, :) = [cvxf, cvy, aP];
%     cylCoord(p+3, :) = [cvxf, -cvy, aP];
    cvy = cvy+dyv;
end;

%plot(cylCoord(:, 1), cylCoord(:, 2), 'o');
    
for p = L*N+NK+1:L*N+NK+D
    r = input('enter r: ');
    xd = input('enter x: ');
    yd = input('enter y: ');
    cylCoord(p, :) = [xd*1e-3, yd*1e-3, r*1e-3];
end;

cxh = xs + 2*aP+aD;
cyh = (NK/2+1)*dyv;
for p = L*N+NK+D+1:2:L*N+NK+D+Nh-1
    cylCoord(p, :) = [cxh, cyh, aP];
    cylCoord(p+1, :) = [cxh, -cyh, aP];
    cxh = cxh + dxh;
end;


%find A matrix
C = zeros((L*N+D+NK+Nh)*(2*M+1), (L*N+D+NK+Nh)*(2*M+1));
F = zeros((L*N+D+NK+Nh)*(2*M+1), 1);
T = zeros((L*N+D+NK+Nh)*(2*M+1), 1);
for p=1:L*N+D+NK+Nh
    df = distance([xs,ystart], [cylCoord(p, 1), cylCoord(p,2)]);
    ff = angle2P([xs, ystart], [cylCoord(p, 1), cylCoord(p,2)]);
    for m = -M:M
        fc = (1i)^(m)*exp(-1i*m*fii)*exp(1i*kE*df*cos(ff-fii));
        F((p-1)*(2*M+1) + m+M+1, 1) = fc;
    end;
    for q=1:L*N+D+NK+Nh
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
                         j2 = besselj(m, kD*cylCoord(p, 3));
                         jt = diff(besselj(m, rk));
                         jd1 = subs(jt, rk, kE*cylCoord(p, 3));
                         jd2 = subs(jt, rk, kD*cylCoord(p, 3));
                         c = -(dielD*kE*h1*jd2-dielE*kD*hd1*j2)/(dielD*kE*j1*jd2-dielE*kD*jd1*j2);
                         t = dielE*kD*(h1*jd1-hd1*j1)/(dielE*kD*jd1*j2-dielD*kE*j1*jd2);
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
ystart = min(cylCoord(:, 2))+5.8*wg;
yfinal = max(cylCoord(:, 2))-5.8*wg;
rangey = abs(ystart)+abs(yfinal);
rangex = rangey*2;
xstart = min(cylCoord(:, 1))+dx;
x = linspace(xstart, xstart+rangex, 70);
y = linspace(ystart, yfinal, 50);
E = zeros(size(y, 2), size(x, 2));
% Ei = zeros(size(y, 2), size(x, 2));
% Es = zeros(size(y, 2), size(x, 2));
% Eint = zeros(size(y, 2), size(x, 2));
figure;
for yp = 1:size(y, 2)
    for xp = 1:size(x, 2)
        if cylContainsXY(cylCoord, x(xp), y(yp))==-1% && (y(yp)<wg/2 && y(yp)>-wg/2)
            dei = distance([xs, ystart], [x(xp), y(yp)]);
            fei = angle2P([xs, ystart], [x(xp), y(yp)]);
            ei = exp(1i*kE*dei*cos(fei-fii));
            es = 0;
            eesc = 0;
            for p=1:L*N+D+NK+Nh
                des = distance([cylCoord(p, 1), cylCoord(p, 2)], [x(xp), y(yp)]);
                fes = angle2P([cylCoord(p, 1), cylCoord(p, 2)], [x(xp), y(yp)]);
                for m = -M:M
                    es = es + A((p-1)*(2*M+1) + m+M+1)*besselh(m, 2, kE*des)*exp(1i*m*fes);
                end;
            end;
            E(yp, xp) = es+ei;
%             Ei(yp, xp) = ei;
%             Es(yp, xp) = es;
        else
            curCyl = cylContainsXY(cylCoord, x(xp), y(yp));
            if(cylCoord(curCyl, 3) == aD)
                des = distance([cylCoord(curCyl, 1), cylCoord(curCyl, 2)], [x(xp), y(yp)]);
                fes = angle2P([cylCoord(curCyl, 1), cylCoord(curCyl, 2)], [x(xp), y(yp)]);
                eint = 0;
                for m = -M:M
                    %binJn = F((curCyl-1)*(2*M+1)+m+M+1)*besselj(m, kE*des)+A((curCyl-1)*(2*M+1)+m+M+1)*besselh(m, 2, kE*des);
                    b = T((curCyl-1)*(2*M+1)+m+M+1, 1)*A((curCyl-1)*(2*M+1)+m+M+1, 1);
                    eint = eint + b*besselj(m, kD*des)*exp(1i*m*fes);
                end;
                E(yp, xp) = eint;
            else
                E(yp, xp) = 0;
            end;
        end;
    end;
end;
imagesc(x, y, abs(real(E)));
colorbar;
hold on;
% axis([xstart xstart+rangex ystart yfinal]);
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







    
    
    