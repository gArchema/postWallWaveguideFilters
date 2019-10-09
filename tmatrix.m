%init
format long;
dielC = 10;
dielout = 1;
dielCin = 0;
a = 0.632e-3;
dx = 2.54e-3;
wg = 8.82e-3;
dy=wg/6;
fii = pi;
M=5;
L=1;
N = 12;
H = 3*4;
V = 3*4;
cylCoord = ones(L*N, 3);
D = input('enter amount of defects: ');
f = 10; %input('enter frequency: ');
freq = f*1e9;
lambda1 =  physconst('LightSpeed')/(freq*sqrt(dielC));
lambdaOut = physconst('LightSpeed')/(freq*sqrt(dielout));
lambda2 =  physconst('LightSpeed')/(freq*sqrt(dielCin));
%k1 = 2*pi/lambda1;
k1 = 2*pi*freq*sqrt(dielC)/physconst('LightSpeed');
k2 = 2*pi/lambda2;
kout = 2*pi/lambdaOut;
xs = 1*dx; ys = 0;
%creating coordinate matrix.
for l = 1:L
    cx = (H/4+1)*dx;
    cy = (l-1)*dy+wg/2;
    for p = 1:N/2
        cylCoord((l-1)*N + p, 1) = cx;
        cylCoord((l-1)*N + p, 2) = cy;
        cylCoord((l-1)*N + p, 3) = a;
        cylCoord((l-1)*N + (p+N/2), 1) = cx;
        cylCoord((l-1)*N + (p+N/2), 2) = -cy;
        cylCoord((l-1)*N + (p+N/2), 3) = a;
        cx = cx+dx;
    end;
end;

chy = dy;
chx = xs-xs;
for p = L*N+1:4:L*N+H-3
    cylCoord(p, :) = [chx, chy, a];
    cylCoord(p+1, :) = [chx, -chy, a];
    cylCoord(p+2, :) = [chx+(H/4+N/2+2)*dx, chy, a];
    cylCoord(p+3, :) = [chx+(H/4+N/2+2)*dx, -chy, a];
    chx = chx + dx;
end;


cvy = dy;
cvxi = xs+H/4*dx-xs;
cvxf = xs+(H/4+N/2+1)*dx-xs;
for p = L*N+H+1:4:L*N+H+V-3
    cylCoord(p, :) = [cvxi, cvy, a];
    cylCoord(p+1, :) = [cvxi, -cvy, a];
    cylCoord(p+2, :) = [cvxf, cvy, a];
    cylCoord(p+3, :) = [cvxf, -cvy, a];
    cvy = cvy+dy;
end;

%plot(cylCoord(:, 1), cylCoord(:, 2), 'o');

% cy = dy/2;
% cxs = cylCoord(1, 1)-dx/2;
% cxf = cylCoord(L*N, 1)+dx/2;
% for p = L*N+1:4:L*N+D-3
%     cylCoord(p, :) = [cxs, cy, a];
%     cylCoord(p+1, :) = [cxs, -cy, a];
%     cylCoord(p+2, :) = [cxf, cy, a];
%     cylCoord(p+3, :) = [cxf, -cy, a];
%     cy = cy+dy;
% end;
    
for p = L*N+H+V+1:L*N+H+V+D
    r = input('enter r: ');
    xd = input('enter x: ');
    yd = input('enter y: ');
    cylCoord(p, :) = [xd*1e-3, yd*1e-3, r*1e-3];
end;

% fileID = fopen('defects.txt','r');
% formatSpec = '%d %f';
% sizeDef = [13, 3];
% Defects = fscanf(fileID,formatSpec,sizeDef);
% D = Defects(1,1);
% fclose(fileID);
% for p = L*N+1:L*N+D
%     d = p-L*N;
%     cylCoord(p, :) = Defects(d+1, :)*1e-3;
% end;


%find A matrix
C = zeros((L*N+H+V+D)*(2*M+1), (L*N+H+V+D)*(2*M+1));
F = zeros((L*N+H+V+D)*(2*M+1), 1);
for p=1:L*N+H+V+D
    df = distance([xs,ys], [cylCoord(p, 1), cylCoord(p,2)]);
    ff = angle2P([xs, ys], [cylCoord(p, 1), cylCoord(p,2)]);
%     if cylCoord(p, 2)<0
%         ff = ff+pi;
%     end;
%     if p>=L*N+H+1
%         fiii = fii+pi/2;
%     else
%         fiii = fii;
%     end;
    for m = -M:M
        fc = (1i)^(m)*exp(-1i*m*fii)*exp(1i*k1*df*cos(ff-fii));
        F((p-1)*(2*M+1) + m+M+1, 1) = fc;
    end;
    for q=1:L*N+H+V+D
        dc = distance([cylCoord(p, 1), cylCoord(p, 2)], [cylCoord(q, 1), cylCoord(q, 2)]);
        fc = angle2P([cylCoord(p, 1), cylCoord(p, 2)], [cylCoord(q, 1), cylCoord(q, 2)]);
        for m = -M:M
            for n = -M:M
                if p~=q
                    c = -besselh(m-n, 2, k1*dc)*exp(-1*(1i)*(m-n)*fc);
                elseif m==n
%                      syms xx real;
%                      beshd = diff(besselj(m,xx) - 1i*bessely(m, xx));
%                      beshd1 = subs(beshd, xx, k1*cylCoord(p, 3));
%                      besjd = diff(besselj(m, xx));
%                      besjd1 = subs(besjd, xx, k1*cylCoord(p, 3));
%                     besjd2 = subs(besjd, xx, k2*cylCoord(p, 3));
%                     c = (k2*besselh(m, 2, k1*cylCoord(p, 3))*besjd1-k1*beshd*besselj(m, k2*cylCoord(p, 3)))/(k2*besselj(m, k1*cylCoord(p, 3))*besjd2-k1*besjd1*besselj(m, k2*cylCoord(p, 3)));
                     c = -besselh(m, 2, k1*cylCoord(p, 3))/besselj(m, k1*cylCoord(p, 3));
%                     c = -beshd1/besjd1;
                    
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
x = linspace(xs-5*dx, xs+(N/2+H/2+5)*dx, 100);
%x = linspace(-20e-3, 20e-3, 70);

y = linspace(ys-wg/2-6*a-(L-1)*4*a-2*dy, ys+wg/2+6*a+(L-1)*4*a+2*dy, 70);
E = zeros(size(y, 2), size(x, 2));
Ei = zeros(size(y, 2), size(x, 2));
Es = zeros(size(y, 2), size(x, 2));
Eesc = zeros(size(y, 2), size(x, 2));
figure;
k = k1;
for yp = 1:size(y, 2)
    for xp = 1:size(x, 2)
%         if abs(y(yp))>wg/2
%             k = kout;
%         else
%             k = k1;
%         end;

        if ~cylContainsXY(cylCoord, x(xp), y(yp))% && (y(yp)<wg/2 && y(yp)>-wg/2)
           
            dei = distance([xs, ys], [x(xp), y(yp)]);
            fei = angle2P([xs, ys], [x(xp), y(yp)]);
            ei = exp(1i*k*dei*cos(fei-fii));
%             for m = -M:M
%                 ei = ei + (1i)^(-m)*besselj(m, k*dei)*exp(1i*m*(fei-fii));
%             end;
            es = 0;
            eesc = 0;
            for p=1:L*N+H+V+D
                esctmp = 0;
%                 if p>=L*N+H+1
%                     fiii = fii-pi/2;
%                 else
%                     fiii = fii;
%                 end;
                des = distance([cylCoord(p, 1), cylCoord(p, 2)], [x(xp), y(yp)]);
                fes = angle2P([cylCoord(p, 1), cylCoord(p, 2)], [x(xp), y(yp)]);
                for m = -M:M
                    es = es + A((p-1)*(2*M+1) + m+M+1)*besselh(m, 2, k*des)*exp(1i*m*fes);
                   % esctmp = esctmp + (1i)^(m)*besselj(m, k*des)*exp(1i*m*cos(fes-fii));
                end;
                esctmp = esctmp*exp(1i*k*des*cos(fes-fii));
                eesc = eesc+esctmp;
            end;
            E(yp, xp) = es+ei;
            Ei(yp, xp) = ei;
            Es(yp, xp) = es;
            Eesc(yp, xp) = eesc;
        else
            E(yp, xp) = 0;
            Ei(yp, xp) = 0;
            Es(yp, xp) = 0;
            Eesc(yp, xp) = 0;
        end;
    end;
end;
imagesc(x, y, abs(real(E)));
colorbar;
figure;
imagesc(x, y, (real(Es)));
colorbar;
% figure;
% imagesc(x, y, (real(Ei)));
% colorbar;
% figure;
% imagesc(x, y, (real(Eesc)));
% colorbar;







    
    
    