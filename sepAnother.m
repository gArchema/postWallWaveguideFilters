%init
xs = 0; ys = 0;
format short;
dielC = 10;
dielCin = 1;
a = 1.5e-3;
dx = 3.3e-3;
dy=4.04e-3;
wg = 8.82e-3;
fii = 0;
M=5;
L=3;
N = 6;
D = input('enter amount of defects: ');
cylCoord = ones(N+D, 3);
f = 10; %input('enter frequency: ');
freq = f*1e9;
lambda1 =  physconst('LightSpeed')/(freq*sqrt(dielC));
lambda2 =  physconst('LightSpeed')/(freq*sqrt(dielCin));
k1 = 2*pi/lambda1;
k2 = 2*pi/lambda2;
nc = sqrt(dielCin/dielC); %an piriqit
eta = nc^2;
%creating coordinate matrix.
for l = 1:L
    cx = 0;
    cy = (l-1)*dy+wg/2;
    for p = 1:N/2
        cylCoord(p, 1) = cx;
        cylCoord(p, 2) = cy;
        cylCoord(p, 3) = a;
        cylCoord(p+N/2, 1) = cx;
        cylCoord(p+N/2, 2) = -cy;
        cylCoord(p+N/2, 3) = a;
        cx = cx+dx;
    end;
end;

for p = N+1:N+D
    r = input('enter r: ');
    xd = input('enter x: ');
    yd = input('enter y: ');
    cylCoord(p, :) = [xd*1e-3, yd*1e-3, r*1e-3];
end;

syms x real;
%creating matrixes
for p = 1:N+D
    t = zeros(2*M+1, 2*M+1);
    aa = zeros(2*M+1, 1);
    ds = distance([xs, ys], [cylCoord(p, 1), cylCoord(p, 2)]);
    fis = angle2P([xs, ys], [cylCoord(p, 1), cylCoord(p, 2)]);
    for m = -M:M
        aa(m+M+1, 1) = (-1i)^m*exp(-1i*k1*ds*cos(fii-fis))*exp(-1i*m*fii);
        for n = -M:M
            if m==n
                besj = diff(besselj(m, x), x);
                besj1 = subs(besj, x, k1*a);
                besj2 = subs(besj, x, k2*a);
                besh = diff(besselj(m, x) - 1i*bessely(m,  x), x);
                besh1 = subs(besh, x, k1*a);
                besh2 = subs(besh, x, k2*a);
                t(m+M+1, m+M+1) = (eta*besselj(m, k1*a)*besj2-nc*besj1*besselj(m, k2*a))/(nc*besh1*besselh(m, 2, k2*a)-eta*besselh(m, 2, k1*a)*besh2);
            end;
        end;
    end;
    T{p,1} = t;
    AA{p, 1} = aa;
    alfa = zeros(2*M+1, 2*M+1);
    for q = 1:N+D
        dbit = distance([cylCoord(p, 1), cylCoord(p, 2)], [cylCoord(q, 1), cylCoord(q, 2)]);
        fibit = angle2P([cylCoord(p, 1), cylCoord(p, 2)], [cylCoord(q, 1), cylCoord(q, 2)]);
        for m = -M:M
            for n = -M:M
                alfa(m+M+1, n+M+1) = besselh(m-n, 2, k1*dbit)*exp(-1i*(m-n)*fibit);
            end;
        end;
        Alfa{p, q} = alfa;
    end;
end;

for p = 1:N+D
    Fi{p, 1} = cell2mat(T(p))*cell2mat(AA(p));
    for q = 1:N+D
        if p == q
            C{p, q} = zeros(2*M+1, 2*M+1);
        else
            C{p, q} = -1*cell2mat(T(p))*cell2mat(Alfa(p, q));
        end;
    end;
end;


A = cell2mat(C)\cell2mat(Fi);


        
%incident field and seccondary field;
x = linspace(xs, cylCoord(N, 1), 75);
y = linspace(ys-wg/2-6*a-(L-1)*2*a, ys+wg/2+6*a+(L-1)*2*a, 75);
E = zeros(size(y, 2), size(x, 2));
Ei = zeros(size(y, 2), size(x, 2));
Es = zeros(size(y, 2), size(x, 2));
for yp = 1:size(y, 2)
    for xp = 1:size(x, 2)
        if ~cylContainsXY(cylCoord, x(xp), y(yp));
            ei = 0;
            dei = distance([xs, ys], [x(xp), y(yp)]);
            fei = angle2P([xs, ys], [x(xp), y(yp)]);
            for m = -M:M
                ei = ei + (1i)^(-m)*besselj(m, k1*dei)*exp(1i*m*(fei-fii));
            end;
            es = 0;
            for p=1:N+D
                des = distance([cylCoord(p, 1), cylCoord(p, 2)], [x(xp), y(yp)]);
                fes = angle2P([cylCoord(p, 1), cylCoord(p, 2)], [x(xp), y(yp)]);
                for m = -M:M
                    es = es + A((p-1)*(2*M+1) + m+M+1)*(1i)^(-m)*besselh(m, 2, k1*des)*exp(1i*m*fes);
                end;
            end;
            E(yp, xp) = es+ei;
            Ei(yp, xp) = ei;
            Es(yp, xp) = es;
        else
            E(yp, xp) = 0;
            Ei(yp, xp) = 0;
            Es(yp, xp) = 0;
        end;
    end;
end;
imagesc(x, y, abs(real(E)));
colorbar;











                
                
                
                
    
    