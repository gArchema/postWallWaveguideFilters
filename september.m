%initialization
format long;
dc = 10; %dielectric constant
xs=0; ys=0; %sourcce coordinates
a = 0.632e-3; %radius
dx = 2.54e-3; %period x axis
dy = dx; %period y axis
wg = 8.82e-3; %waveguide port length
M = 3; %symsum bounds
f = input('enter freq: ');
freq = f*1e9; %frequency
lambda =  physconst('LightSpeed')/(freq*sqrt(dc)); %wavelength in meters
k1 = 2*pi/lambda; %wavenumber in meters too
fii = 0; %incident angle
N = 6; %amount of cylinders
K = input('enter amount of defects! ');

%make coordinate matrix for cylinders

cylCoord = zeros(N+K, 3);
cx = 0; cy = wg/2;
for p = 0:N/2-1
    q = 2*p+1; %for coordinate matrix indexes
%     viscircles([cx, cy], a, 'Color', 'k', 'LineStyle', '-');
%     viscircles([cx, -cy], a, 'Color', 'k', 'LineStyle', '-');
    cylCoord(q+1, :) = [cx, -cy, a];
    cylCoord(q, :) = [cx, cy, a];
    cx = cx + dx;
end;

for p = N+1:N+K
    r = input('enter r: '); r = r*1e-3;
    xd = input('enter x: '); xd = xd*1e-3;
    yd = input('enter y: '); yd = yd*1e-3;
    cylCoord(p, :) = [xd, yd, r];
end;
    
%axis([0 cylCoord(N, 1) cylCoord(N, 2)-5*a cylCoord(N-1, 2)+5*a]);


%secondary field
for posr = 1:N+K
    %phase matrix
    d = distance(cylCoord(posr, 1), cylCoord(posr, 2), xs, ys); %distance from source to cylinder
    theta = angle2P([xs, ys], [cylCoord(posr, 1), cylCoord(posr, 2)]);
    Fi = zeros(2*M+1,1);
    for mm = -M:M
        Ficur = (-1i)^mm*exp(-1i*mm*fii)*exp(1i*k1*d*cos(theta-fii));
        Fi(mm+M+1,1) = Ficur;
    end;
    FiRes{posr, 1} = Fi;
    %C matrix
    for posc = 1:N+K
        ccurnm = zeros(2*M+1);
        for mm = -M:M
            for nn  = -M:M
                ccur = 0;
                if posr~=posc
                    d = distance(cylCoord(posr, 1), cylCoord(posr, 2), cylCoord(posc, 1), cylCoord(posc, 2)); %distance between p and q element
                    theta = angle2P([cylCoord(posr, 1), cylCoord(posr, 2)], [cylCoord(posc, 1), cylCoord(posc, 2)]);
                    ccur = besselh(mm-nn, 2, k1*d)*exp(-1i*(mm-nn)*theta);
                elseif mm == nn
                    ccur = besselh(nn, 2, k1*cylCoord(posr, 3))/besselj(mm, k1*cylCoord(posr, 3));
                else
                    ccur = 0;
                end;
                ccurnm(mm+M+1, nn+M+1) = ccur; 
            end;
        end;
        C{posr, posc} = ccurnm; 
    end;
end;
A = cell2mat(C)\cell2mat(FiRes);




%incident field and seccondary field;
syms m n real;
x = linspace(0, cylCoord(N, 1), 100); y = linspace(cylCoord(N, 2)-5*a, cylCoord(N-1, 2)+5*a, 100);  
E = ones(size(y, 2), size(x, 2));
Eii = ones(size(y, 2), size(x, 2));
Ess = ones(size(y, 2), size(x, 2));
for xpos = 1:size(x, 2)
    for ypos = 1:size(y, 2)
        if cylContainsXY(cylCoord, x(xpos), y(ypos)) == 0
            %incident field for x y coordinate;
            d = distance(x(xpos), y(ypos), xs, ys); %distance from incident field to x y
            theta = angle2P([xs, ys], [x(xpos), y(ypos)]);
            Ei = 0;
            for mm = -M:M
                Ei = Ei + 1i^(-mm)*besselj(mm, k1*d)*exp(-1i*mm*(fii-theta));
            end;
            
            %scattered field for x y coordinate
            %Nth cilynder
            Es = 0;
            for posr = 1:N+K
                d = distance(cylCoord(posr, 1), cylCoord(posr, 2), x(xpos), y(ypos)); %distance from nth cylinder to x y
                theta = angle2P([cylCoord(posr, 1), cylCoord(posr, 2)], [x(xpos), y(ypos)]);               
                for mm = -M:M
                    Acur = A((posr-1)*(2*M+1) + mm+M+1);                
                    Es = Es + Acur*1i^(-mm)*besselh(mm, 2, k1*d)*exp(1i*mm*theta);
                end;
            end;
            E(ypos, xpos) = Ei + Es;
            Eii(ypos, xpos) = Ei;
            Ess(ypos, xpos) = Es;
        else
            E(ypos, xpos) = 0;
            Eii(ypos, xpos) = 0;
            Ess(ypos, xpos) = 0;
        end;
    end;
end;
imagesc(x, y, abs(real(E)));
colorbar;

