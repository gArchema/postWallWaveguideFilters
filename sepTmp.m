% %initialization
% dc = 10; %dielectric constant
% xs=0; ys=0; %sourcce coordinates
% a = 0.632e-3; %radius
% dx = 2.54e-3; %period x axis
% dy = dx; %period y axis
% wg = 8.82e-3; %waveguide port length
% M = 3; %symsum bounds
% L=1;
% f = input('enter freq: ');
% freq = f*1e9; %frequency
% lambda =  physconst('LightSpeed')/(freq*sqrt(dc)); %wavelength in meters
% k1 = 2*pi/lambda; %wavenumber in meters too
% fii = 0; %incident angle
% x = linspace(xs, 20e-3, 100);
% y = linspace(ys-wg/2-6*a-(L-1)*2*a, ys+wg/2+6*a+(L-1)*2*a, 100);  
% Eii = zeros(size(y, 2), size(x, 2));
% 
% for xpos = 1:size(x, 2)
%     for ypos = 1:size(y, 2)
%         d = distance(x(xpos), y(ypos), xs, ys); %distance from incident field to x y
%         theta = angle2P([xs, ys], [x(xpos), y(ypos)]);
%         Ei = 0;
%         for mm = -M:M
%             Ei = Ei + 1i^(-mm)*besselj(mm, k1*d)*exp(-1i*mm*(fii-theta));
%         end;
%         Eii(ypos, xpos) = Ei;
%     end;
% end;
% imagesc(x, y, abs(real(Eii)));
% colorbar;

fileID = fopen('defects.txt','r');
formatSpec = '%d %f';
sizeDef = [3, 13];
Defects = fscanf(fileID,formatSpec,sizeDef);
D = Defects(1,1);
fclose(fileID);
display(Defects);