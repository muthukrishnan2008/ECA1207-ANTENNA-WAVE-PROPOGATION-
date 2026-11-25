%% Exp 19 : Critical frequency, MUF and skip distance
clc; clear; close all;

Nmax = 1e12;                                  % peak electron density (m^-3)
fc = 9*sqrt(Nmax)/1e6;                        % critical frequency (MHz)
fprintf('Critical frequency fc = 9*sqrt(Nmax) = %.3f MHz\n', fc);
hm = 300; ym = 100;                           % layer height and semi-thickness (km)

%% Electron density profile (parabolic layer)
h = linspace(hm-ym, hm+ym, 200);
N = Nmax*(1 - ((h-hm)/ym).^2);
figure; plot(N,h,'LineWidth',1.5); grid on; xlabel('Electron density N (m^{-3})'); ylabel('Height (km)');
title('Parabolic ionospheric layer');

%% Refractive index at the layer peak
f_test = [6 9 12]*1e6;
n_peak = sqrt(1 - 81*Nmax./f_test.^2);        % complex/zero => reflection, real => penetration
for i = 1:3
    if f_test(i)/1e6 < fc, s = 'reflected (n imaginary)'; elseif f_test(i)/1e6 == fc, s = 'critical (n = 0)'; else, s = 'penetrates (n real)'; end
    fprintf('f = %2d MHz : n = %s -> %s\n', f_test(i)/1e6, num2str(n_peak(i),'%.3f'), s);
end

%% MUF: flat earth (secant law) and curved earth
hp = 300; Re = 6371;                          % virtual height, earth radius (km)
d = 100:50:4000;                              % ground range (km)
muf_flat = fc*sqrt(1 + (d/(2*hp)).^2);
psi = d/Re;                                   % central angle (rad)
phi0 = atan(Re*sin(psi/2)./((Re+hp) - Re*cos(psi/2)));   % incidence angle, curved earth
muf_curv = fc./cos(phi0);
err = 100*(muf_flat - muf_curv)./muf_curv;    % flat-earth error (%)

figure; subplot(2,1,1); plot(d,muf_flat,'--',d,muf_curv,'LineWidth',1.5); grid on;
xlabel('Ground range d (km)'); ylabel('MUF (MHz)'); legend('Flat earth','Curved earth','Location','northwest');
subplot(2,1,2); plot(d,err,'LineWidth',1.5); grid on; xlabel('Ground range d (km)'); ylabel('Flat-earth error (%)');

%% Skip distance
f0 = 15;                                      % operating frequency (MHz) > fc
d_skip = 2*hp*sqrt((f0/fc)^2 - 1);
fprintf('Skip distance at f0 = %g MHz (flat earth) = %.1f km\n', f0, d_skip);

%% Table + working frequency (0.85 x MUF)
ds = [500 1000 2000 3000];
fprintf('\n d(km)   MUF flat(MHz)  MUF curved(MHz)  Error(%%)  FOT = 0.85*MUF curved\n');
for x = ds
    mf = interp1(d,muf_flat,x); mc = interp1(d,muf_curv,x);
    fprintf('%5d %12.2f %15.2f %9.2f %16.2f\n', x, mf, mc, 100*(mf-mc)/mc, 0.85*mc);
end
