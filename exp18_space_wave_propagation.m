%% Exp 18 : Space wave propagation - Friis free space vs two-ray ground reflection
clc; clear; close all;

c = 3e8; f = 900e6; lam = c/f;                 % 900 MHz
h1 = 30; h2 = 1.5; Pt = 1; Gt = 1; Gr = 1;     % antenna heights (m), power (W), gains
d  = logspace(1.5,4.5,600);                    % distance (m)

model = @(d,h1,h2) deal( Pt*Gt*Gr*(lam./(4*pi*d)).^2, ...
    Pt*Gt*Gr*(lam./(4*pi*sqrt(d.^2+(h1-h2)^2))).^2 .* ...
    abs(1 - exp(-1j*2*pi*(sqrt(d.^2+(h1+h2)^2)-sqrt(d.^2+(h1-h2)^2))/lam)).^2 );   % Gamma = -1
[Pfs,Ptr] = model(d,h1,h2);

figure; semilogx(d,10*log10(Pfs/Pt),d,10*log10(Ptr/Pt),'LineWidth',1.5); grid on;
xlabel('Distance (m)'); ylabel('Path gain (dB)'); legend('Free space','Two-ray'); title('Space-wave path gain, 900 MHz');
d_b = 4*h1*h2/lam;  xline(d_b,'--k','d_b');

%% Breakpoint and radio horizons
dh_1  = 3.57*(sqrt(h1)+sqrt(h2));              % k = 1   (km)
dh_43 = 4.12*(sqrt(h1)+sqrt(h2));              % k = 4/3 (km)
fprintf('Wavelength = %.4f m\nBreakpoint d_b = 4*h1*h2/lambda = %.1f m\n', lam, d_b);
fprintf('Radio horizon k=1   : %.2f km\nRadio horizon k=4/3 : %.2f km\n', dh_1, dh_43);

%% Table at specified distances
ds = [100 500 1000 2000 5000];
[Pf,Pt2] = model(ds,h1,h2);
fprintf('\n d(m)    FSPL(dB)   Two-ray PL(dB)   Difference(dB)\n');
for i = 1:numel(ds)
    fprintf('%6d %10.2f %14.2f %14.2f\n', ds(i), -10*log10(Pf(i)/Pt), -10*log10(Pt2(i)/Pt), 10*log10(Pf(i)/Pt2(i)));
end

%% Effect of antenna heights (two-ray curves)
figure; hold on; grid on; set(gca,'XScale','log');
hs = [10 1.5; 30 1.5; 30 5; 60 5];
for i = 1:size(hs,1)
    [~,Pr] = model(d,hs(i,1),hs(i,2));
    plot(d,10*log10(Pr/Pt),'LineWidth',1.2,'DisplayName',sprintf('h_1=%g m, h_2=%g m',hs(i,1),hs(i,2)));
end
[Pf0,~] = model(d,h1,h2); plot(d,10*log10(Pf0/Pt),'k--','DisplayName','Free space');
xlabel('Distance (m)'); ylabel('Path gain (dB)'); legend('Location','southwest'); title('Effect of antenna heights');

%% Slope beyond breakpoint (should be ~ 40 dB/decade)
[~,Pr] = model(d,h1,h2); idx = d > 3*d_b;
pf = polyfit(log10(d(idx)),10*log10(Pr(idx)/Pt),1);
fprintf('\nTwo-ray mean slope beyond d_b ~ %.1f dB/decade (free space = -20 dB/decade)\n', pf(1));
