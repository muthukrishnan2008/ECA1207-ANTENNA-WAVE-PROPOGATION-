%% Exp 20 : Ionospheric virtual height and tropospheric refractivity / duct analysis
clc; clear; close all;

%% Part A : virtual height h'(f) of a parabolic layer (ionogram)
fc = 9; hb = 200; ym = 100;                   % MHz, km
f  = linspace(0.5,8.95,300);
hv = hb + (ym/2)*(f/fc).*log((1+f/fc)./(1-f/fc));
figure; plot(f,hv,'LineWidth',1.5); grid on; xlabel('f (MHz)'); ylabel("h' (km)"); title('Ionogram trace (parabolic layer)');
dt = 2e-3; hv_pulse = 3e8*dt/2/1e3;           % h' = c*dt/2 (km)
fprintf("Echo delay %.1f ms -> virtual height h' = %.0f km\n", dt*1e3, hv_pulse);
hf = @(x) hb + (ym/2)*(x/fc).*log((1+x/fc)./(1-x/fc));
exp_h = [221.2 253.6 325.9]; fq = [4 6 8];
fprintf('\n f(MHz)  computed h''(km)  expected  deviation(%%)\n');
for i = 1:3, fprintf('%5g %14.1f %10.1f %10.2f\n', fq(i), hf(fq(i)), exp_h(i), 100*(hf(fq(i))-exp_h(i))/exp_h(i)); end

%% Part B : tropospheric refractivity and M profile
T = 288; P = 1013; e = 10;                    % K, hPa, hPa
N0 = 77.6/T*(P + 4810*e/T);
fprintf('\nSurface refractivity N0 = %.1f N-units\n', N0);
z = 0:5:400;                                  % height (m)
grad = -40*ones(size(z)); grad(z>=100 & z<150) = -300;   % dN/dh (N-units/km)
N = N0 + cumtrapz(z/1e3, grad);               % refractivity profile
M = N + 0.157*z;                              % modified refractivity
trap = grad < -157;                           % trapping layer (dM/dh < 0)
a = 6371; k = 1./(1 + a*grad*1e-6);           % effective earth-radius factor
d_hor = 3.57*sqrt(k(1))*(sqrt(30)+sqrt(1.5)); % radio horizon (km), h1 = 30 m, h2 = 1.5 m
fprintf('k (lower layer, -40 N/km) = %.3f ; radio horizon = %.2f km\n', k(1), d_hor);
fprintf('Trapping layer (dM/dh < 0) from %g m to %g m\n', min(z(trap)), max(z(trap))+5);

figure; plot(M,z,'LineWidth',1.5); hold on; plot(M(trap),z(trap),'r.','MarkerSize',12); grid on;
xlabel('Modified refractivity M'); ylabel('Height z (m)'); title('M-profile with duct (red)'); legend('M profile','Trapping layer');

%% Classification table
gs = [20 -40 -100 -300];                      % example gradients (N-units/km)
fprintf('\n dN/dh      condition               k        dM/dh (M/km)\n');
for g = gs
    if g > 0, c = 'sub-refraction'; elseif g > -79, c = 'normal'; elseif g > -157, c = 'super-refraction'; else, c = 'ducting (trapping)'; end
    kk = 1/(1 + a*g*1e-6); fprintf('%5d   %-20s %8.3f %12d\n', g, c, kk, g+157);
end
fprintf('(k < 0 means the ray curvature exceeds the earth''s: the wave is trapped in the duct.)\n');
