%% Exp 17 : Ground wave and sky wave propagation (Norton attenuation, crossover distance)
clc; clear; close all;

c = 3e8; f = 1e6; lam = c/f; P = 1;            % 1 MHz, P = 1 kW
d = logspace(0,3,2000);                        % distance (km)
sig = [0.01 5];                                % land, sea conductivity (S/m)
hv = 100; rho = 0.5;                           % virtual height (km), ionospheric reflection factor

Eg_fun = @(d,s) 300*sqrt(P)./d .* ((2+0.3*(pi*(d*1e3)/(lam*(18000*s/(f/1e6))))) ./ ...
         (2 + (pi*(d*1e3)/(lam*(18000*s/(f/1e6)))) + 0.6*(pi*(d*1e3)/(lam*(18000*s/(f/1e6)))).^2));
Eg = [Eg_fun(d,sig(1)); Eg_fun(d,sig(2))];     % ground-wave field (mV/m)
slant = @(d) 2*sqrt((d/2).^2 + hv^2);          % sky-wave path (km)
Es = 300*sqrt(P)./slant(d) * rho;              % sky-wave field (mV/m)

figure; loglog(d,Eg(1,:),d,Eg(2,:),d,Es,'--','LineWidth',1.5); grid on;
xlabel('Distance d (km)'); ylabel('Field strength (mV/m)');
legend('Ground wave (land)','Ground wave (sea)','Sky wave'); title('Ground wave and sky wave, f = 1 MHz, P = 1 kW');

%% Crossover distance (where E_g = E_s)
for k = 1:2
    [~,i] = min(abs(log(Eg(k,:)) - log(Es)));
    nm = {'land','sea'};
    fprintf('Crossover distance (%s) = %.1f km\n', nm{k}, d(i));
end

%% Spot values  (d = 50, 100, 200 km)
ds = [50 100 200];
fprintf('\n d(km)   p(land)     A(p)    Eg land   Eg sea    slant s(km)   Es (mV/m)\n');
for dd = ds
    x = 18000*sig(1)/(f/1e6); p = pi*(dd*1e3)/(lam*x); A = (2+0.3*p)/(2+p+0.6*p^2);
    fprintf('%5d %9.3f %9.4f %9.3f %9.3f %12.2f %11.4f\n', dd, p, A, Eg_fun(dd,sig(1)), Eg_fun(dd,sig(2)), slant(dd), 300*sqrt(P)/slant(dd)*rho);
end
fprintf('\nFading zone: where Eg and Es are comparable (near the crossover) the two waves interfere\n');
fprintf('and, because the sky wave varies with the ionosphere (night-time), the signal fades.\n');
