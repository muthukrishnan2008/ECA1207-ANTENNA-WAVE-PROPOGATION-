%% Exp 16 : Phased array antenna (4 elements, d = lambda/2) - MATLAB simulation of the trainer-kit experiment
% Aim  : patterns for beta = 0,30,60,90 (and -60) deg, steering angle, HPBW broadening, scan loss, SLL, grating lobes
clc; clear; close all;

N    = 4;                       % number of elements
d_l  = 0.5;                     % spacing in wavelengths
betas = [0 30 60 90 -60];       % progressive phase shift (deg)
th   = linspace(-90,90,3601);   % angle from broadside (deg)
kd   = 2*pi*d_l;                % k*d (rad)  -> 180 deg for d = lambda/2

AF   = @(th,b) abs( sin(N*(kd*sind(th)-deg2rad(b))/2) ./ (N*sin((kd*sind(th)-deg2rad(b))/2)) );
AF   = @(th,b) fixnan(AF(th,b));                   % remove 0/0 at the peak
elem = @(th) sqrt(max(cosd(th),0));                % element field pattern (power ~ cos(theta)), set to 1 for isotropic

res = zeros(numel(betas),6);
figure('Name','Rectangular patterns'); hold on; grid on;
for k = 1:numel(betas)
    E  = AF(th,betas(k)).*elem(th);
    P  = 20*log10(E+1e-12);                        % dB re array factor maximum
    [pk,ip] = max(P);   th0 = th(ip);
    % HPBW: -3 dB crossings around the peak
    l = find(P(1:ip) <= pk-3, 1, 'last');  r = ip-1+find(P(ip:end) <= pk-3, 1, 'first');
    if isempty(l), l = 1; end, if isempty(r), r = numel(th); end
    hp = th(r) - th(l);
    % side-lobe level: highest local maximum outside the main lobe
    loc = find(P(2:end-1)>P(1:end-2) & P(2:end-1)>=P(3:end))+1;
    loc = loc(abs(th(loc)-th0) > hp);   sll = max([P(loc) -Inf]) - pk;
    th_pred = asind(deg2rad(betas(k))/kd);         % theta0 = asin(beta/kd)
    res(k,:) = [betas(k) th_pred th0 hp sll pk];
    plot(th,P,'LineWidth',1.3,'DisplayName',sprintf('\\beta = %d^\\circ',betas(k)));
end
ylim([-40 2]); xlabel('\theta from broadside (deg)'); ylabel('Relative power (dB)');
title('Phased array N = 4, d = \lambda/2'); legend('Location','south');

%% Polar plot (beta = 60 deg) over 0-360 deg  (as in the trainer-kit table)
ang = (0:1:360); a_ = ang; a_(a_>180) = a_(a_>180)-360;     % 0..360 -> -180..180
th_b = a_; E = AF(th_b,60).*(abs(th_b)<=90).*elem(th_b) + 0;
figure('Name','Polar'); polarplot(deg2rad(ang),max(20*log10(E+1e-6),-40)+40,'LineWidth',1.5);
title('\beta = 60^\circ polar pattern (dB + 40)');

%% Tables
fprintf('\n beta   pred.th0  meas.th0   HPBW    SLL(dB)  peak(dB)  HPBW(th0) pred  scan loss pred\n');
hp0 = res(1,4);
for k = 1:numel(betas)
    fprintf('%5d %9.2f %9.2f %8.2f %8.2f %9.2f %12.2f %12.2f\n', res(k,1), res(k,2), res(k,3), res(k,4), ...
            res(k,5), res(k,6), hp0/cosd(res(k,2)), 10*log10(cosd(res(k,2))));
end
fprintf('\nBeam direction error is within 3 deg if |pred - meas| <= 3.  SLL must be <= -10 dB.\n');
th_max = asind(90/180);                      % largest scan angle (beta = 90 deg)
fprintf('Grating-lobe-free spacing for beta = 90 deg: d_max = lambda/(1+|sin th0|) = %.3f lambda (used %.2f lambda) -> %s\n', ...
        1/(1+abs(sind(th_max))), d_l, ternary(d_l <= 1/(1+abs(sind(th_max))),'OK','grating lobe!'));

function y = fixnan(y), y(isnan(y)) = 1; end
function s = ternary(c,a,b), if c, s = a; else, s = b; end, end
