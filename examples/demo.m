%% Demo: comparing forecasts made at different frequencies
% Zhang, Zu and Li, Journal of Business and Economic Statistics.
%
% Run from the repository root, or add it to the path first:
%   addpath(genpath('path/to/dm-different-frequencies'));

rng(1);

% Simulate two forecast-loss series on a common monthly grid, then observe
% forecast 1 monthly (delta_u = 1) and forecast 2 quarterly (delta_v = 3).
n = 600; rho = 0.3; phi = 0.5;
[e1, e2] = barma11(phi, 0, phi, 0, rho, n);   % correlated ARMA(1,1) errors

show = @(r) fprintf(['  two-sample  = %6.3f  (p = %.3f)\n' ...
                     '  clustered-t = %6.3f  (p = %.3f)\n' ...
                     '  aligned DM  = %6.3f  (p = %.3f)\n'], ...
                     r.dm2s, r.pval_dm2s, r.dmcl, r.pval_dmcl, r.dm, r.pval_dm);

%% (a) Equal expected loss (the null)
u = e1.^2;  u = u(1:1:end);        % monthly losses of forecast 1
v = e2.^2;  v = v(1:3:end);        % quarterly losses of forecast 2
fprintf('\n--- Equal expected loss ---\n');
show(dm_difffreq(u, v, 1, 3));

%% (b) Forecast 1 worse by a constant margin (the alternative)
u2 = e1.^2 + 0.3;  u2 = u2(1:1:end);
fprintf('\n--- Forecast 1 has larger expected loss ---\n');
show(dm_difffreq(u2, v, 1, 3));

%% (c) One-sided p-values (H1: forecast 1 is better), as in the empirical section
fprintf('\n--- One-sided (left tail: forecast 1 better) ---\n');
show(dm_difffreq(u2, v, 1, 3, 'tail', 'left'));
