% Substrate Inhibition Bioprocess Model
% Biomass (X), Substrate (S), and Product (P) over time
clear; clc;
% ------------------------
% Parameters
% ------------------------
mu_max = 0.4;    % 1/h - Maximum specific growth rate
Ks     = 0.5;    % g/L - Saturation constant
Ki     = 2.0;    % g/L - Inhibition constant
Yxs    = 0.5;    % g biomass / g substrate
Ypx    = 0.2;    % g product / g biomass
% Initial Conditions
X0 = 0.1;    % g/L - Initial biomass
S0 = 10.0;   % g/L - Initial substrate
P0 = 0.0;    % g/L - Initial product
y0 = [X0, S0, P0];
% Time Span
tspan = linspace(0, 50, 500);  % 0 to 50 hours, 500 points
% ------------------------
% ODE Solver
% ------------------------
odefun = @(t, y) substrateInhibitionODE(t, y, mu_max, Ks, Ki, Yxs, Ypx);
[t, Y] = ode45(odefun, tspan, y0);
% Extract Results
X = Y(:, 1);
S = Y(:, 2);
P = Y(:, 3);
% ------------------------
% Plotting
% ------------------------
figure('Color', 'w');
plot(t, X, 'g-', 'LineWidth', 2); hold on;
plot(t, S, 'b--', 'LineWidth', 2);
plot(t, P, 'r-.', 'LineWidth', 2);

xlabel('Time (hours)', 'FontSize', 12);
ylabel('Concentration (g/L)', 'FontSize', 12);
title('Bioprocess Dynamics with Substrate Inhibition', 'FontSize', 14);
legend({'Biomass (X)', 'Substrate (S)', 'Product (P)'}, ...
       'Location', 'best', 'FontSize', 11);

grid on;
xlim([0 max(t)]);
ylim([0 max([X; S; P]) * 1.1]);
set(gca, 'FontSize', 11);
box on;

% ------------------------
% ODE Function Definition
% ------------------------
function dydt = substrateInhibitionODE(~, y, mu_max, Ks, Ki, Yxs, Ypx)
    X = y(1);  % Biomass
    S = y(2);  % Substrate

    % Substrate inhibition growth rate
    mu = mu_max * S / (Ks + S + (S^2 / Ki));

    % ODEs
    dXdt = mu * X;
    dSdt = - (1 / Yxs) * mu * X;
    dPdt = Ypx * mu * X;

    dydt = [dXdt; dSdt; dPdt];
end
