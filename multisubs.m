% Multiple Substrate Inhibition Bioprocess Model


clear; clc;

% ------------------------
% Parameters
% ------------------------
mu_max = 0.5;      % 1/h
Ks1 = 0.3;         % g/L (S1 saturation constant)
Ki1 = 3.0;         % g/L (S1 inhibition constant)
Ks2 = 0.4;         % g/L (S2 saturation constant)
Ki2 = 2.0;         % g/L (S2 inhibition constant)

Yxs1 = 0.4;        % g biomass / g substrate S1
Yxs2 = 0.3;        % g biomass / g substrate S2
Ypx  = 0.2;        % g product / g biomass

% Initial Conditions [X, S1, S2, P]
X0 = 0.1; S1_0 = 5.0; S2_0 = 4.0; P0 = 0.0;
y0 = [X0, S1_0, S2_0, P0];

% Time span
tspan = linspace(0, 50, 500);

% ------------------------
% Solve ODE
% ------------------------
odefun = @(t, y) multiSubstrateODE(t, y, mu_max, Ks1, Ki1, Ks2, Ki2, Yxs1, Yxs2, Ypx);
[t, Y] = ode45(odefun, tspan, y0);

% Extract results
X  = Y(:,1);
S1 = Y(:,2);
S2 = Y(:,3);
P  = Y(:,4);

% ------------------------
% Plot
% ------------------------
figure('Color', 'w');
plot(t, X,  'g-',  'LineWidth', 2); hold on;
plot(t, S1, 'b--', 'LineWidth', 2);
plot(t, S2, 'c-.', 'LineWidth', 2);
plot(t, P,  'r:',  'LineWidth', 2);

xlabel('Time (hours)', 'FontSize', 12);
ylabel('Concentration (g/L)', 'FontSize', 12);
title('Bioprocess with Two Substrate Inhibition', 'FontSize', 14);
legend({'Biomass (X)', 'Substrate 1 (S1)', 'Substrate 2 (S2)', 'Product (P)'}, ...
       'Location', 'best', 'FontSize', 11);

grid on;
xlim([0 max(t)]);
ylim([0 max([X; S1; S2; P]) * 1.1]);
set(gca, 'FontSize', 11);
box on;

% ------------------------
% ODE Function
% ------------------------
function dydt = multiSubstrateODE(~, y, mu_max, Ks1, Ki1, Ks2, Ki2, Yxs1, Yxs2, Ypx)
    X  = y(1);
    S1 = y(2);
    S2 = y(3);
    
    % Growth rate with inhibition from both substrates
    mu = mu_max * (S1 / (Ks1 + S1 + (S1^2 / Ki1))) * ...
                   (S2 / (Ks2 + S2 + (S2^2 / Ki2)));

    dXdt  = mu * X;
    dS1dt = -(1 / Yxs1) * mu * X;
    dS2dt = -(1 / Yxs2) * mu * X;
    dPdt  = Ypx * mu * X;

    dydt = [dXdt; dS1dt; dS2dt; dPdt];
end
