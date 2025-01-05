clc; clear;

% C1 FFL (AND gate) ------------------------------------------------------
% Promoter is activated by Transcription factor (Binary dynamics)
% Promoter initiates Gene Y transcription
% Gene Y crosses threshold Kyz, activates transcription of Gene Z
% Z is transcribed ONLY when X = 1 and Y > threshold

alpha = 0.2;
beta = 1;
tSpan = [0, 15];
Kyz = 0.3;
z0 = 0;
y0 = 0;
activationTime = 2;
deactivationTime = 8;

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

t = 0:0.01:max(tSpan);
x = @(t) beta * (t > activationTime & t < deactivationTime);
x_arr = arrayfun(x, t);

dydt = @(t1, y) (x(t1) - alpha * y);
[t1, y] = ode45(dydt, tSpan, y0, options);

y_interpolated = @(t2) interp1(t1, y, t2, 'linear', 'extrap');

dzdt = @(t2, z) (x(t2) * (y_interpolated(t2) > Kyz * max(y)) - alpha * z);
[t2, z] = ode45(dzdt, tSpan, z0, options);

figure;
hold on;
ylim([0 max(max(y), max(z))]);
plot(t, x_arr, 'black', 'DisplayName', 'Promoter');
plot(t1, y, 'red','DisplayName', 'Gene Y');
plot(t2, z, 'blue', 'DisplayName', 'Gene Z');
xline(deactivationTime, 'k--', 'LineWidth', 1.2, 'DisplayName', 'deactivation');
title('C1 Feed Forward Loop (AND Gate)');
xlabel('Time');
ylabel('Dynamics of promoter and gene transcription');
legend('Location', 'northeast');
grid on;
hold off;
