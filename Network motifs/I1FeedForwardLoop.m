clc; clear;

% I1 FFL (AND gate) ------------------------------------------------------
% Promoter is activated by Transcription factor (Binary dynamics)
% Promoter initiates Gene Y and Gene Z transcription
% Gene Y crosses threshold Kyz, stars INHIBITING transcription of Gene Z

alpha = 0.5;
beta = 1;
alpha2 = 0.85;
beta2 = 1;
tSpan = [0, 20];
Kyz = 0.6;
z0 = 0;
y0 = 0;
activationTime = 2;
deactivationTime = 8;

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

t = 0:0.01:max(tSpan);
x = @(t) 1 * (t > activationTime & t < deactivationTime);
x_arr = arrayfun(x, t);

dydt = @(t1, y) (x(t1) - alpha * y);
[t1, y] = ode45(dydt, tSpan, y0, options);

y_interpolated = @(t2) interp1(t1, y, t2, 'linear', 'extrap');

dzdt = @(t2, z) (beta * ((x(t2) == 1) & (y_interpolated(t2) < Kyz * max(y))) - alpha2 * z);
[t2, z] = ode45(dzdt, tSpan, z0, options);



figure;
hold on;

plot(t, x_arr, 'black', 'DisplayName', 'Promoter [Binary]');
plot(t1, y, 'red','DisplayName', 'Gene Y');
plot(t2, z, 'blue', 'DisplayName', 'Gene Z');
yline(Kyz * max(y), 'k--', 'Linewidth', 1.2, 'DisplayName', 'y = Kyz threshold');
title('I1 Feed Forward Loop (AND Gate)');
xlabel('Time');
ylabel('Dynamics of promoter and gene transcription');
legend('Location', 'northeast');
grid on;
hold off;
