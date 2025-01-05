clc; clear;

% Single Input Module ------------------------------------------------------
% Promoter is activated by environmental stimulus (Binary)
% Promoter initiates Gene X when its concentration crosses K1
% Promoter initiates Gene Y when its concentration crosses K2
% Promoter initiates Gene Z when its concentration crosses K3
% K3 > K2 > K1

alpha = 0.2;
beta = 1;

alpha1 = 0.3;
alpha2 = 0.6;
alpha3 = 0.85;

tSpan = [0, 20];
K1 = 0.2;
K2 = 0.4;
K3 = 0.6;
y30 = 0;
y20 = 0;
y10 = 0;
x0 = 0;

activationTime = 2;
deactivationTime = 8;

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

t = 0:0.01:max(tSpan);
stimulus = @(t) 1 * (t > activationTime & t < deactivationTime);
stimulus_arr = arrayfun(stimulus, t);

dxdt = @(t_, x) (stimulus(t_) - alpha * x);
[t_, x] = ode45(dxdt, tSpan, x0, options);

x_interpolated1 = @(t1) interp1(t_, x, t1, 'linear', 'extrap');
x_interpolated2 = @(t2) interp1(t_, x, t2, 'linear', 'extrap');
x_interpolated3 = @(t3) interp1(t_, x, t3, 'linear', 'extrap');

dy1_dt = @(t1, y1) (beta * ((x_interpolated1(t1) > K1 * max(x))) - alpha1 * y1);
[t1, y1] = ode45(dy1_dt, tSpan, y10, options);

dy2_dt = @(t2, y2) (beta * ((x_interpolated2(t2) > K2 * max(x))) - alpha2 * y2);
[t2, y2] = ode45(dy2_dt, tSpan, y20, options);

dy3_dt = @(t3, y3) (beta * ((x_interpolated3(t3) > K3 * max(x))) - alpha3 * y3);
[t3, y3] = ode45(dy3_dt, tSpan, y30, options);

figure;
hold on;
ylim([0 2.5]);
plot(t, stimulus_arr, 'DisplayName', 'Stimulus for X');
plot(t_, x/1.5, 'DisplayName', 'Gene X');
plot(t1, y1/max(y1), 'DisplayName', 'Gene Y1', 'LineWidth', 1.3);
plot(t2, y2/max(y2), 'DisplayName', 'Gene Y2', 'LineWidth', 1.3);
plot(t3, y3/max(y3), 'DisplayName', 'Gene Y3', 'LineWidth', 1.3);
title('Single Input Module - Last In First Out (LIFO)');
xlabel('Time');
ylabel('Dynamics of promoter and gene transcription');
legend('Location', 'northeast');
grid on;
hold off;
