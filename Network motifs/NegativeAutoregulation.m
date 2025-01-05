clc; clear;

% NEGATIVE AUTOREGULATION --------------------------------------------------
% Environmental factor (always present) --> Activates Stimulus (X)
% X ----> X*
% X drives response; and response deactivates itself (Negative autoregulation)

alpha = 0.2;
beta = 1;
tSpan = [0, 15];
K = 0.6;
z0 = 0;
activationTime = 2;
deactivationTime = 8;

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

dydt = @(t1, y) (((beta - alpha * y) * (y < K) * (t1 > activationTime & t1 < deactivationTime) - alpha * y * (y > K))) - alpha * y * (t1 > deactivationTime); 
dzdt = @(t2, z) (beta * (t2 > activationTime & t2 < deactivationTime) - alpha * z) * (t2 > activationTime);

[t1, y] = ode45(dydt, tSpan, z0, options);
[t2, z] = ode45(dzdt, tSpan, z0, options);

t = 0:0.01:max(tSpan);
stimulus = @(t) beta * (t > activationTime & t < deactivationTime);
stimulus_arr = arrayfun(stimulus, t);

z_normalized = z/max(z);

figure;
hold on;
ylim([0 1.5]);
plot(t1, y, 'r', 'DisplayName', 'Response (NAR)');
plot(t2, z_normalized, 'b', 'DisplayName', 'Response (without NAR)');
plot(t, stimulus_arr, 'black', 'DisplayName', 'Stimulus (Binary)');
xlabel('Time');
ylabel('Stimulus and Response (Normalized)');
title('Negative Autoregulation');
legend('Location', 'northeast');
grid on;
hold off;

